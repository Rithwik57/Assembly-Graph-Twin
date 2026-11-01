"""Explainable assembly-graph digital twin demo.

This module intentionally uses only the Python standard library. It builds a
small de Bruijn graph, detects common graph signals, scores them with an
explainable baseline, and stores twin snapshots over time.
"""

from __future__ import annotations

import argparse
from array import array
import gzip
import hashlib
import json
import random
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, Iterable, Iterator, TextIO


DEFAULT_BUBBLE_DEPTH = 3
MAX_BUBBLE_DEPTH = 12
DEFAULT_BUBBLE_SEARCH_NODES = 64
MAX_BUBBLE_SEARCH_NODES = 256
DEFAULT_STRAND_MODE = "bidirected"
STRAND_MODES = {"bidirected", "forward"}
_COMPLEMENT = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")


def reverse_complement(sequence: str) -> str:
    """Return the IUPAC DNA reverse complement of a sequence."""
    return sequence.upper().translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Read:
    read_id: str
    sequence: str
    quality: float = 1.0


@dataclass(slots=True)
class GraphNode:
    sequence: str
    incoming: set[str] | list[str] = field(default_factory=set)
    outgoing: set[str] | list[str] = field(default_factory=set)
    supporting_reads: set[str] | None = field(default_factory=set)
    support_count: int = 0
    _last_support_read: str | None = None

    @property
    def coverage(self) -> int:
        return len(self.supporting_reads) if self.supporting_reads is not None else self.support_count


@dataclass(slots=True)
class GraphEdge:
    source: str
    target: str
    kmer: str
    supporting_reads: set[str] | None = field(default_factory=set)
    quality_support: float = 0.0
    support_count: int = 0
    _last_support_read: str | None = None


@dataclass(slots=True)
class AssemblyGraph:
    k: int
    strand_mode: str = DEFAULT_STRAND_MODE
    nodes: dict[str, GraphNode] = field(default_factory=dict)
    edges: dict[tuple[str, str], GraphEdge] = field(default_factory=dict)
    read_node_paths: dict[str, list[str]] = field(default_factory=dict)
    compact: bool = False

    def add_edge(
        self, source: str, target: str, kmer: str, read_id: str, quality: float = 1.0
    ) -> None:
        if source not in self.nodes:
            self.nodes[source] = GraphNode(
                source,
                incoming=[] if self.compact else set(),
                outgoing=[] if self.compact else set(),
                supporting_reads=None if self.compact else set(),
            )
        if target not in self.nodes:
            self.nodes[target] = GraphNode(
                target,
                incoming=[] if self.compact else set(),
                outgoing=[] if self.compact else set(),
                supporting_reads=None if self.compact else set(),
            )
        key = (source, target)
        edge = self.edges.get(key)
        if edge is None:
            edge = GraphEdge(
                source,
                target,
                "" if self.compact else kmer,
                supporting_reads=None if self.compact else set(),
            )
            self.edges[key] = edge
            if self.compact:
                self.nodes[source].outgoing.append(target)
                self.nodes[target].incoming.append(source)
            else:
                self.nodes[source].outgoing.add(target)
                self.nodes[target].incoming.add(source)

        for node in (self.nodes[source], self.nodes[target]):
            if node.supporting_reads is not None:
                previous_count = len(node.supporting_reads)
                node.supporting_reads.add(read_id)
                if len(node.supporting_reads) != previous_count:
                    node.support_count += 1
            elif node._last_support_read != read_id:
                node.support_count += 1
                node._last_support_read = read_id

        if edge.supporting_reads is not None:
            previous_count = len(edge.supporting_reads)
            edge.supporting_reads.add(read_id)
            is_new_support = len(edge.supporting_reads) != previous_count
        else:
            is_new_support = edge._last_support_read != read_id
            if is_new_support:
                edge.support_count += 1
                edge._last_support_read = read_id
        if edge.supporting_reads is not None:
            edge.support_count = len(edge.supporting_reads)
        if is_new_support:
            edge.quality_support += quality


class _CompactNeighbors:
    __slots__ = ("graph", "offsets", "values", "node_id")

    def __init__(self, graph: "CompactAssemblyGraph", offsets: array, values: array, node_id: int):
        self.graph = graph
        self.offsets = offsets
        self.values = values
        self.node_id = node_id

    def __len__(self) -> int:
        return self.offsets[self.node_id + 1] - self.offsets[self.node_id]

    def __iter__(self) -> Iterator[str]:
        for index in range(self.offsets[self.node_id], self.offsets[self.node_id + 1]):
            yield self.graph._sequences[self.values[index]]


class _CompactNodeView:
    __slots__ = ("graph", "node_id")

    def __init__(self, graph: "CompactAssemblyGraph", node_id: int):
        self.graph = graph
        self.node_id = node_id

    @property
    def sequence(self) -> str:
        return self.graph._sequences[self.node_id]

    @property
    def incoming(self) -> _CompactNeighbors:
        return _CompactNeighbors(self.graph, self.graph._incoming_offsets, self.graph._incoming_values, self.node_id)

    @property
    def outgoing(self) -> _CompactNeighbors:
        return _CompactNeighbors(self.graph, self.graph._outgoing_offsets, self.graph._outgoing_values, self.node_id)

    @property
    def coverage(self) -> int:
        return self.graph._node_support[self.node_id]

    @property
    def support_count(self) -> int:
        return self.coverage

    @property
    def supporting_reads(self) -> None:
        return None


class _CompactNodeMap:
    __slots__ = ("graph",)

    def __init__(self, graph: "CompactAssemblyGraph"):
        self.graph = graph

    def __len__(self) -> int:
        return len(self.graph._sequences)

    def __iter__(self) -> Iterator[str]:
        return iter(self.graph._sequences)

    def __getitem__(self, sequence: str) -> _CompactNodeView:
        return _CompactNodeView(self.graph, self.graph._node_ids[sequence])

    def values(self) -> Iterator[_CompactNodeView]:
        for node_id in range(len(self.graph._sequences)):
            yield _CompactNodeView(self.graph, node_id)

    def items(self) -> Iterator[tuple[str, _CompactNodeView]]:
        for node_id, sequence in enumerate(self.graph._sequences):
            yield sequence, _CompactNodeView(self.graph, node_id)


class _CompactEdgeView:
    __slots__ = ("graph", "edge_id")

    def __init__(self, graph: "CompactAssemblyGraph", edge_id: int):
        self.graph = graph
        self.edge_id = edge_id

    @property
    def source(self) -> str:
        return self.graph._sequences[self.graph._edge_sources[self.edge_id]]

    @property
    def target(self) -> str:
        return self.graph._sequences[self.graph._edge_targets[self.edge_id]]

    @property
    def kmer(self) -> str:
        return ""

    @property
    def supporting_reads(self) -> None:
        return None

    @property
    def support_count(self) -> int:
        return self.graph._edge_support[self.edge_id]

    @property
    def quality_support(self) -> float:
        return self.graph._edge_quality[self.edge_id]


class _CompactEdgeMap:
    __slots__ = ("graph",)

    def __init__(self, graph: "CompactAssemblyGraph"):
        self.graph = graph

    def __len__(self) -> int:
        return len(self.graph._edge_sources)

    def __getitem__(self, key: tuple[str, str]) -> _CompactEdgeView:
        source_id = self.graph._node_ids[key[0]]
        target_id = self.graph._node_ids[key[1]]
        return _CompactEdgeView(self.graph, self.graph._edge_ids[(source_id << 32) | target_id])

    def values(self) -> Iterator[_CompactEdgeView]:
        for edge_id in range(len(self.graph._edge_sources)):
            yield _CompactEdgeView(self.graph, edge_id)

    def items(self) -> Iterator[tuple[tuple[str, str], _CompactEdgeView]]:
        for edge_id, (source_id, target_id) in enumerate(zip(self.graph._edge_sources, self.graph._edge_targets)):
            key = (self.graph._sequences[source_id], self.graph._sequences[target_id])
            yield key, _CompactEdgeView(self.graph, edge_id)


class CompactAssemblyGraph:
    """Integer-indexed graph storage for memory-conscious benchmark runs."""

    __slots__ = (
        "k", "compact", "strand_mode", "read_node_paths", "_node_ids", "_sequences", "_node_support",
        "_node_last_read", "_edge_ids", "_edge_sources", "_edge_targets", "_edge_support",
        "_edge_last_read", "_edge_quality", "_incoming_offsets", "_incoming_values",
        "_outgoing_offsets", "_outgoing_values", "nodes", "edges", "_finalized", "backend", "device",
        "gpu_peak_allocated_bytes", "gpu_peak_reserved_bytes",
    )

    def __init__(self, k: int, strand_mode: str = DEFAULT_STRAND_MODE):
        self.k = k
        self.compact = True
        self.strand_mode = strand_mode
        self.backend = "cpu"
        self.device = "cpu"
        self.gpu_peak_allocated_bytes = 0
        self.gpu_peak_reserved_bytes = 0
        self.read_node_paths: dict[str, list[str]] = {}
        self._node_ids: dict[str, int] = {}
        self._sequences: list[str] = []
        self._node_support = array("I")
        self._node_last_read = array("q")
        self._edge_ids: dict[tuple[int, int], int] = {}
        self._edge_sources = array("I")
        self._edge_targets = array("I")
        self._edge_support = array("I")
        self._edge_last_read = array("q")
        self._edge_quality = array("d")
        self._incoming_offsets = array("I", [0])
        self._incoming_values = array("I")
        self._outgoing_offsets = array("I", [0])
        self._outgoing_values = array("I")
        self.nodes = _CompactNodeMap(self)
        self.edges = _CompactEdgeMap(self)
        self._finalized = False

    def _node_id(self, sequence: str) -> int:
        node_id = self._node_ids.get(sequence)
        if node_id is None:
            node_id = len(self._sequences)
            self._node_ids[sequence] = node_id
            self._sequences.append(sequence)
            self._node_support.append(0)
            self._node_last_read.append(-1)
        return node_id

    def add_read_edge(self, source: str, target: str, read_index: int, quality: float) -> None:
        source_id = self._node_id(source)
        target_id = self._node_id(target)
        for node_id in (source_id, target_id):
            if self._node_last_read[node_id] != read_index:
                self._node_support[node_id] += 1
                self._node_last_read[node_id] = read_index
        if source_id >= 2**32 or target_id >= 2**32:
            raise ValueError("compact graph exceeds its 32-bit node index capacity")
        pair = (source_id << 32) | target_id
        edge_id = self._edge_ids.get(pair)
        if edge_id is None:
            edge_id = len(self._edge_sources)
            self._edge_ids[pair] = edge_id
            self._edge_sources.append(source_id)
            self._edge_targets.append(target_id)
            self._edge_support.append(0)
            self._edge_last_read.append(-1)
            self._edge_quality.append(0.0)
        if self._edge_last_read[edge_id] != read_index:
            self._edge_support[edge_id] += 1
            self._edge_last_read[edge_id] = read_index
            self._edge_quality[edge_id] += quality

    def add_aggregated_edge(
        self, source: str, target: str, support: int, quality_support: float,
        source_coverage: int, target_coverage: int,
    ) -> None:
        """Insert one edge with per-read support already aggregated upstream."""
        if self._finalized:
            raise RuntimeError("cannot add edges after compact graph finalization")
        source_id = self._node_ids.get(source)
        if source_id is None:
            source_id = len(self._sequences)
            self._node_ids[source] = source_id
            self._sequences.append(source)
            self._node_support.append(source_coverage)
            self._node_last_read.append(-1)
        target_id = self._node_ids.get(target)
        if target_id is None:
            target_id = len(self._sequences)
            self._node_ids[target] = target_id
            self._sequences.append(target)
            self._node_support.append(target_coverage)
            self._node_last_read.append(-1)
        if source_id >= 2**32 or target_id >= 2**32:
            raise ValueError("compact graph exceeds its 32-bit node index capacity")
        pair = (source_id << 32) | target_id
        edge_id = self._edge_ids.get(pair)
        if edge_id is None:
            edge_id = len(self._edge_sources)
            self._edge_ids[pair] = edge_id
            self._edge_sources.append(source_id)
            self._edge_targets.append(target_id)
            self._edge_support.append(support)
            self._edge_last_read.append(-1)
            self._edge_quality.append(quality_support)
        else:
            self._edge_support[edge_id] += support
            self._edge_quality[edge_id] += quality_support

    def finalize(self) -> None:
        if self._finalized:
            return
        node_count = len(self._sequences)
        edge_count = len(self._edge_sources)
        incoming_counts = array("I", [0]) * node_count
        outgoing_counts = array("I", [0]) * node_count
        for source_id, target_id in zip(self._edge_sources, self._edge_targets):
            outgoing_counts[source_id] += 1
            incoming_counts[target_id] += 1
        self._incoming_offsets = array("I", [0]) * (node_count + 1)
        self._outgoing_offsets = array("I", [0]) * (node_count + 1)
        for node_id in range(node_count):
            self._incoming_offsets[node_id + 1] = self._incoming_offsets[node_id] + incoming_counts[node_id]
            self._outgoing_offsets[node_id + 1] = self._outgoing_offsets[node_id] + outgoing_counts[node_id]
        self._incoming_values = array("I", [0]) * edge_count
        self._outgoing_values = array("I", [0]) * edge_count
        incoming_cursor = array("I", self._incoming_offsets[:-1])
        outgoing_cursor = array("I", self._outgoing_offsets[:-1])
        for source_id, target_id in zip(self._edge_sources, self._edge_targets):
            out_position = outgoing_cursor[source_id]
            self._outgoing_values[out_position] = target_id
            outgoing_cursor[source_id] += 1
            in_position = incoming_cursor[target_id]
            self._incoming_values[in_position] = source_id
            incoming_cursor[target_id] += 1
        self._finalized = True


def extract_kmers(sequence: str, k: int) -> list[str]:
    if k < 2:
        raise ValueError("k must be at least 2")
    if len(sequence) < k:
        return []
    return [sequence[index : index + k] for index in range(len(sequence) - k + 1)]


def build_graph(
    reads: Iterable[Read], k: int, *, compact: bool = False,
    strand_mode: str = DEFAULT_STRAND_MODE,
    backend: str = "cpu", device: str = "auto",
    progress: Callable[[str], None] | None = None,
) -> AssemblyGraph | CompactAssemblyGraph:
    if k < 2:
        raise ValueError("k must be at least 2")
    if strand_mode not in STRAND_MODES:
        raise ValueError(f"strand_mode must be one of {sorted(STRAND_MODES)}")
    if backend not in {"cpu", "torch"}:
        raise ValueError("backend must be cpu or torch")
    if backend == "torch":
        if not compact:
            raise ValueError("the torch backend currently requires compact=True")
        from gpu_backend import build_compact_graph_torch

        graph, actual_device = build_compact_graph_torch(
            reads, k, strand_mode=strand_mode, device=device, progress=progress
        )
        graph.backend = "torch"
        graph.device = actual_device
        return graph
    if compact:
        graph = CompactAssemblyGraph(k, strand_mode=strand_mode)
    else:
        graph = AssemblyGraph(k=k, strand_mode=strand_mode)
    compact_read_ids: set[str] = set()
    for read_index, read in enumerate(reads):
        if compact:
            if read.read_id in compact_read_ids:
                raise ValueError("compact graph mode requires unique read IDs")
            compact_read_ids.add(read.read_id)
        sequence = read.sequence.upper()
        oriented_sequences = [sequence]
        if strand_mode == "bidirected":
            reverse = reverse_complement(sequence)
            if reverse != sequence:
                oriented_sequences.append(reverse)
        if compact:
            transitions: set[tuple[str, str]] = set()
            for oriented in oriented_sequences:
                for index in range(max(0, len(oriented) - k + 1)):
                    transitions.add((oriented[index : index + k - 1], oriented[index + 1 : index + k]))
            for source, target in sorted(transitions):
                graph.add_read_edge(source, target, read_index, read.quality)
            continue
        kmers = extract_kmers(sequence, k)
        if kmers:
            graph.read_node_paths[read.read_id] = [kmers[0][:-1]] + [
                kmer[1:] for kmer in kmers
            ]
        else:
            graph.read_node_paths[read.read_id] = []
        observed_edges: set[tuple[str, str, str]] = set()
        for oriented in oriented_sequences:
            for kmer in extract_kmers(oriented, k):
                observed_edges.add((kmer[:-1], kmer[1:], kmer))
        for source, target, kmer in sorted(observed_edges):
            graph.add_edge(source, target, kmer, read.read_id, read.quality)
    if compact:
        graph.finalize()
    return graph


def graph_metrics(graph: AssemblyGraph | CompactAssemblyGraph) -> dict[str, int | str]:
    if isinstance(graph, CompactAssemblyGraph):
        component_count = len(_compact_component_summaries(graph))
        if graph.strand_mode == "bidirected":
            canonical_node_count = sum(
                sequence <= reverse_complement(sequence) for sequence in graph._sequences
            )
            canonical_edge_count = sum(
                (graph._sequences[source], graph._sequences[target])
                <= (reverse_complement(graph._sequences[target]), reverse_complement(graph._sequences[source]))
                for source, target in zip(graph._edge_sources, graph._edge_targets)
            )
        else:
            canonical_node_count = len(graph.nodes)
            canonical_edge_count = len(graph.edges)
    else:
        components = connected_components(graph)
        if graph.strand_mode == "bidirected":
            components = _strand_collapsed_components(graph, components)
            canonical_node_count = sum(
                sequence <= reverse_complement(sequence) for sequence in graph.nodes
            )
            canonical_edge_count = sum(
                (source, target)
                <= (reverse_complement(target), reverse_complement(source))
                for source, target in graph.edges
            )
        else:
            canonical_node_count = len(graph.nodes)
            canonical_edge_count = len(graph.edges)
        component_count = len(components)
    branching_keys = {
        min(node.sequence, reverse_complement(node.sequence))
        if graph.strand_mode == "bidirected" else node.sequence
        for node in graph.nodes.values()
        if len(node.incoming) > 1 or len(node.outgoing) > 1
    }
    return {
        "k": graph.k,
        "strand_mode": graph.strand_mode,
        "node_count": len(graph.nodes),
        "edge_count": len(graph.edges),
        "strand_collapsed_node_count": canonical_node_count,
        "strand_collapsed_edge_count": canonical_edge_count,
        "component_count": component_count,
        "branching_node_count": len(branching_keys),
    }


def _compact_component_summaries(graph: CompactAssemblyGraph) -> list[dict[str, object]]:
    """Summarize weak components without materializing Python sets of every node."""
    node_count = len(graph._sequences)
    component_of = array("i", [-1]) * node_count
    summaries: list[dict[str, object]] = []
    for start_id in range(node_count):
        if component_of[start_id] >= 0:
            continue
        component_id = len(summaries)
        component_of[start_id] = component_id
        stack = [start_id]
        count = 0
        smallest = graph._sequences[start_id]
        smallest_reverse = reverse_complement(smallest)
        while stack:
            node_id = stack.pop()
            count += 1
            sequence = graph._sequences[node_id]
            if sequence < smallest:
                smallest = sequence
            reverse = reverse_complement(sequence)
            if reverse < smallest_reverse:
                smallest_reverse = reverse
            for values, offsets in (
                (graph._incoming_values, graph._incoming_offsets),
                (graph._outgoing_values, graph._outgoing_offsets),
            ):
                for index in range(offsets[node_id], offsets[node_id + 1]):
                    neighbor_id = values[index]
                    if component_of[neighbor_id] < 0:
                        component_of[neighbor_id] = component_id
                        stack.append(neighbor_id)
        summaries.append({
            "component": component_id + 1,
            "min_node": smallest,
            "min_reverse_complement_node": smallest_reverse,
            "node_count": count,
            "edge_count": 0,
        })
    for source_id, target_id in zip(graph._edge_sources, graph._edge_targets):
        component_id = component_of[source_id]
        summaries[component_id]["edge_count"] += 1
    if graph.strand_mode == "bidirected":
        summaries = [
            row for row in summaries
            if row["min_node"] <= row["min_reverse_complement_node"]
        ]
    return sorted(summaries, key=lambda row: row["min_node"])


def connected_components(graph: AssemblyGraph) -> list[set[str]]:
    remaining = set(graph.nodes)
    components: list[set[str]] = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        component = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            node = graph.nodes[current]
            for neighbor in (*node.incoming, *node.outgoing):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return sorted(components, key=lambda component: min(component))


def _strand_collapsed_components(
    graph: AssemblyGraph, components: list[set[str]] | None = None
) -> list[set[str]]:
    components = components if components is not None else connected_components(graph)
    if graph.strand_mode != "bidirected":
        return components
    representatives: dict[str, set[str]] = {}
    for component in components:
        forward_min = min(component)
        reverse_min = min(reverse_complement(node) for node in component)
        representatives.setdefault(min(forward_min, reverse_min), component)
    return sorted(representatives.values(), key=lambda component: min(component))


def detect_disconnected_components(graph: AssemblyGraph) -> list[dict[str, object]]:
    components = _strand_collapsed_components(graph)
    if len(components) < 2:
        return []

    diagnostics = []
    for index, component in enumerate(components, start=1):
        component_edges = [
            edge
            for edge in graph.edges.values()
            if edge.source in component and edge.target in component
        ]
        diagnostics.append(
            {
                "component": index,
                "nodes": sorted(component),
                "node_count": len(component),
                "edge_count": len(component_edges),
                "component_count": len(components),
                "reason": "weakly disconnected graph component may indicate a coverage gap or separate fragment",
            }
        )
    return diagnostics


def detect_tips(graph: AssemblyGraph) -> list[dict[str, object]]:
    tips: dict[str, dict[str, object]] = {}
    for node in graph.nodes.values():
        terminal = not node.incoming or not node.outgoing
        if terminal and node.coverage <= 1:
            key = (
                min(node.sequence, reverse_complement(node.sequence))
                if graph.strand_mode == "bidirected" else node.sequence
            )
            tips.setdefault(
                key,
                {
                    "node": key,
                    "coverage": node.coverage,
                    "reason": "terminal node has support from at most one read",
                }
            )
    return sorted(tips.values(), key=lambda item: item["node"])


def detect_repeat_candidates(graph: AssemblyGraph) -> list[dict[str, object]]:
    candidates: dict[str, dict[str, object]] = {}
    for node in graph.nodes.values():
        high_degree = len(node.incoming) > 1 or len(node.outgoing) > 1
        high_support = node.coverage >= 3
        if high_degree or high_support:
            key = (
                min(node.sequence, reverse_complement(node.sequence))
                if graph.strand_mode == "bidirected" else node.sequence
            )
            candidates.setdefault(
                key,
                {
                    "node": key,
                    "coverage": node.coverage,
                    "in_degree": len(node.incoming),
                    "out_degree": len(node.outgoing),
                    "reason": "high support or multiple graph connections",
                }
            )
    return sorted(candidates.values(), key=lambda item: item["node"])


def reachable_within(
    graph: AssemblyGraph,
    start: str,
    depth: int,
    blocked: set[str] | None = None,
) -> set[str]:
    reached: set[str] = set()
    visited = {start}
    blocked_nodes = blocked or set()
    frontier = {start}
    for _ in range(depth):
        next_frontier = set()
        for node_name in frontier:
            next_frontier.update(
                neighbor
                for neighbor in graph.nodes[node_name].outgoing
                if neighbor not in blocked_nodes and neighbor not in visited
            )
        visited.update(next_frontier)
        reached.update(next_frontier)
        frontier = next_frontier
        if not frontier:
            break
    return reached


def shortest_paths_within(
    graph: AssemblyGraph,
    start: str,
    depth: int,
    max_nodes: int = DEFAULT_BUBBLE_SEARCH_NODES,
    blocked: set[str] | None = None,
) -> tuple[dict[str, int], dict[str, str], bool, str | None]:
    """Return bounded paths, predecessors, and whether node expansion was limited."""
    distances = {start: 0}
    predecessors: dict[str, str] = {}
    frontier = [start]
    blocked_nodes = blocked or set()
    limited = False
    limit_reason: str | None = None
    for distance in range(1, depth + 1):
        next_frontier = []
        for current in frontier:
            for neighbor in sorted(graph.nodes[current].outgoing):
                if neighbor in blocked_nodes or neighbor in distances:
                    continue
                if len(distances) >= max_nodes:
                    limited = True
                    limit_reason = "node_cap"
                    break
                distances[neighbor] = distance
                predecessors[neighbor] = current
                next_frontier.append(neighbor)
            if limited:
                break
        frontier = next_frontier
        if not frontier:
            break
    if frontier and distance >= depth:
        limited = True
        if limit_reason is None:
            limit_reason = "depth"
    return distances, predecessors, limited, limit_reason


def reconstruct_path(
    start: str, end: str, predecessors: dict[str, str]
) -> list[str]:
    path = [end]
    while path[-1] != start:
        path.append(predecessors[path[-1]])
    return list(reversed(path))


def reads_supporting_path(graph: AssemblyGraph, path: list[str]) -> list[str]:
    """Find reads whose observed node walk contains the complete path contiguously."""
    if len(path) < 2:
        return []
    width = len(path)
    reverse_path = [reverse_complement(node) for node in reversed(path)]
    return sorted(
        read_id
        for read_id, read_path in graph.read_node_paths.items()
        if any(
            read_path[index : index + width] == path
            or (graph.strand_mode == "bidirected" and read_path[index : index + width] == reverse_path)
            for index in range(len(read_path) - width + 1)
        )
    )


def _validate_bubble_depth(depth: int) -> None:
    if not 1 <= depth <= MAX_BUBBLE_DEPTH:
        raise ValueError(
            f"bubble depth must be between 1 and {MAX_BUBBLE_DEPTH}"
        )


def _validate_bubble_search_nodes(max_nodes: int) -> None:
    if not 2 <= max_nodes <= MAX_BUBBLE_SEARCH_NODES:
        raise ValueError(
            f"bubble search node limit must be between 2 and {MAX_BUBBLE_SEARCH_NODES}"
        )


def _candidate_id(start: str, rejoin: str, paths: dict[str, list[str]]) -> str:
    forward_topology = {
        "start": start,
        "rejoin": rejoin,
        "paths": sorted(tuple(path) for path in paths.values()),
    }
    reverse_topology = {
        "start": reverse_complement(rejoin),
        "rejoin": reverse_complement(start),
        "paths": sorted(
            tuple(reverse_complement(node) for node in reversed(path))
            for path in paths.values()
        ),
    }
    encoded = min(
        json.dumps(topology, separators=(",", ":"), sort_keys=True)
        for topology in (forward_topology, reverse_topology)
    )
    return "bubble-" + hashlib.sha256(encoded.encode("ascii")).hexdigest()[:16]


def _path_edges(path: list[str]) -> set[tuple[str, str]]:
    return set(zip(path, path[1:]))


def _classify_relationship(
    first: dict[str, object], second: dict[str, object]
) -> tuple[str, list[str], list[str]]:
    first_nodes = set(first["paths"][first["branches"][0]])
    second_nodes = set(second["paths"][second["branches"][0]])
    first_internal = set().union(*(set(path[1:-1]) for path in first["paths"].values()))
    second_internal = set().union(*(set(path[1:-1]) for path in second["paths"].values()))
    shared_nodes = sorted(first_internal & second_internal)
    shared_edges = sorted(
        "->".join(edge)
        for edge in set().union(*(_path_edges(path) for path in first["paths"].values()))
        & set().union(*(_path_edges(path) for path in second["paths"].values()))
    )
    first_contains_second = second["start"] in first_internal and second["rejoin"] in first_nodes
    second_contains_first = first["start"] in second_internal and first["rejoin"] in second_nodes
    if first_contains_second or second_contains_first:
        relationship_type = "nested"
    else:
        relationship_type = "overlapping"
    return relationship_type, shared_nodes, shared_edges


def _add_candidate_relationships(bubbles: list[dict[str, object]]) -> None:
    for bubble in bubbles:
        bubble["relationships"] = []
    for index, first in enumerate(bubbles):
        for second in bubbles[index + 1 :]:
            relationship_type, shared_nodes, shared_edges = _classify_relationship(first, second)
            if not shared_nodes and not shared_edges:
                continue
            first["relationships"].append(
                {
                    "candidate_id": second["candidate_id"],
                    "type": relationship_type,
                    "shared_internal_nodes": shared_nodes,
                    "shared_edges": shared_edges,
                }
            )
            second["relationships"].append(
                {
                    "candidate_id": first["candidate_id"],
                    "type": relationship_type,
                    "shared_internal_nodes": shared_nodes,
                    "shared_edges": shared_edges,
                }
            )
    for bubble in bubbles:
        bubble["relationships"].sort(key=lambda relation: relation["candidate_id"])


def detect_bubbles_with_limits(
    graph: AssemblyGraph,
    depth: int = DEFAULT_BUBBLE_DEPTH,
    max_nodes: int = DEFAULT_BUBBLE_SEARCH_NODES,
    progress: Callable[[str], None] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    _validate_bubble_depth(depth)
    _validate_bubble_search_nodes(max_nodes)
    candidate_pairs: list[tuple[str, str, str, list[str], list[str], bool, bool]] = []
    limited_searches: list[dict[str, object]] = []
    # Only divergence nodes can start a bubble. Filtering before sorting avoids
    # materializing/sorting every node in large sparse graphs.
    branch_nodes = sorted(
        (node for node in graph.nodes.values() if len(node.outgoing) > 1),
        key=lambda item: item.sequence,
    )
    for branch_index, node in enumerate(branch_nodes, start=1):
        if progress is not None and branch_index % 100 == 0:
            progress(f"Bubble search: checked {branch_index:,}/{len(branch_nodes):,} divergence nodes")
        successors = sorted(successor for successor in node.outgoing if successor != node.sequence)
        if len(successors) < 2:
            continue
        for index, left in enumerate(successors):
            for right in successors[index + 1 :]:
                left_distances, left_predecessors, left_limited, left_limit_reason = shortest_paths_within(
                    graph, left, depth=depth, max_nodes=max_nodes, blocked={node.sequence, right}
                )
                right_distances, right_predecessors, right_limited, right_limit_reason = shortest_paths_within(
                    graph, right, depth=depth, max_nodes=max_nodes, blocked={node.sequence, left}
                )
                common = set(left_distances) & set(right_distances)
                common.discard(left)
                common.discard(right)
                ordered_common = sorted(
                    common,
                    key=lambda candidate: (
                        max(left_distances[candidate], right_distances[candidate]),
                        left_distances[candidate] + right_distances[candidate],
                        candidate,
                    ),
                )
                for rejoin in ordered_common:
                    left_path = [node.sequence] + reconstruct_path(
                        left, rejoin, left_predecessors
                    )
                    right_path = [node.sequence] + reconstruct_path(
                        right, rejoin, right_predecessors
                    )
                    left_internal = set(left_path[1:-1])
                    right_internal = set(right_path[1:-1])
                    if left_internal & right_internal:
                        continue
                    candidate_pairs.append(
                        (
                            node.sequence,
                            left,
                            right,
                            left_path,
                            right_path,
                            left_limited,
                            right_limited,
                        )
                    )
                    break
                else:
                    if left_limited or right_limited:
                        limited_searches.append(
                            {
                                "start": node.sequence,
                                "branches": [left, right],
                                "search_depth": depth,
                                "max_nodes": max_nodes,
                                "branch_search_limited": {
                                    left: left_limited,
                                    right: right_limited,
                                },
                                "branch_limit_reason": {
                                    left: left_limit_reason,
                                    right: right_limit_reason,
                                },
                                "reason": "bounded search stopped before a complete reconvergence was found; absence of a bubble is not established",
                            }
                        )

    if progress is not None:
        progress(f"Bubble search complete: checked {len(branch_nodes):,} divergence nodes")

    grouped: dict[tuple[str, str], dict[str, object]] = {}
    for start, left, right, left_path, right_path, left_limited, right_limited in candidate_pairs:
        rejoin = left_path[-1]
        key = (start, rejoin)
        bubble = grouped.setdefault(
            key,
            {
                "start": start,
                "branches": [],
                "rejoin": rejoin,
                "paths": {},
                "path_read_support": {},
                "branch_support": {},
                "search_depth": depth,
                "max_nodes": max_nodes,
                "branch_search_limited": {},
                "reason": "two or more outgoing paths reconverge nearby",
            },
        )
        paths = bubble["paths"]
        branches = bubble["branches"]
        path_read_support = bubble["path_read_support"]
        branch_support = bubble["branch_support"]
        for branch, path in ((left, left_path), (right, right_path)):
            if branch in paths:
                continue
            internal = set(path[1:-1])
            if any(internal & set(existing[1:-1]) for existing in paths.values()):
                continue
            paths[branch] = path
            branches.append(branch)
            path_read_support[branch] = reads_supporting_path(graph, path)
            branch_support[branch] = graph.edges[(start, branch)].support_count
            bubble["branch_search_limited"][branch] = (
                left_limited if branch == left else right_limited
            )

    bubbles = [bubble for bubble in grouped.values() if len(bubble["branches"]) >= 2]
    for bubble in bubbles:
        bubble["branches"] = sorted(bubble["branches"])
        bubble["paths"] = {
            branch: bubble["paths"][branch] for branch in bubble["branches"]
        }
        bubble["path_read_support"] = {
            branch: bubble["path_read_support"][branch]
            for branch in bubble["branches"]
        }
        bubble["branch_support"] = {
            branch: bubble["branch_support"][branch]
            for branch in bubble["branches"]
        }
        bubble["branch_search_limited"] = {
            branch: bubble["branch_search_limited"][branch]
            for branch in bubble["branches"]
        }
        bubble["candidate_id"] = _candidate_id(
            bubble["start"], bubble["rejoin"], bubble["paths"]
        )
    bubbles = sorted(
        bubbles,
        key=lambda item: (item["start"], item["rejoin"], item["candidate_id"]),
    )
    # The directed double-cover contains both a path and its reverse-complement
    # view. Keep one canonical candidate record for that physical structure.
    unique_bubbles: dict[str, dict[str, object]] = {}
    for bubble in bubbles:
        unique_bubbles.setdefault(str(bubble["candidate_id"]), bubble)
    bubbles = list(unique_bubbles.values())
    _add_candidate_relationships(bubbles)
    limited_searches = sorted(
        limited_searches,
        key=lambda item: (item["start"], item["branches"]),
    )
    return bubbles, limited_searches


def detect_bubbles(
    graph: AssemblyGraph,
    depth: int = DEFAULT_BUBBLE_DEPTH,
    max_nodes: int = DEFAULT_BUBBLE_SEARCH_NODES,
) -> list[dict[str, object]]:
    """Return complete bubble candidates while preserving the original API."""
    bubbles, _ = detect_bubbles_with_limits(graph, depth=depth, max_nodes=max_nodes)
    return bubbles


def classify_diagnostics(
    graph: AssemblyGraph,
    bubble_depth: int = DEFAULT_BUBBLE_DEPTH,
    max_bubble_search_nodes: int = DEFAULT_BUBBLE_SEARCH_NODES,
) -> list[dict[str, object]]:
    tips = detect_tips(graph)
    repeats = detect_repeat_candidates(graph)
    bubbles, limited_searches = detect_bubbles_with_limits(
        graph, depth=bubble_depth, max_nodes=max_bubble_search_nodes
    )
    disconnected_components = detect_disconnected_components(graph)
    results = []

    for tip in tips:
        results.append(
            {
                "type": "tip",
                "label": "error",
                "heuristic_score": 0.86,
                "evidence": tip,
                "reason": "Low-support terminal graph structure is suspicious for a sequencing error.",
            }
        )

    for repeat in repeats:
        heuristic_score = (
            0.72 if repeat["in_degree"] > 1 or repeat["out_degree"] > 1 else 0.65
        )
        results.append(
            {
                "type": "repeat_candidate",
                "label": "repeat_or_branch",
                "heuristic_score": heuristic_score,
                "evidence": repeat,
                "reason": "High support or branching can indicate repeated sequence or unresolved structure.",
            }
        )

    for bubble in bubbles:
        minimum_branch_support = min(bubble["branch_support"].values())
        heuristic_score = min(0.95, 0.5 + 0.1 * minimum_branch_support)
        results.append(
            {
                "type": "bubble",
                "label": "true_variant_or_unresolved",
                "heuristic_score": heuristic_score,
                "evidence": bubble,
                "reason": "Reconverging paths may represent a real variant, haplotypes, or insufficient evidence; the score is a structural heuristic, not a probability.",
            }
        )

    for limited_search in limited_searches:
        results.append(
            {
                "type": "bubble_search_limited",
                "label": "search_limited",
                "heuristic_score": 0.0,
                "evidence": limited_search,
                "reason": "The bounded graph search was truncated; no complete bubble conclusion is made.",
            }
        )

    for component in disconnected_components:
        results.append(
            {
                "type": "disconnected_component",
                "label": "coverage_gap_or_separate_fragment",
                "heuristic_score": 0.5,
                "evidence": component,
                "reason": "Disconnected components are reported as a structural signal; this does not establish a coverage gap or biological separation.",
            }
        )

    return results


def read_fingerprint(reads: list[Read]) -> str:
    content = "\n".join(
        f"{read.read_id}\t{read.sequence.upper()}\t{read.quality:.8f}"
        for read in sorted(reads, key=lambda item: item.read_id)
    )
    return hashlib.sha256(content.encode("ascii")).hexdigest()


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def _read_id(header: str, marker: str, source: Path, line_number: int) -> str:
    if not header.startswith(marker):
        raise ValueError(f"{source}:{line_number}: expected a {marker!r} record header")
    fields = header[1:].strip().split()
    if not fields:
        raise ValueError(f"{source}:{line_number}: record header has no identifier")
    return fields[0]


def _validate_sequence(sequence: str, source: Path, read_id: str) -> str:
    normalized = sequence.upper()
    invalid = sorted(set(normalized) - set("ACGTN"))
    if invalid:
        raise ValueError(
            f"{source}: read {read_id!r} contains unsupported bases: {''.join(invalid)}"
        )
    return normalized


def iter_parse_fasta(stream: TextIO, source: Path) -> Iterable[Read]:
    current_id: str | None = None
    sequence_parts: list[str] = []
    header_line = 0
    records = 0
    for line_number, raw_line in enumerate(stream, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequence = _validate_sequence("".join(sequence_parts), source, current_id)
                if not sequence:
                    raise ValueError(f"{source}:{header_line}: read {current_id!r} is empty")
                records += 1
                yield Read(current_id, sequence)
            current_id = _read_id(line, ">", source, line_number)
            sequence_parts = []
            header_line = line_number
        else:
            if current_id is None:
                raise ValueError(f"{source}:{line_number}: sequence appears before a FASTA header")
            sequence_parts.append(line)
    if current_id is not None:
        sequence = _validate_sequence("".join(sequence_parts), source, current_id)
        if not sequence:
            raise ValueError(f"{source}:{header_line}: read {current_id!r} is empty")
        records += 1
        yield Read(current_id, sequence)
    if not records:
        raise ValueError(f"{source}: no FASTA reads found")


def parse_fasta(stream: TextIO, source: Path) -> list[Read]:
    return list(iter_parse_fasta(stream, source))


def iter_parse_fastq(stream: TextIO, source: Path) -> Iterable[Read]:
    """Parse standard four-line FASTQ records (Phred+33 qualities)."""
    line_number = 0
    records = 0
    while True:
        header = stream.readline()
        if not header:
            break
        line_number += 1
        if not header.strip():
            continue
        read_id = _read_id(header.strip(), "@", source, line_number)
        sequence_line = stream.readline()
        plus_line = stream.readline()
        quality_line = stream.readline()
        if not sequence_line or not plus_line or not quality_line:
            raise ValueError(f"{source}:{line_number}: incomplete four-line FASTQ record")
        line_number += 3
        sequence = _validate_sequence(sequence_line.strip(), source, read_id)
        quality_text = quality_line.rstrip("\r\n")
        if not sequence:
            raise ValueError(f"{source}:{line_number - 2}: read {read_id!r} is empty")
        if not plus_line.startswith("+"):
            raise ValueError(f"{source}:{line_number - 1}: expected '+' separator")
        if len(quality_text) != len(sequence):
            raise ValueError(
                f"{source}:{line_number}: quality length does not match read {read_id!r}"
            )
        phred_scores = [ord(char) - 33 for char in quality_text]
        if any(score < 0 for score in phred_scores):
            raise ValueError(f"{source}:{line_number}: invalid Phred+33 quality character")
        mean_error_probability = sum(10 ** (-score / 10) for score in phred_scores) / len(phred_scores)
        records += 1
        yield Read(read_id, sequence, quality=1.0 - mean_error_probability)
    if not records:
        raise ValueError(f"{source}: no FASTQ reads found")


def parse_fastq(stream: TextIO, source: Path) -> list[Read]:
    return list(iter_parse_fastq(stream, source))


def iter_load_reads(path: Path) -> Iterable[Read]:
    """Load FASTA or four-line FASTQ, optionally gzip-compressed."""
    with _open_text(path) as stream:
        first_content = ""
        for line in stream:
            if line.strip():
                first_content = line
                break
        if not first_content:
            raise ValueError(f"{path}: input is empty")
        stream.seek(0)
        if first_content.startswith(">"):
            yield from iter_parse_fasta(stream, path)
            return
        if first_content.startswith("@"):
            yield from iter_parse_fastq(stream, path)
            return
        raise ValueError(f"{path}: expected FASTA ('>') or FASTQ ('@') input")


def load_reads(path: Path) -> list[Read]:
    return list(iter_load_reads(path))


def serialize_graph(graph: AssemblyGraph) -> dict[str, object]:
    return {
        "strand_mode": graph.strand_mode,
        "k": graph.k,
        "nodes": [
            {
                "sequence": node.sequence,
                "coverage": node.coverage,
                "incoming": sorted(node.incoming),
                "outgoing": sorted(node.outgoing),
                "supporting_reads": sorted(node.supporting_reads or ()),
                "support_count": node.support_count,
            }
            for node in sorted(graph.nodes.values(), key=lambda item: item.sequence)
        ],
        "edges": [
            {
                "source": edge.source,
                "target": edge.target,
                "kmer": edge.kmer,
                "supporting_reads": sorted(edge.supporting_reads or ()),
                "support_count": edge.support_count,
                "quality_support": round(edge.quality_support, 6),
            }
            for edge in sorted(graph.edges.values(), key=lambda item: (item.source, item.target))
        ],
    }


class AssemblyGraphTwin:
    def __init__(
        self,
        k: int = 3,
        bubble_depth: int = DEFAULT_BUBBLE_DEPTH,
        max_bubble_search_nodes: int = DEFAULT_BUBBLE_SEARCH_NODES,
        strand_mode: str = DEFAULT_STRAND_MODE,
    ) -> None:
        if k < 2:
            raise ValueError("k must be at least 2")
        _validate_bubble_depth(bubble_depth)
        _validate_bubble_search_nodes(max_bubble_search_nodes)
        if strand_mode not in STRAND_MODES:
            raise ValueError(f"strand_mode must be one of {sorted(STRAND_MODES)}")
        self.k = k
        self.strand_mode = strand_mode
        self.bubble_depth = bubble_depth
        self.max_bubble_search_nodes = max_bubble_search_nodes
        self.snapshots: list[dict[str, object]] = []
        self._read_ids: set[str] = set()

    def update(
        self,
        reads: list[Read],
        label: str,
        source: str = "synthetic",
        timestamp_utc: str | None = None,
    ) -> dict[str, object]:
        read_ids = [read.read_id for read in reads]
        if len(read_ids) != len(set(read_ids)):
            raise ValueError("read IDs must be unique within a snapshot")
        if any(not 0.0 <= read.quality <= 1.0 for read in reads):
            raise ValueError("read quality must be an accuracy value between 0 and 1")
        graph = build_graph(reads, self.k, strand_mode=self.strand_mode)
        long_read_count = sum(len(read.sequence) >= 8 for read in reads)
        diagnostics = classify_diagnostics(
            graph,
            bubble_depth=self.bubble_depth,
            max_bubble_search_nodes=self.max_bubble_search_nodes,
        )
        previous_graph = self.snapshots[-1]["graph"] if self.snapshots else {"nodes": [], "edges": []}
        previous_nodes = {node["sequence"] for node in previous_graph["nodes"]}
        previous_edges = {(edge["source"], edge["target"]): edge for edge in previous_graph["edges"]}
        current_graph = serialize_graph(graph)
        graph_sha256 = hashlib.sha256(
            json.dumps(current_graph, sort_keys=True, separators=(",", ":")).encode("ascii")
        ).hexdigest()
        current_edges = {(edge["source"], edge["target"]): edge for edge in current_graph["edges"]}
        current_read_ids = set(read_ids)
        added_read_ids = sorted(current_read_ids - self._read_ids)
        edge_support_increases = []
        for edge_key in sorted(set(current_edges) & set(previous_edges)):
            new_support = sorted(set(current_edges[edge_key]["supporting_reads"]) - set(previous_edges[edge_key]["supporting_reads"]))
            if new_support:
                edge_support_increases.append({"source": edge_key[0], "target": edge_key[1], "new_read_ids": new_support})
        snapshot = {
            "snapshot_id": len(self.snapshots) + 1,
            "label": label,
            "source": source,
            "timestamp_utc": timestamp_utc or datetime.now(timezone.utc).isoformat(),
            "input_sha256": read_fingerprint(reads),
            "batch_sha256": read_fingerprint([read for read in reads if read.read_id in set(added_read_ids)]),
            "graph_sha256": graph_sha256,
            "read_count": len(reads),
            "batch_read_count": len(added_read_ids),
            "state_delta": {
                "added_read_ids": added_read_ids,
                "added_nodes": sorted(set(graph.nodes) - previous_nodes),
                "added_edges": [list(edge) for edge in sorted(set(current_edges) - set(previous_edges))],
                "edge_support_increases": edge_support_increases,
            },
            "long_read_count": long_read_count,
            "bubble_depth": self.bubble_depth,
            "max_bubble_search_nodes": self.max_bubble_search_nodes,
            "strand_mode": self.strand_mode,
            "metrics": graph_metrics(graph),
            "diagnostics": diagnostics,
            "graph": current_graph,
        }
        self.snapshots.append(snapshot)
        self._read_ids = current_read_ids
        return snapshot


def initial_reads() -> list[Read]:
    return [
        Read("main_1", "AAAC"),
        Read("main_2", "AACG"),
        Read("main_3", "ACGT"),
        Read("main_4", "CGTT"),
        Read("main_5", "GTTA"),
        Read("main_6", "TTAA"),
        Read("variant_1", "AACC"),
        Read("variant_2", "ACCT"),
        Read("variant_3", "CCTT"),
        Read("variant_4", "CTTA"),
        Read("repeat_support_1", "CGGA"),
        Read("repeat_support_2", "CGGA"),
        Read("error_tip", "AACA", quality=0.35),
    ]


def bridge_reads() -> list[Read]:
    return [
        Read("bridge_reference", "AAACGTTAA", quality=0.98),
        Read("bridge_variant", "AAACCTTAA", quality=0.98),
    ]


def summarize_snapshot(snapshot: dict[str, object]) -> str:
    metrics = snapshot["metrics"]
    diagnostics = snapshot["diagnostics"]
    counts: dict[str, int] = defaultdict(int)
    for diagnostic in diagnostics:
        counts[str(diagnostic["label"])] += 1
    return (
        f"Snapshot {snapshot['snapshot_id']} ({snapshot['label']}): "
        f"reads={snapshot['read_count']}, nodes={metrics['node_count']}, "
        f"edges={metrics['edge_count']}, components={metrics['component_count']}, "
        f"diagnostics={dict(counts)}"
    )


def run_demo(
    output: Path | None = None,
    bubble_depth: int = DEFAULT_BUBBLE_DEPTH,
) -> dict[str, object]:
    twin = AssemblyGraphTwin(k=3, bubble_depth=bubble_depth)
    first_reads = initial_reads()
    first = twin.update(first_reads, "initial reads")
    second = twin.update(first_reads + bridge_reads(), "after bridge reads")
    first_bubble = next(
        (item for item in first["diagnostics"] if item["type"] == "bubble"),
        None,
    )
    second_bubble = next(
        (item for item in second["diagnostics"] if item["type"] == "bubble"),
        None,
    )
    snapshots = [
        {key: value for key, value in item.items() if key != "_read_ids"}
        for item in twin.snapshots
    ]
    example_bubble = (second_bubble or first_bubble or {}).get("evidence")
    report = {
        "demonstration": {
            "type": "synthetic_positive_control",
            "represents_real_sample_data": False,
            "bubble_count_by_snapshot": {
                "initial_reads": sum(item["type"] == "bubble" for item in first["diagnostics"]),
                "after_bridge_reads": sum(item["type"] == "bubble" for item in second["diagnostics"]),
            },
            "example_bubble": example_bubble,
        },
        "snapshots": snapshots,
    }

    print("Assembly graph twin demo (synthetic reads)")
    print(summarize_snapshot(first))
    print(summarize_snapshot(second))
    if first_bubble and second_bubble:
        print(
            "Bubble heuristic scores: "
            f"{first_bubble['heuristic_score']:.2f}, {second_bubble['heuristic_score']:.2f}"
        )
    if example_bubble:
        print("Controlled synthetic bubble (not a result from the real-read benchmark):")
        print(
            f"- candidate={example_bubble['candidate_id']} "
            f"divergence={example_bubble['start']} rejoin={example_bubble['rejoin']}"
        )
        for branch in example_bubble["branches"]:
            path = example_bubble["paths"][branch]
            path_sequence = path[0] + "".join(node[-1] for node in path[1:])
            read_support = len(example_bubble["path_read_support"].get(branch, []))
            print(
                f"  path={path_sequence} branch_edge_support={example_bubble['branch_support'][branch]} "
                f"complete_path_reads={read_support}"
            )
    else:
        print("No complete synthetic bubble was found at this search depth.")
    print("Decision examples:")
    for diagnostic in first["diagnostics"][:5]:
        print(
            f"- {diagnostic['type']}: {diagnostic['label']} "
            f"(heuristic_score={diagnostic['heuristic_score']:.2f}) - {diagnostic['reason']}"
        )
    print("The twin preserves graph evidence and records the update as snapshot 2.")

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"Wrote report to {output}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    demo_parser = subparsers.add_parser("demo", help="run the synthetic twin demo")
    demo_parser.add_argument("--output", type=Path, help="write the JSON report to this path")
    demo_parser.add_argument(
        "--bubble-depth",
        type=int,
        default=DEFAULT_BUBBLE_DEPTH,
        help=f"bubble search depth (1-{MAX_BUBBLE_DEPTH}, default: {DEFAULT_BUBBLE_DEPTH})",
    )
    analyze_parser = subparsers.add_parser(
        "analyze", help="analyze one or more FASTA/FASTQ read files as successive batches"
    )
    analyze_parser.add_argument("inputs", nargs="+", type=Path, help="read files (.gz supported)")
    analyze_parser.add_argument("--k", type=int, default=31, help="k-mer length (default: 31)")
    analyze_parser.add_argument("--strand-mode", choices=sorted(STRAND_MODES), default=DEFAULT_STRAND_MODE)
    analyze_parser.add_argument(
        "--bubble-depth",
        type=int,
        default=DEFAULT_BUBBLE_DEPTH,
        help=f"bubble search depth (1-{MAX_BUBBLE_DEPTH}, default: {DEFAULT_BUBBLE_DEPTH})",
    )
    analyze_parser.add_argument("--output", type=Path, help="write the JSON report to this path")
    analyze_parser.add_argument("--model-checkpoint", type=Path, help="optional trained path scorer JSON checkpoint")
    replay_parser = subparsers.add_parser(
        "replay", help="deterministically replay ordered read files as offline batches"
    )
    replay_parser.add_argument("inputs", nargs="+", type=Path, help="ordered read batches (.gz supported)")
    replay_parser.add_argument("--k", type=int, default=31, help="k-mer length (default: 31)")
    replay_parser.add_argument("--strand-mode", choices=sorted(STRAND_MODES), default=DEFAULT_STRAND_MODE)
    replay_parser.add_argument("--bubble-depth", type=int, default=DEFAULT_BUBBLE_DEPTH)
    replay_parser.add_argument("--output", type=Path, required=True, help="write replay JSON report")
    replay_parser.add_argument("--model-checkpoint", type=Path, help="optional trained path scorer JSON checkpoint")
    benchmark_parser = subparsers.add_parser(
        "benchmark", help="evaluate heuristic bubble choices against exact reference paths"
    )
    benchmark_parser.add_argument("--reads", nargs="+", required=True, type=Path)
    benchmark_parser.add_argument("--reference", required=True, type=Path)
    benchmark_parser.add_argument("--k", type=int, default=31)
    benchmark_parser.add_argument("--seed", type=int, default=1)
    benchmark_parser.add_argument("--max-reads", type=int)
    benchmark_parser.add_argument(
        "--bubble-depth", type=int, default=DEFAULT_BUBBLE_DEPTH,
        help=f"bubble search depth (1-{MAX_BUBBLE_DEPTH}, default: {DEFAULT_BUBBLE_DEPTH})",
    )
    benchmark_parser.add_argument(
        "--bubble-search-node-cap", type=int, default=DEFAULT_BUBBLE_SEARCH_NODES,
        help=f"per-branch search cap (2-{MAX_BUBBLE_SEARCH_NODES}, default: {DEFAULT_BUBBLE_SEARCH_NODES})",
    )
    benchmark_parser.add_argument("--min-reference-identity", type=float, default=0.90)
    benchmark_parser.add_argument("--strand-mode", choices=sorted(STRAND_MODES), default=DEFAULT_STRAND_MODE)
    benchmark_parser.add_argument("--graph-backend", choices=("cpu", "torch"), default="cpu",
                                  help="graph construction backend; torch is optional")
    benchmark_parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto",
                                  help="tensor device for --graph-backend torch")
    benchmark_parser.add_argument(
        "--no-python-memory-tracing", action="store_true",
        help="skip tracemalloc overhead for faster large-input runs (peak metric will be null)",
    )
    benchmark_parser.add_argument("--model-checkpoint", type=Path, help="optional checkpoint for held-out NN evaluation")
    benchmark_parser.add_argument("--dataset-id", help="genome/sample identifier for split auditing")
    benchmark_parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "demo":
        try:
            run_demo(args.output, bubble_depth=args.bubble_depth)
        except ValueError as error:
            parser.error(str(error))
    elif args.command == "analyze":
        try:
            twin = AssemblyGraphTwin(k=args.k, bubble_depth=args.bubble_depth, strand_mode=args.strand_mode)
            accumulated_reads: list[Read] = []
            for input_path in args.inputs:
                batch = load_reads(input_path)
                accumulated_reads.extend(batch)
                snapshot = twin.update(
                    accumulated_reads,
                    label=f"after {input_path.name}",
                    source=str(input_path),
                )
                if args.model_checkpoint:
                    from path_model import resolve_snapshot
                    snapshot["model_resolution"] = resolve_snapshot(snapshot, args.model_checkpoint)
                print(summarize_snapshot(snapshot))
            report = {"snapshots": twin.snapshots}
            if args.output is not None:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
                print(f"Wrote report to {args.output}")
        except (OSError, ValueError, RuntimeError) as error:
            parser.error(str(error))
    elif args.command == "replay":
        try:
            twin = AssemblyGraphTwin(k=args.k, bubble_depth=args.bubble_depth, strand_mode=args.strand_mode)
            accumulated_reads: list[Read] = []
            for batch_index, input_path in enumerate(args.inputs, start=1):
                batch = load_reads(input_path)
                accumulated_reads.extend(batch)
                timestamp = (datetime(2000, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=batch_index)).isoformat()
                snapshot = twin.update(
                    accumulated_reads,
                    label=f"replay batch {batch_index}: {input_path.name}",
                    source=str(input_path),
                    timestamp_utc=timestamp,
                )
                if args.model_checkpoint:
                    from path_model import resolve_snapshot
                    snapshot["model_resolution"] = resolve_snapshot(snapshot, args.model_checkpoint)
                print(summarize_snapshot(snapshot))
            report = {
                "schema_version": "1.0",
                "replay": {"mode": "offline_deterministic", "batch_order": [str(path) for path in args.inputs], "k": args.k, "bubble_depth": args.bubble_depth, "strand_mode": args.strand_mode},
                "snapshots": twin.snapshots,
            }
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            print(f"Wrote replay report to {args.output}")
        except (OSError, ValueError, RuntimeError) as error:
            parser.error(str(error))
    elif args.command == "benchmark":
        try:
            from benchmark import run_benchmark

            report = run_benchmark(
                args.reads, args.reference, args.k, args.seed, args.max_reads,
                args.model_checkpoint, args.dataset_id, args.bubble_depth,
                args.bubble_search_node_cap,
                args.min_reference_identity,
                args.strand_mode,
                args.graph_backend,
                args.device,
                progress=lambda message: print(message, flush=True),
                trace_python_memory=not args.no_python_memory_tracing,
            )
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            print(json.dumps(report["metrics"], sort_keys=True))
            print(json.dumps({
                "graph_backend": report["protocol"]["graph_backend"],
                "graph_device": report["protocol"]["graph_device"],
            }, sort_keys=True))
            summary = report["inputs"]["diagnostic_summary"]
            print(json.dumps({key: value for key, value in summary.items() if key != "examples"}, sort_keys=True))
            print(json.dumps({"reference_validation": report["reference_validation"]}, sort_keys=True))
            if summary["examples"]["repeat_or_branch_nodes"]:
                example = summary["examples"]["repeat_or_branch_nodes"][0]
                print(
                    "Example branch/repeat node: "
                    f"{example['node']} (coverage={example['coverage']}, "
                    f"in={example['in_degree']}, out={example['out_degree']})"
                )
            if summary["examples"]["bubble_search_limited"]:
                example = summary["examples"]["bubble_search_limited"][0]
                print(
                    "Example unresolved search: "
                    f"start={example['start']}, limit_reasons={example['branch_limit_reason']}"
                )
            print(f"Wrote benchmark report to {args.output}")
        except (OSError, ValueError, RuntimeError) as error:
            parser.error(str(error))


if __name__ == "__main__":
    main()
