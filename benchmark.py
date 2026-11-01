"""Deterministic, reference-based graph candidate benchmark utilities.

This is a graph-level diagnostic benchmark, not a complete genome assembler
evaluation. Reference-compatible paths are proxy labels, not variant truth.
"""

from __future__ import annotations

import hashlib
import heapq
import json
import math
import time
import tracemalloc
from pathlib import Path
from typing import Any, Callable

from assembly_graph_twin import (
    AssemblyGraph,
    CompactAssemblyGraph,
    DEFAULT_BUBBLE_DEPTH,
    DEFAULT_BUBBLE_SEARCH_NODES,
    MAX_BUBBLE_SEARCH_NODES,
    Read,
    build_graph,
    detect_bubbles_with_limits,
    detect_repeat_candidates,
    detect_tips,
    graph_metrics,
    iter_load_reads,
    parse_fasta,
    read_fingerprint,
    reverse_complement,
    _open_text,
    _compact_component_summaries,
)


def _path_sequence(path: list[str]) -> str:
    if not path:
        return ""
    return path[0] + "".join(node[-1] for node in path[1:])


def _sample_rank(seed: int, read_id: str) -> int:
    return int.from_bytes(
        hashlib.sha256(f"{seed}\0{read_id}".encode("utf-8")).digest(), "big"
    )


def reference_records(path: Path) -> dict[str, str]:
    with _open_text(path) as stream:
        reads = parse_fasta(stream, path)
    return {read.read_id: read.sequence for read in reads}


def _reference_anchor_index(
    references: dict[str, str], seed_length: int = 15, reference_stride: int = 12,
    max_occurrences: int = 24,
) -> dict[str, dict[str, list[int]]]:
    """Build a bounded-memory sampled exact-seed index for read placement."""
    indexes: dict[str, dict[str, list[int]]] = {}
    for reference_id, raw_sequence in references.items():
        sequence = raw_sequence.upper()
        positions: dict[str, list[int]] = {}
        for offset in range(0, max(0, len(sequence) - seed_length + 1), reference_stride):
            seed = sequence[offset : offset + seed_length]
            if "N" in seed:
                continue
            occurrences = positions.get(seed)
            if occurrences is None:
                positions[seed] = [offset]
            elif len(occurrences) < max_occurrences:
                occurrences.append(offset)
        indexes[reference_id] = positions
    return indexes


def _anchor_placement_candidates(
    sequence: str,
    indexes: dict[str, dict[str, list[int]]],
    seed_length: int,
) -> list[dict[str, Any]]:
    """Find read placements from clusters of exact seeds on either strand.

    This is a lightweight placement proxy, not a gapped read aligner. Returned
    intervals therefore describe estimated whole-read spans; exact-seed support
    and the method are included in the report so the estimate is auditable.
    """
    sequence = sequence.upper()
    if len(sequence) < seed_length:
        return []
    minimum_anchors = 1 if len(sequence) < seed_length * 2 else 2
    cluster_width = max(3, min(16, len(sequence) // 12))
    candidates: list[dict[str, Any]] = []
    for strand, oriented in (
        ("+", sequence),
        ("-", reverse_complement(sequence)),
    ):
        grouped: dict[str, list[tuple[int, int, int]]] = {}
        query_stride = 2 if seed_length < 8 else 8
        for query_offset in range(0, len(oriented) - seed_length + 1, query_stride):
            seed = oriented[query_offset : query_offset + seed_length]
            if "N" in seed:
                continue
            for reference_id, index in indexes.items():
                for reference_offset in index.get(seed, ()):
                    grouped.setdefault(reference_id, []).append(
                        (reference_offset - query_offset, query_offset, reference_offset)
                    )
        for reference_id, anchors in grouped.items():
            anchors.sort()
            right = 0
            for left in range(len(anchors)):
                right = max(right, left)
                while right + 1 < len(anchors) and anchors[right + 1][0] - anchors[left][0] <= cluster_width:
                    right += 1
                cluster = anchors[left : right + 1]
                query_offsets = sorted({anchor[1] for anchor in cluster})
                if len(query_offsets) < minimum_anchors:
                    continue
                if len(query_offsets) > 1 and query_offsets[-1] - query_offsets[0] < seed_length:
                    continue
                diagonal_values = sorted(anchor[0] for anchor in cluster)
                diagonal = diagonal_values[len(diagonal_values) // 2]
                start = diagonal
                # Nearby sliding windows describe the same placement; keep one.
                if any(
                    item["reference"] == reference_id
                    and item["strand"] == strand
                    and abs(item["start"] - start) <= cluster_width
                    for item in candidates
                ):
                    continue
                candidates.append({
                    "reference": reference_id,
                    "strand": strand,
                    "start": start,
                    "anchor_count": len(query_offsets),
                    "anchor_span": query_offsets[-1] - query_offsets[0] + seed_length,
                })
    return sorted(
        candidates,
        key=lambda item: (-item["anchor_count"], -item["anchor_span"], item["reference"], item["strand"], item["start"]),
    )


def reference_backed_graph_validation(
    reads: list[Read], references: dict[str, str], graph: CompactAssemblyGraph,
    k: int, strand_mode: str,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Summarize reference-anchored read and graph-edge support.

    Read spans are estimated from sampled exact-seed placements. Ambiguous or
    weak placements are excluded from coverage and edge-support counts.
    """
    shortest_reference = min(
        (len(sequence) for sequence in references.values() if sequence),
        default=15,
    )
    # Read-to-reference anchors should have enough specificity regardless of
    # graph k. Tying this seed to k made small-k experiments index millions of
    # repetitive short seeds and could make reference validation dominate runs.
    seed_length = min(15, shortest_reference)
    # Sparse reference indexing limits Python object count on Colab. Queries
    # still sample multiple separated seeds; this remains an estimate, not a
    # production aligner.
    reference_stride = 1 if seed_length < 8 else 12
    index_started = time.perf_counter()
    if progress is not None:
        progress(
            f"Reference index: building sparse {seed_length}-mer index "
            f"(reference stride {reference_stride})"
        )
    indexes = _reference_anchor_index(references, seed_length, reference_stride)
    if progress is not None:
        indexed_seeds = sum(len(index) for index in indexes.values())
        progress(
            f"Reference index complete: {indexed_seeds:,} distinct anchors "
            f"in {time.perf_counter() - index_started:.1f}s"
        )
    reference_coverage = {
        reference_id: bytearray(len(sequence))
        for reference_id, sequence in references.items()
    }
    edge_supported = bytearray(len(graph._edge_sources))
    placed_reads = 0
    ambiguous_reads = 0
    placed_bases = 0
    placed_anchor_count = 0
    for read_index, read in enumerate(reads, start=1):
        if progress is not None and read_index % 10 == 0:
            progress(f"Reference placement: processed {read_index:,}/{len(reads):,} selected reads")
        placements = _anchor_placement_candidates(read.sequence, indexes, seed_length)
        minimum_anchors = 1 if len(read.sequence) < seed_length * 2 else 2
        eligible = [item for item in placements if item["anchor_count"] >= minimum_anchors]
        if not eligible:
            continue
        best = eligible[0]
        equally_supported = [
            item for item in eligible
            if item["anchor_count"] == best["anchor_count"]
            and item["anchor_span"] == best["anchor_span"]
        ]
        distinct_loci = {
            (item["reference"], item["strand"], item["start"])
            for item in equally_supported
        }
        if len(distinct_loci) > 1:
            ambiguous_reads += 1
            continue
        reference_id = best["reference"]
        reference = references[reference_id]
        start = max(0, min(len(reference), best["start"]))
        end = max(start, min(len(reference), start + len(read.sequence)))
        if end <= start:
            continue
        reference_coverage[reference_id][start:end] = b"\1" * (end - start)
        placed_reads += 1
        placed_bases += len(read.sequence)
        placed_anchor_count += best["anchor_count"]

        sequence = read.sequence.upper()
        orientations = (sequence, reverse_complement(sequence)) if strand_mode == "bidirected" else (sequence,)
        for oriented in orientations:
            for offset in range(max(0, len(oriented) - k + 1)):
                kmer = oriented[offset : offset + k]
                if any(base not in "ACGT" for base in kmer):
                    continue
                source_id = graph._node_ids.get(kmer[:-1])
                target_id = graph._node_ids.get(kmer[1:])
                if source_id is None or target_id is None:
                    continue
                edge_id = graph._edge_ids.get((source_id << 32) | target_id)
                if edge_id is not None:
                    edge_supported[edge_id] = 1

    reference_length = sum(len(sequence) for sequence in references.values())
    covered_bases = sum(sum(mask) for mask in reference_coverage.values())
    graph_edge_count = len(graph._edge_sources)
    supported_edge_count = sum(edge_supported)
    return {
        "method": f"sampled exact {seed_length}-mer anchor placement on both read strands; estimated whole-read spans",
        "is_production_alignment": False,
        "anchor_seed_length": seed_length,
        "reference_index_stride": reference_stride,
        "minimum_anchor_count": "1 for reads shorter than two seeds; otherwise 2 separated anchors",
        "selected_read_count": len(reads),
        "reference_placed_read_count": placed_reads,
        "ambiguous_anchor_placement_read_count": ambiguous_reads,
        "unplaced_read_count": len(reads) - placed_reads - ambiguous_reads,
        "reference_placed_read_fraction": placed_reads / len(reads) if reads else None,
        "reference_placed_read_bases": placed_bases,
        "reference_anchor_count_on_placed_reads": placed_anchor_count,
        "estimated_reference_covered_bases": covered_bases,
        "total_reference_bases": reference_length,
        "estimated_reference_breadth": covered_bases / reference_length if reference_length else None,
        "estimated_coverage_by_reference": {
            reference_id: {
                "covered_bases": sum(reference_coverage[reference_id]),
                "reference_bases": len(sequence),
                "breadth": sum(reference_coverage[reference_id]) / len(sequence) if sequence else None,
            }
            for reference_id, sequence in sorted(references.items())
        },
        "graph_edge_count": graph_edge_count,
        "graph_edges_supported_by_reference_placed_reads": supported_edge_count,
        "graph_edge_support_fraction": supported_edge_count / graph_edge_count if graph_edge_count else None,
        "interpretation": (
            "Read placement and breadth are estimates from exact-seed anchor clusters; they are not gapped alignments. "
            "An edge is counted when a uniquely anchor-placed read traverses it."
        ),
    }


def _edit_distance_with_limit(first: str, second: str, limit: int) -> int | None:
    if abs(len(first) - len(second)) > limit:
        return None
    previous = list(range(len(second) + 1))
    for row_index, first_base in enumerate(first, start=1):
        current = [row_index] + [limit + 1] * len(second)
        low = max(1, row_index - limit)
        high = min(len(second), row_index + limit)
        for column in range(low, high + 1):
            substitution = previous[column - 1] + (first_base != second[column - 1])
            current[column] = min(
                previous[column] + 1,
                current[column - 1] + 1,
                substitution,
            )
        if min(current) > limit:
            return None
        previous = current
    distance = previous[-1]
    return distance if distance <= limit else None


def reference_supported_alignments(
    sequence: str,
    references: dict[str, str],
    min_identity: float = 0.90,
    max_seed_occurrences: int = 128,
    max_candidate_starts: int = 2048,
) -> tuple[list[dict[str, Any]], bool]:
    """Seed-and-extend approximate alignment of a candidate path to references.

    This intentionally small standard-library aligner is suitable for labeling
    short graph paths. It is not a substitute for a production read aligner.
    """
    if not 0.5 <= min_identity <= 1.0:
        raise ValueError("minimum reference identity must be between 0.5 and 1.0")
    query = sequence.upper()
    if not query:
        return [], False
    reverse = query.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    max_edits = max(1, int(len(query) * (1.0 - min_identity)))
    seed_length = min(15, max(1, len(query) // 2))
    step = max(1, seed_length // 2)
    offsets = sorted(set(range(0, max(1, len(query) - seed_length + 1), step)) | {len(query) - seed_length})
    best: dict[tuple[str, str], dict[str, Any]] = {}
    truncated = False
    for strand, oriented_query in (("+", query), ("-", reverse)):
        for reference_id, reference in references.items():
            starts: set[int] = set()
            for offset in offsets:
                seed = oriented_query[offset : offset + seed_length]
                cursor = 0
                occurrences = 0
                while True:
                    position = reference.find(seed, cursor)
                    if position < 0:
                        break
                    occurrences += 1
                    if occurrences > max_seed_occurrences:
                        truncated = True
                        break
                    predicted_start = position - offset
                    for adjustment in range(-max_edits, max_edits + 1):
                        candidate_start = predicted_start + adjustment
                        if 0 <= candidate_start < len(reference):
                            starts.add(candidate_start)
                            if len(starts) >= max_candidate_starts:
                                truncated = True
                                break
                    if len(starts) >= max_candidate_starts:
                        break
                    cursor = position + 1
                if len(starts) >= max_candidate_starts:
                    break
            for candidate_start in starts:
                for length_delta in range(-max_edits, max_edits + 1):
                    target_length = len(oriented_query) + length_delta
                    end = candidate_start + target_length
                    if target_length <= 0 or end > len(reference):
                        continue
                    target = reference[candidate_start:end]
                    distance = _edit_distance_with_limit(oriented_query, target, max_edits)
                    if distance is None:
                        continue
                    identity = 1.0 - distance / max(len(oriented_query), len(target))
                    if identity < min_identity:
                        continue
                    key = (reference_id, strand)
                    hit = {
                        "reference": reference_id,
                        "strand": strand,
                        "start": candidate_start,
                        "end": end,
                        "edit_distance": distance,
                        "identity": round(identity, 6),
                    }
                    if key not in best or (distance, candidate_start, end) < (
                        best[key]["edit_distance"], best[key]["start"], best[key]["end"]
                    ):
                        best[key] = hit
    return sorted(best.values(), key=lambda hit: (hit["reference"], hit["strand"], hit["start"])), truncated


def _heuristic_choice(bubble: dict[str, Any]) -> tuple[str | None, str]:
    ranked = sorted(
        bubble["branches"],
        key=lambda branch: (
            -len(bubble["path_read_support"][branch]),
            -bubble["branch_support"][branch],
            branch,
        ),
    )
    if len(ranked) < 2:
        return None, "insufficient_branches"
    first, second = ranked[:2]
    first_score = (len(bubble["path_read_support"][first]), bubble["branch_support"][first])
    second_score = (len(bubble["path_read_support"][second]), bubble["branch_support"][second])
    if first_score == second_score:
        return None, "tied_support"
    return first, "highest_complete_path_then_edge_support"


def path_features(
    bubble: dict[str, Any], branch: str,
    graph_data: dict[str, Any] | AssemblyGraph | CompactAssemblyGraph,
) -> list[float]:
    path = bubble["paths"][branch]
    if isinstance(graph_data, (AssemblyGraph, CompactAssemblyGraph)):
        path_edges = [graph_data.edges[(path[index], path[index + 1])] for index in range(len(path) - 1)]
        edge_support = [edge.support_count for edge in path_edges]
        quality_support = [edge.quality_support for edge in path_edges]
        internal_coverage = [graph_data.nodes[name].coverage for name in path[1:-1]]
    else:
        nodes = {node["sequence"]: node for node in graph_data["nodes"]}
        edges = {(edge["source"], edge["target"]): edge for edge in graph_data["edges"]}
        path_edges = [edges[(path[index], path[index + 1])] for index in range(len(path) - 1)]
        edge_support = [len(edge["supporting_reads"]) for edge in path_edges]
        quality_support = [edge["quality_support"] for edge in path_edges]
        internal_coverage = [nodes[name]["coverage"] for name in path[1:-1]]
    return [
        math.log1p(max(0, len(path) - 1)),
        math.log1p(bubble["branch_support"][branch]),
        math.log1p(len(bubble["path_read_support"][branch])),
        math.log1p(sum(edge_support) / len(edge_support)) if edge_support else 0.0,
        sum(quality_support) / len(quality_support) if quality_support else 0.0,
        sum(internal_coverage) / max(1, len(internal_coverage)),
    ]


def _run_benchmark_impl(
    read_paths: list[Path], reference_path: Path, k: int, seed: int, max_reads: int | None,
    model_checkpoint: Path | None, dataset_id: str | None, started: float,
    bubble_depth: int, bubble_search_node_cap: int,
    min_reference_identity: float, strand_mode: str,
    graph_backend: str, device: str,
    progress: Callable[[str], None] | None,
) -> dict[str, Any]:
    if max_reads is not None and max_reads < 1:
        raise ValueError("--max-reads must be at least 1")
    if not 0.5 <= min_reference_identity <= 1.0:
        raise ValueError("minimum reference identity must be between 0.5 and 1.0")
    selected_heap: list[tuple[int, str, Read]] = []
    seen_ids: set[str] = set()
    total_read_count = 0
    read_selection_started = time.perf_counter()
    for path in read_paths:
        for read in iter_load_reads(path):
            if read.read_id in seen_ids:
                raise ValueError("read IDs must be unique across benchmark input files")
            seen_ids.add(read.read_id)
            total_read_count += 1
            rank = _sample_rank(seed, read.read_id)
            if max_reads is None:
                heapq.heappush(selected_heap, (-rank, read.read_id, read))
            elif len(selected_heap) < max_reads:
                heapq.heappush(selected_heap, (-rank, read.read_id, read))
            else:
                current_worst_rank = -selected_heap[0][0]
                if rank < current_worst_rank:
                    heapq.heapreplace(selected_heap, (-rank, read.read_id, read))
            if progress is not None and total_read_count % 1000 == 0:
                progress(f"Read selection: scanned {total_read_count:,} reads")
    reads = [entry[2] for entry in sorted(selected_heap, key=lambda entry: -entry[0])]
    read_selection_seconds = time.perf_counter() - read_selection_started
    if progress is not None:
        progress(
            f"Read selection complete: retained {len(reads):,} of {total_read_count:,} reads "
            f"in {read_selection_seconds:.1f}s"
        )
    references = reference_records(reference_path)
    if not references:
        raise ValueError("reference FASTA contains no sequences")

    graph_started = time.perf_counter()
    graph = build_graph(
        reads, k, compact=True, strand_mode=strand_mode,
        backend=graph_backend, device=device, progress=progress,
    )
    graph_construction_seconds = time.perf_counter() - graph_started
    if progress is not None:
        progress(
            f"Graph construction complete on {graph.device} ({graph.backend}): "
            f"{len(graph.nodes):,} nodes, {len(graph.edges):,} oriented edges "
            f"in {graph_construction_seconds:.1f}s"
        )
    diagnostics_started = time.perf_counter()
    if progress is not None:
        progress("Running graph diagnostics on CPU")
    tips = detect_tips(graph)
    repeats = detect_repeat_candidates(graph)
    bubbles, limited_searches = detect_bubbles_with_limits(
        graph, depth=bubble_depth, max_nodes=bubble_search_node_cap, progress=progress
    )
    components = _compact_component_summaries(graph)
    diagnostic_summary = {
        "tip_count": len(tips),
        "repeat_or_branch_count": len(repeats),
        "branching_repeat_count": sum(
            item["in_degree"] > 1 or item["out_degree"] > 1 for item in repeats
        ),
        "high_support_repeat_count": sum(item["coverage"] >= 3 for item in repeats),
        "bubble_candidate_count": len(bubbles),
        "bubble_search_limited_count": len(limited_searches),
        "bubble_search_limit_reasons": {
            reason: sum(
                value == reason
                for search in limited_searches
                for value in search.get("branch_limit_reason", {}).values()
            )
            for reason in ("depth", "node_cap")
        },
        "disconnected_component_count": len(components),
        "examples": {
            "tips": tips[:5],
            "repeat_or_branch_nodes": repeats[:10],
            "largest_components": sorted(
                components,
                key=lambda item: (-int(item["node_count"]), str(item["min_node"])),
            )[:10],
            "bubble_search_limited": limited_searches[:10],
        },
        "interpretation": (
            "No bounded bubble candidate was found. Structural diagnostics are evidence about this sampled graph, "
            "not confirmed variants, repeats, errors, or coverage gaps."
            if not bubbles else
            "Structural diagnostics are evidence about this sampled graph, not confirmed variants, repeats, errors, or coverage gaps."
        ),
    }
    diagnostic_seconds = time.perf_counter() - diagnostics_started
    if progress is not None:
        progress(
            f"Graph diagnostics complete: {len(bubbles):,} bubbles, {len(repeats):,} branch/repeat "
            f"signals in {diagnostic_seconds:.1f}s"
        )
    for bubble in bubbles:
        for branch, path in bubble["paths"].items():
            branch_sequence = _path_sequence(path)
            bubble["path_read_support"][branch] = [
                read.read_id for read in reads
                if branch_sequence in read.sequence.upper()
                or reverse_complement(branch_sequence) in read.sequence.upper()
            ]
    metrics = graph_metrics(graph)
    reference_started = time.perf_counter()
    reference_validation = reference_backed_graph_validation(
        reads, references, graph, k, strand_mode, progress=progress
    )
    reference_validation_seconds = time.perf_counter() - reference_started
    if progress is not None:
        progress(
            f"Reference placement complete: {reference_validation['reference_placed_read_count']:,} "
            f"reads placed (estimated) in {reference_validation_seconds:.1f}s"
        )
    reads_fingerprint = read_fingerprint(reads)
    protocol_path = Path(__file__).parent / "config" / "benchmark_protocol.json"
    protocol_config = json.loads(protocol_path.read_text(encoding="utf-8")) if protocol_path.exists() else {}
    protocol_sha256 = hashlib.sha256(protocol_path.read_bytes()).hexdigest() if protocol_path.exists() else None
    rows: list[dict[str, Any]] = []
    training_examples: list[dict[str, Any]] = []
    evaluated = 0
    correct = 0
    resolved = 0
    false_resolutions = 0
    model_evaluated = model_resolved = model_correct = model_false_resolutions = 0
    checkpoint = None
    model_scorer = None
    if model_checkpoint is not None:
        checkpoint = json.loads(model_checkpoint.read_text(encoding="utf-8"))
        if checkpoint.get("feature_count") != 6:
            raise ValueError("model checkpoint has an incompatible feature schema")
        evaluation_group = dataset_id or reference_path.stem
        trained_groups = set(checkpoint.get("train_dataset_ids", [])) | set(checkpoint.get("validation_dataset_ids", []))
        if evaluation_group in trained_groups:
            raise ValueError("held-out evaluation dataset_id overlaps the model's training or validation groups")
        from path_model import make_scorer
        model_scorer = make_scorer(checkpoint)
    path_alignment_started = time.perf_counter()
    for bubble in bubbles:
        compatible = []
        branch_reference_hits = {}
        branch_reference_alignments = {}
        for branch, path in bubble["paths"].items():
            sequence = _path_sequence(path)
            alignments, alignment_truncated = reference_supported_alignments(
                sequence, references, min_identity=min_reference_identity
            )
            branch_reference_alignments[branch] = {
                "hits": alignments,
                "search_truncated": alignment_truncated,
            }
            branch_reference_hits[branch] = sorted({hit["reference"] for hit in alignments})
            if alignments:
                compatible.append(branch)
        choice, rationale = _heuristic_choice(bubble)
        branch_features = {
            branch: path_features(bubble, branch, graph)
            for branch in bubble["branches"]
        }
        model_selection = None
        if checkpoint is not None:
            probabilities = {branch: model_scorer(features) for branch, features in branch_features.items()}
            ranked_model = sorted(probabilities, key=lambda branch: (-probabilities[branch], branch))
            model_choice = ranked_model[0] if ranked_model and probabilities[ranked_model[0]] >= checkpoint["auto_resolve_threshold"] else None
            model_selection = {
                "selected_branch": model_choice,
                "decision": "auto_resolved" if model_choice else "abstain",
                "branch_probabilities": probabilities,
                "threshold": checkpoint["auto_resolve_threshold"],
                "model_version": checkpoint.get("model_version"),
            }
        else:
            model_choice = None
        alignment_was_truncated = any(
            value["search_truncated"] for value in branch_reference_alignments.values()
        )
        truth = compatible[0] if len(compatible) == 1 and not alignment_was_truncated else None
        if alignment_was_truncated:
            status = "alignment_search_truncated"
        else:
            status = "unique_reference_path" if truth else (
                "multiple_reference_paths" if compatible else "no_reference_supported_path"
            )
        if truth:
            evaluated += 1
            if choice is not None:
                resolved += 1
                if choice == truth:
                    correct += 1
                else:
                    false_resolutions += 1
            if model_selection is not None:
                model_evaluated += 1
                if model_choice is not None:
                    model_resolved += 1
                    if model_choice == truth:
                        model_correct += 1
                    else:
                        model_false_resolutions += 1
            for branch, path in bubble["paths"].items():
                features = path_features(bubble, branch, graph)
                training_examples.append({
                    "candidate_id": bubble["candidate_id"],
                    "branch": branch,
                    "features": features,
                    "label": int(branch == truth),
                })
        rows.append({
            "candidate_id": bubble["candidate_id"],
            "start": bubble["start"],
            "rejoin": bubble["rejoin"],
            "reference_label_status": status,
            "reference_compatible_branches": compatible,
            "branch_reference_hits": branch_reference_hits,
            "branch_reference_alignments": branch_reference_alignments,
            "reference_alignment_search_truncated": alignment_was_truncated,
            "heuristic_selected_branch": choice,
            "heuristic_decision": rationale,
            "heuristic_correct": (choice == truth) if truth and choice else None,
            "paths": bubble["paths"],
            "path_read_support": {key: len(value) for key, value in bubble["path_read_support"].items()},
            "branch_support": bubble["branch_support"],
            "branch_features": branch_features,
            "model_selection": model_selection,
        })
    path_alignment_seconds = time.perf_counter() - path_alignment_started
    if progress is not None:
        progress(f"Path labeling complete: {len(rows):,} candidates in {path_alignment_seconds:.1f}s")
    report = {
        "schema_version": "1.0",
        "benchmark": "reference_supported_path_compatibility",
        "reference_validation": reference_validation,
        "protocol": {
            "protocol_version": protocol_config.get("protocol_version", "unknown"),
            "protocol_sha256": protocol_sha256,
            "seed": seed,
            "k": k,
            "max_reads": max_reads,
            "bubble_depth": bubble_depth,
            "bubble_search_node_cap": bubble_search_node_cap,
            "strand_mode": strand_mode,
            "graph_backend": graph.backend,
            "graph_device": graph.device,
            "gpu_peak_allocated_bytes": graph.gpu_peak_allocated_bytes,
            "gpu_peak_reserved_bytes": graph.gpu_peak_reserved_bytes,
            "accelerated_stages": (
                ["kmer_window_encoding", "per_read_edge_deduplication", "edge_support_aggregation", "node_support_aggregation", "compact_adjacency_csr"]
                if graph.backend == "torch" else []
            ),
            "strand_representation": (
                "bidirected double-cover with reverse-complement twin transitions"
                if strand_mode == "bidirected" else
                "forward-only directed graph preserving input read orientation"
            ),
            "reference_label_method": "seed-and-extend edit-distance alignment against the reference and its reverse complement",
            "read_reference_mapping_method": reference_validation["method"],
            "minimum_reference_identity": min_reference_identity,
            "read_selection": "sha256 bottom-k selection by seed and read ID; subsets are nested across max_reads values" if max_reads is not None else "all input reads sorted by deterministic read-ID hash",
            "sampling_algorithm": "sha256_bottom_k_by_read_id_v1",
            "reference_matching": "bounded seed-and-extend path alignment; branch labels require the configured minimum identity",
            "label_semantics": "reference_compatibility_proxy",
            "limitations": [
                "Reference compatibility is a proxy label, not biological variant truth; approximate alignment can accept near matches.",
                "A path absent from the reference is not necessarily an error or incorrect path.",
                "Multiple reference-compatible paths are left ambiguous; one reference does not establish variant truth.",
                "This report evaluates graph bubble selection only and does not emit or benchmark contig assemblies.",
                "The bidirected double-cover stores each observed transition and its reverse-complement twin; it is not a compact canonical-k-mer graph.",
                "Coverage is per physical read per oriented node/edge, not a count of distinct canonical k-molecules.",
                "Read placement, placed-base totals, and reference breadth are anchor-derived estimates, not gapped alignments or base-level consensus.",
            ],
        },
        "inputs": {
            "read_files": [str(path) for path in read_paths],
            "reference_file": str(reference_path),
            "dataset_id": dataset_id or reference_path.stem,
            "read_fingerprint_sha256": reads_fingerprint,
            "reference_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
            "model_checkpoint_sha256": hashlib.sha256(model_checkpoint.read_bytes()).hexdigest() if model_checkpoint else None,
            "read_count": len(reads),
            "total_input_read_count": total_read_count,
            "selected_read_count": len(reads),
            "graph_metrics": metrics,
            "diagnostic_summary": diagnostic_summary,
        },
        "metrics": {
            "read_selection_seconds": round(read_selection_seconds, 6),
            "graph_construction_seconds": round(graph_construction_seconds, 6),
            "gpu_peak_allocated_bytes": graph.gpu_peak_allocated_bytes,
            "gpu_peak_reserved_bytes": graph.gpu_peak_reserved_bytes,
            "graph_diagnostics_seconds": round(diagnostic_seconds, 6),
            "reference_validation_seconds": round(reference_validation_seconds, 6),
            "path_alignment_seconds": round(path_alignment_seconds, 6),
            "candidate_count": len(rows),
            "uniquely_labeled_candidate_count": evaluated,
            "resolved_labeled_candidate_count": resolved,
            "abstention_rate_on_labeled": (evaluated - resolved) / evaluated if evaluated else None,
            "accuracy_on_resolved_labeled": correct / resolved if resolved else None,
            "false_resolutions": false_resolutions,
            "coverage_on_labeled": resolved / evaluated if evaluated else None,
            "runtime_seconds": round(time.perf_counter() - started, 6),
            "python_tracemalloc_peak_bytes": None,
            "model": ({
                "evaluated_labeled_candidates": model_evaluated,
                "resolved_labeled_candidates": model_resolved,
                "abstention_rate_on_labeled": (model_evaluated - model_resolved) / model_evaluated if model_evaluated else None,
                "accuracy_on_resolved_labeled": model_correct / model_resolved if model_resolved else None,
                "false_resolutions": model_false_resolutions,
            } if model_checkpoint is not None else None),
        },
        "candidates": rows,
        "training_examples": training_examples,
    }
    return report


def benchmark_evidence_fingerprint(report: dict[str, Any]) -> str:
    """Hash deterministic benchmark evidence, excluding timing and memory metrics."""
    stable = {
        "protocol": report["protocol"],
        "inputs": report["inputs"],
        "metrics": {
            key: value for key, value in report["metrics"].items()
            if key not in {
                "runtime_seconds", "python_tracemalloc_peak_bytes",
                "read_selection_seconds", "graph_construction_seconds",
                "graph_diagnostics_seconds", "reference_validation_seconds",
                "path_alignment_seconds", "gpu_peak_allocated_bytes",
                "gpu_peak_reserved_bytes",
            }
        },
        "candidates": report["candidates"],
        "training_examples": report["training_examples"],
    }
    # Runtime and memory are expected to vary; all graph evidence and labels
    # should remain byte-stable for fixed inputs and options.
    return hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def run_benchmark(
    read_paths: list[Path], reference_path: Path, k: int, seed: int, max_reads: int | None,
    model_checkpoint: Path | None = None, dataset_id: str | None = None,
    bubble_depth: int = DEFAULT_BUBBLE_DEPTH,
    bubble_search_node_cap: int = DEFAULT_BUBBLE_SEARCH_NODES,
    min_reference_identity: float = 0.90,
    strand_mode: str = "bidirected",
    graph_backend: str = "cpu",
    device: str = "auto",
    progress: Callable[[str], None] | None = None,
    trace_python_memory: bool = True,
) -> dict[str, Any]:
    if trace_python_memory and tracemalloc.is_tracing():
        raise RuntimeError("benchmark manages Python memory tracing; stop external tracemalloc first")
    started = time.perf_counter()
    if trace_python_memory:
        tracemalloc.start()
    try:
        report = _run_benchmark_impl(
            read_paths, reference_path, k, seed, max_reads, model_checkpoint, dataset_id, started,
            bubble_depth, bubble_search_node_cap, min_reference_identity, strand_mode,
            graph_backend, device, progress,
        )
        if trace_python_memory:
            _, peak_bytes = tracemalloc.get_traced_memory()
            report["metrics"]["python_tracemalloc_peak_bytes"] = peak_bytes
        return report
    finally:
        if trace_python_memory:
            tracemalloc.stop()


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", nargs="+", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--max-reads", type=int)
    parser.add_argument("--bubble-depth", type=int, default=DEFAULT_BUBBLE_DEPTH)
    parser.add_argument("--bubble-search-node-cap", type=int, default=DEFAULT_BUBBLE_SEARCH_NODES)
    parser.add_argument("--min-reference-identity", type=float, default=0.90)
    parser.add_argument("--strand-mode", choices=("bidirected", "forward"), default="bidirected")
    parser.add_argument("--graph-backend", choices=("cpu", "torch"), default="cpu")
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument(
        "--no-python-memory-tracing", action="store_true",
        help="skip tracemalloc overhead for faster large-input runs (peak metric will be null)",
    )
    parser.add_argument("--dataset-id", help="genome/sample identifier used to audit train/validation separation")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        report = run_benchmark(
            args.reads, args.reference, args.k, args.seed, args.max_reads,
            dataset_id=args.dataset_id, bubble_depth=args.bubble_depth,
            bubble_search_node_cap=args.bubble_search_node_cap,
            min_reference_identity=args.min_reference_identity,
            strand_mode=args.strand_mode,
            graph_backend=args.graph_backend,
            device=args.device,
            progress=lambda message: print(message, flush=True),
            trace_python_memory=not args.no_python_memory_tracing,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
