"""Optional PyTorch k-mer extraction for compact benchmark graph construction.

Only the embarrassingly parallel sequence-window encoding and per-read edge
deduplication run on the requested tensor device. Compact graph assembly,
diagnostics, reference placement, and path labeling remain CPU operations.
"""

from __future__ import annotations

from array import array
from typing import Callable, Iterable

from assembly_graph_twin import CompactAssemblyGraph, Read, reverse_complement


def _torch_module():
    try:
        import torch
    except ImportError as error:
        raise RuntimeError(
            "The PyTorch graph backend is optional. Install PyTorch in this environment, "
            "or use --graph-backend cpu."
        ) from error
    return torch


def resolve_device(device: str) -> str:
    torch = _torch_module()
    if device == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA was requested, but this Python runtime cannot access a CUDA GPU. "
            "In Colab choose Runtime > Change runtime type > T4 GPU, then reconnect."
        )
    if device not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be auto, cpu, or cuda")
    return device


def _base_tensor(sequence: str, torch, device: str):
    # A/C/G/T map to 0/1/2/3. Ambiguous bases map to 4 and invalidate windows.
    encode = {"A": 0, "C": 1, "G": 2, "T": 3}
    values = [encode.get(base, 4) for base in sequence.upper()]
    return torch.tensor(values, dtype=torch.int64, device=device)


def _encoded_kmers(bases, k: int, weights, torch):
    if bases.numel() < k:
        return torch.empty(0, dtype=torch.int64, device=bases.device)
    windows = bases.unfold(0, k, 1)
    valid = (windows < 4).all(dim=1)
    # k <= 31 makes every two-bit encoded k-mer fit in a signed int64.
    codes = (windows * weights).sum(dim=1)
    return codes[valid]


def _decode_kmer(code: int, k: int) -> str:
    alphabet = "ACGT"
    chars = ["A"] * k
    for position in range(k - 1, -1, -1):
        chars[position] = alphabet[code & 3]
        code >>= 2
    return "".join(chars)


def _finalize_adjacency_on_device(graph: CompactAssemblyGraph, torch, device: str) -> None:
    """Build compact incoming/outgoing CSR arrays with tensor operations."""
    import numpy as np

    node_count = len(graph._sequences)
    edge_count = len(graph._edge_sources)
    source_np = np.frombuffer(graph._edge_sources, dtype=np.uint32)
    target_np = np.frombuffer(graph._edge_targets, dtype=np.uint32)
    sources = torch.as_tensor(source_np.astype(np.int64), device=device)
    targets = torch.as_tensor(target_np.astype(np.int64), device=device)

    out_counts = torch.bincount(sources, minlength=node_count)
    in_counts = torch.bincount(targets, minlength=node_count)
    zero = torch.zeros(1, dtype=torch.int64, device=device)
    out_offsets = torch.cat((zero, torch.cumsum(out_counts, dim=0)))
    in_offsets = torch.cat((zero, torch.cumsum(in_counts, dim=0)))
    out_order = torch.argsort(sources, stable=True)
    in_order = torch.argsort(targets, stable=True)

    def to_uint_array(values) -> array:
        packed = array("I")
        packed.frombytes(values.to("cpu").numpy().astype(np.uint32, copy=False).tobytes())
        return packed

    graph._outgoing_offsets = to_uint_array(out_offsets)
    graph._incoming_offsets = to_uint_array(in_offsets)
    graph._outgoing_values = to_uint_array(targets[out_order])
    graph._incoming_values = to_uint_array(sources[in_order])
    graph._finalized = True


def build_compact_graph_torch(
    reads: Iterable[Read], k: int, *, strand_mode: str, device: str = "auto",
    progress: Callable[[str], None] | None = None,
) -> tuple[CompactAssemblyGraph, str]:
    """Extract distinct per-read k-mer edges with PyTorch, then build CPU graph storage."""
    if k > 31:
        raise ValueError("the PyTorch encoded-kmer backend supports k values up to 31")
    if strand_mode not in {"bidirected", "forward"}:
        raise ValueError("strand_mode must be bidirected or forward")
    torch = _torch_module()
    actual_device = resolve_device(device)
    if actual_device == "cuda":
        torch.cuda.reset_peak_memory_stats()
    graph = CompactAssemblyGraph(k, strand_mode=strand_mode)
    weights = torch.tensor(
        [1 << (2 * (k - index - 1)) for index in range(k)],
        dtype=torch.int64,
        device=actual_device,
    )
    seen_read_ids: set[str] = set()
    read_edge_batches = []
    read_qualities: list[float] = []
    ambiguous_reads: list[tuple[int, Read]] = []

    for read_index, read in enumerate(reads):
        if read.read_id in seen_read_ids:
            raise ValueError("compact graph mode requires unique read IDs")
        seen_read_ids.add(read.read_id)
        read_qualities.append(read.quality)
        sequence = read.sequence.upper()
        if any(base not in "ACGT" for base in sequence):
            # Preserve exact CPU/IUPAC semantics; add these reads after GPU
            # aggregation so their support composes correctly with ACGT reads.
            ambiguous_reads.append((read_index, read))
            continue
        bases = _base_tensor(read.sequence, torch, actual_device)
        forward = _encoded_kmers(bases, k, weights, torch)
        if strand_mode == "bidirected" and bases.numel():
            reversed_bases = torch.flip(bases, dims=(0,))
            reverse_complement_bases = torch.where(reversed_bases < 4, 3 - reversed_bases, 4)
            reverse = _encoded_kmers(reverse_complement_bases, k, weights, torch)
            per_read_edges = torch.unique(torch.cat((forward, reverse)), sorted=True)
        else:
            per_read_edges = torch.unique(forward, sorted=True)

        # Retain only distinct (read, edge) pairs on device, then transfer once
        # after all reads. This avoids synchronizing the GPU once per read.
        if per_read_edges.numel():
            read_ids = torch.full_like(per_read_edges, read_index)
            read_edge_batches.append(torch.stack((read_ids, per_read_edges), dim=1))
        if progress is not None and (read_index + 1) % 10 == 0:
            progress(f"CUDA k-mer extraction: processed {read_index + 1:,} reads")

    if read_edge_batches:
        if progress is not None:
            progress(f"CUDA support aggregation: reducing {len(read_edge_batches):,} per-read edge batches")
        all_read_edges = torch.cat(read_edge_batches, dim=0)
        edge_codes, edge_inverse, edge_support = torch.unique(
            all_read_edges[:, 1], sorted=True, return_inverse=True, return_counts=True
        )
        read_indices = all_read_edges[:, 0]
        qualities = torch.tensor(read_qualities, dtype=torch.float64, device=actual_device)
        edge_quality = torch.zeros(edge_codes.numel(), dtype=torch.float64, device=actual_device)
        edge_quality.index_add_(0, edge_inverse, qualities[read_indices])

        node_mask = (1 << (2 * (k - 1))) - 1
        source_codes = all_read_edges[:, 1] >> 2
        target_codes = all_read_edges[:, 1] & node_mask
        node_pairs = torch.cat((
            torch.stack((read_indices, source_codes), dim=1),
            torch.stack((read_indices, target_codes), dim=1),
        ), dim=0)
        node_pairs = torch.unique(node_pairs, dim=0, sorted=True)
        node_codes, node_support = torch.unique(
            node_pairs[:, 1], sorted=True, return_counts=True
        )

        # Transfer aggregate counts once. CPU graph materialization now visits
        # each distinct edge once, instead of once for every supporting read.
        import numpy as np

        edge_codes_cpu = edge_codes.cpu().numpy().astype(np.int64, copy=False)
        edge_support_cpu = edge_support.cpu().numpy().astype(np.uint32, copy=False)
        edge_quality_cpu = edge_quality.cpu().numpy().astype(np.float64, copy=False)
        node_codes_cpu = node_codes.cpu().numpy().astype(np.int64, copy=False)
        node_support_cpu = node_support.cpu().numpy().astype(np.uint32, copy=False)
        source_codes_cpu = edge_codes_cpu >> 2
        target_codes_cpu = edge_codes_cpu & node_mask
        source_ids_cpu = np.searchsorted(node_codes_cpu, source_codes_cpu).astype(np.uint32, copy=False)
        target_ids_cpu = np.searchsorted(node_codes_cpu, target_codes_cpu).astype(np.uint32, copy=False)
        if progress is not None:
            progress(f"CUDA aggregation complete: materializing {len(edge_codes_cpu):,} distinct graph edges on CPU")

        if len(node_codes_cpu) >= 2**32:
            raise ValueError("compact graph exceeds its 32-bit node index capacity")
        node_total = len(node_codes_cpu)
        for node_id, code in enumerate(node_codes_cpu):
            sequence = _decode_kmer(int(code), k - 1)
            graph._node_ids[sequence] = node_id
            graph._sequences.append(sequence)
            if progress is not None and (node_id + 1) % 100_000 == 0:
                progress(f"CPU graph-key materialization: {node_id + 1:,}/{node_total:,} nodes")
        graph._node_support = array("I")
        graph._node_support.frombytes(node_support_cpu.tobytes())
        graph._node_last_read = array("q", [-1]) * len(node_codes_cpu)

        graph._edge_sources = array("I")
        graph._edge_sources.frombytes(source_ids_cpu.tobytes())
        graph._edge_targets = array("I")
        graph._edge_targets.frombytes(target_ids_cpu.tobytes())
        graph._edge_support = array("I")
        graph._edge_support.frombytes(edge_support_cpu.tobytes())
        graph._edge_last_read = array("q", [-1]) * len(edge_codes_cpu)
        graph._edge_quality = array("d")
        graph._edge_quality.frombytes(edge_quality_cpu.tobytes())
        graph._edge_ids = {}
        edge_total = len(source_ids_cpu)
        for edge_id, (source_id, target_id) in enumerate(zip(source_ids_cpu, target_ids_cpu)):
            graph._edge_ids[(int(source_id) << 32) | int(target_id)] = edge_id
            if progress is not None and (edge_id + 1) % 100_000 == 0:
                progress(f"CPU edge-index materialization: {edge_id + 1:,}/{edge_total:,} edges")

    # Ambiguous reads are usually a small fraction of the input and retain the
    # existing CPU handling for IUPAC bases.
    for read_index, read in ambiguous_reads:
        sequence = read.sequence.upper()
        oriented_sequences = [sequence]
        if strand_mode == "bidirected":
            reverse = reverse_complement(sequence)
            if reverse != sequence:
                oriented_sequences.append(reverse)
        transitions = set()
        for oriented in oriented_sequences:
            for offset in range(max(0, len(oriented) - k + 1)):
                kmer = oriented[offset : offset + k]
                transitions.add((kmer[:-1], kmer[1:]))
        for source, target in sorted(transitions):
            graph.add_read_edge(source, target, read_index, read.quality)

    if progress is not None:
        progress(f"Building compact adjacency on {actual_device}")
    _finalize_adjacency_on_device(graph, torch, actual_device)
    if actual_device == "cuda":
        torch.cuda.synchronize()
        graph.gpu_peak_allocated_bytes = torch.cuda.max_memory_allocated()
        graph.gpu_peak_reserved_bytes = torch.cuda.max_memory_reserved()
    return graph, actual_device
