# Assembly Graph Twin

A small, dependency-free genome assembly graph exploration tool. It builds a directed de Bruijn graph from DNA reads, records read support and quality metadata, reports structural signals, and stores cumulative snapshots. Its rules are educational heuristics, not a production assembler or a biological classifier.

## Run the synthetic demo

From the repository root:

```text
python assembly_graph_twin/assembly_graph_twin.py demo
```

Write the JSON snapshots to a file:

```text
python assembly_graph_twin/assembly_graph_twin.py demo --output assembly_graph_twin/demo_report.json
```

## Analyze read files

The `analyze` command accepts FASTA or standard four-line FASTQ files. Gzip-compressed files are supported when their name ends in `.gz`. Multiple inputs are treated as successive batches; each snapshot contains the cumulative graph after that batch.

```text
python assembly_graph_twin/assembly_graph_twin.py analyze reads.fastq --k 31 --output report.json
python assembly_graph_twin/assembly_graph_twin.py analyze batch_1.fasta batch_2.fasta.gz --k 31 --output report.json
python assembly_graph_twin/assembly_graph_twin.py analyze reads.fastq --k 31 --bubble-depth 6 --output report.json
```

The default `k` is 31. Choose a value compatible with the read lengths; `k` must be at least 2. Input bases are limited to `A`, `C`, `G`, `T`, and `N`. FASTQ qualities are interpreted as Phred+33 and summarized as mean per-base accuracy. FASTA reads have no quality data and use the default accuracy value of 1.0. Read orientation is preserved as supplied; reverse-complement normalization is not implemented.

Bubble search depth defaults to 3 and accepts values from 1 through 12. Each branch search also stops after 64 discovered nodes. These conservative limits keep this educational tool bounded; they cannot be raised through the CLI. A complete bubble records its effective depth and per-branch limit status. If a search reaches a safeguard before finding a complete candidate, the report emits a `bubble_search_limited` diagnostic instead of treating the candidate as absent.

Read identifiers must be unique across all batches supplied to one `analyze` command. Repeated identifiers are rejected instead of silently combining their support.

Reports contain read counts, a content fingerprint, graph topology, supporting read IDs, edge quality-support summaries, and diagnostic evidence. Raw read sequences and per-base quality strings are not copied into the report. Snapshot timestamps vary by run, while graph traversal and serialized ordering are deterministic for the same ordered inputs.

## Tests

```text
python -m unittest discover -s assembly_graph_twin -p "test_*.py" -v
```

## What it does and does not claim

- Builds a small directed de Bruijn graph and tracks read-level node and edge support.
- Reports low-support tips, branch/repeat candidates, bounded reconverging paths, and disconnected components.
- Defines a bubble as two or more outgoing branches from one divergence that reach the same nearest valid reconvergence within three hops per branch. Branches in one divergence/rejoin group are emitted in one record; nested or overlapping structures with a different divergence or rejoin remain separate records.
- Each bubble has a stable `candidate_id` derived from its divergence, rejoin, and canonical branch paths. Support counts, timestamps, and read order are excluded from the ID.
- Bubble evidence contains reciprocal `relationships` for candidates that share internal nodes or edges. Relationships are classified as `nested` when one candidate contains the other's divergence and rejoin, otherwise `overlapping`; unrelated candidates have an empty relationship list.
- Stores cumulative graph snapshots when input batches are added.
- Uses explicit deterministic rules and heuristic scores; it has no trained AI model and scores are not probabilities.
- Does not produce contigs, resolve biological variants, or establish that a graph structure is an error, repeat, or coverage gap.
- Does not yet model reverse-complement bidirected graphs, paired reads, adapters, trimming, ambiguity codes beyond `N`, or human-scale data efficiently.
- Bubble candidates beyond the configured depth or 64-node branch-search safeguard are not reported as complete. A `bubble_search_limited` diagnostic records the branches, effective depth, node cap, and which branch searches were truncated. Cycles are traversed once per search; a cycle is evidence only when a separate branch also reaches a valid reconvergence. Ambiguous overlapping paths are preserved as structural evidence rather than assigned a biological interpretation.

See [`assembly_graph_twin_reference.md`](assembly_graph_twin_reference.md) for the project context and model assumptions.
