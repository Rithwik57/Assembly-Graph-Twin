# Explainable Assembly Graph Twin

## Purpose

This prototype demonstrates how an evolving graph model can support exploration of genome assembly evidence. It is an educational tool, not a production assembler, biological classifier, or clinical system. It currently uses deterministic rules; it contains no trained AI model.

Genome sequencing produces unordered DNA reads. Assembly connects reads through overlaps, but sequencing errors, repeats, low coverage, and real variants create ambiguous graph structures. The lectures describe these problems through overlap graphs, greedy assembly, de Bruijn graphs, unitigs, bubbles, and long-read evidence.

## Core idea

The prototype maintains a sequence of graph snapshots. Each snapshot records:

- Input read count, source label, and content fingerprint (raw sequences are not copied into the report)
- Graph nodes, edges, and read support
- Coverage and local graph structure
- Suspected tips, repeats, and bubbles
- Explainable classifications and transparent heuristic scores
- Every graph update as a timestamped snapshot

Rule-based diagnostics annotate structural signals without deleting graph paths. An ambiguous branch remains represented in the graph. Labels such as `true_variant_or_unresolved` are descriptive categories, not biological conclusions.

## Graph model

For a chosen k-mer size `k`:

- Each `(k-1)`-mer is a graph node.
- Each k-mer creates a directed edge from its prefix to its suffix.
- Edges retain the read IDs that support them.
- Node coverage is the number of distinct supporting reads.
- FASTQ quality is summarized as mean expected per-base accuracy and accumulated on supported edges; it is evidence metadata, not a calibrated error model. FASTA reads use the default quality value of 1.0 because they have no quality string.

This is the transformed de Bruijn representation discussed in `lec-31.md`.

## Diagnostics

The prototype detects:

- Low-support terminal nodes: possible sequencing-error tips.
- High-support or high-degree nodes: possible repeats or branch points.
- Branches that reconverge within a bounded search: possible bubbles or variants. A bubble is defined here as two or more outgoing branches from one divergence that reach the same nearest valid reconvergence, with internally disjoint branch paths. The default search depth is 3 and the permitted maximum is 12. The reported rejoin minimizes the larger branch distance, then total distance, then lexical order. Path-level read support is listed when one input read spans the complete branch path.
- Disconnected graph components: an explicit structural diagnostic with component nodes and edge counts.

The current diagnostics are intentionally transparent deterministic heuristics. `heuristic_score` values are not calibrated probabilities and should not be used for biological decisions. Bubble evidence includes branch-edge support, complete-path read support where available, and the bounded search depth. FASTQ quality summaries are recorded on edges but do not change the heuristic labels or scores.

Bubble detection excludes the divergence node and the alternate branch during each bounded search. This rejects obvious self-loop and cross-branch false positives, but it does not prove that a reconverging structure is a biological variant. Disconnected-component diagnostics likewise report structure without establishing a coverage gap or biological separation. A later version can replace the scorer with a random forest, gradient-boosted model, or graph neural network only after defining and validating appropriate evidence without changing the twin snapshot concept.

Multiple compatible branches with the same divergence and rejoin are grouped into one bubble record. Each record receives a stable `candidate_id` derived from its start node, rejoin node, and lexically canonical branch paths; support counts, timestamps, and input order do not affect it. Nested or overlapping structures with different divergence or rejoin nodes are retained as separate records, so one graph can contain both an outer and an inner candidate. Every related pair receives reciprocal relationship evidence with the related ID, a `nested` or `overlapping` type, and shared internal nodes and edges. A candidate is `nested` when one candidate contains the other's divergence and rejoin; other shared path material is `overlapping`; unrelated candidates have no relationship entries.

Each branch search defaults to a 64-node expansion cap; benchmark runs can raise it up to the hard ceiling of 256. When depth or node expansion stops a search before a complete candidate is found, the tool emits a `bubble_search_limited` diagnostic with the effective depth, cap, affected branches, and whether depth or the node cap stopped expansion. This is an uncertainty signal, not evidence that no bubble exists. Complete candidates include the same depth and per-branch limit flags in their evidence. Cycles are traversed without revisiting nodes and are reported only when they still lead to a valid reconvergence with another branch. Branch-edge support counts reads supporting the first edge only, while path-level support requires one read to contain the complete node path.

## Demo dataset

The synthetic reads include:

1. A main sequence path.
2. An alternate path that creates a bubble.
3. A repeated read that increases support for a repeat candidate.
4. A short low-support error read that creates a suspicious tip.
5. An additional bridge-read batch used to demonstrate a twin update.

The synthetic data is deterministic, small, and safe to inspect during a presentation.

## Two-day prototype plan

### Day 1

1. Load synthetic reads.
2. Build the k-mer graph.
3. Calculate coverage, degrees, and supporting reads.
4. Detect tips, repeats, bounded bubbles, and disconnected components.
5. Produce an explainable classification for each diagnostic.

### Day 2

1. Store an initial twin snapshot.
2. Add a bridge-read batch.
3. Rebuild the graph and store a second snapshot.
4. Compare graph metrics and structural heuristic scores.
5. Export a JSON report and run the automated tests.

## Why this is a digital twin

A static graph shows one assembly state. Successive read batches create cumulative snapshots with graph evidence, input fingerprints, and timestamps, making it possible to ask, "What changed after this batch arrived?" Updates are file-driven and offline; the prototype is not synchronized with a live physical source.

The demo's bridge-read batch simulates an update. The `analyze` command can load FASTA or four-line FASTQ files, including gzip-compressed inputs, as successive batches.

## Suggested five-minute demo

1. Explain that reads are unordered fragments.
2. Run the initial dataset with `k=3`.
3. Show graph nodes, edges, coverage, and detected diagnostics.
4. Point out the low-support tip and the branch/reconvergence bubble.
5. Explain the heuristic score and evidence behind each label.
6. Add the bridge-read batch.
7. Run the update and show the second snapshot.
8. Emphasize that uncertainty is preserved rather than hidden.

## Limitations

- The graph builder is small and educational.
- The diagnostic layer is a set of explainable rules, not a trained biological model.
- The data is synthetic and does not represent a real genome.
- Input orientation is preserved; reverse-complement normalization and bidirected graph semantics are not implemented.
- FASTQ quality is reduced to one mean accuracy value per read; FASTA assumes quality 1.0.
- No contigs are produced and no biological validation is performed; no claim is made about human-genome-scale performance.
- Real deployment would require validated assemblers, calibrated models, provenance controls, privacy, and independent biological validation.

## Research-prototype extensions

The executable now also offers a deterministic `replay` command for ordered FASTA/FASTQ batches. Replay uses logical timestamps beginning at 2000-01-01 UTC and records added read IDs, graph nodes, edges, and newly supported edges in each state delta. The source files remain the event log; snapshots do not copy raw read sequences. Duplicate read IDs across cumulative batches are rejected because support is tracked by read ID.

The `benchmark` command compares the support-based heuristic with reference-supported bubble paths. With `--max-reads`, it selects reads using SHA-256 bottom-k ranking by seed and read ID, making fixed-seed sample sizes nested and independent of input order. The benchmark uses an integer-indexed graph with packed adjacency, support, and quality arrays. Regular analysis and replay continue to keep detailed read evidence. Its report includes diagnostic counts and examples for tips, branch/repeat nodes, disconnected components, bounded bubble searches, and detected candidates. It reports accuracy, coverage, runtime, and Python traced peak allocations (not total process RSS) for candidates with exactly one branch aligned to the reference. Path labeling uses a bounded seed-and-extend edit-distance search on both reference orientations with a configurable identity cutoff; ambiguous, unlabeled, and truncated cases remain explicit. This establishes reference compatibility only, not variant truth or assembly quality. Frozen defaults, subset sizes, and pending evaluation gates are in `assembly_graph_twin/benchmark_protocol.json`.

An optional `path_model.py` module trains a small PyTorch path scorer on graph-derived features from reference-labeled candidates. It calibrates a temperature and selects an automatic-resolution threshold on separate validation reports, then abstains if the validation precision target is not reached. Applying a model adds a proposal to each snapshot while preserving the original graph and alternatives. This is a compact neural baseline, not a validated graph neural network or biological decision system. PyTorch is isolated in `requirements-ml.txt`; core analysis remains standard-library-only.

See `assembly_graph_twin/README.md` for commands, benchmark limitations, data-split guidance, and the current roadmap boundary. A credible benchmark still requires pinned public read/reference files, disjoint held-out genomes or samples, and comparison with an established assembler and reference-based assembly evaluation. No real public dataset or trained checkpoint is bundled.
