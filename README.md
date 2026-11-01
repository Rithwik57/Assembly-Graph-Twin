# Assembly Graph Twin

This project builds an explainable de Bruijn graph from FASTA/FASTQ reads, tracks graph evidence across offline batches, and benchmarks bounded bubble-path choices against reference-compatible paths. It is a research prototype, not a complete assembler.

## Project layout

- `assembly_graph_twin.py`, `benchmark.py`, `path_model.py`: command-line application and implementation modules.
- `tests/`: unit and behavior checks.
- `fixtures/`: small deterministic synthetic fixtures used for positive-control checks.
- `config/`: benchmark protocol and configuration.
- `data/reads/`, `data/references/`: locally downloaded datasets and references; they are not required to run the synthetic demo.
- `reports/demo/`, `reports/checks/`: current generated reports.
- `reports/legacy_benchmarks/`: earlier benchmark reports made with the forward-only graph. Keep these separate from the strand-aware results.

Run commands from this project folder, or prefix paths with `assembly_graph_twin/` from the repository root.

## Strand-aware de Bruijn graph

Graph construction now defaults to `--strand-mode bidirected`. It represents the bidirected graph as a directed double-cover: for each observed transition `u -> v`, it adds the valid reverse-complement twin `rc(v) -> rc(u)`. It never reverses an isolated edge without complementing and reversing its endpoints. The mirror transition carries the same physical read support, so it is not counted as a second read. Reverse-complement views of a bubble share one canonical candidate ID and are reported once.

This follows the standard double-stranded/bidirected de Bruijn graph model used in assembly work, including BCALM 2 and Cuttlefish 2, while keeping the prototype’s existing path traversal straightforward. The representation is not a compact canonical-k-mer graph and can use more memory than a canonical node/port implementation. Reports show both oriented node/edge counts and counts collapsed across reverse-complement pairs; component and complete-candidate summaries collapse mirror views. Bounded-search-limit events are still counted per oriented search, so mirror searches may both appear. Support counts refer to oriented graph transitions; do not sum a transition and its reverse-complement twin as independent molecules. [BCALM 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC4908363/) describes bidirected graph compaction for double-stranded DNA; [Cuttlefish 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9454175/) uses canonical k-mer forms in a bidirected edge-centric graph.

Use `--strand-mode forward` only to reproduce old forward-only behavior. New benchmark comparisons should use the default `bidirected` mode and record it with the run.

## Run the synthetic demo

```powershell
python assembly_graph_twin.py demo --output reports/demo/demo_report.json
```

The demo creates two snapshots from synthetic reads, then prints graph summaries, the synthetic bubble's divergence/rejoin and paths, and branch support. This is a controlled positive example, not a result from a real dataset. The benchmark continues to report zero bubbles when none are found in the supplied reads.

## Check the benchmark pipeline

```powershell
python check_benchmark_fixture.py --output reports/checks/known_bubble_fixture_report.json
```

The positive control has known branches and a tiny reference. The check repeats the same seeded run and fails if no candidate or unique reference-compatible path is found, or if candidate evidence is nondeterministic. These are proxy labels on a synthetic fixture, not biological measurements.

Run unit checks with:

```powershell
python -m unittest discover -s tests -p "test_*.py" -v
```

## Analyze and replay reads

```powershell
python assembly_graph_twin.py analyze data/reads/SRR31691735_1.fastq.gz --k 31 --output reports/demo/sa_analysis.json
python assembly_graph_twin.py replay batch_001.fastq.gz batch_002.fastq.gz --k 31 --output reports/demo/replay.json
```

FASTA/FASTQ gzip input is supported. Read identifiers must be unique across cumulative batches. FASTQ quality is descriptive metadata. Replay is offline and deterministic for unchanged input files and settings.

## Reference benchmark

The frozen defaults and validation gates are in `config/benchmark_protocol.json`. Candidate labels are **reference-compatibility proxies**, not variant truth. The lightweight seed-and-extend path aligner accepts paths at the configured identity threshold and leaves ambiguous or unsupported candidates unlabeled. It is not a production read aligner.

Example using the locally saved *S. aureus* sample and its reference:

```powershell
python assembly_graph_twin.py benchmark --reads data/reads/SRR31691735_1.fastq.gz --reference data/references/sa_reference.fa --k 31 --max-reads 100 --seed 1 --dataset-id sa-atcc43300 --output reports/benchmarks/sa-atcc43300-100.json
```

With `--max-reads`, SHA-256 bottom-k sampling gives nested, input-order-independent subsets for a fixed seed. Reports include tips, branch/repeat signals, disconnected components, bounded-search limits, candidate evidence, labels, runtime, and Python traced peak memory. They also include a reference-backed read/edge validation summary: reads are placed using sampled exact-seed anchors on both strands, uniquely placed read bases and estimated reference breadth are reported, and graph edges traversed by uniquely placed reads are counted.

The read placements and breadth are explicitly **estimates**, not gapped alignments or base-level consensus. The index samples reference seeds every 12 bases and queries reads at multiple offsets to limit Python memory and runtime; it is not a production mapper. A zero-bubble run can still report reference placement, breadth, and edge support. A zero-candidate run has undefined bubble accuracy and coverage, not a zero-error score. Increase read limits cautiously and watch memory.

Benchmark accuracy measures path choice among candidates with exactly one reference-compatible branch. It is a graph-level baseline, not an assembly evaluation; the tool does not emit contigs. Held-out species, comparison with a pinned assembler, and reference-based assembly evaluation remain pending gates. Multiple subsets from one genome do not establish species-level generalization.

## Optional PyTorch k-mer backend and Google Colab

The default graph implementation remains standard-library CPU code. An optional PyTorch backend runs k-mer window encoding, per-read edge deduplication, read/node-support aggregation, and compact adjacency construction on CUDA, then transfers each distinct graph edge once into the compact CPU graph. This avoids sending every read-edge observation through Python one by one. Use `--graph-backend torch --device cuda` to request CUDA explicitly; it errors if CUDA is unavailable. `--device auto` selects CUDA when available and otherwise uses PyTorch on CPU.

K-mer encoding, per-read deduplication, edge/node support aggregation, and adjacency CSR construction run on the GPU. String-keyed graph materialization, bubble searches, diagnostics, reference placement, and path labeling remain CPU work, so the overall speedup depends on which stage dominates a run. Progress messages identify the active phase; the report records time spent in selection, graph construction, diagnostics, and reference validation, plus peak CUDA memory. The notebook disables Python allocation tracing for runtime; use `--no-python-memory-tracing` in other large runs if tracing overhead is too high (the Python peak-memory metric is then null). The reference anchor length is independent of graph `k` (up to 15 bases), and the index is sparse to avoid excessive Python allocations. This does not promise exponential acceleration, and GPU memory does not replace system RAM. The PyTorch dependency stays optional; default CLI behavior is unchanged.

Example:

```powershell
python assembly_graph_twin.py benchmark --reads data/reads/SRR31691735_1.fastq.gz --reference data/references/sa_reference.fa --k 15 --max-reads 50 --seed 1 --dataset-id sa-atcc43300-k15-gpu --graph-backend torch --device cuda --output reports/benchmarks/sa-atcc43300-k15-gpu.json
```

The report records `graph_backend`, `graph_device`, and `graph_construction_seconds` so CPU and CUDA runs can be compared.

The Colab notebook is a hosted runner for the existing CLI; it does not duplicate or replace the graph implementation:

The Colab notebook is a hosted runner for the existing CLI; it does not duplicate or replace the graph implementation:

```text
notebooks/colab_benchmark.ipynb
```

Copy the `assembly_graph_twin` project folder to Google Drive, open the notebook in Colab, and run its setup and path-configuration cells. The defaults expect the project at `MyDrive/assembly_graph_twin` and use the bundled *S. aureus* read/reference paths. Change `PROJECT_DIR`, `READS_PATH`, or `REFERENCE_PATH` in the notebook if your Drive layout differs. The notebook saves reports under `reports/benchmarks/colab/`.

Begin with the 50-read smoke run, inspect runtime and peak memory, and only then run the optional 250-read cell. Select **Runtime > Change runtime type > T4 GPU** before running the notebook. It requests the PyTorch CUDA backend and displays the detected GPU. The notebook invokes `assembly_graph_twin.py benchmark` through a subprocess, so the command-line workflow remains available and unchanged.

## Optional PyTorch path scorer

The standard graph workflow has no non-standard dependencies. To install the optional scorer:

```powershell
python -m pip install -r requirements-ml.txt
```

The model is a small multilayer perceptron over graph-derived path features, not a full message-passing GNN or a validated assembler. Training and calibration require enough uniquely labeled candidates from disjoint genomes. Until that evidence exists, model development remains deferred.
