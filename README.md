# Assembly Graph Twin

Assembly Graph Twin is a beginner-friendly research prototype for turning DNA read files into an explainable de Bruijn graph, tracking graph evidence, and checking graph paths against a reference genome. It is not a production genome assembler and does not produce finished contigs.

This guide walks from a fresh clone through a small local demo, downloading real reads and a matching reference, then running the Colab notebook on a T4 GPU.

## 1. What you need

- Windows 10/11, macOS, or Linux for cloning and local setup.
- Git installed.
- Python 3.10 or newer installed and available as `python`.
- A Google account with Google Drive and Google Colab access for the T4 run.

The local demo and standard graph tools use only Python's standard library; no `pip install` is needed. The notebook uses the PyTorch already available in Colab. The downloaded dataset files are intentionally not stored in GitHub because read files can be large.

## 2. Clone the project

Open Command Prompt or PowerShell and run:

```cmd
git clone https://github.com/Rithwik57/Assembly-Graph-Twin.git
cd Assembly-Graph-Twin
python --version
```

The project folder contains `assembly_graph_twin.py`, `benchmark.py`, `gpu_backend.py`, `notebooks/`, `data/`, and `reports/`. If `python --version` is older than 3.10, install a newer Python release and reopen the terminal.

## 3. Run the built-in demo first

The demo uses synthetic reads, so it works before downloading any data:

```cmd
python assembly_graph_twin.py demo --output reports\demo\demo_report.json
```

It prints two graph snapshots and an example synthetic bubble. This confirms that Python can run the project; it is not a result from biological data.

## 4. Download real reads and a matching reference

Use one read dataset and the reference for the same organism/strain. The project does **not** include either file after cloning. Downloaded files belong in `data/reads/` and `data/references/`; the existing `.gitignore` keeps them local and out of future Git commits.

### Recommended first real-data run: compact *S. aureus* sample

The Colab notebook is preconfigured for run `SRR31691735` and its *Staphylococcus aureus* reference. In PowerShell, from the cloned project folder, download the FASTQ link from ENA and the reference FASTA from NCBI:

```powershell
$api = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR31691735&result=read_run&fields=run_accession,fastq_ftp&format=tsv"
$report = Invoke-RestMethod -Uri $api
$row = ($report -split "`n" | Select-Object -Skip 1 -First 1).Trim()
$readUrl = (($row -split "`t")[1] -split ";")[0]
$readUrl = $readUrl -replace "^ftp://", "https://"
if ($readUrl -notmatch "^https?://") { $readUrl = "https://$readUrl" }
curl.exe -L --fail --retry 3 $readUrl -o "data\reads\SRR31691735_1.fastq.gz"

curl.exe -L --fail --retry 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP178839.1&rettype=fasta&retmode=text" -o "data\references\sa_reference.fa"
```

ENA lists `SRR31691735` as an Oxford Nanopore *S. aureus* run; the reference accession used by this project is `CP178839.1`. [ENA run record](https://www.ncbi.nlm.nih.gov/sra/SRX27053554%5Baccn%5D), [NCBI reference sequence](https://www.ncbi.nlm.nih.gov/nuccore/CP178839.1).

Check that both files exist:

```cmd
dir data\reads\SRR31691735_1.fastq.gz
dir data\references\sa_reference.fa
```

### Optional master-plan anchor: *E. coli* K-12 MG1655

For the planned *E. coli* benchmark, download the Loman Lab ONT reads and matching reference instead. The read FASTA is about 1.5 GB, so allow time and disk space:

```cmd
curl.exe -L --fail --retry 3 "https://s3.climb.ac.uk/nanopore/E_coli_K12_1D_R9.2_SpotON_2.pass.fasta" -o "data\reads\ecoli_reads.fasta"
curl.exe -L --fail --retry 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" -o "data\references\ecoli_reference.fa"
```

The *E. coli* read file and reference must be used together. For the notebook, change `READS_PATH`, `REFERENCE_PATH`, and `DATASET_ID` in its configuration cell to these `ecoli_*` names. [Loman Lab read dataset](https://lab.loman.net/2016/07/30/nanopore-r9-data-release/), [NCBI reference NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3).

## 5. Optional local benchmark check

From the project folder, the smaller *S. aureus* subset can be checked on your computer:

```cmd
python assembly_graph_twin.py benchmark --reads data\reads\SRR31691735_1.fastq.gz --reference data\references\sa_reference.fa --k 15 --max-reads 50 --seed 1 --dataset-id sa-atcc43300-smoke --output reports\benchmarks\sa-smoke.json
```

`--max-reads 50` controls how many reads are used to build the graph. To select a reproducible subset, the program still scans the input file, so this command can take a while even with 50 selected reads. The benchmark reports diagnostics and reference-compatibility proxy labels; it does not claim biological variant truth.

## 6. Put the project and downloaded data in Google Drive

The notebook runs the project code from Google Drive and expects the data files to be in the project’s `data/` folders.

1. In Google Drive, create or open `MyDrive`.
2. Upload the entire cloned `Assembly-Graph-Twin` folder, including the downloaded files under `data/reads/` and `data/references/`.
3. For the default *S. aureus* notebook settings, rename the uploaded project folder to `assembly_graph_twin`, or plan to edit `PROJECT_DIR` in the notebook.
4. In Drive, open `assembly_graph_twin/notebooks/colab_benchmark.ipynb`. Choose **Open with → Google Colaboratory** if prompted.

The `.gitignore` only affects Git. It does not prevent you from uploading your local read/reference files to Drive. The *E. coli* FASTA is much larger than the *S. aureus* FASTQ; uploading or processing it can take longer.

## 7. Select the T4 GPU and configure the notebook

In Colab, select **Runtime → Change runtime type → T4 GPU**, save, and reconnect if Colab asks. Then run the notebook cells from the top, in order:

1. **Mount Google Drive**: approve the Drive access prompt.
2. **Configure project and input paths**: the defaults are:
   - `PROJECT_DIR = /content/drive/MyDrive/assembly_graph_twin`
   - `READS_PATH = .../data/reads/SRR31691735_1.fastq.gz`
   - `REFERENCE_PATH = .../data/references/sa_reference.fa`
   - `GRAPH_BACKEND = "torch"` and `DEVICE = "cuda"`
3. **Preflight**: confirm the project files and protocol version all show `OK`; it also prints the detected GPU name. If your Drive folder is named differently, change `PROJECT_DIR` before running this cell. If you use the *E. coli* files, change the two input paths and `DATASET_ID` here too.
4. **Run the benchmark**: leave the initial `MAX_READS = 50` smoke run in place. The cell streams progress messages as it scans reads and advances through graph construction, diagnostics, and reference checks.

At the end, look for `CUDA confirmed`, `graph_backend: torch`, and `graph_device: cuda`. The report should show nonzero peak CUDA memory. Those confirm that CUDA operations ran on the T4. Increase to the notebook's separate 250-read run only after the smoke run finishes and you have reviewed its runtime and memory.

### What the T4 accelerates

CUDA performs k-mer encoding, per-read edge deduplication, edge/node support aggregation, and adjacency-index construction. Python graph-key creation, bubble searches, diagnostics, reference placement, and path labeling still run on the CPU. The report gives stage timings and CUDA memory use. This is a mixed CPU/GPU pipeline; a T4 does not guarantee every run will be faster, and GPU memory does not replace system RAM.

## Other commands

Analyze one read file:

```cmd
python assembly_graph_twin.py analyze data\reads\SRR31691735_1.fastq.gz --k 31 --output reports\demo\sa_analysis.json
```

Replay ordered read batches (each file should contain distinct read IDs):

```cmd
python assembly_graph_twin.py replay batch_001.fastq.gz batch_002.fastq.gz --k 31 --output reports\demo\replay.json
```

## How to read benchmark results

- A bubble is a graph structure with diverging paths that reconverge; it is not automatically a biological variant.
- Candidate labels mean a path is compatible with the supplied reference under the prototype's approximate aligner. They are not validated variant truth.
- A zero-candidate result means no complete bubble was reported under the selected settings; it does not mean the genome has no variation.
- The program is a research prototype and does not output a finished assembly.
- The optional path model is not required to run the demo, analysis, replay, or benchmark workflows.

## Graph and strand model

Graph construction defaults to `--strand-mode bidirected`. It represents reverse-complement transitions as a directed double-cover: for an observed transition `u -> v`, it adds the valid twin `rc(v) -> rc(u)`. It does not reverse an isolated edge without also reversing and complementing its endpoints. The twin carries the same physical read support and is not counted as another molecule.

This keeps traversal straightforward but is not a compact canonical-k-mer or node-side-port graph. The representation can use more memory than a canonical implementation. Reports include oriented and reverse-complement-collapsed counts where available. Do not sum support for a transition and its twin as if they were independent reads. Use `--strand-mode forward` only to reproduce earlier forward-only behavior.

## Optional path model

The standard workflows need no third-party packages. The optional path scorer has separate dependencies listed in `requirements-ml.txt`:

```cmd
python -m pip install -r requirements-ml.txt
```

It is a small multilayer perceptron over graph-derived path features, not a full message-passing GNN or a validated assembler. Training and calibration require enough uniquely labeled candidates from disjoint genomes; model evaluation remains a later research milestone.
