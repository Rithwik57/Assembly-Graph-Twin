# Practical Pangenome Graph Construction and Visualization: A Hands-On Masterclass

Pangenome graphs represent a major shift in genomics, replacing single linear reference genomes with rich network structures that capture the genetic diversity of entire species. This masterclass provides a complete breakdown of constructing, parsing, and visualizing a pangenome graph for a highly complex genomic region using industry-standard tools **AGC**, **Minigraph-Cactus**, and **Bandage**.

---

## 1. The Target Biological Region: Why the C4 Genes on Chromosome 6?

### What It Is
The **C4 (Complement Component 4)** genes are located on human **Chromosome 6** and encode proteins responsible for the human immune response.

### Why It Matters
The C4 gene region is **highly polymorphic**—meaning its DNA sequence, structural arrangement, and gene copy numbers vary drastically from person to person. Because a single linear reference genome cannot capture these structural variations, the C4 region is a classic target for pangenomic analysis.

```
Individual A:  [ Conserved Tail ] ────── Path Option 1 (C4A) ────── [ Conserved Tail ]
Individual B:  [ Conserved Tail ] ────── Path Option 2 (C4B) ────── [ Conserved Tail ]
Individual C:  [ Conserved Tail ] ── Path Option 3 (C4A + C4B) ── [ Conserved Tail ]
```

### Real-World Analogy
Imagine a navigation app mapping a city where 90% of drivers take the same main highway, but at one complex junction (the C4 region), every driver chooses a different combination of overpasses, underpasses, and side streets. A single line drawn on a map fails to show all available routes; you need a multi-lane, branching network map to represent all driver paths accurately.

> 💡 **In Simple Terms:** The C4 region controls immune defense and varies so drastically between individuals that comparing people requires a multi-route pangenome graph rather than a single reference sequence.

---

## 2. Overview of the Software Toolchain

Constructing and visualizing a pangenome graph requires three specialized software tools:

```
1. AGC (Archive Extraction) ──► 2. Minigraph-Cactus (Graph Construction) ──► 3. Bandage (Network Visualization)
```

### 1. AGC (Assembled Genome Compressor)
* **What it is:** A specialized bioinformatic archive utility designed to compress, store, and manage multiple assembled genome sequences.
* **Why it matters:** Raw FASTA files for dozens of human haplotypes require large amounts of disk space. AGC stores them in a compressed `.agc` format and extracts individual FASTA files on demand.
* **Real-World Analogy:** AGC is like a `.zip` archive manager built specifically for DNA sequence files.

### 2. Minigraph-Cactus
* **What it is:** An advanced pangenome graph construction pipeline that aligns multiple sequence haplotypes relative to a chosen reference sequence.
* **Why it matters:** It converts flat, disconnected FASTA files into a unified **GFA (Graphical Fragment Assembly)** graph file.
* **Real-World Analogy:** Minigraph-Cactus is like an architect taking individual road maps from 50 drivers and stitching them into a single master map that highlights shared highways and unique detours.

### 3. Bandage (Bioinformatic Application for Navigating De Novo Assembly Graphs Easily)
* **What it is:** A graphical user interface (GUI) desktop application designed to render GFA text files as interactive visual networks.
* **Why it matters:** Raw GFA files contain thousands of lines of text that are impossible to inspect by eye. Bandage turns text data into a visual map of nodes and edges.
* **Real-World Analogy:** Bandage is like an image viewer that turns raw HTML layout code into a visual web page.

---

## 3. Step-by-Step Practical Workflow & Execution

### Step 1: Environment Setup & Binary Installation
To run this pipeline, pre-compiled binaries for **AGC** and **Minigraph-Cactus** are downloaded directly from their GitHub releases and saved into the system's executable path (`/bin`). 

Because Minigraph-Cactus requires specialized dependencies, it is executed inside a dedicated **Python virtual environment**. In high-performance research environments, these jobs are submitted to remote computing clusters using job schedulers such as **PBS (Portable Batch System)**.

```
Download Binaries (GitHub) ──► Place in /bin ──► Activate Virtual Environment ──► Run via PBS Job Scheduler
```

---

### Step 2: Extracting Sequence Files with AGC
The demonstration dataset is retrieved from a public **Zenodo** repository associated with the draft human pangenome reference paper. AGC unpacks the compressed `.agc` archive to extract individual FASTA files for each C4 gene haplotype into the working directory.

---

### Step 3: Preparing the Manifest & Running Minigraph-Cactus
Minigraph-Cactus requires two primary inputs:

1. **Tab-Delimited Manifest File:** A mapping text file where each row contains a unique sequence name and its corresponding FASTA file path, separated by a tab character (`\t`).
2. **Reference Sequence Backbone:** A designated sequence along which all other haplotypes are aligned. In this demonstration, the **C4 CHM13** haplotype is selected as the reference backbone.

```
Sample_ID_1    /path/to/sample1.fasta
Sample_ID_2    /path/to/sample2.fasta
Sample_ID_3    /path/to/sample3.fasta
```

Running the pipeline on the localized C4 gene region takes approximately **3 to 4 minutes** on a remote server, producing an output **GFA graph file**.

---

### Step 4: Deconstructing the GFA File Format
The generated GFA file stores the graph network using three primary line types:

```
S   node_101   ATCGGCTAGCT...      <── Segment Line (Node ID + Sequence)
L   node_101   +   node_102   +    <── Link Line (Edge between Nodes)
W   GRCh38     0   node_101,node_102... <── Walk Line (Haplotype Trajectory)
```

* **`S` Lines (Segments / Nodes):** Define discrete DNA sequence blocks, listing the segment ID and its corresponding nucleotide string.
* **`L` Lines (Links / Edges):** Define directed connections (edges) between pairs of segments.
* **`W` Lines (Walks / Paths):** Store the exact ordered list of nodes traversed by a specific individual haplotype (such as GRCh38 or CHM13) through the graph.

> 💡 **In Simple Terms:** A GFA file is a text blueprint where `S` lines list the rooms (nodes), `L` lines list the hallways between them (edges), and `W` lines list the guided tours through the building (haplotype paths).

---

### Step 5: Visualizing and Path-Tracing in Bandage

#### 1. Graph Rendering
Loading the GFA file into Bandage and clicking **"Draw Graph"** renders the visual network. 

#### 2. Observed Graph Topology
The resulting graph structure displays **two long linear tails** connected by a **complex central network of loops**. The long tails represent conserved regions shared identically across individuals, while the central loops represent polymorphic C4 gene variants.

```
[ Long Conserved Tail ] ─────── ( Central Polymorphic Loops ) ─────── [ Long Conserved Tail ]
```

#### 3. Highlighting a Specific Haplotype Path (GRCh38)
To trace how a specific individual's genome (e.g., the GRCh38 reference haplotype) navigates through the graph:
1. Extract the `W` (walk) line corresponding to **GRCh38** from the GFA file.
2. Run a Python script to convert bracketed node formatting (e.g., `(node1)(node2)`) into a comma-separated format (`node1,node2`).
3. Paste the comma-separated node list into Bandage's node search text box.
4. Apply a distinct highlight color to trace the exact trajectory taken by the GRCh38 haplotype through the graph.

```
Raw Walk String: (101)(102)(105) ──(Python Script)──► Formatted String: 101,102,105 ──► Highlight in Bandage
```

#### 4. Comparison & Tool Limitations
* **Validation:** The visualized path covers both conserved tails and specific central loops, matching benchmark outputs from the published Human Pangenome Reference Consortium paper.
* **Tool Limitation:** Bandage is designed for localized gene-level graphs (like C4). It cannot render massive whole-genome pangenome graphs due to memory and visual layout limits.

---

## 4. Mental / Conceptual Map

```
                     PANGENOME GRAPH CONSTRUCTION & VISUALIZATION
                                          │
       ┌──────────────────────────────────┴──────────────────────────────────┐
       ▼                                                                     ▼
1. TARGET DATASET                                                      2. PIPELINE EXECUTION
   • C4 Genes on Chromosome 6                                         • AGC: Unpack FASTA files
   • Immune response function                                         • Manifest File: Sample \t Path
   • Highly polymorphic region                                        • Reference: C4 CHM13 backbone
   • Dataset sourced from Zenodo                                      • Minigraph-Cactus run (~3-4 min)
                                                                             │
       ┌─────────────────────────────────────────────────────────────────────┘
       ▼
3. OUTPUT PARSING & VISUALIZATION
   • GFA File Format: S (Nodes), L (Edges), W (Haplotype Walks)
   • Bandage GUI: "Draw Graph" rendering
   • Observed Topology: Conserved tails flanking central loops
   • Path Tracing: Convert walk brackets to commas via Python script
   • Highlight GRCh38 trajectory in Bandage
   • Bandage Limitation: Localized subregions only (cannot render full genomes)
```

---

## 5. Key Terms Glossary

| Term | Simple Meaning |
| :--- | :--- |
| **C4 Genes** | Genes on human chromosome 6 that control immune response and vary greatly between people. |
| **Polymorphic Region** | A genomic region whose sequence, structure, or copy number differs significantly across individuals. |
| **AGC** | Assembled Genome Compressor; a tool used to compress and extract genome FASTA files. |
| **Minigraph-Cactus** | A pipeline that aligns multiple sequence haplotypes to build a pangenome graph in GFA format. |
| **GFA Format** | Graphical Fragment Assembly format; a text file storing nodes (`S`), edges (`L`), and paths (`W`). |
| **Segment (`S` Line)** | A line in a GFA file representing a graph node and its DNA sequence. |
| **Link (`L` Line)** | A line in a GFA file representing an edge between two graph nodes. |
| **Walk (`W` Line)** | A line in a GFA file tracing a specific haplotype's path through graph nodes. |
| **Bandage** | A GUI desktop tool used to visually render and interact with GFA assembly graphs. |
| **PBS** | Portable Batch System; a job scheduler used to execute software on remote computing servers. |

---

## 6. The 5 Most Important Things to Remember

1. **C4 Genes Exemplify Polymorphism:** Located on chromosome 6, C4 genes control immune response and serve as an ideal target for pangenomics due to extreme structural variation.
2. **Minigraph-Cactus Requires a Manifest and Reference:** The construction tool takes a tab-separated file mapping samples to FASTA paths and aligns them against a reference backbone (e.g., C4 CHM13).
3. **GFA Files Store Graph Structure:** GFA text files store graph networks using `S` lines for segment nodes, `L` lines for link edges, and `W` lines for haplotype walks.
4. **Haplotype Tracing Requires Formatting:** To highlight a specific path (e.g., GRCh38) in Bandage, bracketed walk strings from the GFA file must be converted into comma-separated lists using Python.
5. **Bandage Has Scale Limitations:** Bandage is effective for visualizing localized subregion graphs (like C4) but cannot render large whole-genome graphs due to memory constraints.

---

## 7. Common Misconceptions & Confusing Points

* ❌ **Misconception 1: Minigraph-Cactus builds pangenome graphs without needing a reference sequence.**
  * *Correction:* Minigraph-Cactus explicitly requires a designated reference sequence (e.g., C4 CHM13) along which all other sample sequences are aligned.
* ❌ **Misconception 2: Bandage can be used to visualize an entire human pangenome graph at once.**
  * *Correction:* Bandage cannot handle whole-genome pangenome graphs due to memory and visual layout constraints; it is intended for localized regions like C4.
* ❌ **Misconception 3: `S` lines in a GFA file represent connections between sequences.**
  * *Correction:* `S` lines represent **Segments** (nodes) containing the actual DNA sequence text. Connections between nodes are stored in `L` (Link) lines.
* ❌ **Misconception 4: Walk strings from GFA files can be pasted directly into Bandage without modification.**
  * *Correction:* Walk strings contain bracketed formatting `(node1)(node2)` that must be converted into a comma-separated format `node1,node2` before Bandage can highlight the path.

---

## 8. 10 Self-Assessment Questions

### Questions

1. **[Basic Recall]** On which human chromosome are the C4 genes located, and what is their primary biological function?
2. **[Basic Recall]** What three primary line types define the network structure in a GFA file?
3. **[Tool Function]** What is the role of AGC in the pangenome graph construction workflow?
4. **[Procedural]** What two column values are required in the input manifest file for Minigraph-Cactus, and how must they be separated?
5. **[Conceptual]** Which specific haplotype sequence was selected as the reference backbone for aligning the C4 gene sequences in this demonstration?
6. **[Data Formatting]** Why must a Python script be used when preparing a GFA walk string for visualization in Bandage?
7. **[Graph Interpretation]** Describe the visual graph topology observed when loading the C4 GFA file into Bandage. What do the tails and loops represent biologically?
8. **[Tool Limitation]** Why was a localized region like the C4 gene chosen for this demonstration instead of an entire human chromosome?
9. **[System Administration]** What job scheduling system was used to execute the pipeline on the remote server, and what environment was created to manage software requirements?
10. **[Synthesis]** Trace the complete data flow of a DNA sequence from a compressed archive to a color-highlighted path on a visual graph.

---

### Separate Answer Key

1. **Answer:** Human Chromosome 6; they encode proteins responsible for the human immune response.
2. **Answer:** `S` (Segment) lines, `L` (Link) lines, and `W` (Walk) lines.
3. **Answer:** AGC compresses and manages assembled genome files, allowing fast extraction of specific haplotype FASTA files on demand.
4. **Answer:** A sample/sequence identifier and its corresponding FASTA file path, separated by a tab character (`\t`).
5. **Answer:** The C4 CHM13 haplotype.
6. **Answer:** Bandage requires a comma-separated list of node IDs (`node1,node2`) to highlight a path, whereas raw GFA walk strings store nodes in bracketed format `(node1)(node2)`.
7. **Answer:** The topology consists of two long linear tails flanking central loops. The tails represent conserved sequence regions shared across individuals, while the loops represent polymorphic C4 gene variants.
8. **Answer:** Bandage cannot render massive whole-genome pangenome graphs due to memory and rendering limits; localized subregion graphs like C4 are small enough to render smoothly.
9. **Answer:** The PBS (Portable Batch System) job scheduler was used on the server, and a Python virtual environment was created to manage dependencies.
10. **Answer:** 
    1. AGC extracts the FASTA sequence file from a `.agc` archive.
    2. Minigraph-Cactus aligns the FASTA sequence against the CHM13 reference to build a `.gfa` graph file.
    3. The `.gfa` file defines sequence nodes (`S`), edges (`L`), and haplotype walks (`W`).
    4. Bandage loads the GFA file and renders the node-and-edge network.
    5. A Python script reformats the `W` walk string into a comma-separated node list.
    6. Bandage searches the node list and applies a custom highlight color to display the haplotype's path.
