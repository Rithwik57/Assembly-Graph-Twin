# Masterclass: Analysis of the Greedy Algorithm for Genome Assembly

Genome assembly is the computational challenge of piecing together millions of short, random DNA fragments (called **reads**) to reconstruct an organism's original, complete genome. 

While theoretical formulations like the **Shortest Common Superstring (SCS)** model assembly using the principle of **parsimony** (finding the shortest single sequence that contains every read as a substring), exact solutions to SCS are **NP-complete** and computationally intractable for large biological datasets. 

To solve this problem in practice, computer scientists rely on the **Greedy Algorithm for Genome Assembly**. This masterclass provides an exhaustive, step-by-step breakdown of how the greedy algorithm works, its computational efficiency, its mathematical worst-case limits, and why it can fail on biological genomes with repeated sequences.

---

## 1. How the Greedy Algorithm Works (Step-by-Step)

### What It Is
The **Greedy Algorithm** is a heuristic approach that constructs a genome sequence by making locally optimal choices at each step—specifically, by repeatedly joining the pair of reads that share the longest remaining **suffix-prefix overlap**.

### Graph Representation
To execute the algorithm, the set of input reads \\(F\\) is represented as an **Overlap Graph**:
* **Vertices (\\(V\\)):** Each read \\(f_i \in F\\) becomes a vertex labeled \\(f_i\\).
* **Directed Edges (\\(E\\)):** A directed arrow from \\(f_1 \to f_2\\) indicates that read \\(f_1\\) precedes read \\(f_2\\) in the assembly.
* **Degrees:** 
  * **Out-degree (\\(\Delta_{\text{out}}\\)):** The number of arrows leaving a vertex.
  * **In-degree (\\(\Delta_{\text{in}}\\)):** The number of arrows entering a vertex.

```
Initial State:           Greedy Step: Pick Max Overlap (e.g., Length 2)
[ Read 1 ]   [ Read 2 ]    [ Read 1 ] ──────────────► [ Read 2 ]
  (Δin=0,      (Δin=0,       (Δin=0,   Overlap = 2    (Δin=1,
  Δout=0)      Δout=0)       Δout=1)                  Δout=0)
```

### The Step-by-Step Process
1. **Initialization:** Place all \\(|F|\\) reads as unconnected vertices in a graph with **0 directed edges**.
2. **Find Maximum Overlap:** Scan all pairs of vertices \\((f_1, f_2)\\) that satisfy two degree constraints:
   * Read \\(f_1\\) has no outgoing edge (\\(\Delta_{\text{out}}(f_1) = 0\\)).
   * Read \\(f_2\\) has no incoming edge (\\(\Delta_{\text{in}}(f_2) = 0\\)).
3. **Add Directed Edge:** Identify the pair with the **maximum suffix-prefix overlap length** and draw a directed edge \\(f_1 \to f_2\\).
4. **Repeat:** Continue adding edges between valid pairs until all vertices are merged into a **single connected component**.
5. **Chain Formation & Output:** Because every vertex is constrained to at most one incoming and one outgoing edge, the resulting graph forms a single linear **chain**. The algorithm spells out the final superstring by merging the overlapping prefixes and suffixes along this path.

---

### Worked Example from the Lecture

Consider an input set \\(F\\) containing 5 short reads:

```
Reads:  f1 = "AAA",  f2 = "AAB",  f3 = "ABB",  f4 = "BBB",  f5 = "BBA"
```

```
Step 1: The algorithm evaluates all pairs and finds a maximum overlap length of 2 
        (e.g., suffix "AA" of "AAA" matches prefix "AA" of "AAB").
        Edge added:  AAA ──(2)──► AAB

Step 2: Next max overlap of 2: suffix "AB" of "AAB" matches prefix "AB" of "ABB".
        Edge added:  AAB ──(2)──► ABB

Step 3: Continue connecting remaining valid vertices with overlap length 2:
        ABB ──(2)──► BBB ──(2)──► BBA

Final Chain: AAA ──► AAB ──► ABB ──► BBB ──► BBA
```

To **spell** the output string, we write the first read in full and append only the non-overlapping tail characters of each subsequent read:

\\[\text{Spelled String} = \text{AAA} + \text{B} + \text{B} + \text{B} + \text{A} = \text{"AAABBBBA"} \quad\\]

Every input read is a contiguous substring of `"AAABBBBA"`:
* `"AAA"` \\(\to\\) indices 1–3
* `"AAB"` \\(\to\\) indices 2–4
* `"ABB"` \\(\to\\) indices 3–5
* `"BBB"` \\(\to\\) indices 4–6
* `"BBA"` \\(\to\\) indices 5–7

> 💡 **In Simple Terms:** Imagine a stack of overlapping word cards. At every step, you find the two cards that share the longest matching sequence of letters at their ends, tape them together, and repeat. When all cards are taped into a single line, you read the sentence from left to right without repeating the overlapping letters.

---

## 2. Question 1: Algorithmic Runtime Complexity Analysis

To analyze how fast the greedy algorithm runs, we define two mathematical parameters for a read set \\(F\\):
* **\\(|F|\\):** The total **count** of reads (number of strings).
* **\\(\|F\|\\):** The **sum of the lengths** of all strings in \\(F\\).

*(Example: For 5 reads of length 3 each, \\(|F| = 5\\) and \\(\|F\| = 15\\)).*

```
              ALGORITHMIC WORKFLOW & TIME COMPLEXITY
  ┌─────────────────────────────────────────────────────────────┐
  │ Step 1: Pre-compute All-Versus-All Overlaps                 │
  │ Tool: Generalized Suffix Trees                              │
  │ Time: O(||F|| + |F|²)                                       │
  └──────────────────────────────┬──────────────────────────────┘
                                 ▼
  ┌─────────────────────────────────────────────────────────────┐
  │ Step 2: Sort All |F|² Overlap Lengths                       │
  │ Method: Descending order sort + Degree tracking             │
  │ Time: O(|F|² log |F|)                                       │
  └──────────────────────────────┬──────────────────────────────┘
                                 ▼
  ┌─────────────────────────────────────────────────────────────┐
  │ Step 3: Edge Insertion & Sequence Spelling                  │
  │ Method: Linear graph updates along the sorted array          │
  │ Combined Time Complexity: O(||F|| + |F|² log |F|)           │
  └─────────────────────────────────────────────────────────────┘
```

---

### Step-by-Step Complexity Derivation

#### Step 1: Pre-computing Overlaps using Suffix Trees
At every iteration, the algorithm needs the longest suffix-prefix overlap between valid pairs. Rather than re-calculating overlaps on the fly, it **pre-computes** the overlaps for all \\(|F| \times (|F| - 1)\\) pairs.

Using a **Generalized Suffix Tree** (a data structure that stores all suffixes of a set of strings), all-pairs suffix-prefix overlaps can be computed efficiently:

\\[\text{Pre-computation Time} = O(\|F\| + |F|^2) \quad\\]

* \\(\|F\|\\) represents the time required to construct the generalized suffix tree.
* \\(|F|^2\\) represents the output array size storing the overlap integer for every pair of reads.

#### Step 2: Sorting Overlaps
Once the \\(|F|^2\\) overlap integers are computed, they are sorted in **non-increasing (descending) order** so that the algorithm can process larger overlaps first.

Sorting an array of size \\(|F|^2\\) using standard comparison sorting (like mergesort) takes:

\\[O(|F|^2 \log(|F|^2)) = O(2 \cdot |F|^2 \log |F|) = O(|F|^2 \log |F|) \quad\\]

#### Step 3: Graph Construction & Bookkeeping
As the algorithm iterates through the sorted overlap array, it uses a **Boolean array** to track the in-degree and out-degree of each vertex in \\(O(1)\\) time per lookup. Valid edges are added until a single chain is built.

---

### Combined Total Complexity
Combining the pre-computation and sorting phases yields:

\\[\text{Total Time Complexity} = O(\|F\| + |F|^2 \log |F|) \quad\\]

Because this running time is a **polynomial function** of the input size, the greedy algorithm is computationally efficient and scales well compared to exact exponential-time algorithms.

> 💡 **In Simple Terms:** Instead of trying every possible combination of reads (which would take longer than the age of the universe), the greedy algorithm builds a lookup table of all overlap lengths, sorts them from largest to smallest, and glues the reads together in order. This turns an impossible problem into one a modern computer can solve quickly.

---

## 3. Question 2: Theoretical Bounds & Worst-Case Analysis

While the greedy algorithm is guaranteed to produce a valid superstring, it is **not guaranteed to produce the shortest possible superstring**. To evaluate its theoretical quality, computer scientists analyze its **approximation ratio** (\\(c\\)):

\\[\text{Length}(\text{Greedy Output } S') \le c \times \text{Length}(\text{Shortest Superstring } S) \quad\\]

---

### History of Theoretical Bounds
Over the past three decades, computer scientists have proven progressively tighter upper bounds on the approximation factor \\(c\\):
* **1991 (Blum et al.):** Proved \\(c \le 4\\) (Greedy output is at most 4 times the length of the optimal SCS).
* **Subsequent Refinement:** Tightened to \\(c \le 3.5\\).
* **Recent Benchmark (Mucha, 2013):** Proved \\(c \le 3.425\\).

---

### The Hard / Worst-Case Example for Greedy

To test the lower limit of \\(c\\), consider a constructed "bad example" with 3 input reads using exponent notation (\\(a^k\\) represents character \\(a\\) repeated \\(k\\) times):

\\[\text{Read 1: } S_1 = a b^k, \quad \text{Read 2: } S_2 = b^k c, \quad \text{Read 3: } S_3 = b^{k+1} \quad\\]

```
                   WORST-CASE ASSEMBLY COMPARISON
                   
1. Greedy Algorithm Execution:
   • Overlap(S1, S2) = bᵏ (Length k)  ──► Greedy connects S1 ──► S2 first
   • Remaining choice: Connect S3 at the end
   • Greedy Path:   (a bᵏ) ──► (bᵏ c) ──► (bᵏ⁺¹)
   • Greedy Output:  a bᵏ c bᵏ⁺¹                                        
   • Greedy Length: 1 + k + 1 + (k + 1) = 2k + 3                       

2. Optimal Shortest Common Superstring (SCS):
   • Optimal String: a bᵏ⁺¹ c                                           
   • Verification:
     - S1 (a bᵏ)   is inside "a bᵏ⁺¹ c"
     - S2 (bᵏ c)   is inside "a bᵏ⁺¹ c"
     - S3 (bᵏ⁺¹)   is inside "a bᵏ⁺¹ c"
   • Optimal Length: 1 + (k + 1) + 1 = k + 3                           
```

---

### Calculating the Approximation Ratio

Comparing the length of the Greedy output to the optimal SCS output:

\\[\text{Ratio} = \frac{\text{Length}(\text{Greedy})}{\text{Length}(\text{Optimal SCS})} = \frac{2k + 3}{k + 3} \quad\\]

As \\(k\\) becomes extremely large (\\(k \to \infty\\)):

\\[\lim_{k \to \infty} \frac{2k + 3}{k + 3} = 2 \quad\\]

This hard example proves that the greedy algorithm can produce an output **twice as long** as the true shortest common superstring.

---

### The Open Conjecture
Because the worst-case example shows a factor of 2, but the best mathematical proof only proves an upper bound of \\(3.425\\), theoretical computer scientists have formulated a famous **open conjecture**:

\\[\text{Conjecture: The Greedy Algorithm has an approximation factor of exactly } c = 2 \quad\\]

Proving or disproving whether \\(c = 2\\) remains an unsolved open problem in computer science.

> 💡 **In Simple Terms:** In the worst possible scenario, a greedy algorithm can fall into a trap and build a sequence that is twice as long as the ideal answer. Computer scientists know it will never be more than 3.425 times as long, and strongly suspect the real limit is capped at 2.

---

## 4. Question 3: Biological Reality & The Repeat Collapse Problem

Even if the greedy algorithm produces a superstring close to the theoretical minimum length, **can a biologist trust the output?** Does the shortest superstring match the true biological genome?.

### The English Text Counterexample
To test this, consider a hypothetical target genome represented as English text containing 3 identical copies of the repeated word `"LONG"`:

\\[\text{True Target Genome } G = \text{"A_LONG_LONG_LONG_TIME"} \quad\\]

Suppose we sample 10 short reads of fixed length \\(L = 6\\) from this genome:

\\[\text{Reads: } \{\text{"NG_LON"}, \text{"_LONG_"}, \text{"ONG_TI"}, \dots\} \quad\\]

If we run the greedy algorithm on these length-6 reads, it produces:

\\[\text{Greedy Output } S' = \text{"A_LONG_TIME"} \quad\\]

```
Target Genome:   A _ L O N G _ L O N G _ L O N G _ T I M E  (3 copies of "LONG")
                                 │
                         [ Greedy Assembly ]
                                 ▼
Greedy Output:   A _ L O N G _ T I M E                       (Only 1 copy of "LONG"!)
```

---

### Why Did the Algorithm Fail?

The output `"A_LONG_TIME"` is a **valid superstring**—every single read of length 6 is present inside it. However, **it is biologically incorrect** because it contains fewer repeat copies than the original genome.

This reveals a fundamental flaw in using Shortest Common Superstring for biology:
1. SCS relies on the **parsimony principle**, which aggressively minimizes total sequence length.
2. Because repeated regions share identical overlapping sequences, minimizing total length **collapses multiple repeat copies into a single copy**.
3. As a result, the greedy algorithm returns a sequence that is too short, leading biologists to misinterpret the true structure of the genome.

---

### The Fix: Increasing Read Length
If we increase the read length from \\(L = 6\\) to \\(L = 8\\), individual reads become long enough to extend past the repeat boundaries into unique flanking sequences. With length-8 reads, the greedy algorithm preserves all 3 repeat copies and successfully reconstructs `"A_LONG_LONG_LONG_TIME"`.

> 💡 **In Simple Terms:** The greedy algorithm is obsessed with making the shortest sequence possible. When it sees identical repeating patterns (like `LONG_LONG_LONG`), it merges them all into one `LONG` to save space. To stop the algorithm from cheating, we must use longer reads that span across the repeat boundaries.

---

## 5. Mathematical Proof: When is Greedy Guaranteed to Reconstruct the True Genome?

Under what formal conditions can a biologist mathematically rely on the greedy algorithm to reconstruct the exact biological genome \\(G\\)?.

---

### Formal Definitions & Assumptions

1. **Error-Free Reads:** Reads contain zero sequencing errors.
2. **Substring-Free Read Set:** No read is completely contained within another read.
3. **Full Coverage:** Reads cover every single base of the target genome \\(G\\).
4. **Genome Substring Notation (\\(S_T^L\\)):** A substring of genome \\(G\\) starting at base coordinate \\(T\\) with length \\(L\\).
5. **Genome Repeat:** A pair of identical substrings \\(S_{T_1}^L = S_{T_2}^L\\) where \\(T_1 \neq T_2\\).
6. **Bridged Substring:** Substring \\(S_T^L\\) is **bridged** if there exists at least one read \\(R\\) that completely covers \\(S_T^L\\) and extends by **at least 1 character to the left** and **at least 1 character to the right**.

```
Genome G:   ... A [  S_T^L  ] C ...
Read R:       │ ┴───────────┴ │       <── Read covers S_T^L PLUS 
              └─ Left   Right ┘           at least 1 base on both sides!
```

7. **Bridged Repeat:** A repeat is **bridged** if **at least one** of its copies in the genome is bridged by a read.

---

### Theorem / Claim
Given a target genome \\(G\\) and an error-free read set \\(F\\), the Greedy Algorithm is **guaranteed to reconstruct the exact genome sequence \\(G\\)** if **every repeat in \\(G\\) is bridged**.

---

### Step-by-Step Proof by Contradiction

Suppose the greedy algorithm fails to reconstruct \\(G\\). Let the connection of read \\(R_i \to R_j\\) with overlap length \\(L\\) be the **very first incorrect edge** added by the algorithm.

There are two mutually exclusive cases for why \\(R_i \to R_j\\) could be wrong:

```
                                PROOF BY CONTRADICTION
                                          │
            ┌─────────────────────────────┴─────────────────────────────┐
            ▼                                                           ▼
  CASE 1: Rj IS the true successor                          CASE 2: Rj is NOT the true successor
  • Greedy picked overlap L, but the                        • Rj is sampled from a distant region.
    true overlap was shorter.                      • Overlap L implies a repeat in G.
  • Implies an unbridged repeat of                          • Repeat is bridged by some read Rk.
    length L in G.                                 • Overlap(Ri, Rk) > L.
  • CONTRADICTION: Assumed all repeats                      • Greedy would have picked Rk over Rj!
    are bridged!                                       • CONTRADICTION!
            │                                                           │
            └─────────────────────────────┬─────────────────────────────┘
                                          ▼
                         GREEDY GUARANTEED TO RECONSTRUCT G
```

#### Case 1: \\(R_j\\) IS the true successor of \\(R_i\\) in \\(G\\), but Greedy miscalculated the overlap length
* The algorithm joined \\(R_i \to R_j\\) using an overlap of length \\(L\\), whereas the true genomic overlap was shorter.
* This shared sequence of length \\(L\\) implies the existence of a repeat of length \\(L\\) in the genome.
* Because \\(R_j\\) is the immediate successor of \\(R_i\\), no read exists that extends across this repeat with flanking characters.
* This implies the repeat is **unbridged**, which directly **contradicts** our assumption that *every* repeat in \\(G\\) is bridged.

#### Case 2: \\(R_j\\) is NOT the true successor of \\(R_i\\) in \\(G\\)
* \\(R_j\\) actually originated from a distant region of the genome.
* The suffix-prefix overlap of length \\(L\\) between \\(R_i\\) and \\(R_j\\) implies the presence of a repeat of length \\(L\\) in \\(G\\).
* By our assumption, this repeat MUST be bridged by some read \\(R_k\\) in the dataset.
* Because \\(R_k\\) bridges this repeat (extending beyond it on both sides), the overlap between \\(R_i\\) and \\(R_k\\) must be **strictly greater than \\(L\\)** (\\(\text{Overlap}(R_i, R_k) > L\\)).
* Because the greedy algorithm always selects the largest remaining overlap, it would have selected \\(R_k\\) over \\(R_j\\).
* Joining \\(R_i \to R_j\\) with overlap \\(L\\) **contradicts** the fundamental operation of the greedy algorithm.

---

### Conclusion of Proof
Because both possible failure cases lead to a mathematical contradiction, the greedy algorithm cannot make a mistake. Thus, **if every repeat is bridged, the greedy algorithm is guaranteed to reconstruct the correct genome sequence \\(G\\)**.

> 💡 **In Simple Terms:** If your sequencing reads are long enough to "bridge" across every repeated region in the DNA (like stepping stones that reach all the way across a river), the greedy algorithm is mathematically proven to assemble the genome perfectly without making a single mistake.

---

## Conceptual Map

```
                     GREEDY GENOME ASSEMBLY ANALYSIS
                                     │
         ┌───────────────────────────┼───────────────────────────┐
         ▼                           ▼                           ▼
1. TIME COMPLEXITY          2. APPROXIMATION BOUNDS      3. BIOLOGICAL VALIDITY
   • Pre-compute Overlaps:     • Compare Greedy (S') to     • Parsimony Fallacy:
     O(||F|| + |F|²) via         Optimal SCS (S)         SCS minimizes length,
     Suffix Trees     • Proven Upper Bound:          collapsing repeats
   • Sort Overlaps:              c ≤ 3.425 (Mucha)     • Text Counterexample:
     O(|F|² log |F|)  • Hard Example: Ratio 2      "A_LONG_TIME" loses
   • Total Time:                 for (2k+3)/(k+3)    repeat copies
     O(||F|| + |F|² log |F|)   • Open Conjecture:           • The Solution:
     (Polynomial)           c = 2                   Longer reads that bridge
                                                              repeats
                                                                 │
                                                                 ▼
                                                        4. BRIDGED REPEAT PROOF
                                                           • Assumptions: Error-free,
                                                             full coverage, substring-
                                                             free, bridged repeats
                                                           • Proof by Contradiction:
                                                             Case 1 & 2 fail
                                                           • Result: Guaranteed
                                                             exact assembly
```

---

## Key Terms Glossary

| Term | Simple Definition |
| :--- | :--- |
| **Greedy Algorithm** | An algorithm that makes the locally optimal choice at each step (picking the largest overlap) to find a global solution. |
| **Overlap Graph** | A network where vertices are reads and directed edges represent suffix-prefix overlap connections. |
| **In-Degree / Out-Degree** | The number of incoming (\\(\Delta_{\text{in}}\\)) or outgoing (\\(\Delta_{\text{out}}\\)) directed edges attached to a vertex. |
| **Generalized Suffix Tree** | A fast tree data structure used to find all suffix-prefix overlaps across a set of strings in linear time. |
| **Approximation Ratio (\\(c\\))** | A factor measuring how much longer an approximate algorithm's output is compared to the optimal shortest answer. |
| **Open Conjecture** | A mathematical proposition that is strongly suspected to be true but has not yet been rigorously proven or disproven. |
| **Repeat Collapse** | A assembly error where multiple identical copies of a repeated sequence are merged into a single copy. |
| **Bridged Substring** | A genomic region covered by a read that extends by at least 1 base past both the left and right boundaries. |
| **Bridged Repeat** | A repetitive sequence where at least one of its genomic copies is covered by a bridging read. |
| **Polynomial Time** | An algorithm running time bounded by \\(O(n^k)\\), making it computationally efficient and scalable. |

---

## The 5 Most Important Things to Remember

1. **Greedy Builds a Linear Chain:** The algorithm repeatedly joins the valid read pair with the maximum suffix-prefix overlap, maintaining \\(\Delta_{\text{in}} \le 1\\) and \\(\Delta_{\text{out}} \le 1\\) to form a single path.
2. **Polynomial Time Efficiency:** Using generalized suffix trees and sorting, the greedy algorithm runs in \\(O(\|F\| + |F|^2 \log |F|)\\) time.
3. **The \\(c = 2\\) Open Conjecture:** The greedy algorithm is proven to be at most \\(3.425\\) times longer than the optimal SCS, but worst-case examples suggest the true limit is \\(c = 2\\).
4. **SCS Collapses Repeats:** Because Shortest Common Superstring minimizes total length, it naturally collapses biological repeats into fewer copies, producing incorrect assemblies.
5. **Bridged Repeats Guarantee Exact Assembly:** If every repeat in a genome is bridged by at least one read, the greedy algorithm is mathematically guaranteed to reconstruct the exact genome.

---

## Common Misconceptions & Confusing Points

* ❌ **Misconception 1: The greedy algorithm always finds the Shortest Common Superstring.**
  * *Correction:* Greedy finds a valid superstring in polynomial time, but local choices can prevent it from finding the global minimum length.
* ❌ **Misconception 2: A shorter superstring is always biologically better.**
  * *Correction:* In biology, shorter superstrings often mean that repetitive DNA regions have been incorrectly collapsed and lost.
* ❌ **Misconception 3: To bridge a repeat, EVERY copy of the repeat must be covered by a bridging read.**
  * *Correction:* A repeat is classified as bridged if **at least one** of its copies in the genome is spanned by a bridging read.
* ❌ **Misconception 4: Sorting overlaps is the most time-consuming part of pre-computation.**
  * *Correction:* Constructing the suffix tree takes \\(O(\|F\| + |F|^2)\\), while sorting takes \\(O(|F|^2 \log |F|)\\). Sorting dominates the asymptotic complexity for large read counts, but both steps are polynomial.

---

## 10 Self-Assessment Questions

### Questions

1. **[Basic Recall]** What two graph degree conditions must be satisfied before the greedy algorithm can connect read \\(f_1 \to f_2\\)?
2. **[Basic Recall]** What is the total combined time complexity of the greedy algorithm for a read set \\(F\\)?
3. **[Calculation]** If a read set contains \\(|F| = 100\\) reads with a total combined length of \\(\|F\| = 10,000\\) bases, how many pairwise overlap entries are stored in the pre-computed overlap array?
4. **[Conceptual]** Why does the greedy algorithm produce a single linear chain of reads rather than a branched graph?
5. **[Conceptual]** In the worst-case example (\\(a b^k\\), \\(b^k c\\), \\(b^{k+1}\\)), why does the greedy algorithm choose to overlap \\(a b^k\\) and \\(b^k c\\) first?
6. **[Theoretical]** What is the currently proven upper bound on the approximation ratio \\(c\\) for the greedy algorithm, and what is the conjectured bound?
7. **[Biological Application]** Why did the greedy algorithm assemble `"A_LONG_LONG_LONG_TIME"` into `"A_LONG_TIME"`?
8. **[Conceptual]** What does it mean for a genomic substring \\(S_T^L\\) to be "bridged" by a read?
9. **[Proof Logic]** In Case 2 of the proof by contradiction, why would the existence of a bridging read \\(R_k\\) prevent the greedy algorithm from connecting \\(R_i \to R_j\\)?
10. **[Synthesis]** How does increasing the read length \\(L\\) during sequencing help resolve both theoretical repeat collapse and practical assembly errors?

---

### Separate Answer Key

1. **Answer:** Read \\(f_1\\) must have no outgoing edge (\\(\Delta_{\text{out}}(f_1) = 0\\)), and read \\(f_2\\) must have no incoming edge (\\(\Delta_{\text{in}}(f_2) = 0\\)).
2. **Answer:** \\(O(\|F\| + |F|^2 \log |F|)\\).
3. **Answer:** \\(|F| \times (|F| - 1) = 100 \times 99 = 9,900\\) pairwise overlap entries (or \\(|F|^2 = 10,000\\) entries in an \\(|F| \times |F|\\) matrix).
4. **Answer:** Because the algorithm restricts connections to vertices with \\(\Delta_{\text{in}} = 0\\) and \\(\Delta_{\text{out}} = 0\\), every node can have at most 1 incoming and 1 outgoing edge, preventing branching.
5. **Answer:** Because \\(a b^k\\) and \\(b^k c\\) share an overlap of length \\(k\\) (the sequence \\(b^k\\)), which is the maximum available overlap length in the read set.
6. **Answer:** The best proven upper bound is \\(c \le 3.425\\) (Mucha, 2013), and the conjectured bound is \\(c = 2\\).
7. **Answer:** Because the SCS formulation optimizes for minimal string length (parsimony), causing the algorithm to collapse 3 identical repeat copies of `"LONG"` into a single copy.
8. **Answer:** A read completely covers the substring \\(S_T^L\\) and extends by at least 1 base past the left boundary and at least 1 base past the right boundary.
9. **Answer:** If \\(R_k\\) bridges the repeat, the overlap between \\(R_i\\) and \\(R_k\\) is strictly greater than \\(L\\). Since greedy selects the largest remaining overlap, it would pick \\(R_k\\) over \\(R_j\\).
10. **Answer:** Longer reads are more likely to extend past repeat boundaries into unique flanking regions, bridging repeats and allowing the greedy algorithm to correctly reconstruct the true genome without collapsing repeat copies.
