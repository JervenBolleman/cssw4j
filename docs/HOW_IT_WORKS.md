# How cssw4j works

This document explains what the code does and why, starting from the plain Smith-Waterman algorithm and adding one
idea at a time. It follows the design of the SSW library (see [Credits](#credits)), which cssw4j ports to Java.

Contents:

1. [The problem: local alignment](#1-the-problem-local-alignment)
2. [Plain Smith-Waterman](#2-plain-smith-waterman)
3. [Gaps that cost more to open than to extend](#3-gaps-that-cost-more-to-open-than-to-extend)
4. [The three steps of SSW](#4-the-three-steps-of-ssw)
5. [Doing many cells at once: SIMD](#5-doing-many-cells-at-once-simd)
6. [The striped layout](#6-the-striped-layout)
7. [The inner loop](#7-the-inner-loop)
8. [The lazy-F loop](#8-the-lazy-f-loop)
9. [8 bits first, 16 bits if needed](#9-8-bits-first-16-bits-if-needed)
10. [Where does the best alignment end?](#10-where-does-the-best-alignment-end)
11. [The second best score](#11-the-second-best-score)
12. [Where does it start? The reverse pass](#12-where-does-it-start-the-reverse-pass)
13. [The path: banded alignment and traceback](#13-the-path-banded-alignment-and-traceback)
14. [CIGAR strings](#14-cigar-strings)
15. [What is different in the Java version](#15-what-is-different-in-the-java-version)
16. [Limits](#16-limits)
17. [Map of the code](#17-map-of-the-code)
18. [Credits](#credits)

---

## 1. The problem: local alignment

We have two sequences: a **query** (also called the read) and a **reference** (the target). We want the part of the
query and the part of the reference that are most alike. The parts may be anywhere in the two sequences, and they may
contain gaps (letters present in one sequence but not in the other).

```
reference:  ...A C G T T A C...
                | |   | | | |
query:         A C - T T A C
                   ^
                   a gap: the reference has a G that the query does not
```

"Most alike" is measured with a score:

* each pair of letters gets a score from a **substitution matrix** (for proteins usually BLOSUM50 or BLOSUM62;
  for DNA simply "+2 if equal, -2 if not"),
* each gap costs a penalty.

The best **local alignment** is the one with the highest total score. Smith and Waterman showed in 1981 how to find it
exactly.

## 2. Plain Smith-Waterman

Put the reference along the top and the query down the side, and fill a table `H`. Cell `H[j][i]` holds the best score
of any alignment that **ends** at query letter `j` and reference letter `i`.

Each cell looks at three neighbours:

```
             reference letter i-1   reference letter i
                     |                    |
query letter j-1 -- [ diag ]           [ up   ]
                            \             |
                             \            v
query letter j   -- [ left ] -------> [ H[j][i] ]
```

* from `diag`: extend the alignment with the pair (query j, reference i), adding the matrix score,
* from `up`: a gap in the reference (an extra query letter, an **insertion**),
* from `left`: a gap in the query (an extra reference letter, a **deletion**),
* or start fresh at 0. This "never below zero" rule is what makes the alignment *local*.

The answer is the largest number anywhere in the table.

Worked example: reference `ACGTTAC`, query `ACTTAC`, match +2, mismatch -2, gap open 3, gap extend 1 (the DNA defaults
of the C library's `ssw_test`):

```
            A   C   G   T   T   A   C      <- reference
        0   0   0   0   0   0   0   0
    A   0  [2]  0   0   0   0   2   0
    C   0   0  [4]-[1]  0   0   0   4
    T   0   0   1   2  [3]  2   0   1
    T   0   0   0   0   4  [5]  2   1
    A   0   2   0   0   1   2  [7]  4
    C   0   0   4   1   0   1   4  [9]     <- best score 9
    ^
  query
```

Following the bracketed cells back from 9 gives the alignment: `AC`, one deletion (the reference `G`), then `TTAC`.
Six matches (6 × 2 = 12) minus one gap opening (3) = 9. In CIGAR notation this is `2M1D4M` (see section 14).

For a query of length `m` and a reference of length `n` the table has `m × n` cells. A 400 letter protein against a
1,000 letter protein is 400,000 cells; a query against a whole database is billions. All the work in SSW is about
filling this table fast.

## 3. Gaps that cost more to open than to extend

Real gaps tend to be a few letters long, not many single-letter gaps. So opening a gap costs more than making an
existing one longer. In SSW a gap of length `k` costs:

```
gapOpen + (k - 1) * gapExtend
```

(the first letter of a gap pays `gapOpen`, every next letter pays `gapExtend`).

To handle this, Gotoh (1982) keeps two extra tables next to `H`:

* `E[j][i]`: best score ending at this cell **inside a horizontal gap** (a deletion),
* `F[j][i]`: best score ending at this cell **inside a vertical gap** (an insertion).

```
E[j][i] = max( E[j][i-1] - gapExtend,      keep the gap going
               H[j][i-1] - gapOpen )       open a new gap

F[j][i] = max( F[j-1][i] - gapExtend,
               H[j-1][i] - gapOpen )

H[j][i] = max( 0,
               H[j-1][i-1] + score(reference[i], query[j]),
               E[j][i],
               F[j][i] )
```

Everything below is a faster way to compute exactly these three formulas.

## 4. The three steps of SSW

SSW splits the job into three steps, from cheap to expensive, and only does the expensive ones when asked:

```
 query + reference
        |
        v
 +--------------------------+   fast, vectorized (sections 5 to 11)
 | 1. forward pass          |   -> best score, where it ENDS,
 |    over the whole table  |      second best score
 +--------------------------+
        |
        v
 +--------------------------+   same kernel, but on a small table: the
 | 2. reverse pass          |   reversed query prefix against the reversed
 |                          |   reference prefix (section 12)
 +--------------------------+   -> where the best alignment STARTS
        |
        v
 +--------------------------+   only the region between start and end,
 | 3. banded alignment +    |   only cells near the diagonal (section 13)
 |    traceback             |   -> the path (CIGAR)
 +--------------------------+
```

The forward pass does not store the whole table; it keeps only one column at a time. That is why it is fast and uses
little memory, and also why it cannot tell where the alignment starts or what path it takes. Steps 2 and 3 recover
that.

In the code: `Aligner.align` runs the three steps; `swByte` and `swWord` are the kernels for steps 1 and 2;
`bandedSw` and `bandRow` do step 3.

## 5. Doing many cells at once: SIMD

Modern CPUs have **vector registers**: one register holds many small numbers (called **lanes**), and one instruction
works on all lanes at once.

```
            lane: 0   1   2   3   ...  31
register a:     [ 5 | 1 | 7 | 0 | ... | 2 ]
register b:     [ 3 | 4 | 2 | 9 | ... | 8 ]
max(a, b):      [ 5 | 4 | 7 | 9 | ... | 8 ]      one instruction, 32 results
```

A 256 bit register (AVX2) holds 32 bytes or 16 shorts. If we can compute 32 cells of the table with one instruction,
we can go up to 32 times faster.

The catch: cells depend on their neighbours. We may only compute cells together if none of them needs another one's
result. In a column, each cell needs the cell above it (through `F`), so we cannot simply put a column of neighbours
into one register. The striped layout solves this.

The C library uses SSE2 (128 bits: 16 bytes or 8 shorts). cssw4j uses the Java Vector API and takes the widest
registers the CPU has (`Lanes.java`): 32 byte lanes on AVX2, 64 on AVX-512.

## 6. The striped layout

Michael Farrar (2007) proposed to cut the query into `lanes` equal pieces and give each piece its own lane.
For a query of length `m` and `L` lanes, each piece is `segLen = ceil(m / L)` letters long.

Example: query of 10 letters `q0 .. q9`, 4 lanes, so `segLen = 3`. The positions past the end (10 and 11) are padding.

```
query positions 0..9 cut into 4 pieces of 3:

  piece for lane 0: q0 q1 q2
  piece for lane 1: q3 q4 q5
  piece for lane 2: q6 q7 q8
  piece for lane 3: q9 -- --      (padding)

stored as segLen = 3 vectors, one position from every piece per vector:

              lane 0   lane 1   lane 2   lane 3
  vector 0:  [  q0   |   q3   |   q6   |   q9   ]
  vector 1:  [  q1   |   q4   |   q7   |   --   ]
  vector 2:  [  q2   |   q5   |   q8   |   --   ]
```

Within one vector, the lanes are far apart in the query (`segLen` positions), so they do not depend on each other.
Walking from vector 0 to vector 2 walks *down* each piece at the same time. When we pass the end of a piece, the next
query position is at the top of the *next* lane:

```
lane 0 bottom (q2)  --->  lane 1 top (q3)
lane 1 bottom (q5)  --->  lane 2 top (q6)
```

This hand-over is done by **shifting** the vector one lane up (lane k moves to lane k + 1, lane 0 gets "zero").

### The query profile

For every column we need `score(reference[i], query[j])` for all `j`. Looking these up one by one would be slow, so
before aligning we build a **query profile** (`Profile.java`): for each possible reference letter, the scores against
the whole query, already in striped order.

```
profile row for reference letter 'A':
  vector 0:  [ s(A,q0) | s(A,q3) | s(A,q6) | s(A,q9) ]
  vector 1:  [ s(A,q1) | s(A,q4) | s(A,q7) |    0    ]
  vector 2:  [ s(A,q2) | s(A,q5) | s(A,q8) |    0    ]
profile row for 'C': ...
profile row for 'G': ...
```

During alignment, column `i` simply reads the profile row of letter `reference[i]`: one vector load per segment.
The profile is built once per query and reused for every reference (`Profile.of`).

## 7. The inner loop

The forward pass walks over the reference, one column at a time. For each column it walks over the `segLen` vectors.
All arrays below hold `segLen` vectors:

* `hLoad`: `H` of the previous column,
* `hStore`: `H` of this column (being written),
* `pvE`: `E` values carried from column to column.

For one column, in words:

```
vH = last vector of the previous column, shifted one lane up     (the diagonal for vector 0)
vF = zero

for each vector j of the column:
    vH = vH + profile[reference[i]][j]          diagonal + score
    vH = max(vH, E[j], vF)                      best of the three ways in (and 0)
    store vH as this column's H[j]
    remember the column maximum

    E[j] = max(E[j] - gapExtend, vH - gapOpen)  for the next column
    vF   = max(vF  - gapExtend, vH - gapOpen)   for the next vector down
    vH   = previous column's H[j]               becomes the diagonal for vector j + 1
```

Picture of what one step uses, for lane 1 of the example (query positions q3, q4, q5):

```
                previous column   this column
vector 0 (q3):       H  ------\      H
                               \
vector 1 (q4):       H          \->  H  <- E from the left (pvE)
                                     ^
                                     F from above (vF, from vector 0)
```

The arithmetic **saturates**: with 8 bit unsigned numbers, `3 - 5` gives `0` (not a wrap-around), and `250 + 10`
gives `255`. Saturating at 0 is exactly the "never below zero" rule of Smith-Waterman, for free.

## 8. The lazy-F loop

There is one problem left. In the inner loop, `vF` moves *down* the vectors of one lane. But a vertical gap can also
cross from the bottom of one lane's piece to the top of the next lane's piece:

```
lane 0 piece: q0 q1 q2  |  a gap that starts at q1 should continue
lane 1 piece: q3 q4 q5  |  into q3, q4 ... but q3 was computed in vector 0,
                        v  before we knew about it
```

Farrar's trick: assume this rarely matters, and fix it up afterwards (**lazily**). After the inner loop, shift `vF`
one lane up (so the bottom of lane 0 feeds the top of lane 1) and run down the vectors again, only updating `H` where
`F` is larger:

```
repeat up to L times:
    vF = vF shifted one lane up
    for each vector j:
        H[j] = max(H[j], vF)
        vF   = vF - gapExtend
        if in no lane vF > H[j] - gapOpen:   stop, nothing more to fix
```

In practice this stops after one or two vectors, because a gap that has already lost `gapExtend` many times rarely
beats a fresh `H`. The early stop is what makes the method fast. (Section 16 explains the one case where it stops
too early.)

The SSW version (taken from SWPS3) does not update `E` in this loop. That means an insertion may not be directly
followed by a deletion. With normal gap penalties this practically never changes the best score.

## 9. 8 bits first, 16 bits if needed

A byte lane can only count up to 255. Using bytes gives 32 lanes per AVX2 register; shorts (16 bits) give only 16.
So SSW first tries with bytes (`swByte`):

* If the best score stays below the limit, done, at full speed.
* If any cell reaches the limit (overflow), the score is reported as 255 and `Aligner.align` runs the whole pass again
  with 16 bit lanes (`swWord`), which go up to 32,767.

Short or distant proteins usually fit in bytes. Long, similar proteins (for example human against mouse orthologs) do
not and take the 16 bit path.

Substitution scores can be negative, but unsigned bytes cannot. The C library adds a **bias** (the most negative score
in the matrix, e.g. 5 for BLOSUM50) to every profile value and subtracts it again after each add. The Java version
uses a different trick for the same problem, see section 15.

## 10. Where does the best alignment end?

During the forward pass two things are kept per column:

* the column maximum, saved in `maxColumn[i]` (used for the second best score),
* if the overall best score went up in this column, a copy of the whole column (`hMax`) and the column number `i`.

After the pass, the end on the reference is that column, and the end on the query is the first query position in the
saved column holding the best score. Because the column is in striped order, position `p` is found at vector
`p % segLen`, lane `p / segLen`.

If two places have the same best score, SSW keeps the first column (smallest reference position) and, within it, the
smallest query position.

## 11. The second best score

It is often useful to know whether the query also matches somewhere else (for example a repeat). SSW gives a cheap
estimate: the largest column maximum *outside* a window of `maskLen` columns around the best end.

```
column:      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
maxColumn:   2  3  5  4  3  9 12 30 41 38 20  7  5 16 18  9  3
                                     ^
                                     best end (column 8)
with maskLen = 4 (real use needs at least 15; 4 keeps the picture small), columns 4..12 are masked out:
             2  3  5  4  .  .  .  .  .  .  .  .  .  16 18  9  3
second best = 18 (column 14)
```

`maskLen` must be at least 15, otherwise no second best score is reported. `ssw_test` uses half the query length.

## 12. Where does it start? The reverse pass

The best alignment ends at (query `e_q`, reference `e_r`). If we read both sequences **backwards** from that point,
the same alignment now *ends* where it originally *started*. So we run the same kernel again on:

* the query prefix `query[0 .. e_q]`, reversed (a new, small profile),
* the reference `reference[0 .. e_r]`, walked from right to left,

and stop as soon as a column reaches the known best score (the `terminate` argument).

```
forward:   reference  . . . . [ s . . . . . e ] . . .
                                            ^ best end found

reverse:   reference  . . . . [ s . . . . . e ]
                                ^ <---------- walk backwards from e
                                start = where the score is reached
```

The end of the reversed alignment gives the start of the original one.

## 13. The path: banded alignment and traceback

Now we know the box the alignment lives in: query `[start_q .. e_q]` × reference `[start_r .. e_r]`. To find the
path (which letters are matched, where the gaps are) we fill the table again for that box only, and this time we keep
a **direction** for every cell, so we can walk back from the end.

To save more work, only cells near the diagonal are filled (a **band**). The band starts as wide as the difference in
length between the two parts, plus one. If the best score is not reached inside the band, it is doubled and tried
again.

```
          reference part
       +--------------------+
       |\\\\                |
 query |  \\\\\             |   only the cells marked \ are computed:
 part  |     \\\\\          |   a band of width w around the diagonal
       |         \\\\\      |
       |             \\\\\  |
       |                 \\\|
       +--------------------+
```

The band is computed one query row at a time (`bandRow`), with the same E/F/H formulas as section 3. Each cell
stores one byte with three pieces of information:

```
bit:     4        3        2 1 0
      +--------+--------+---------+
      | F open | E open |   H     |
      +--------+--------+---------+
H: where the best value came from: 1 = diagonal, 2/3 = from E, 4/5 = from F
E open / F open: 1 if that gap was just opened here, 0 if it was extended
```

The traceback starts at the best cell and follows the directions back to the start, writing down `M` (diagonal),
`I` (insertion, moving up) and `D` (deletion, moving left).

As a final check, the path is re-scored. If it does not give the best score from step 1, the alignment is tried once
more with a full width band. If that still fails, the alignment is returned without a path and `flag = 1` (this also
happens in the C library, see section 16).

## 14. CIGAR strings

A **CIGAR** string describes the path as runs of operations:

```
reference:  A C G T T A C
query:      A C - T T A C
            M M D M M M M      ->  2M 1D 4M  ->  "2M1D4M"
```

* `M`: a pair of letters (match or mismatch),
* `I`: a letter only in the query,
* `D`: a letter only in the reference.

Like in the BAM file format, each run is packed in one `int`: the length in the high 28 bits and the operation in the
low 4 bits (`Cigar.java`).

`Cigar.markMismatch` (the C `mark_mismatch`) rewrites `M` into `=` (equal letters) and `X` (different letters), adds
soft clips (`S`) for the parts of the query outside the alignment, and counts the edit distance (mismatches plus
inserted and deleted letters). For the example above: `2=1D4=`, edit distance 1.

## 15. What is different in the Java version

cssw4j follows the C code closely, and its tests check that the results are the same (section 17). The differences
are in *how* things are computed, not *what* is computed.

### Wider vectors

The C code always uses 128 bits. cssw4j asks the JVM for the best vector size of the CPU (`Lanes.java`), so on AVX2 it
does 32 byte cells per instruction instead of 16. This can be changed with `-Dcssw4j.vectorBits=128`.

### Excess-128 bytes instead of a bias

The C code keeps scores as unsigned bytes and uses unsigned saturating instructions. The Java Vector API has unsigned
operations, but in JDK 25 and 27 some of them are broken or crash the JIT compiler (see `JDK_BUG_REPORTS.md`).

So the 8 bit kernel stores every score `s` as the signed byte `s - 128`:

```
real score:      0     1     2   ...   127   128  ...  255
stored byte:  -128  -127  -126  ...    -1     0   ...  127
```

Comparing and saturating stored bytes with the ordinary *signed* instructions gives exactly the same results as the
unsigned instructions on the real scores. The bonus: the profile can hold the plain (possibly negative) matrix scores,
so the "add bias, subtract bias" step of the C code disappears.

### Cheaper lane shift

Shifting a vector one lane up (the `_mm_slli_si128` of the C code) has no direct counterpart in the Vector API. The
obvious `slice()` compiles to about 30 instructions. `Lanes.shiftBytes` instead treats the vector as 32 bit numbers:
each number is shifted left by 8 bits and the byte that falls off the top of the previous number is moved in, which
is about 7 instructions.

```
 int lanes:     [ d c b a ][ h g f e ]           (bytes, little end first)
 shift by 8:    [ c b a 0 ][ g f e 0 ]
 carry in:      [ c b a Z ][ g f e d ]           Z = the "zero" fill, d = top byte of the previous int
```

### Other small things

* Gap penalty vectors are loaded from a small array instead of being created with `broadcast`, because the JIT would
  otherwise rebuild them inside the hot loop.
* An `Aligner` keeps its work arrays between calls, so aligning does not allocate memory except for the result. One
  `Aligner` per thread; a `Profile` can be shared.
* The traceback direction is one byte per cell (the C code uses three), and the directions are computed without
  `if`/`else`, which avoids mispredicted branches.
* A path that runs outside the alignment box is treated as a failed path instead of being read out of bounds.

### When results can differ from C

With wider vectors, only the **second best score** (section 11) can differ. The padding positions at the end of the
striped query (section 6) are included in the column maxima, and the amount of padding depends on the number of lanes.
With 128 bit vectors (16 lanes, like the C code) every field is identical to C.

## 16. Limits

These come from the original algorithm and are the same in the C library:

* **Gap open must be larger than gap extend.** The early stop of the lazy-F loop (section 8) assumes that a fresh gap
  costs more than continuing one. With `gapOpen <= gapExtend` it can stop too early, and the score can come out lower
  than the true best score.
* **The path is sometimes missing** (`flag == 1`). For some alignments the banded step does not find a path with the
  exact score. With BLOSUM50 and very cheap gaps (3/1) this happens for about 1 in 10 protein pairs; with typical gap
  penalties (10/2, 11/1) it is rare.
* **Scores above 32,767** do not fit in the 16 bit kernel.

## 17. Map of the code

```
src/main/java/swiss/sib/cssw4j/
  ScoringMatrix.java   substitution matrices (BLOSUM50, BLOSUM62, DNA, parser) and letter encoding
  Profile.java         the striped query profile (section 6), 8 and 16 bit
  Lanes.java           vector size choice and the lane shift (sections 5 and 15)
  Aligner.java         align(): the three steps (section 4)
                         swByte / swWord   forward and reverse pass (sections 7 to 12)
                         bandedSw / bandRow  banded alignment and traceback (section 13)
  Alignment.java       the result (scores, positions, CIGAR, flag)
  Cigar.java           CIGAR packing, printing and markMismatch (section 14)

src/test/
  resources/dataset/   proteins from UniProtKB/Swiss-Prot, generated edge cases, DNA demo data
  resources/expected/  results of the original C library on that data
  c/ssw_reference.c    the program that produced those results
  c/ssw_bench.c        C speed benchmark
  java/.../CReferenceTest.java    compares every field with the C results
  java/.../ScalarOracleTest.java  compares scores with a plain, slow Smith-Waterman
  java/.../Benchmark.java         Java speed benchmark
```

## Credits

cssw4j is a port. The algorithm, the structure of the code and most of the details are the work of others:

* **Mengyao Zhao** wrote the SSW library (Complete Striped Smith-Waterman Library) and open-sourced it under the MIT
  license, together with **Wan-Ping Lee** (C++ wrapper), **Yongan Zhao** (Python wrapper), **Daniel Cameron**
  (Java JNI wrapper), **Nan Xiao** (R package) and the other SSW developers.
  Source: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
  Paper: Zhao M, Lee W-P, Garrison EP, Marth GT (2013). *SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in
  Genomic Applications.* PLoS ONE 8(12): e82138. https://doi.org/10.1371/journal.pone.0082138
* **Michael Farrar** invented the striped layout and the lazy-F loop, and released his implementation (swsse2) under
  the BSD license. The core loop of SSW is based on it.
  Farrar M (2007). *Striped Smith-Waterman speeds database searches six times over other SIMD implementations.*
  Bioinformatics 23(2): 156-161.
* **SWPS3** (Adam Szalkowski, Christian Ledergerber, Philipp Krähenbühl and Christophe Dessimoz, ETH Zürich), MIT
  licensed, is the source of the SSW version of the lazy-F loop.
  Szalkowski A et al. (2008). *SWPS3 - fast multi-threaded vectorized Smith-Waterman for IBM Cell/B.E. and x86/SSE2.*
  BMC Research Notes 1: 107.
* The underlying methods: Smith TF, Waterman MS (1981). *Identification of common molecular subsequences.*
  J Mol Biol 147: 195-197, and Gotoh O (1982). *An improved algorithm for matching biological sequences.*
  J Mol Biol 162: 705-708.

The test data contains protein sequences from UniProtKB/Swiss-Prot (https://www.uniprot.org), and BLOSUM62 and the DNA
demo files as distributed with the SSW library.

cssw4j is released under the MIT license, the same license as the SSW library; see [LICENSE](../LICENSE).
