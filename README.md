# cssw4j

A pure Java port of the [Complete Striped Smith-Waterman (SSW) library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
using the Java Vector API (`jdk.incubator.vector`). No JNI, no native code.

How the code works, step by step from plain Smith-Waterman: [docs/HOW_IT_WORKS.md](docs/HOW_IT_WORKS.md).

It reproduces `ssw_init` / `ssw_align` / `mark_mismatch`: optimal score, sub-optimal score, begin and end positions,
CIGAR (BAM packed) and the `=`/`X`/soft-clip CIGAR with edit distance.

## Usage

```java
ScoringMatrix m = ScoringMatrix.BLOSUM50;          // or ScoringMatrix.blosum62(), nucleotide(2, 2), parse(reader)
Profile query = Profile.of("MVHLTPEEKSAVTALWGKVNVDEVGG", m); // build once per query, immutable, shareable
Aligner aligner = new Aligner();                     // owns reusable scratch memory: one per thread
byte[] target = m.encode("MVHLTGEEKAAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFF");

Alignment a = aligner.align(query, target, 10, 2);  // gap open 10, extend 2, full alignment
a.score1(); a.refBegin1(); a.refEnd1(); a.cigarString();
Cigar.Marked x = Cigar.markMismatch(a, target, m.encode("..."), query.length());
```

The full C signature is available as
`align(profile, ref, refLen, gapOpen, gapExtend, flag, filters, filterd, maskLen)` with `Aligner.FLAG_*` constants.

Requirements: JDK 25 or later, run with `--add-modules jdk.incubator.vector`. Tested on JDK 25 and JDK 27.0.1.

## Design

* Same algorithm as the C code: Farrar striped query profile, 8 bit kernel with fallback to 16 bit on overflow,
  reverse pass to find the begin, banded global alignment for the CIGAR.
* Vector width is the CPU's preferred shape (32 byte lanes on AVX2, 64 on AVX-512) instead of SSE2's 16.
  `-Dcssw4j.vectorBits=128` selects the C library's width.
* The 8 bit kernel keeps scores in excess-128 form (score − 128 as a signed byte), so it only needs signed saturating
  add/sub and signed max. The profile holds raw scores, which saves the per-cell bias subtraction of the C code.
* No allocation in steady state besides the result: an `Aligner` reuses its buffers between calls.
* Lane shifts (`_mm_slli_si128`) are done on 32 bit lanes (`vpermd` + shifts): a byte `slice()` compiles to ~30
  instructions on AVX2.
* Gap penalty vectors are loaded from a small array instead of `broadcast()`, because C2 otherwise re-broadcasts them
  for every use in the hot loop.
* The banded traceback stores one packed direction byte per cell and computes directions without branches.

## Correctness tests

`src/test/resources/dataset` holds the test dataset:

* `uniprot_human.fasta`: 24 reviewed human proteins from UniProtKB/Swiss-Prot (103 to 1,863 aa; globins,
  kinases, p53, BRCA1, CFTR, EGFR, ...).
* `uniprot_homologs.fasta`: 48 reviewed homologs and unrelated proteins from mouse, zebrafish, yeast and E. coli.
* `synthetic_queries.fasta` / `synthetic_targets.fasta`: generated deterministically by
  `src/test/scripts/make_synthetic.py`. They cover lengths around SIMD lane boundaries (1 to 129), random sequences,
  mutated copies with indels, low complexity repeats, ambiguity codes `B Z X J U O *`, lower case, chimeras, and a
  near copy of EGFR that forces the 16 bit kernel.
* `dna_*.fa`: the DNA demo data of the C library.

`src/test/resources/expected/*.tsv.gz` are the C library's results for all 55 × 55 protein pairs, under BLOSUM50 with
gaps 3/1 and 10/2 and BLOSUM62 with 11/1, plus the 100 DNA pairs. They are produced by `src/test/c/ssw_reference.c`;
regenerate them with `src/test/scripts/generate_expected.sh` (set `SSW_SRC` to the C library's `src` directory).

* `CReferenceTest` compares every field. `mvn test` runs it twice: with the preferred vector width, and with 128 bits
  where all fields, including the sub-optimal score, are bit-identical to C. With wider vectors only the sub-optimal
  score and its position can differ: the C code includes the padding cells of the last segment in its per-column
  maxima, and the padding depends on the lane count.
* `ScalarOracleTest` checks scores against a plain scalar Gotoh implementation on 2,000 random and mutated protein
  pairs, and checks that every CIGAR re-scores to the reported score.

## Performance

`src/test/c/ssw_bench.c` and `swiss.sib.cssw4j.Benchmark` align the protein dataset all against all (BLOSUM50, gaps
10/2, 20 rounds). Results on an i7-6700HQ (AVX2), C built with `-O3 -march=native`, Java on Temurin 25:

| GCUPS                          | C SSW | cssw4j |
|--------------------------------|------:|-------:|
| whole dataset, scores only     |  4.67 |   4.90 |
| whole dataset, with CIGAR      |  1.95 |   1.90 |
| UniProt queries, scores only   |  4.62 |   5.20 |
| one 1,000 × 5,000 pair         |  6.09 |   8.92 |

Short queries are slower than C (16 aa: 0.60 vs 0.91 GCUPS; 64 aa: 1.69 vs 2.26), because with 32 lanes the
per-column fixed costs are spread over very little work. 128 bit vectors are faster below roughly 200 residues.

```
mvn test-compile
java --add-modules jdk.incubator.vector -cp target/classes:target/test-classes swiss.sib.cssw4j.Benchmark 20 10 2
```

## Known limitations (inherited from the C library)

* Scores are exact only for `gapOpen > gapExtend`. With linear or inverted gap penalties the lazy-F early exit can stop
  too early: on random proteins with gaps 3/3, both C and this port report a lower than optimal score for about 30%
  of pairs.
* The banded traceback sometimes fails (`flag == 1`, no CIGAR). The C library hits this on 10% of the BLOSUM50 3/1
  protein pairs in the dataset; the port fails on exactly the same pairs.
* Scores above 32,767 saturate in the 16 bit kernel.

## Native library (GraalVM Native Image)

cssw4j can also be built as a C shared library, `libcssw4j.so`, with GraalVM Native Image. The C API
(`swiss.sib.cssw4j.capi.CApi`) mirrors the SSW functions: `cssw4j_profile_new`, `cssw4j_aligner_new`, `cssw4j_align`,
and the matching `*_free` functions. Native Image writes the header `libcssw4j.h` next to the library.

Build (needs a GraalVM JDK; the Java JDK 25 build of the library is unchanged):

    GRAALVM_HOME=/path/to/graalvm mvn -Pnative package     # writes target/native/libcssw4j.so and libcssw4j.h
    # or, without Maven: src/native/build-native.sh target/classes target/native

Extra Native Image options can be given in `NATIVE_IMAGE_ARGS`, e.g. `NATIVE_IMAGE_ARGS=-O3`.

Check against the C reference results (the same data as the JVM tests):

    src/test/scripts/test_native.sh

Status, measured on the same machine as the Performance section (Oracle GraalVM 25.4.4, Native Image default `-O2`):

* **Correct:** every alignment of the protein and DNA reference sets matches the C library in score, positions, flag
  and CIGAR (3,025 + 100 pairs).
* **Slow:** the native kernel runs at about 0.02 GCUPS, against about 10 GCUPS for the same code on the JVM (long
  1,000 × 5,000 pair, scores only). `-O3` doubles the native speed (0.03 GCUPS); GraalVM CE gives the same result as
  Oracle GraalVM. Simple vector loops are fast natively (about 1.3 ns per vector operation, about the same as the JIT),
  so the loss is specific to the large alignment kernel, which Native Image does not compile to SIMD code.
  The cause has not been pinned down further.
* So for now use the library through the JVM when speed matters. The native build is useful for C/C++ callers that
  need correct results and can accept the speed.

## JDK notes

* `ByteVector.reduceLanes(UMAX)` is wrong on JDK 25 (it returns 128 for vectors whose lanes are all below 128), and on
  JDK 27 (27-ea and 27.0.1) C2 crashes (`ShouldNotReachHere` in `ReductionNode::opcode`) when compiling a loop with a
  loop-carried `ByteVector.lanewise(UMAX, ...)`. The excess-128 kernel avoids all unsigned operators.
* `VectorOperators.SUADD` is declared `Binary` in JDK 25 and `Associative` in JDK 27, so code compiled against one and
  run on the other fails with `NoSuchFieldError`. cssw4j does not use it.

## Credits and license

cssw4j is a port of the SSW library by **Mengyao Zhao** and the other SSW developers (Wan-Ping Lee, Yongan Zhao,
Daniel Cameron, Nan Xiao and contributors), who open-sourced it under the MIT license:
https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library. Please cite their paper when you use this code:

> Zhao M, Lee W-P, Garrison EP, Marth GT (2013). SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in Genomic
> Applications. PLoS ONE 8(12): e82138. https://doi.org/10.1371/journal.pone.0082138

The striped algorithm and lazy-F loop are by Michael Farrar (swsse2, BSD license); the SSW variant of the lazy-F loop
comes from SWPS3 (ETH Zürich, MIT license). See [docs/HOW_IT_WORKS.md](docs/HOW_IT_WORKS.md#credits) for full
references.

cssw4j is released under the MIT license, the same license as the SSW library. See [LICENSE](LICENSE), which also
reproduces the BSD notice of Michael Farrar carried by the SSW source.
