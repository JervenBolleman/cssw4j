# TornadoVM kernel (experimental)

A port of the scoring step to [TornadoVM](https://www.tornadovm.org/), for running Smith-Waterman on OpenCL, CUDA or
SPIR-V devices. It is separate from the main cssw4j build: TornadoVM needs its own JDK setup, and this module does not
depend on the Vector API code.

## What it computes

For every pair of a query and a target: the optimal local alignment **score** and its **end positions** on the target
and the query, with the same tie-breaking as the SSW C library. It does not compute begin positions or CIGARs; run the
CPU aligner (`Aligner`) on the hits that matter.

## Why it is not a port of the striped kernel

The striped layout packs one alignment into the lanes of one vector register. A GPU has thousands of threads and no
shared register file to stripe across, so the natural port is inter-sequence parallelism, as in CUDASW++: **one
work-item per query-target pair**, each running a plain scalar affine-gap (Gotoh) recurrence column by column.
The recurrence needs two rows (H and E, one value per query position) per pair; these live in scratch arrays of
`(longest query + 1) * pairs` ints each.

`Layout` chooses where each pair keeps its row in the scratch arrays:

* `INTERLEAVED`: the values of one query position for all pairs are adjacent (coalesced loads on a GPU). Default.
* `CONTIGUOUS`: the values of one pair are adjacent (better for a CPU device).

## Files

* `SwKernel`: the kernel, `scorePairs`. Ordinary Java with one `@Parallel` loop, so it also runs sequentially on the JVM.
* `TornadoScorer`: builds the task graph once and runs it.
* `Matrix`, `Sequences`: NCBI matrix parser and sequence packing (same letter encoding as the C library).
* `Check`: compares kernel results with `src/test/resources/expected/protein_blosum62_o11_e1.tsv.gz`.

## Build and run

The SDKs from sdkman: `7.3.1-jdk22plus-opencl` runs on JDK 22 or newer, `7.3.1-jdk21-opencl` runs on JDK 21 only (with
preview features; the scripts detect this from `$TORNADOVM_HOME/etc/tornado.jdk`).

    export TORNADOVM_HOME=$HOME/.sdkman/candidates/tornadovm/7.3.1-jdk22plus-opencl
    export TORNADO_JDK=/usr/lib/jvm/java-25-temurin-jdk
    tornado/build.sh
    tornado/run-check.sh jvm          # the kernel as plain Java, no device needed
    tornado/run-check.sh tornado      # through TornadoVM on the default device
    LAYOUT=contiguous tornado/run-check.sh tornado

`$TORNADOVM_HOME/bin/tornado --devices` lists the devices. A device needs an OpenCL platform; on a machine without a GPU
driver the `pocl` package provides a CPU device.

## Results

Checked on the 55 x 55 protein pairs of the test dataset (BLOSUM62, gap open 11, gap extend 1): score and both end
positions are identical to the C library for all 3,025 pairs, as plain Java and through TornadoVM, with both layouts.

Timings on the development machine (Intel i7-6700HQ, 8 threads). **No GPU was available:** the NVIDIA driver was not
loaded, so the only OpenCL device was pocl on the CPU. These numbers say nothing about GPU speed.

| 3,025 pairs, scores and end positions | time |
|---|---:|
| cssw4j `Aligner` (striped Vector API, CPU, scores only) | about 0.07 s |
| `SwKernel` as sequential Java, interleaved / contiguous | about 2.4 s / 1.6 s |
| `SwKernel` through TornadoVM on pocl, interleaved / contiguous | about 3.3 s / 2.9 s |

On a CPU device this kernel is about 40 times slower than the striped CPU code, and TornadoVM on pocl is not faster than
a single Java thread. The kernel is meant for a GPU with many pairs per launch; measure it there before relying on it.

## Limits

* Scores equal the exact Smith-Waterman scores for gap open greater than gap extend, as in the C library.
* Scratch memory grows with `pairs * longest query`: 3,025 pairs of up to 1,863 residues use about 45 MB.
* Only scores and end positions; no begin positions or CIGAR.
