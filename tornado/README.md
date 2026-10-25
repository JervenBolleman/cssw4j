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
driver the `pocl` package provides a CPU device. Use `--jvm=-Dcssw4j.device=driver:index` to pick one (see Results).

## Results

Checked on the 55 x 55 protein pairs of the test dataset (BLOSUM62, gap open 11, gap extend 1): score and both end
positions are identical to the C library for all 3,025 pairs, as plain Java, through TornadoVM on the Quadro M1000M
(OpenCL) and on the pocl CPU device, with both layouts.

Also checked through the TornadoVM CUDA backend on the Quadro (CUDA 12.8 libraries, see Build and run), both layouts.

Timings: 3,025 pairs = 292 million DP cells. Development machine: Intel i7-6700HQ, NVIDIA Quadro M1000M (Maxwell,
2 GB, driver 580.178.04), OpenCL from NVIDIA and pocl. Times are the median of runs 2 to 4 of `Check` (kernel plus
transfers of the results to the host; the first run also includes compilation and is not counted).

| Implementation | Device | time | GCUPS |
|---|---|---:|---:|
| `Aligner`, striped Vector API (scores only, BLOSUM50 10/2) | CPU, 1 thread | about 0.065 s | 4.49 |
| `SwKernel` as sequential Java (scores and end positions) | CPU, 1 thread | 1.6 s (contiguous) | 0.18 |
| `SwKernel` via TornadoVM OpenCL | pocl, CPU | 2.76 s (contiguous) | 0.106 |
| `SwKernel` via TornadoVM OpenCL | pocl, CPU | 3.49 s (interleaved) | 0.084 |
| `SwKernel` via TornadoVM OpenCL | Quadro M1000M | 5.23 s (interleaved) | 0.056 |
| `SwKernel` via TornadoVM OpenCL | Quadro M1000M | 6.07 s (contiguous) | 0.048 |
| `SwKernel` via TornadoVM CUDA | Quadro M1000M | 4.55 s (interleaved) | 0.064 |
| `SwKernel` via TornadoVM CUDA | Quadro M1000M | 5.58 s (contiguous) | 0.052 |

The striped CPU code is about 40 times faster than the best TornadoVM run, and the Quadro is slower than the pocl CPU
device for this kernel. The kernel does a scalar recurrence and reads and writes its H and E rows in global memory for
every cell. Why the Quadro is slower was not profiled; the old Maxwell chip and that memory traffic are the likely
causes. The interleaved layout, which suits the GPU in principle, is faster on the Quadro and slower on the CPU.

Selecting a device: `tornado --devices` lists `driver:index`. Pass `--jvm=-Dcssw4j.device=0:0` (the Quadro, the default) or
`--jvm=-Dcssw4j.device=0:1` (the pocl CPU) to `run-check.sh`, with the matching `TORNADOVM_HOME`. The CUDA SDK uses
device `0:0`. `Check` prints the backend class it actually runs on (`CUDABackendImpl` or `OCLBackendImpl`).

The OpenCL and CUDA rows are the same kernel on the same GPU. The CUDA backend is a little faster than OpenCL on this
card, but both are slower than the pocl CPU device for this kernel, so the GPU is not the right target for it yet.

CUDA on this Maxwell GPU (sm_50): CUDA 13 dropped sm_50, so CUDA 13.4 NVRTC refuses to compile for it. CUDA 12.8 is the
last toolkit that targets it. Its NVRTC works with the 580 driver, but its `nvcc` does not accept the host compiler here
(GCC 16 is newer than CUDA 12.8 supports), so TornadoVM must use the NVRTC path. Set:

    export CUDA_HOME=/usr/local/cuda-12.8
    export LD_LIBRARY_PATH=/usr/local/cuda-12.8/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
    TORNADOVM_HOME=$HOME/.sdkman/candidates/tornadovm/7.3.1-jdk22plus-cuda tornado/run-check.sh tornado

## JMH benchmark (CUDA and OpenCL)

`bench.sh` runs a JMH benchmark (`src/bench`) of the same task graph: transfers, kernel and the copy of the results to
the host, on the 3,025 protein pairs. JMH runs non-forked, in the launching JVM, because the TornadoVM launcher supplies
module and export flags that a forked JVM would not get; JMH prints a warning for that, and the numbers are still valid
for this comparison.

    export TORNADO_JDK=/usr/lib/jvm/java-25-temurin-jdk
    export TORNADOVM_HOME=$HOME/.sdkman/candidates/tornadovm/7.3.1-jdk22plus-cuda
    export CUDA_HOME=/usr/local/cuda-12.8
    export LD_LIBRARY_PATH=/usr/local/cuda-12.8/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
    DEVICE=0:0 tornado/bench.sh          # DEVICE is driver:index, from tornado --devices

The throughput is reported in ops/ms, where one invocation counts as 292,393,080 cells. Dividing by 10^6 gives GCUPS.

| Backend and device | `layout` | throughput (ops/ms) | GCUPS |
|---|---|---:|---:|
| CUDA, Quadro M1000M | INTERLEAVED | 64,095 ± 421 | 0.064 |
| CUDA, Quadro M1000M | CONTIGUOUS | 51,560 ± 22 | 0.052 |
| OpenCL, Quadro M1000M | INTERLEAVED | 55,797 ± 298 | 0.056 |
| OpenCL, Quadro M1000M | CONTIGUOUS | 47,530 ± 418 | 0.048 |
| OpenCL, pocl CPU | INTERLEAVED | 90,150 ± 2,302 | 0.090 |
| OpenCL, pocl CPU | CONTIGUOUS | 113,619 ± 1,923 | 0.114 |

Five measurement iterations of 10 s after two warm-up iterations of 5 s each. The CUDA rows agree with the separate
`run-check.sh` timings above (4.55 s and 5.58 s per run); the OpenCL rows agree with those timings too.
The OpenCL runs use `TORNADOVM_HOME` of the `7.3.1-jdk22plus-opencl` SDK and `DEVICE=0:0` (Quadro) or `DEVICE=0:1`
(pocl CPU); the CUDA environment variables are not needed for them.

## Why the GPU is slower: investigation

Each step has its own JMH benchmark in `src/bench` (run with `BENCH_FILTER=<class> tornado/bench.sh`) and its own
commit. All GPU numbers: CUDA backend, Quadro M1000M. Every benchmark checks its results against the plain Java row
kernel in setup.

| # | Benchmark | Question | Result |
|---|---|---|---|
| 1 | `TransferBench` | Transfers or launch overhead? | No: 0.07 ms of a 4,566 ms call. |
| 2 | `ImbalanceBench` | Does one pair set the time? | Yes: the largest pair alone takes 4,202 ms of 4,564. |
| 3 | `BalancedScalingBench` | Speed with balanced work? | 1.12 GCUPS from about 4,096 equal pairs, 18x the dataset. |
| 4 | `CpuBalancedBench` (main project) | Striped CPU on the same work? | 5.3 GCUPS on 1 thread; about 23.5 on 8 while cool. |
| 5 | `StripBench` | Is the kernel bound by H/E memory traffic? | Mostly not: 4x less traffic gives +10% (1.23 GCUPS). |
| 6 | `MatrixGatherBench` | Cost of the matrix lookups? | 56% of the time: 2.82 GCUPS without them. |
| 7 | `LocalMatrixBench` | Matrix in local memory? | 2.85 GCUPS (2.3x); dataset 0.119 GCUPS. |

**Is it the data dependencies?** Partly, in this sense: inside a pair every cell depends on its left, upper and diagonal
neighbours, so one work-item computes its pair serially, one cell after another. A single GPU thread is slow (about
1.2 us per cell for the row kernel, against a few ns per cell on a CPU core). The GPU is only fast when many pairs run at
the same time and none is much larger than the others. The test dataset breaks that: one pair has 35 times the mean
cells and takes 92% of the run (experiment 2).

**Does a different benchmark setup give better results?** Yes. With balanced pairs (many pairs of equal length) the
same kernel is 18 times faster than on the test dataset (experiment 3). The test dataset measures the latency of one
GPU thread, not the throughput of the GPU.

**What limits the balanced case?** Not the row memory traffic (experiment 5), but the substitution-matrix lookups:
neighbouring work-items read different, scattered matrix entries from global memory (experiment 6). Keeping the matrix
in local memory removes that cost (experiment 7, 2.85 GCUPS).

**And against the CPU?** Even the best GPU version (2.85 GCUPS) is slower than one CPU thread running the striped
Vector API aligner on the same pairs (5.3 GCUPS, experiment 4). The striped kernel scores up to 32 cells per instruction
with 8-bit lanes, the GPU kernel one cell per thread with 32-bit values, and the Quadro M1000M is a small 2015 laptop
GPU (5 multiprocessors). A larger GPU, 8- or 16-bit packed scores, and splitting long pairs over many threads
(anti-diagonal parallelism, as CUDASW++ does for long sequences) are the next steps; none of them is measured here.

**Measurement notes.** The 8-thread CPU number is not stable on this laptop: under full load the package reaches 85 C
within about 5 s and the clock is cut to 0.25-1.8 GHz from outside the CPU (its own throttle counter does not change).
One TornadoVM 7.3.1 code-generation problem was found and worked around in the strip kernel: a value loaded from a
`byte[]` is computed and stored with the byte type, so `row + letter` was truncated. `reproducers/` has a minimal
reproducer and the report for the TornadoVM team. A second problem suspected in the experiment 5 commit (a conditional
load compiled into an unconditional one) did not reproduce on its own and is not a known bug.

## Limits

* Scores equal the exact Smith-Waterman scores for gap open greater than gap extend, as in the C library.
* Scratch memory grows with `pairs * longest query`: 3,025 pairs of up to 1,863 residues use about 45 MB.
* Only scores and end positions; no begin positions or CIGAR.
