# TornadoVM reproducer: a byte[] value stored into an int[] is written as one byte

## Summary

When a kernel stores into an `int[]` element a value that comes from a `byte[]` load (directly or after int
arithmetic), the generated kernel writes only one byte into the int element. When all inputs of the arithmetic come
from byte loads, the arithmetic itself is also done in the 8-bit type. As a result, the
upper three bytes of the int are never written: negative values lose their sign, values above 127 are cut modulo 256,
and where the device buffer was not zero (pocl) the upper bytes contain garbage.

The same code is generated for the CUDA and OpenCL backends, so the problem looks like it is in the shared compiler
part (store width or the type of the value derived from the byte load), not in one backend.

## Environment

* TornadoVM 7.3.1, SDKs `7.3.1-jdk22plus-cuda` and `7.3.1-jdk22plus-opencl` (from sdkman)
* JDK: Temurin 25.0.5+8 (Linux x86_64, Fedora 44, kernel 7.2.11)
* Devices:
  * NVIDIA Quadro M1000M (Maxwell, sm_50), driver 580.178.04; CUDA backend with CUDA 12.8 NVRTC, and OpenCL
  * pocl 7.1 (OpenCL 3.0, LLVM 21.1.8) CPU device, Intel i7-6700HQ

## Reproduce

    export JAVA_HOME=/path/to/jdk-25
    export TORNADOVM_HOME=/path/to/tornadovm-7.3.1-jdk22plus-opencl   # or the -cuda SDK (then also CUDA_HOME etc.)
    ./run.sh 0:0     # device as driver:index from `tornado --devices`

`ByteToInt.java` runs four one-line kernels over all 256 byte values (-128 .. 127) and compares each with the same
method run as plain Java.

| Kernel | Java |
|---|---|
| `copy` | `out[i] = a[i];` |
| `addConstant` | `out[i] = a[i] + 1000;` |
| `multiplyConstant` | `out[i] = a[i] * 24;` |
| `withIndex` | `out[i] = (i % 24) * 24 + a[i];` |

## Expected

Four `PASS` lines (for example `copy` stores -128 for `a[0] = -128`).

## Actual

CUDA backend, Quadro M1000M (OpenCL on the same GPU prints the same):

    FAIL copy             out[i] = a[i]: 128 of 256 wrong, e.g. a[0] = -128: expected -128, got 128
    FAIL addConstant      out[i] = a[i] + 1000: 256 of 256 wrong, e.g. a[0] = -128: expected 872, got 104
    FAIL multiplyConstant out[i] = a[i] * 24: 245 of 256 wrong, e.g. a[0] = -128: expected -3072, got 0
    FAIL withIndex        out[i] = (i % 24) * 24 + a[i]: 156 of 256 wrong, e.g. a[0] = -128: expected -128, got 128

OpenCL backend, pocl CPU device (the upper bytes keep whatever was in the device buffer):

    FAIL copy             out[i] = a[i]: 174 of 256 wrong, e.g. a[0] = -128: expected -128, got 128
    FAIL addConstant      out[i] = a[i] + 1000: 256 of 256 wrong, e.g. a[0] = -128: expected 872, got 560246888
    FAIL multiplyConstant out[i] = a[i] * 24: 256 of 256 wrong, e.g. a[0] = -128: expected -3072, got 1076438016
    FAIL withIndex        out[i] = (i % 24) * 24 + a[i]: 158 of 256 wrong, e.g. a[0] = -128: expected -128, got 128

## Generated code (`tornado --printKernel`)

OpenCL, `copy`: the int[] element is written with a `char` store.

```c
ch_8  =  *((__global char *) ul_7);      // a[i]
l_9  =  l_5 << 2;                        // i * 4: the out element is an int
l_10  =  l_9 + 16L;
ul_11  =  ul_1 + l_10;
*((__global char *) ul_11)  =  ch_8;     // one byte written into out[i]
```

OpenCL, `withIndex`: the arithmetic is done in int, but the store is still one byte.

```c
ch_8  =  *((__global char *) ul_7);      // a[i]
i_12  =  i_4 % 24;
i_13  =  i_12 << 4;
i_14  =  i_12 << 3;
i_15  =  i_13 + i_14;
i_16  =  i_15 + ch_8;                    // int, correct
*((__global char *) ul_11)  =  i_16;     // one byte written into out[i]
```

CUDA, `multiplyConstant` (`out[i] = a[i] * 24`): the arithmetic is done in `signed char` and stored as one byte.

```c
signed char ch_8, ch_12, ch_13, ch_14;
ch_8  =  *(( signed char *) ul_7);       // a[i]
ch_12  =  ch_8 << 4;
ch_13  =  ch_8 << 3;
ch_14  =  ch_12 + ch_13;                 // a[i] * 24 in 8 bits
*(( signed char *) ul_11)  =  ch_14;     // one byte written into out[i]
```

## Impact

Any kernel that stores into an int array a value derived from byte data gives wrong results. An index computed only
from byte data is affected too: for `out[i] = table[a[i] * 24 + b[i]]` (a substitution-matrix lookup with `a`, `b`
byte arrays) the generated code computes the index in 8 bits and widens it with `& 4294967295UL`, so a negative 8-bit
value becomes an index of about 4 billion. On pocl that crashed the JVM with SIGSEGV; on the Quadro the kernel returned
zeros. Whether a given expression is affected depends on the surrounding code: in our Smith-Waterman kernel the same
`matrix[target * n + query]` lookup was correct in one kernel and truncated in another.

## Workaround

Load the byte data into an `int[]` on the host and pass that to the kernel, so no value in the kernel has the byte
type. Masking with `& 0xFF` does not help: the compiler removes the mask.

## Not a bug (checked)

We had also seen a `CUDA_ERROR_ILLEGAL_ADDRESS` in a kernel that used conditional loads of byte data
(`j + 1 < n ? q[i + 1] : 0`), and suspected that the load in the untaken arm is executed unconditionally. A standalone
reproducer with indices far out of bounds in the untaken arm (`k < n ? data[k] : -1`, int data) passed on all three
backends, so we do not report that. After replacing those loads with clamped indices the crash was gone but the results
were still wrong because of the byte problem above; what exactly caused the crash was not determined.
