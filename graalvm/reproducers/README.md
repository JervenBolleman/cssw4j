# GraalVM reproducer: Vector API code falls back to the Java implementation when a vector-typed callee is not inlined

For the [GraalVM issue tracker](https://github.com/oracle/graal/issues) (sections as in its bug template).

## Describe the issue

A Vector API loop compiled by Graal runs every vector operation through its Java fallback when the loop passes a
vector to a method that is not inlined. This affects `ByteVector.fromArray`, `lanewise`, `max` and the other
operations of the whole loop, not only the call itself. The result is about 140x slower than with the method inlined,
and about 30x slower than C2 with the same method kept out of line.

C2 handles the same situation by boxing the vector at the call. The surrounding operations still use their
intrinsics, so the loop is only about 5x slower.

We hit this without any flags, in a protein database search built on
[cssw4j](../..), a Java Vector API port of the SSW Smith-Waterman library:

* One search over UniProtKB/Swiss-Prot took over 500 s on Oracle GraalVM 25.4.4, instead of 7.5 s on Temurin 25.
* Graal's inliner left a 73-byte helper, `Lanes.shiftBytes(ByteVector)`, out of line: `-Djdk.graal.Log=Inlining`
  reports `not-within-budget`. This happened only once the calling loop had grown (an E-value computation and a second
  alignment call for hits).
* Stack samples were all in `ByteVector.bOpTemplate` and other `*Template` fallbacks.
* `-Djdk.graal.TypicalGraphSize=20000` or `-Djdk.graal.TuneInlinerExploration=1` restore the normal speed.

Whether the inliner hits this depends on the size of the surrounding code, so the reproducer forces the same inlining
decision with `-XX:CompileCommand=dontinline` instead.

## Steps to reproduce the issue

`VectorInlining.java` is a single file with no dependencies:

* `kernel` is the 8-bit striped Smith-Waterman column loop of cssw4j, reduced to computing the score.
* `scan` calls `kernel` for 20,000 random targets.
* `shiftBytes(ByteVector)` shifts the lanes of a vector by one byte. `kernel` calls it about twice per target column.

```
./run.sh                    # uses java from PATH, or: JAVA=/path/to/bin/java ./run.sh
```

which runs

```
java --add-modules jdk.incubator.vector VectorInlining.java
java -XX:CompileCommand=quiet -XX:CompileCommand=dontinline,VectorInlining::shiftBytes \
     --add-modules jdk.incubator.vector VectorInlining.java
```

## Describe GraalVM and your environment

* GraalVM: Oracle GraalVM 25.4.4 (`Java HotSpot(TM) 64-Bit Server VM 25.0.4.1.1+1-LTS-jvmci-25.4-b23`) and GraalVM
  Community 25.4.4.1 (`OpenJDK 64-Bit Server VM 25.0.4.1.1+1-jvmci-25.4-b23`), both installed with sdkman
* JDK version: 25.0.4
* OS: Fedora 44, Linux 7.2.11 x86_64
* Architecture: Intel Core i7-6700HQ (AVX2, `Species[byte, 32, S_256_BIT]`)

## More details

Time of one `scan` round once warmed up (round 7 of 8):

| JVM | `shiftBytes` inlined (default) | `dontinline shiftBytes` |
|---|---:|---:|
| Oracle GraalVM 25.4.4 | 187 ms | 26,655 ms |
| GraalVM CE 25.4.4.1 | 160 ms | 21,524 ms |
| Temurin 25.0.4 (C2) | 155 ms | 816 ms |

`-XX:+PrintCompilation` (Oracle GraalVM, `dontinline`) shows that the slowdown starts when `scan` is compiled with
`kernel` inlined into it. While `kernel` runs as its own tier 4 compilation, round 1 takes 658 ms. That is close to
C2's 816 ms, so `kernel` compiled on its own does use the intrinsics.

```
1541 2864       4       VectorInlining::kernel (656 bytes)
round  1:     658 ms (checksum 487669)
3071 2900 %     4       VectorInlining::scan @ 14 (55 bytes)
3350 2901       4       VectorInlining::scan (55 bytes)
round  2:   13279 ms (checksum 487669)
round  3:   25635 ms (checksum 487669)
```

A stack sample in the slow phase (`jstack`):

```
"main" #3 ... runnable
   java.lang.Thread.State: RUNNABLE
	at jdk.incubator.vector.ByteVector.fromArray0Template(jdk.incubator.vector@25.0.4.1.1/ByteVector.java:3832)
	at jdk.incubator.vector.Byte256Vector.fromArray0(jdk.incubator.vector@25.0.4.1.1/Byte256Vector.java:1020)
	at jdk.incubator.vector.ByteVector.fromArray(jdk.incubator.vector@25.0.4.1.1/ByteVector.java:3032)
	at VectorInlining.kernel(VectorInlining.java:102)
	at VectorInlining.scan(VectorInlining.java:133)
	at VectorInlining.main(VectorInlining.java:173)
```

In the database search, the samples showed `ByteVector.bOpTemplate` under `ByteVector.lambda$binaryOperations$34`
from the same kind of loop.

## Expected behavior

When a method that takes or returns a vector is not inlined, the vector should be boxed at that call (as C2 does).
The rest of the loop should keep its intrinsics, so the cost is a few times slower, not 140 times. Alternatively, the
inliner could treat small vector-typed helpers as must-inline, so that a growing caller cannot push them out of the
budget.

## Workarounds

* Use C2 (`-XX:-UseJVMCICompiler`, or a JDK without Graal).
* Raise the inlining budget, for example `-Djdk.graal.TypicalGraphSize=20000`.
* Keep vector-typed helpers out of hot loops, writing the code inline instead.
