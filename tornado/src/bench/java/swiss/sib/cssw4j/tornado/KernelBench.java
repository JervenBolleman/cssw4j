package swiss.sib.cssw4j.tornado;

import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OperationsPerInvocation;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * JMH benchmark of the scoring kernel on a TornadoVM device, for the 55 x 55 protein pairs of the test dataset
 * (BLOSUM62, gap open 11, gap extend 1).
 * <p>
 * Each invocation runs the whole task graph: transfers, kernel and the copy of the results back to the host. One
 * invocation is {@link #CELLS} dynamic-programming cells, so the reported throughput in ops/ms divided by 10^6 is
 * GCUPS (giga cell updates per second).
 * <p>
 * The device is the one named by {@code -Dcssw4j.device} (default {@code 0:0}); the TornadoVM SDK decides the backend.
 * The benchmark runs in the launching JVM (no fork), because the TornadoVM launcher supplies the module and export
 * flags that a forked JVM would not get; see {@link KernelBenchRunner}.
 */
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@OperationsPerInvocation(KernelBench.CELLS)
@Warmup(iterations = 2, time = 5)
@Measurement(iterations = 5, time = 10)
@Fork(0)
@State(Scope.Benchmark)
public class KernelBench {
	/** Dynamic-programming cells per invocation: the sum of query length times target length over all pairs. */
	public static final int CELLS = 292_393_080;

	/** The test dataset; {@code -Dcssw4j.data} overrides the default, which is relative to the tornado directory. */
	private static final Path DATA = Path.of(System.getProperty("cssw4j.data", "../src/test/resources/dataset"));
	private static final Path MATRIX = DATA.resolve("blosum62.txt");

	@Param({ "INTERLEAVED", "CONTIGUOUS" })
	public String layout;

	private TornadoScorer scorer;

	@Setup(Level.Trial)
	public void setup() throws Exception {
		Matrix matrix = Matrix.read(MATRIX);
		List<Path> queries = List.of(DATA.resolve("uniprot_human.fasta"), DATA.resolve("synthetic_queries.fasta"));
		List<Path> targets = List.of(DATA.resolve("uniprot_homologs.fasta"), DATA.resolve("synthetic_targets.fasta"));
		Sequences q = Sequences.read(matrix, queries);
		Sequences t = Sequences.read(matrix, targets);
		scorer = new TornadoScorer(matrix, q, t, 11, 1, Layout.valueOf(layout));
		scorer.run(); // compiles the kernel and warms up the device
		if (scorer.scores[0] < 0) {
			throw new IllegalStateException("unexpected result");
		}
	}

	@TearDown(Level.Trial)
	public void tearDown() {
		scorer.close();
	}

	@Benchmark
	public void score(Blackhole bh) {
		scorer.run();
		bh.consume(scorer.scores);
		bh.consume(scorer.refEnds);
		bh.consume(scorer.readEnds);
	}
}
