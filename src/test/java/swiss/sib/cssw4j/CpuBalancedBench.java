package swiss.sib.cssw4j;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
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
 * Experiment 4 of the GPU investigation (see tornado/README.md): the striped Vector API {@link Aligner} on the
 * balanced workload of the TornadoVM benchmarks, for a fair CPU-against-GPU comparison.
 * <p>
 * 128 random queries against 128 random targets, all 300 residues long, generated exactly as
 * {@code swiss.sib.cssw4j.tornado.Workloads.random} does (same letters, seeds 11 and 12), BLOSUM62, gap open 11, gap
 * extend 1. Like the GPU kernel it computes scores and end positions only ({@code flag} 0). The pairs are split over
 * {@code threads} worker threads, each with its own {@link Aligner}. Throughput in ops/ms with one op per cell, so
 * divided by 10^6 it is GCUPS.
 */
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3, time = 2)
@Measurement(iterations = 5, time = 3)
@Fork(value = 2, jvmArgs = { "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class CpuBalancedBench {
	static final int COUNT = 128, LENGTH = 300;
	static final int CELLS = COUNT * COUNT * LENGTH * LENGTH;

	@Param({ "1", "8" })
	public int threads;

	private List<byte[]> queries, targets;
	private ExecutorService pool;
	private ThreadLocal<Aligner> aligners;
	private final ScoringMatrix matrix = ScoringMatrix.blosum62();

	@Setup(Level.Trial)
	public void setup() {
		queries = random(COUNT, LENGTH, 11);
		targets = random(COUNT, LENGTH, 12);
		pool = Executors.newFixedThreadPool(threads);
		aligners = ThreadLocal.withInitial(Aligner::new);
	}

	@TearDown(Level.Trial)
	public void tearDown() {
		pool.shutdown();
	}

	/** Same generator as the TornadoVM module, so both benchmarks align identical sequences. */
	private List<byte[]> random(int count, int length, long seed) {
		String letters = "ARNDCQEGHILKMFPSTWYV";
		Random rnd = new Random(seed);
		List<byte[]> out = new ArrayList<>();
		for (int i = 0; i < count; i++) {
			StringBuilder sb = new StringBuilder(length);
			for (int k = 0; k < length; k++) {
				sb.append(letters.charAt(rnd.nextInt(letters.length())));
			}
			out.add(matrix.encode(sb));
		}
		return out;
	}

	@Benchmark
	@OperationsPerInvocation(CELLS)
	public void balanced(Blackhole bh) throws Exception {
		List<Future<Long>> parts = new ArrayList<>();
		for (byte[] q : queries) {
			parts.add(pool.submit(() -> {
				Profile profile = Profile.of(q, matrix);
				Aligner aligner = aligners.get();
				long checksum = 0;
				for (byte[] t : targets) {
					Alignment a = aligner.align(profile, t, t.length, 11, 1, 0, 0, 0, q.length / 2);
					checksum += a.score1() + a.refEnd1() + a.readEnd1();
				}
				return checksum;
			}));
		}
		long checksum = 0;
		for (Future<Long> part : parts) {
			checksum += part.get();
		}
		bh.consume(checksum);
	}
}
