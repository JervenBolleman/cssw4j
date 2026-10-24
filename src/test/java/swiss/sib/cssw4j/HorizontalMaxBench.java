package swiss.sib.cssw4j;

import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

/**
 * Strategy: take the maximum of all lanes of a byte vector, once per column, in three ways.
 * <p>
 * Why: the striped kernel needs the column maximum after every column (for the second best score and the best end
 * position). The horizontal reduction is the only step that works across lanes, so its cost is paid once per column.
 * Variants tried:
 * <ol>
 * <li>{@code reduceLanes(MAX)}: the library operation, which the port uses.</li>
 * <li>a shuffle tree: fold the vector onto itself with a permute at distance 16, 8, 4, 2 and 1 and take a lane-wise
 * max each time.</li>
 * <li>a scalar scan: store the vector to an array and loop over it.</li>
 * </ol>
 * Result in the port: the first variant measured fastest in a separate microbenchmark (about 2.9 ns per vector,
 * against 5.6 ns for the shuffle tree), which is why the port uses it. The benchmark below repeats that comparison
 * with the scalar scan added.
 * <p>
 * The setup checks that all three return the same maximum for every vector.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(value = 2, jvmArgs = { "--add-modules", "jdk.incubator.vector" })
@State(Scope.Thread)
public class HorizontalMaxBench {
	static final VectorSpecies<Byte> BYTE = VectorSpecies.of(byte.class, Lanes.SHAPE);
	static final int VECTORS = 256;

	static final VectorShuffle<Byte>[] FOLD = foldShuffles();

	byte[][] data;
	byte[] scratch;

	@SuppressWarnings("unchecked")
	static VectorShuffle<Byte>[] foldShuffles() {
		int lanes = BYTE.length();
		int steps = Integer.numberOfTrailingZeros(lanes);
		VectorShuffle<Byte>[] s = new VectorShuffle[steps];
		for (int k = 0; k < steps; k++) {
			final int distance = lanes >> (k + 1);
			s[k] = VectorShuffle.fromOp(BYTE, i -> i ^ distance);
		}
		return s;
	}

	@Setup(Level.Trial)
	public void setup() {
		Random rnd = new Random(11);
		data = new byte[VECTORS][BYTE.length()];
		for (byte[] v : data) {
			rnd.nextBytes(v);
		}
		scratch = new byte[BYTE.length()];
		for (byte[] v : data) {
			ByteVector x = ByteVector.fromArray(BYTE, v, 0);
			int a = reduce(x), b = tree(x), c = scalar(x);
			if (a != b || a != c) {
				throw new IllegalStateException("maxima differ: " + a + " " + b + " " + c);
			}
		}
	}

	static int reduce(ByteVector x) {
		return x.reduceLanes(VectorOperators.MAX);
	}

	static int tree(ByteVector x) {
		for (VectorShuffle<Byte> s : FOLD) {
			x = x.lanewise(VectorOperators.MAX, x.rearrange(s));
		}
		return x.lane(0);
	}

	int scalar(ByteVector x) {
		x.intoArray(scratch, 0);
		int m = scratch[0];
		for (int i = 1; i < scratch.length; i++) {
			m = Math.max(m, scratch[i]);
		}
		return m;
	}

	@Benchmark
	public void reduceLanes(Blackhole bh) {
		int m = Integer.MIN_VALUE;
		for (byte[] v : data) {
			m = Math.max(m, reduce(ByteVector.fromArray(BYTE, v, 0)));
		}
		bh.consume(m);
	}

	@Benchmark
	public void shuffleTree(Blackhole bh) {
		int m = Integer.MIN_VALUE;
		for (byte[] v : data) {
			m = Math.max(m, tree(ByteVector.fromArray(BYTE, v, 0)));
		}
		bh.consume(m);
	}

	@Benchmark
	public void scalarScan(Blackhole bh) {
		int m = Integer.MIN_VALUE;
		for (byte[] v : data) {
			m = Math.max(m, scalar(ByteVector.fromArray(BYTE, v, 0)));
		}
		bh.consume(m);
	}
}
