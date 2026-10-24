package swiss.sib.cssw4j;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Strategy: build the striped query profile by first striping the query letters and then looking up one row of the
 * matrix per reference letter, instead of computing each entry from the query directly.
 * <p>
 * Why: the profile is built once per query, and for a query that is aligned against many references its cost is
 * small. But for short queries it is not negligible, and the first version had a conditional inside the innermost
 * loop ({@code j >= readLen ? pad : mat[...]}) plus an index computation per entry. Variants:
 * <ol>
 * <li>{@code direct}: for each reference letter, for each vector and lane, compute the striped index, test for the
 * padding and look up the score.</li>
 * <li>{@code stripeThenLookup}: stripe the query letters once into an array (padding included), then for each
 * reference letter map the striped letters through one matrix row. This is what {@code Profile.buildByte} does.</li>
 * </ol>
 * Result: stripe-then-lookup is 26 to 29% faster (2.9 against 3.9 us for 100 residues, 23 against 33 us for 1000),
 * even though it allocates a temporary array on every call. The profile was about 5% of the full-CIGAR run in the
 * port's profiler, so the end-to-end effect is small.
 * <p>
 * The setup checks that both builders produce the same bytes.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(2)
@State(Scope.Thread)
public class ProfileBuildBench {
	static final int LANES = 32;
	static final int ALPHABET = 24;

	@Param({ "100", "1000" })
	int queryLength;

	byte[] query, matrix, direct, stripeThenLookup;

	@Setup(Level.Trial)
	public void setup() {
		Random rnd = new Random(9);
		query = new byte[queryLength];
		for (int i = 0; i < query.length; i++) {
			query[i] = (byte) rnd.nextInt(ALPHABET);
		}
		matrix = ScoringMatrix.BLOSUM50.matrix().clone();
		int stride = (queryLength + LANES - 1) / LANES * LANES;
		direct = new byte[ALPHABET * stride];
		stripeThenLookup = new byte[ALPHABET * stride];
		buildDirect(query, matrix, direct);
		buildStripeThenLookup(query, matrix, stripeThenLookup);
		if (!Arrays.equals(direct, stripeThenLookup)) {
			throw new IllegalStateException("profiles differ");
		}
	}

	static void buildDirect(byte[] read, byte[] mat, byte[] out) {
		final int readLen = read.length;
		final int segLen = (readLen + LANES - 1) / LANES;
		int k = 0;
		for (int nt = 0; nt < ALPHABET; nt++) {
			for (int i = 0; i < segLen; i++) {
				for (int seg = 0, j = i; seg < LANES; seg++, j += segLen) {
					out[k++] = (byte) (j >= readLen ? 0 : mat[nt * ALPHABET + read[j]]);
				}
			}
		}
	}

	static void buildStripeThenLookup(byte[] read, byte[] mat, byte[] out) {
		final int readLen = read.length;
		final int segLen = (readLen + LANES - 1) / LANES;
		final int stride = segLen * LANES;
		final byte[] striped = new byte[stride];
		Arrays.fill(striped, (byte) ALPHABET); // index of the padding entry
		for (int j = 0, k = 0, seg = 0; j < readLen; j++) {
			striped[seg * LANES + k] = read[j];
			if (++seg == segLen) {
				seg = 0;
				k++;
			}
		}
		final byte[] row = new byte[ALPHABET + 1];
		int k = 0;
		for (int nt = 0; nt < ALPHABET; nt++) {
			System.arraycopy(mat, nt * ALPHABET, row, 0, ALPHABET);
			row[ALPHABET] = 0;
			for (int p = 0; p < stride; p++) {
				out[k++] = row[striped[p]];
			}
		}
	}

	@Benchmark
	public void direct(Blackhole bh) {
		buildDirect(query, matrix, direct);
		bh.consume(direct);
	}

	@Benchmark
	public void stripeThenLookup(Blackhole bh) {
		buildStripeThenLookup(query, matrix, stripeThenLookup);
		bh.consume(stripeThenLookup);
	}
}
