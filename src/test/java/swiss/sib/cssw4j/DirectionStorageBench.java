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
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Strategy: store the traceback directions of a cell in one byte instead of three.
 * <p>
 * Why: the first version of the banded alignment kept three direction bytes per cell: one for the E (horizontal gap)
 * choice, one for the F (vertical gap) choice and one for H. Three stores per cell means three times the memory traffic
 * for the direction array, which is the largest array the banded step touches. Packing the three choices into one byte
 * (the layout in {@link Aligner}) cuts the stores to one per cell. Variants:
 * <ol>
 * <li>{@code threeBytes}: three byte stores per cell at {@code 3j}, {@code 3j + 1}, {@code 3j + 2}.</li>
 * <li>{@code onePackedByte}: one store per cell.</li>
 * </ol>
 * Result: in isolation the packed byte is 7 to 9% slower than the three stores (about 2.2 against 2.0 us for 1024
 * cells, 139 against 129 us for 65536). The port saw no measurable end-to-end change (about 1.5 GCUPS on the protein
 * set, full CIGAR). The packed layout is kept for its footprint: one third of the bytes, which matters once the
 * direction array no longer fits in cache. The sizes here do not reach that point.
 * Both versions do the same loads and arithmetic, and each cell depends on the previous one, so only the number of
 * stores differs. (Earlier drafts had no dependency and measured the vectorizer, not the stores.)
 * <p>
 * The setup checks that a later read of the packed byte returns the same values as the three stored bytes.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(2)
@State(Scope.Thread)
public class DirectionStorageBench {
	/** Cells per call: one band row of a long query times a few rows. */
	@Param({ "1024", "65536" })
	int cells;

	byte[] threeBytes, onePacked;
	byte[] choices;

	@Setup(Level.Trial)
	public void setup() {
		Random rnd = new Random(5);
		choices = new byte[cells * 3];
		for (int i = 0; i < choices.length; i++) {
			choices[i] = (byte) rnd.nextInt(8);
		}
		threeBytes = new byte[cells * 3];
		onePacked = new byte[cells];
		writeThree(choices, threeBytes, cells);
		writePacked(choices, onePacked, cells);
		for (int j = 0; j < cells; j++) {
			int packed = onePacked[j];
			int three = threeBytes[3 * j] | threeBytes[3 * j + 1] << 3 | threeBytes[3 * j + 2] << 4;
			if (packed != three) {
				throw new IllegalStateException("cell " + j + ": " + packed + " != " + three);
			}
		}
	}

	/**
	 * Three separate stores, as in the first version of {@code bandRow}. Each cell's choices depend on the previous
	 * cell (a carry), as in the real recurrence, so the loop cannot be turned into a plain copy by the vectorizer.
	 */
	static void writeThree(byte[] choices, byte[] out, int cells) {
		int carry = 0;
		for (int j = 0; j < cells; j++) {
			int h = (choices[3 * j] + carry) & 7;
			int e = choices[3 * j + 1] ^ carry;
			int f = (choices[3 * j + 2] + h) & 7;
			carry = (h + e + f) & 7;
			out[3 * j] = (byte) h;
			out[3 * j + 1] = (byte) e;
			out[3 * j + 2] = (byte) f;
		}
	}

	/** One store per cell, the packed layout used now. Same arithmetic and carry as {@link #writeThree}. */
	static void writePacked(byte[] choices, byte[] out, int cells) {
		int carry = 0;
		for (int j = 0; j < cells; j++) {
			int h = (choices[3 * j] + carry) & 7;
			int e = choices[3 * j + 1] ^ carry;
			int f = (choices[3 * j + 2] + h) & 7;
			carry = (h + e + f) & 7;
			out[j] = (byte) (h | e << 3 | f << 4);
		}
	}

	@Benchmark
	public void threeBytes(Blackhole bh) {
		writeThree(choices, threeBytes, cells);
		bh.consume(threeBytes);
	}

	@Benchmark
	public void onePackedByte(Blackhole bh) {
		writePacked(choices, onePacked, cells);
		bh.consume(onePacked);
	}
}
