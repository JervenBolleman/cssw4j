package swiss.sib.cssw4j;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

/**
 * Checks the vector kernels against a plain scalar Gotoh Smith-Waterman, independently of the C library, on random
 * and mutated protein sequences. Also checks that every returned cigar re-scores to the reported score.
 * <p>
 * Only gap penalties with {@code gapOpen > gapExtend} are used: for {@code gapOpen <= gapExtend} the lazy-F early exit
 * of SSW (inherited from Farrar's algorithm, and present in the C library) can stop propagating vertical gaps too
 * early and under-report the score; how often depends on the vector width.
 */
class ScalarOracleTest {
	private static final String AA = "ARNDCQEGHILKMFPSTWYV";

	@Test
	void randomProteins() {
		Random rnd = new Random(42);
		Aligner aligner = new Aligner();
		int[][] gaps = { { 11, 1 }, { 10, 2 }, { 3, 1 }, { 10, 3 }, { 4, 3 } };
		int noPath = 0;
		for (int iter = 0; iter < 2000; iter++) {
			ScoringMatrix matrix = iter % 2 == 0 ? ScoringMatrix.BLOSUM50 : ScoringMatrix.blosum62();
			int[] gap = gaps[iter % gaps.length];
			String read = random(rnd, 1 + rnd.nextInt(iter % 10 == 0 ? 1200 : 300));
			// Half of the references contain a mutated copy of part of the read, to get high scoring alignments.
			String ref = rnd.nextBoolean() ? random(rnd, 1 + rnd.nextInt(400))
					: random(rnd, rnd.nextInt(50)) + mutate(rnd, read.substring(rnd.nextInt(read.length())))
							+ random(rnd, rnd.nextInt(50));
			byte[] r = matrix.encode(read);
			byte[] t = matrix.encode(ref);
			Alignment a = aligner.align(Profile.of(r, matrix), t, gap[0], gap[1]);
			int expected = gotoh(t, r, matrix, gap[0], gap[1]);
			String msg = "iteration " + iter + " read=" + read + " ref=" + ref;
			assertEquals(expected, a.score1(), msg);
			if (a.flag() == 1) {
				// The banded traceback of the C library gives up on some alignments; the port does the same.
				noPath++;
			} else if (expected > 0) {
				assertNotNull(a.cigar(), msg);
				assertEquals(expected, rescore(a, t, r, matrix, gap[0], gap[1]), msg);
			}
		}
		System.out.println("random proteins: " + noPath + " of 2000 alignments without path (flag 1)");
		// The C library fails on ~10% of protein pairs with cheap gaps (see the BLOSUM50 3/1 reference data).
		assertTrue(noPath < 200, noPath + " alignments without path");
	}

	/** Profiles restricted to one kernel must agree with the default one. */
	@Test
	void scoreSizes() {
		Random rnd = new Random(7);
		Aligner aligner = new Aligner();
		ScoringMatrix m = ScoringMatrix.BLOSUM50;
		for (int iter = 0; iter < 300; iter++) {
			String read = random(rnd, 1 + rnd.nextInt(60));
			String ref = random(rnd, 20) + mutate(rnd, read) + random(rnd, 20);
			byte[] r = m.encode(read);
			byte[] t = m.encode(ref);
			Alignment both = aligner.align(Profile.of(r, r.length, m, Profile.ScoreSize.BOTH), t, 10, 1);
			Alignment word = aligner.align(Profile.of(r, r.length, m, Profile.ScoreSize.WORD), t, 10, 1);
			assertEquals(both.score1(), word.score1());
			assertEquals(both.refEnd1(), word.refEnd1());
			assertEquals(both.readEnd1(), word.readEnd1());
			assertEquals(both.refBegin1(), word.refBegin1());
			assertEquals(both.readBegin1(), word.readBegin1());
			if (both.score1() + 5 < 255) { // 5 = BLOSUM50 bias, the 8 bit kernel saturates at 255 - bias
				Alignment b = aligner.align(Profile.of(r, r.length, m, Profile.ScoreSize.BYTE), t, 10, 1);
				assertEquals(both.toString(), b.toString());
			}
		}
	}

	static int gotoh(byte[] ref, byte[] read, ScoringMatrix m, int gapO, int gapE) {
		int n = read.length;
		int[] h = new int[n + 1], e = new int[n + 1];
		int best = 0;
		for (byte rc : ref) {
			int diag = 0, f = 0;
			for (int j = 1; j <= n; j++) {
				e[j] = Math.max(e[j] - gapE, h[j] - gapO);
				f = Math.max(f - gapE, h[j - 1] - gapO);
				int v = Math.max(0, Math.max(diag + m.score(rc, read[j - 1]), Math.max(e[j], f)));
				diag = h[j];
				h[j] = v;
				best = Math.max(best, v);
			}
		}
		return best;
	}

	static int rescore(Alignment a, byte[] ref, byte[] read, ScoringMatrix m, int gapO, int gapE) {
		int s = 0, rp = a.refBegin1(), qp = a.readBegin1();
		for (int c : a.cigar()) {
			int len = Cigar.length(c);
			switch (Cigar.op(c)) {
				case Cigar.M -> {
					for (int k = 0; k < len; k++) {
						s += m.score(ref[rp++], read[qp++]);
					}
				}
				case Cigar.I -> {
					s -= gapO + (len - 1) * gapE;
					qp += len;
				}
				case Cigar.D -> {
					s -= gapO + (len - 1) * gapE;
					rp += len;
				}
				default -> throw new AssertionError(a);
			}
		}
		assertEquals(a.refEnd1() + 1, rp, "cigar spans the reference");
		assertEquals(a.readEnd1() + 1, qp, "cigar spans the read");
		return s;
	}

	static String random(Random rnd, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++) {
			sb.append(AA.charAt(rnd.nextInt(AA.length())));
		}
		return sb.toString();
	}

	static String mutate(Random rnd, String s) {
		StringBuilder sb = new StringBuilder(s.length());
		for (int i = 0; i < s.length(); i++) {
			double p = rnd.nextDouble();
			if (p < 0.02) {
				i += rnd.nextInt(4); // deletion
			} else if (p < 0.04) {
				sb.append(random(rnd, 1 + rnd.nextInt(4))).append(s.charAt(i)); // insertion
			} else if (p < 0.2) {
				sb.append(AA.charAt(rnd.nextInt(AA.length())));
			} else {
				sb.append(s.charAt(i));
			}
		}
		return sb.toString();
	}
}
