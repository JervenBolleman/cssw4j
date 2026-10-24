package swiss.sib.cssw4j;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * Compares every field of every alignment with the output of the original C library (see
 * {@code src/test/scripts/generate_expected.sh}).
 * <p>
 * With 16 byte lanes (as SSE2) all fields must be identical. With wider vectors the sub-optimal score can legitimately
 * differ: the C code includes the padding cells of the last segment in its per column maxima, and the amount of padding
 * depends on the lane count. The optimal alignment (score, begin, end, cigar) must still be identical.
 */
class CReferenceTest {
	static final String[] PROTEIN_QUERIES = { "uniprot_human.fasta", "synthetic_queries.fasta" };
	static final String[] PROTEIN_TARGETS = { "uniprot_homologs.fasta", "synthetic_targets.fasta" };

	record Case(String file, Supplier<ScoringMatrix> matrix, int gapO, int gapE, String[] queries, String[] targets) {
		@Override
		public String toString() {
			return file;
		}
	}

	static Stream<Arguments> cases() {
		return Stream.of(
				new Case("protein_blosum50_o3_e1", () -> ScoringMatrix.BLOSUM50, 3, 1, PROTEIN_QUERIES, PROTEIN_TARGETS),
				new Case("protein_blosum50_o10_e2", () -> ScoringMatrix.BLOSUM50, 10, 2, PROTEIN_QUERIES,
						PROTEIN_TARGETS),
				new Case("protein_blosum62_o11_e1", ScoringMatrix::blosum62, 11, 1, PROTEIN_QUERIES, PROTEIN_TARGETS),
				new Case("dna_m2_x2_o3_e1", () -> ScoringMatrix.nucleotide(2, 2), 3, 1,
						new String[] { "dna_54mer_hap1_1.100.fa" }, new String[] { "dna_1k.fa" }))
				.map(Arguments::of);
	}

	@ParameterizedTest(name = "{0}")
	@MethodSource("cases")
	void matchesCLibrary(Case c) throws IOException {
		ScoringMatrix matrix = c.matrix().get();
		List<Fasta> queries = Fasta.resource(c.queries());
		List<Fasta> targets = Fasta.resource(c.targets());
		Map<String, byte[]> encodedTargets = new HashMap<>();
		for (Fasta t : targets) {
			encodedTargets.put(t.id(), matrix.encode(t.seq()));
		}
		Map<String, Fasta> queryById = new HashMap<>();
		for (Fasta q : queries) {
			queryById.put(q.id(), q);
		}

		boolean sameLanesAsC = Lanes.BYTE.length() == 16;
		Aligner aligner = new Aligner();
		int pairs = 0, score2Differs = 0, wordPath = 0;
		List<String> failures = new ArrayList<>();
		String currentQuery = null;
		Profile profile = null;
		byte[] read = null;
		for (String line : expected(c.file())) {
			String[] f = line.split("\t");
			if (!f[0].equals(currentQuery)) {
				currentQuery = f[0];
				read = matrix.encode(queryById.get(currentQuery).seq());
				profile = Profile.of(read, matrix);
			}
			byte[] ref = encodedTargets.get(f[1]);
			Alignment a = aligner.align(profile, ref, ref.length, c.gapO(), c.gapE(), Aligner.FLAG_ALWAYS, 0, 0,
					read.length / 2);
			pairs++;
			if (a.score1() >= 255) {
				wordPath++;
			}
			String marked = a.cigar() == null ? "-1\t*" : markedString(Cigar.markMismatch(a, ref, read, read.length));
			String optimal = String.join("\t", f[0], f[1], f[2], f[4], f[5], f[6], f[7], f[9], f[10], f[11], f[12]);
			String actualOptimal = String.join("\t", f[0], f[1], str(a.score1()), str(a.refBegin1()),
					str(a.refEnd1()), str(a.readBegin1()), str(a.readEnd1()), str(a.flag()), a.cigarString(), marked);
			if (!optimal.equals(actualOptimal)) {
				failures.add("expected " + optimal + "\n   actual " + actualOptimal);
			}
			if (!f[3].equals(str(a.score2())) || !f[8].equals(str(a.refEnd2()))) {
				score2Differs++;
				if (sameLanesAsC) {
					failures.add(f[0] + " vs " + f[1] + ": score2/refEnd2 expected " + f[3] + "/" + f[8] + " actual "
							+ a.score2() + "/" + a.refEnd2());
				}
			}
		}
		System.out.printf("%s: %d pairs, %d via 16 bit kernel, %d byte lanes, %d sub-optimal differences%n", c.file(),
				pairs, wordPath, Lanes.BYTE.length(), score2Differs);
		assertTrue(pairs > 0);
		assertEquals(List.of(), failures.subList(0, Math.min(10, failures.size())),
				failures.size() + " of " + pairs + " alignments differ");
	}

	private static String markedString(Cigar.Marked m) {
		return m.mismatches() + "\t" + Cigar.toString(m.cigar());
	}

	private static String str(int i) {
		return Integer.toString(i);
	}

	static List<String> expected(String name) throws IOException {
		try (InputStream is = CReferenceTest.class.getResourceAsStream("/expected/" + name + ".tsv.gz");
				BufferedReader in = new BufferedReader(
						new InputStreamReader(new GZIPInputStream(is), StandardCharsets.US_ASCII))) {
			return in.lines().toList();
		}
	}
}
