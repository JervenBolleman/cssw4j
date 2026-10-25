package swiss.sib.cssw4j.tornado;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/** Inputs shared by the TornadoVM benchmarks: BLOSUM62, gap open 11, gap extend 1, and the test dataset. */
final class Workloads {
	static final int GAP_OPEN = 11, GAP_EXTEND = 1;
	/** The test dataset; {@code -Dcssw4j.data} overrides the default, which is relative to the tornado directory. */
	static final Path DATA = Path.of(System.getProperty("cssw4j.data", "../src/test/resources/dataset"));

	private Workloads() {
	}

	static Matrix blosum62() throws Exception {
		return Matrix.read(DATA.resolve("blosum62.txt"));
	}

	/** The 55 queries of the test dataset (human UniProt proteins and synthetic edge cases). */
	static Sequences datasetQueries(Matrix matrix) throws Exception {
		return Sequences.read(matrix, List.of(DATA.resolve("uniprot_human.fasta"), DATA.resolve("synthetic_queries.fasta")));
	}

	/** The 55 targets of the test dataset (homologs and synthetic edge cases). */
	static Sequences datasetTargets(Matrix matrix) throws Exception {
		return Sequences.read(matrix,
				List.of(DATA.resolve("uniprot_homologs.fasta"), DATA.resolve("synthetic_targets.fasta")));
	}

	/**
	 * The single largest pair of the test dataset (BRCA1 human, 1,863 residues, against BRCA1 mouse, 1,812): 3.4 million
	 * cells, 35 times the mean pair. Returns the query and the target, each as a one-sequence set.
	 */
	static Sequences[] largestPair(Matrix matrix) throws Exception {
		List<String[]> queries = fasta(DATA.resolve("uniprot_human.fasta"), DATA.resolve("synthetic_queries.fasta"));
		List<String[]> targets = fasta(DATA.resolve("uniprot_homologs.fasta"), DATA.resolve("synthetic_targets.fasta"));
		String[] bestQ = null, bestT = null;
		long best = -1;
		for (String[] q : queries) {
			for (String[] t : targets) {
				long cells = (long) q[1].length() * t[1].length();
				if (cells > best) {
					best = cells;
					bestQ = q;
					bestT = t;
				}
			}
		}
		return new Sequences[] { Sequences.of(matrix, List.of(bestQ[0]), List.of(bestQ[1])),
				Sequences.of(matrix, List.of(bestT[0]), List.of(bestT[1])) };
	}

	/** Records of FASTA files as {id, sequence}. */
	private static List<String[]> fasta(Path... files) throws Exception {
		List<String[]> out = new ArrayList<>();
		for (Path file : files) {
			String id = null;
			StringBuilder seq = new StringBuilder();
			for (String line : Files.readAllLines(file)) {
				if (line.startsWith(">")) {
					if (id != null) {
						out.add(new String[] { id, seq.toString() });
					}
					id = line.substring(1).split("[ \t]", 2)[0];
					seq.setLength(0);
				} else {
					seq.append(line.strip());
				}
			}
			if (id != null) {
				out.add(new String[] { id, seq.toString() });
			}
		}
		return out;
	}

	/**
	 * {@code count} random protein sequences of {@code length} residues over the 20 standard amino acids: a balanced
	 * workload in which every pair has the same number of cells. The seed makes it reproducible.
	 */
	static Sequences random(Matrix matrix, String prefix, int count, int length, long seed) {
		String letters = "ARNDCQEGHILKMFPSTWYV";
		Random rnd = new Random(seed);
		List<String> ids = new ArrayList<>(), seqs = new ArrayList<>();
		for (int i = 0; i < count; i++) {
			StringBuilder sb = new StringBuilder(length);
			for (int k = 0; k < length; k++) {
				sb.append(letters.charAt(rnd.nextInt(letters.length())));
			}
			ids.add(prefix + i);
			seqs.add(sb.toString());
		}
		return Sequences.of(matrix, ids, seqs);
	}

	/** Dynamic-programming cells of all query-target pairs. */
	static long cells(Sequences queries, Sequences targets) {
		long cells = 0;
		for (int q = 0; q < queries.count(); q++) {
			for (int t = 0; t < targets.count(); t++) {
				cells += (long) (queries.offsets[q + 1] - queries.offsets[q]) * (targets.offsets[t + 1] - targets.offsets[t]);
			}
		}
		return cells;
	}
}
