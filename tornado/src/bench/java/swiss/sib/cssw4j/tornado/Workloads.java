package swiss.sib.cssw4j.tornado;

import java.nio.file.Path;
import java.util.List;

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
