package swiss.sib.cssw4j;

import java.util.List;

/**
 * Throughput benchmark on the protein test dataset, the counterpart of {@code src/test/c/ssw_bench.c}. Aligns all
 * queries against all targets with BLOSUM50, after a warm-up, and reports giga cell updates per second (GCUPS, query
 * length times target length per alignment, as for the C benchmark).
 *
 * <pre>
 * java --add-modules jdk.incubator.vector -cp target/classes:target/test-classes swiss.sib.cssw4j.Benchmark [rounds] [gapO] [gapE]
 * </pre>
 */
public final class Benchmark {
	private Benchmark() {
	}

	public static void main(String[] args) {
		int rounds = args.length > 0 ? Integer.parseInt(args[0]) : 20;
		int gapO = args.length > 1 ? Integer.parseInt(args[1]) : 10;
		int gapE = args.length > 2 ? Integer.parseInt(args[2]) : 2;
		ScoringMatrix m = ScoringMatrix.BLOSUM50;
		List<byte[]> queries = Fasta.resource(CReferenceTest.PROTEIN_QUERIES).stream().map(f -> m.encode(f.seq()))
				.toList();
		List<byte[]> targets = Fasta.resource(CReferenceTest.PROTEIN_TARGETS).stream().map(f -> m.encode(f.seq()))
				.toList();
		System.out.printf("%d queries x %d targets, %d byte lanes%n", queries.size(), targets.size(),
				Lanes.BYTE.length());
		for (int flag : new int[] { 0, Aligner.FLAG_ALWAYS }) {
			run(queries, targets, gapO, gapE, flag, Math.max(1, rounds / 4)); // warm-up
			run(queries, targets, gapO, gapE, flag, rounds);
		}
	}

	private static void run(List<byte[]> queries, List<byte[]> targets, int gapO, int gapE, int flag, int rounds) {
		Aligner aligner = new Aligner();
		long cells = 0, checksum = 0;
		long t0 = System.nanoTime();
		for (int r = 0; r < rounds; r++) {
			for (byte[] q : queries) {
				Profile p = Profile.of(q, ScoringMatrix.BLOSUM50);
				for (byte[] t : targets) {
					Alignment a = aligner.align(p, t, t.length, gapO, gapE, flag, 0, 0, q.length / 2);
					checksum += a.score1();
					cells += (long) q.length * t.length;
				}
			}
		}
		double s = (System.nanoTime() - t0) / 1e9;
		System.out.printf("Java   flag=%d gap=%d/%d: %8.3f s, %6.2f GCUPS (checksum %d)%n", flag, gapO, gapE, s,
				cells / s / 1e9, checksum);
	}
}
