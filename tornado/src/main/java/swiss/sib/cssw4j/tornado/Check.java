package swiss.sib.cssw4j.tornado;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Compares the TornadoVM kernel with the results of the original C library.
 *
 * <pre>
 * Check jvm|tornado matrix.txt queries.fa[,more.fa] targets.fa[,more.fa] gapOpen gapExtend expected.tsv.gz [interleaved|contiguous]
 * </pre>
 *
 * Mode {@code jvm} runs the kernel as plain Java on the CPU and needs no device; mode {@code tornado} runs it
 * through TornadoVM on the selected device. The compared fields are score, target end and query end of every pair,
 * in the order of the reference file (all targets for the first query, then the second query, ...).
 */
public final class Check {
	private Check() {
	}

	public static void main(String[] args) throws Exception {
		if (args.length != 7 && args.length != 8) {
			System.err.println("usage: Check jvm|tornado matrix.txt queries.fa[,..] targets.fa[,..] gapOpen gapExtend expected.tsv.gz [interleaved|contiguous]");
			System.exit(2);
		}
		boolean tornado = args[0].equals("tornado");
		Matrix matrix = Matrix.read(Path.of(args[1]));
		Sequences queries = Sequences.read(matrix, paths(args[2]));
		Sequences targets = Sequences.read(matrix, paths(args[3]));
		int gapOpen = Integer.parseInt(args[4]), gapExtend = Integer.parseInt(args[5]);
		Layout layout = args.length == 8 ? Layout.parse(args[7]) : Layout.INTERLEAVED;
		int pairs = queries.count() * targets.count();
		System.out.printf("%d queries x %d targets = %d pairs, matrix %dx%d, %s scratch layout%n", queries.count(), targets.count(),
				pairs, matrix.size(), matrix.size(), layout);

		int[] scores, refEnds, readEnds;
		if (tornado) {
			System.out.println("device: " + TornadoScorer.device().getDescription() + " | backend " + uk.ac.manchester.tornado.api.runtime.TornadoRuntimeProvider.getTornadoRuntime().getBackend(0).getClass().getName());
			try (TornadoScorer scorer = new TornadoScorer(matrix, queries, targets, gapOpen, gapExtend, layout)) {
				long t0 = System.nanoTime();
				scorer.run();
				long t1 = System.nanoTime();
				System.out.printf("first run (includes kernel compilation): %.1f ms%n", (t1 - t0) / 1e6);
				for (int i = 0; i < 3; i++) {
					t0 = System.nanoTime();
					scorer.run();
					t1 = System.nanoTime();
					System.out.printf("run %d: %.1f ms%n", i + 2, (t1 - t0) / 1e6);
				}
				scores = scorer.scores;
				refEnds = scorer.refEnds;
				readEnds = scorer.readEnds;
			}
		} else {
			scores = new int[pairs];
			refEnds = new int[pairs];
			readEnds = new int[pairs];
			int[] steps = layout.steps(pairs, queries.maxLength);
			int[] params = { gapOpen, gapExtend, matrix.size(), steps[0], steps[1] };
			int[] h = new int[(queries.maxLength + 1) * pairs], e = new int[h.length];
			long t0 = System.nanoTime();
			SwKernel.scorePairs(params, matrix.scores(), queries.residues, queries.offsets, targets.residues,
					targets.offsets, h, e, scores, refEnds, readEnds);
			System.out.printf("sequential Java: %.1f ms%n", (System.nanoTime() - t0) / 1e6);
		}

		List<String[]> expected = new ArrayList<>();
		try (BufferedReader in = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(Files.newInputStream(Path.of(args[6]))), StandardCharsets.US_ASCII))) {
			for (String line = in.readLine(); line != null; line = in.readLine()) {
				expected.add(line.split("\t"));
			}
		}
		if (expected.size() != pairs) {
			System.out.printf("FAIL: expected file has %d pairs, ran %d%n", expected.size(), pairs);
			System.exit(1);
		}
		int bad = 0;
		for (int p = 0; p < pairs; p++) {
			String[] f = expected.get(p);
			boolean same = f[0].equals(queries.ids.get(p / targets.count())) && f[1].equals(targets.ids.get(p % targets.count()))
					&& Integer.parseInt(f[2]) == scores[p] && Integer.parseInt(f[5]) == refEnds[p]
					&& Integer.parseInt(f[7]) == readEnds[p];
			if (!same && bad++ < 5) {
				System.out.printf("  pair %d (%s vs %s): expected score %s refEnd %s readEnd %s, got %d %d %d%n", p, f[0], f[1],
						f[2], f[5], f[7], scores[p], refEnds[p], readEnds[p]);
			}
		}
		System.out.println(bad == 0 ? "PASS: all " + pairs + " pairs match the C reference (score, target end, query end)"
				: "FAIL: " + bad + " of " + pairs + " pairs differ");
		System.exit(bad == 0 ? 0 : 1);
	}

	private static List<Path> paths(String list) {
		List<Path> out = new ArrayList<>();
		for (String s : list.split(",")) {
			out.add(Path.of(s));
		}
		return out;
	}
}
