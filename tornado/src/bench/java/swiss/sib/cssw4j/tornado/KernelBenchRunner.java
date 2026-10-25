package swiss.sib.cssw4j.tornado;

import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

/**
 * Runs {@link KernelBench} with its options built in code, so no JMH command line is needed.
 * <p>
 * Usage, from the tornado directory with the TornadoVM SDK and a JDK 22+ (the CUDA SDK needs CUDA_HOME and
 * LD_LIBRARY_PATH as in README.md):
 *
 * <pre>
 * tornado --classpath target/classes:target/bench:JMH_JARS swiss.sib.cssw4j.tornado.KernelBenchRunner [-Dcssw4j.device=0:0]
 * </pre>
 *
 * The first argument, if given, is a regular expression selecting benchmarks.
 */
public final class KernelBenchRunner {
	private KernelBenchRunner() {
	}

	public static void main(String[] args) throws Exception {
		Options options = new OptionsBuilder() //
				.include(args.length > 0 ? args[0] : KernelBench.class.getSimpleName()) //
				.build();
		new Runner(options).run();
	}
}
