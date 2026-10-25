import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;

/**
 * Reproducer: a value loaded from a byte[] and stored into an int[] (alone or after int arithmetic) is computed and
 * stored as a byte by TornadoVM 7.3.1. The generated kernel does the arithmetic in the char/signed-char type and writes
 * one byte into the int element, so values outside -128..127 are truncated and negative values lose their sign.
 * Seen with the CUDA, OpenCL (NVIDIA) and OpenCL (pocl) backends.
 */
public class ByteToInt {

	/** out[i] = a[i]: plain widening. */
	public static void copy(byte[] a, int[] out) {
		for (@Parallel int i = 0; i < out.length; i++) {
			out[i] = a[i];
		}
	}

	/** out[i] = a[i] + 1000. */
	public static void addConstant(byte[] a, int[] out) {
		for (@Parallel int i = 0; i < out.length; i++) {
			out[i] = a[i] + 1000;
		}
	}

	/** out[i] = a[i] * 24. */
	public static void multiplyConstant(byte[] a, int[] out) {
		for (@Parallel int i = 0; i < out.length; i++) {
			out[i] = a[i] * 24;
		}
	}

	/** Control: the byte mixed with a value derived from the loop index. This one is correct. */
	public static void withIndex(byte[] a, int[] out) {
		for (@Parallel int i = 0; i < out.length; i++) {
			out[i] = (i % 24) * 24 + a[i];
		}
	}

	interface Kernel {
		void apply(byte[] a, int[] out);
	}

	public static void main(String[] args) {
		int n = 256;
		byte[] a = new byte[n];
		for (int i = 0; i < n; i++) {
			a[i] = (byte) (i - 128); // every byte value, -128 .. 127
		}
		run("copy             out[i] = a[i]", a, ByteToInt::copy, (t, x, y) -> t.task("t", ByteToInt::copy, x, y));
		run("addConstant      out[i] = a[i] + 1000", a, ByteToInt::addConstant,
				(t, x, y) -> t.task("t", ByteToInt::addConstant, x, y));
		run("multiplyConstant out[i] = a[i] * 24", a, ByteToInt::multiplyConstant,
				(t, x, y) -> t.task("t", ByteToInt::multiplyConstant, x, y));
		run("withIndex        out[i] = (i % 24) * 24 + a[i]", a, ByteToInt::withIndex,
				(t, x, y) -> t.task("t", ByteToInt::withIndex, x, y));
	}

	interface AddTask {
		void add(TaskGraph graph, byte[] a, int[] out);
	}

	static void run(String name, byte[] a, Kernel java, AddTask task) {
		int[] out = new int[a.length], expected = new int[a.length];
		TaskGraph graph = new TaskGraph("s").transferToDevice(DataTransferMode.FIRST_EXECUTION, a);
		task.add(graph, a, out);
		graph.transferToHost(DataTransferMode.EVERY_EXECUTION, (Object) out);
		new TornadoExecutionPlan(graph.snapshot()).execute();
		java.apply(a, expected);
		int bad = 0, first = -1;
		for (int i = 0; i < out.length; i++) {
			if (out[i] != expected[i] && bad++ == 0) {
				first = i;
			}
		}
		System.out.println(bad == 0 ? "PASS " + name
				: "FAIL " + name + ": " + bad + " of " + out.length + " wrong, e.g. a[" + first + "] = " + a[first]
						+ ": expected " + expected[first] + ", got " + out[first]);
	}
}
