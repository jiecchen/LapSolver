package lapsolver.solvers;


import lapsolver.Graph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.*;

/**
 * ConcurrentSolver runs multiple solvers in parallel, and returns the result
 * of the fastest one. You should NOT pass ConcurrentSolver another instance
 * of ConcurrentSolver
 */
public class ConcurrentSolver extends Solver {
    private final int nThreads;
    private long timeout;
    private ArrayList<Solver> solvers;

    public ConcurrentSolver(Solver[] solvers) {
        this(solvers, Integer.MAX_VALUE);
    }

    public ConcurrentSolver(Solver[] solvers, long timeout) {
        this.timeout = timeout;
        this.solvers = new ArrayList<>(Arrays.asList(solvers));
        nThreads = Math.min(Runtime.getRuntime().availableProcessors(), solvers.length);
    }

    @Override
    public void init(final Graph graph, double[] d) {
        this.graph = graph; this.d = d;
    }

    @Override
    public double[] solve(final double[] b) {
        ExecutorService executor = Executors.newFixedThreadPool(nThreads);

        List<Callable<double[]>> runnableSolvers = new ArrayList<>();
        HashSet<String> names = new HashSet<>();

        for (final Solver solver : solvers) {
            int i = 2;
            String className = solver.getClass().getName();
            String name = className;
            while (names.contains(name)) {
                name = className + i++;
            }
            names.add(name);

            runnableSolvers.add(new Task<double[]>(name) {
                @Override
                public double[] run() {
                    solver.init(graph);
                    return solver.solve(b);
                }
            });
        }

        double[] result = null;
        try {
            result = executor.invokeAny(runnableSolvers, timeout, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            System.err.println("[error] ConcurrentSolver: execution was interrupted");
            e.printStackTrace();
        } catch (ExecutionException e) {
            System.err.println("[error] ConcurrentSolver: " + e.getCause().getMessage());
            e.printStackTrace();
        } catch (TimeoutException e) {
            System.err.println("[error] ConcurrentSolver: none of the solvers finished in time");
        } finally {
            executor.shutdownNow();
        }
        return result;
    }

    private abstract class Task<V> implements Callable<V> {
        private String name;

        protected Task(String name) {
            this.name = name;
        }

        @Override
        public V call() throws Exception {
            V result = run();
            if (Thread.interrupted()) {
                return null;
            }
            System.out.println("[info] Task " + name + " finished");
            return result;
        }

        public abstract V run();
    }
}
