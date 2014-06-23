package lapsolver.solvers;


import lapsolver.Graph;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * ConcurrentSolver runs multiple solvers in parallel, and returns the result
 * of the fastest one. You should NOT pass ConcurrentSolver another instance
 * of ConcurrentSolver
 */
public class ConcurrentSolver implements Solver {
    private final int nThreads;
    private long timeout;
    private Solver[] solvers;

    public ConcurrentSolver(Solver[] solvers) {
        this(solvers, Integer.MAX_VALUE);
    }

    public ConcurrentSolver(Solver[] solvers, long timeout) {
        this.timeout = timeout;
        this.solvers = solvers;
        nThreads = Math.min(Runtime.getRuntime().availableProcessors(), solvers.length);
    }

    @Override
    public void init(final Graph graph) {
        for (Solver solver : solvers)
            solver.init(graph);
    }

    @Override
    public double[] solve(final double[] b) {
        ExecutorService executor = Executors.newFixedThreadPool(nThreads);

        List<Callable<double[]>> runnableSolvers = new ArrayList<>();

        for (final Solver solver : solvers) {
            if (solver.getClass() == this.getClass())
                throw new RuntimeException("Undefined behavior for recursive invocation");

            runnableSolvers.add(new Callable<double[]>() {
                @Override
                public double[] call() throws Exception {
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
            e.printStackTrace();
        } finally {
            executor.shutdown();
        }
        return result;
    }
}
