package lapsolver.lsst;

import lapsolver.Graph;
import lapsolver.Tree;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.*;


public class SpanningTreeRunner implements SpanningTreeStrategy {
    private int nThreads;
    private SpanningTreeStrategy strategy;

    public SpanningTreeRunner(SpanningTreeStrategy strategy) {
        this(strategy, Runtime.getRuntime().availableProcessors());
    }

    public SpanningTreeRunner(SpanningTreeStrategy strategy, int nThreads) {
        this.strategy = strategy;
        this.nThreads = nThreads;
    }

    @Override
    public Tree getTree(final Graph input) {
        ExecutorService executorService = Executors.newFixedThreadPool(nThreads);
        List<Callable<Tree>> tasks = new LinkedList<>();
        for (int i = 0; i < nThreads; i++) {
            tasks.add(new Callable<Tree>() {
                @Override
                public Tree call() throws Exception {
                    return strategy.getTree(input);
                }
            });
        }
        Tree result = null;
        try {
            result = executorService.invokeAny(tasks);
        } catch (InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException("Error: Spanning tree construction was interrupted.", e);
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new RuntimeException("Error: Spanning tree construction failed unexpectedly.", e);
        } finally {
            executorService.shutdown();
        }
        return result;
    }
}
