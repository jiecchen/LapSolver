CC = g++
CFLAGS= -O3 -Wall
#BOOST_DIR = /opt/local/include/
BOOST_DIR = /Users/spielman/progs/boost/boost_1_55_0

# make file for graphcluster. Must set BOOST_DIR in bash to where boost directory
clean:
	-rm graphcluster graphclusterpr graphclusterwpr lowstretch lowstretch2 subgraphtest calcstretchc treeimprove treeimprove2 *.o
calcer: calcer.o graph.o
	gcc $(CFLAGS) -o calcer calcer.o graph.o
spgraph: spgraph.o graph.o
	gcc $(CFLAGS) -o spgraph spgraph.o graph.o
cyclepoprandtree: cyclepoprandtree.o graph.o cyclepoprandtreemodule.o
	gcc $(CFLAGS) -o cyclepoprandtree cyclepoprandtree.o graph.o cyclepoprandtreemodule.o
cyclepoprandedgeprob: cyclepoprandedgeprob.o cyclepoprandtreemodule.o graph.o
	gcc $(CFLAGS) -o cyclepoprandedgeprob cyclepoprandedgeprob.o graph.o cyclepoprandtreemodule.o
randtree: randtree.o graph.o components.o weightedselection.o
	gcc $(CFLAGS) -o randtree randtree.o graph.o components.o weightedselection.o
graphdist: graphdist.o graph.o components.o weightedselection.o
	gcc $(CFLAGS) -o graphdist graphdist.o graph.o components.o weightedselection.o
randedgeprob: randedgeprob.o graph.o components.o
	gcc $(CFLAGS) -o randedgeprob randedgeprob.o graph.o components.o
doubledijkstra: doubledijkstra.o graph.o priorityq.o
	gcc $(CFLAGS) -o doubledijkstra doubledijkstra.o graph.o priorityq.o
dijkstratree: dijkstratree.o graph.o priorityq.o
	gcc $(CFLAGS) -o dijkstratree dijkstratree.o graph.o priorityq.o
treeimprove: treeimprove.o graph.o treeimprovetree.o
	gcc $(CFLAGS) -o treeimprove treeimprove.o graph.o treeimprovetree.o
treeimprove2: treeimprove2.o graph.o treeimprovetree.o priorityq.o
	gcc $(CFLAGS) -lm -o treeimprove2 treeimprove2.o graph.o treeimprovetree.o priorityq.o
treeimprove3: treeimprove3.o graph.o treeimprovetree.o
	gcc $(CFLAGS) -lm -o treeimprove3 treeimprove3.o graph.o treeimprovetree.o
treeimprove4: treeimprove4.o graph.o treeimprovetree.o
	gcc $(CFLAGS) -lm -o treeimprove4 treeimprove4.o graph.o treeimprovetree.o
randimprove: graph.o randimprove.o components.o weightedselection.o calcstretchctree.o
	gcc $(CFLAGS) -o randimprove randimprove.o graph.o components.o weightedselection.o calcstretchctree.o
calcstretchc: calcstretchc.o graph.o
	gcc $(CFLAGS) -o calcstretchc calcstretchc.o graph.o
calcstretchc2: calcstretchc2.o graph.o calcstretchctree.o
	gcc $(CFLAGS) -o calcstretchc2 calcstretchc2.o graph.o calcstretchctree.o
subgraphtest: subgraphtest.cpp
	$(CC) $(CFLAGS) -o subgraphtest subgraphtest.cpp  -L$(BOOST_DIR)
lowstretch2: lowstretch2.cpp
	$(CC) $(CFLAGS) -o lowstretch2 lowstretch2.cpp  -I$(BOOST_DIR)
lowstretch: lowstretch.cpp
	$(CC) $(CFLAGS) -o lowstretch lowstretch.cpp -I$(BOOST_DIR)
#	$(CC) $(CFLAGS) -o lowstretch lowstretch.cpp -L/Users/spielman/progs/boost/boost_1_55_0
graphcluster: graphcluster.o dijkstra.o
	$(CC) -o graphcluster graphcluster.o dijkstra.o -L$(BOOST_DIR)
graphclusterpr: graphclusterpr.o dijkstra.o 
	$(CC) -o graphclusterpr graphclusterpr.o dijkstra.o graph.o -L$(BOOST_DIR)
graphclusterwpr: graphclusterwpr.o dijkstra.o weightedPageRank.o
	$(CC) -o graphclusterwpr graphclusterwpr.o dijkstra.o graph.o weightedPageRank.o -L$(BOOST_DIR)
spgraph.o: spgraph.c graph.h
	gcc $(CFLAGS) -c spgraph.c
calcer.o: calcer.c graph.h
	gcc $(CFLAGS) -c calcer.c
cyclepoprandtree.o: cyclepoprandtree.c graph.h
	gcc $(CFLAGS) -c cyclepoprandtree.c
cyclepoprandtreemodule.o: cyclepoprandtreemodule.c graph.h weightedselection.h
	gcc $(CFLAGS) -c cyclepoprandtreemodule.c
cyclepoprandedgeprob.o: cyclepoprandedgeprob.c graph.h cyclepoprandtreemodule.h
	gcc $(CFLAGS) -c cyclepoprandedgeprob.c
randtree.o: randtree.c graph.h components.h weightedselection.h
	gcc $(CFLAGS) -c randtree.c
graphdist.o: graphdist.c graph.h components.h weightedselection.h
	gcc $(CFLAGS) -c graphdist.c
randedgeprob.o: randedgeprob.c graph.h  components.h
	gcc $(CFLAGS) -c randedgeprob.c
components.o: components.c components.h
	gcc $(CFLAGS) -c components.c
weightedselection.o: weightedselection.c weightedselection.h
	gcc $(CFLAGS) -c weightedselection.c
doubledijkstra.o: doubledijkstra.c graph.h priorityq.h
	gcc $(CFLAGS) -c doubledijkstra.c
treeimprovetree.o: treeimprovetree.c graph.h
	gcc $(CFLAGS) -c treeimprovetree.c
treeimprove.o: treeimprove.c treeimprovetree.h graph.h
	gcc $(CFLAGS) -c treeimprove.c
treeimprove2.o: treeimprove2.c treeimprovetree.h graph.h
	gcc $(CFLAGS) -c treeimprove2.c
treeimprove3.o: treeimprove3.c treeimprovetree.h graph.h
	gcc $(CFLAGS) -c treeimprove3.c
treeimprove4.o: treeimprove4.c treeimprovetree.h graph.h
	gcc $(CFLAGS) -c treeimprove4.c
randimprove.o: randimprove.c components.h weightedselection.h calcstretchctree.h
	gcc $(CLFAGS) -c randimprove.c
calcstretchc.o: calcstretchc.c
	gcc $(CFLAGS) -c calcstretchc.c
calcstretchc2.o: calcstretchc2.c calcstretchctree.h
	gcc $(CFLAGS) -c calcstretchc2.c
calcstretchctree.o: calcstretchctree.h calcstretchctree.c
	gcc $(CFLAGS) -c calcstretchctree.c
graph.o: graph.c graph.h
	gcc $(CFLAGS) -c graph.c
dijkstratree.o: dijkstratree.c graph.h priorityq.h
	$gcc $(CFLAGS) -c dijkstratree.c
priorityq.o: priorityq.c
	$gcc $(CFLAGS) -c priorityq.c

graphcluster.o: graphcluster.cpp # graphcluster.hpp
	$(CC) $(CFLAGS) -c graphcluster.cpp

graphclusterpr.o: graphclusterpr.cpp #graphclusterpr.hpp
	$(CC) $(CFLAGS) -c graphclusterpr.cpp

graphclusterwpr.o: graphclusterwpr.cpp profiler.hpp #graphclusterwpr.hpp 
	$(CC) $(CFLAGS) -c graphclusterwpr.cpp

dijkstra.o: dijkstra.cpp
	$(CC) $(CFLAGS) -c dijkstra.cpp

weightedPageRank.o: weightedPageRank.cpp profiler.hpp weightedPageRank.hpp
	$(CC) $(CFLAGS) -c weightedPageRank.cpp

graphcluster.o graphclusterpr.o dijkstra.o dijkstra2.o weightedPageRank.o: graph.hpp dijkstra.hpp
