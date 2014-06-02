#import Queue, bisect

def myComp(x, y):
    if x[0] < y[0]:
        return -1
    elif x[0]==y[0]:
        return 0
    else:
        return 1
    
class PriorityQueue():
    def __init__(self):
        self.pq = []
    def empty(self):
        if len(self.pq) <= 0:
            return True
        return False
    def put(self, item):
        self.pq.append(item)
        self.pq.sort(myComp)
    def get(self):
        return self.pq.pop(0)

class graph:
    def __init__(self):
        self.adj = []  # list of tuples: first element is index of neighbor, 2nd is weight of edge
        self.edges = [] # list of edges (tuples like:(i, j, v)) in the order they are read in from the file
    def fromTXTijv(self, fileName):
        "load graph from ijv file"
        fp = open(fileName, "rt")
        lines = fp.readlines()
        line1 = (lines.pop(0)).split()
        numNodes = int(line1[0])
        self.adj = [[] for i in xrange(numNodes)]
        for line in lines:
            splitLine = line.split()
            i = int(splitLine[0]) - 1
            j = int(splitLine[1]) - 1
            v = 1.0/float(splitLine[2])
            self.adj[i].append((j, v))
            self.adj[j].append((i, v))
            self.edges.append((i, j, v))

    def fromTXTijv2(self, fp):
        "load graph ijv file pointer"
        numNodes = int(((fp.readline()).split())[0])
        #print "numNodes", numNodes
        self.adj = [[] for i in range(numNodes)]
        lines = [fp.readline() for i in range(numNodes - 1)]
        for line in lines:
            splitLine = line.split()
            i = int(splitLine[0]) - 1
            #print "i", i
            j = int(splitLine[1]) - 1
            v = 1.0/float(splitLine[2])
            self.adj[i].append((j, v))
            self.adj[j].append((i, v))
            self.edges.append((i, j, v))
    def calcBoundary(self, dists, t):
        "figures out boundary for radius of t"
        bdry = 0
        insideNodes = []
        seen = [False for i in range(len(self.adj))]
        for i, d in enumerate(dists):
            if d <= t:
                insideNodes.append(i)
                seen[i]= True
        for node in insideNodes:
            for nbr in self.adj[node]:
                if not seen[nbr]:
                    bdry += 1
                    seen[nbr] = True
        return bdry
    def calcConducts(self, orderedNodes):
        "conduct = tot border edges / tot edges with at least one node in set"
        out = []
        inside = [False for i in range(len(self.adj))]
        totWithin = 0.0
        totBorder = 0.0
        for node in orderedNodes:
            for (nbr, wt) in self.adj[node]:
                if inside[nbr]:
                    totBorder -= wt
                else:
                    totBorder += wt
                    totWithin += wt
            inside[node] = True
            out.append(totBorder/totWithin)
        return out
    def calcClustConducts(self, clusters):
        "conduct = tot border edges / tot edges with at least one node in set::cluster"
        out = []
        for clust in clusters:
            inside = [False for i in range(len(self.adj))]
            totWithin = 0.0
            totBorder = 0.0
            for node in clust:
                for (nbr, wt) in self.adj[node]:
                    if inside[nbr]:
                        totBorder -= wt
                    else:
                        totBorder += wt
                        totWithin += wt
                inside[node] = True
            out.append(totBorder/totWithin)
        return out
        
    def fromParray(self, fp):
        "looad graph from a text Parray"
        numNodes = int(fp.readline())
        for i in xrange(int(numNodes)):
            self.adj.append([])
        lines = [fp.readline() for i in xrange(numNodes)]
        for i, line in enumerate(lines):
            if i == int(line):
                continue
            self.adj[i].append((int(line), 1))   # Essentially we create an unweighted graph
            self.adj[int(line)].append((i, 1))
            self.edges.append((i, int(line), 1))
        #print self.adj
    def calcDists(self, node):
        """node is the index of nodeA. This function then returns the distance between this node and all others in graph
        this is your standard measure of distance by doing breadth first search. I assume that we are doing unweighted dists"""
        outDists = [0 for i in range(len(self.adj))]
        q = []
        seen = [False for i in range(len(self.adj))]
        q.append((node, 0))
        seen[node] = True
        while len(q) > 0:
            front = q.pop(0)
            curNode = front[0]
            curDist = front[1]
            for nbr, weight in self.adj[curNode]:
                if not seen[nbr]:
                    q.append((nbr, curDist + 1))
                    seen[nbr] = True
                    outDists[nbr] = curDist + 1
        graphSize = len(self.adj)
        for i in range(len(outDists)):
            outDists[i] /= float(graphSize)
        return outDists
    def calcDijkstraDists(self, node):
        """This is the real 'normal' dists function"""
        #print "calling calcDijkstraDists"
        outDists = [0.0 for i in range(len(self.adj))]
        seen = [False for i in range(len(self.adj))]
        pq = PriorityQueue()
        pq.put((0.0, node))
        #seen[node] = True
        while not pq.empty():
            #print "begin while loop"
            curDist, curNode = pq.get()
            if seen[curNode]:
                continue
            seen[curNode] = True
            outDists[curNode] = curDist
            #print "curNode", curNode, "curDist", curDist
            for nbr, wt in self.adj[curNode]:
                if seen[nbr]:
                    continue
                pq.put((wt + curDist, nbr))
        return outDists
    def isConnected(self, nodes):
        visited = [False for i in range(len(self.adj))]
        inCluster = [False for i in range(len(self.adj))]
        for nd in nodes:
            inCluster[nd] = True
        q = [nodes[0]]
        while(len(q) > 0):
            curNode = q.pop(0)
            visited[curNode] = True
            for nbr, wt in self.adj[curNode]:
                if inCluster[nbr] and not visited[nbr]:
                    q.append(nbr)
        for nd in nodes:
            if not visited[nd]:
                return False
        return True
    def calcConnectTime(self, node):
        "we see how many edges we have to add before node is connected to each of the other nodes"
        #print "node", node
        n = len(self.adj)
        components = [[i] for i in range(n)]# list of all components (components are just lists of nodes)
        node2comp = [i for i in range(n)]   # map from node idx to its component index
        dists = [0 for i in range(n)]       # the results we will return (in number of edges added before nodes connected)
        #print "self.edges", self.edges
        for edgeIdx, (i, j, v) in enumerate(self.edges):
            comp1 = node2comp[i]
            comp2 = node2comp[j]
            if len(components[comp1]) < len(components[comp2]):
                tmp1 = comp1
                comp1 = comp2
                comp2 = tmp1
            compConnected = None
            if node2comp[node] == comp1:
                compConnected = comp2
                #print "Xxx"
            elif node2comp[node] == comp2:
                compConnected = comp1
                #print "!!!"
            #print "i", i, "j", j, "v", v, "comp1", comp1, "comp2", comp2, "compConnected", compConnected, "node2comp[node]", node2comp[node]
            if compConnected is not None:
                #print "found compConnected:", compConnected
                for nd in components[compConnected]:
                    dists[nd] = edgeIdx + 1
            for nd in components[comp2]:
                node2comp[nd] = comp1
            components[comp1].extend(components[comp2])
            components[comp2] = []
            #print "node2comp", node2comp
            #print "comps", components
            #print "dists", dists
        #print "dists",dists
        #1/0
        return dists

    
def getDists(fileName, node):
    graphs = []
    fp = open(fileName)
    numGraphs = int(fp.readline())
    for i in range(numGraphs):
        g = graph()
        g.fromTXTijv2(fp)
        graphs.append(g)
    fp.close()
    dists = []
    timeDists = []
    for g in graphs:
        dists.append(g.calcDists(node))
        timeDists.append(g.calcConnectTime(node))
    avgDists = [0.0 for i in range(len(dists[0]))]
    avgTimeDists = [0.0 for i in range(len(timeDists[0]))]
    for dlist in dists:
        for i, d in enumerate(dlist):
            avgDists[i] += d/float(numGraphs)
    for dlist in timeDists:
        for i, d in enumerate(dlist):
            avgTimeDists[i] += d/float(numGraphs)
    return (avgDists, avgTimeDists)


def getOrderedNodes(distList):
    listToSort = []
    for i, val in enumerate(distList):
        listToSort.append((val, i))
    listToSort.sort(myComp)
    return [el[1] for el in listToSort]

def vertPrint(lists):
    "all lists must be same length"
    if len(lists) <= 0:
        return
    size = len(lists[0])
    for i in range(size):
        for elList in lists:
            print '%10f' %elList[i],
        print "\n",

DIST_METHOD = "CONNECT_TIME"

def getAvgDists(realGraph, graphs, node):
    if DIST_METHOD == "CONNECT_TIME":
        dists = []
        for g in graphs:
            dists.append(g.calcConnectTime(node))
        avgDists = [0.0 for i in range(len(dists[0]))]
        for dlist in dists:
            for i, d in enumerate(dlist):
                avgDists[i] += d/float(len(graphs))
        return avgDists
    elif DIST_METHOD == "DIJKSTRA":
        return realGraph.calcDijkstraDists(node)
    else:
        print "bad DIST_METHOD", 1/0

def clusterAround(realGraph, graphs, kNodes):
    "list of nodes to cluster around"
    k = len(kNodes)
    listOfDists = []
    clustersOut = [[] for i in range(k)]
    for nd in kNodes:
        listOfDists.append(getAvgDists(realGraph, graphs, nd))
       #now go through all nodes and assign it to the closest cluster
    n = len(listOfDists[0])
    #print listOfDists
    for nd in range(n):
        choices = []
        for idx, dists in enumerate(listOfDists):
            choices.append((dists[nd], idx))
        choices.sort(myComp)
        clustersOut[choices[0][1]].append(nd)
    return clustersOut

def findCentroids(realGraph, graphs, clusters):
    centersOut = []
    for clust in clusters:
        sumDists = []
        for nd in clust:
            dists = getAvgDists(realGraph, graphs, nd)
            distSum = 0.0
            for n in clust:
                distSum += dists[n]
            sumDists.append(distSum)
        choices = []
        for i, d in enumerate(sumDists):
            choices.append((d, i))
        choices.sort(myComp)
        centersOut.append(clust[choices[0][1]])
    return centersOut
    
def clustTest(graphFileName, treeFileName, kNodes, numIters):
    graphs = []
    fp = open(treeFileName)
    numGraphs = int(fp.readline())
    for i in range(numGraphs):
        g = graph()
        g.fromTXTijv2(fp)
        graphs.append(g)
    fp.close()
    g = graph()
    g.fromTXTijv(graphFileName)
    clusters = clusterAround(g, graphs, kNodes)
    print "clusters", clusters
    print "conducts", g.calcClustConducts(clusters)
    for clust in clusters:
        if not g.isConnected(clust):
            print "!!!UNCONNECTED CLUST", 1/0
    for iteration in range(numIters):
        centroids = findCentroids(g, graphs, clusters)
        print "centroids", centroids
        clusters = clusterAround(g, graphs, centroids)
        for clust in clusters:
            if not g.isConnected(clust):
                print "!!!UNCONNECTED CLUST", 1/0
        print "clusters", clusters
        print "conducts", g.calcClustConducts(clusters)
def runTest(graphFileName, treeFileName, node):
    (avgDists, avgTimeDists) = getDists(treeFileName, node)
    #avgDists = map(lambda x: 1.0/(x + .01), avgDists)
    #avgTimeDists = map(lambda x: 1.0/(x + .01), avgTimeDists)

    g = graph()
    g.fromTXTijv(graphFileName)
    avgDistsConducts = g.calcConducts(getOrderedNodes(avgDists))
    avgTimeDistsConducts = g.calcConducts(getOrderedNodes(avgTimeDists))
    dijkstraDists = g.calcDijkstraDists(node)
    vertPrint([range(len(avgDists)), avgDists, avgTimeDists, dijkstraDists])    
    dijkstraDistsConducts = g.calcConducts(getOrderedNodes(dijkstraDists))
    print "alpha dist\t beta dist \t dijkstra dist"
    vertPrint([avgDistsConducts, avgTimeDistsConducts, dijkstraDistsConducts])
    #print avgDists
    #print avgDistsConducts
    #print avgTimeDists
    #print avgTimeDistsConducts
    
    

"""
def distsFromParrays(fileName, node):
    graphs = []
    fp = open(fileName)
    numGraphs = int(fp.readline())
    for i in range(numGraphs):
        g = graph()
        g.fromParray(fp)
        graphs.append(g)
    fp.close()
    dists = []
    for g in graphs:
        dists.append(g.calcDists(node))
    avgDists = [0.0 for i in range(len(dists[0]))]
    for dlist in dists:
        for i, d in enumerate(dlist):
            avgDists[i] += d/float(numGraphs)
    return avgDists
"""
