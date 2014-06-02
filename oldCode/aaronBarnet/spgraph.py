#!/usr/bin/env python
"""\
 spgraph.py - A set of utilities to create spgraphs
"""
import Polynomial
import random
import os
random.seed(0)
#uidCtr = 0

GRAPH2DERIV = False

def findPolyThresh(poly):
	"find value that makes poly evaluate to .5, binary search"
	# do 20 iterations
	low = 0.0
	high = 1.0
	for i in xrange(20):
		middle = (high - low)/2.0
		if poly(middle) > .5:
			high = middle
		else:
			low = middle
	return (high - low)/2.0

class Edge:
	"""has entries for 
		nbr node ref 'nbr'
		wt of edge"""
	def __init__(self, node, wt, poly=None, polyVal=None):
		"poly: Polynomial object   polyVal: specific float used to evalue poly at specific value"
		self.nbr = node
		self.wt = wt
		if poly is None:
			self.poly = Polynomial.Polynomial('x')
		else:
			self.poly = poly
		self.polyVal = polyVal
		
def removeEdge(edges, node, number=None):
	"search through list and remove Edge that has nbr of node, number is max num removes"
	#print "edges:",
	#for edge in edges:
	#	print edge.nbr,
	#print "\n",
	removeList = []
	for edge in edges:
		if edge.nbr is node:
			removeList.append(edge)
	if number is not None:
		for i in xrange(number):
			edges.remove(removeList[i])
	else:
		for edge in removeList:
			edges.remove(edge)
def countEdges(edges, node):
	"count how many edges connect to node"
	ctr = 0
	for edge in edges:
		if edge.nbr is node:
			ctr += 1
	return ctr
			
def funk(a, b):
	print a, b
class Node:
	"""has list of edges, each edge is a list of form [nbr, wt]"""
	def __init__(self):
		self.edges = []
		self.name = "unnamed"
		#self.uid = random.random()#uidCtr
		#uidCtr += 1
	def adjoin(self, edges):
		self.edges.extend(edges)
		for edge in edges:
			nbr = edge.nbr
			nbr.edges.append(Edge(self, edge.wt, edge.poly, edge.polyVal))
		return self
	def merge(self, n1, n2):
		"""combine edges from n1 and n2 into self
		make other nodes refer to self instead of n1 or n2
		Assumes that there is no overlap"""
		#print "n1:", n1, "n2:", n2
		self.edges.extend(n1.edges)
		self.edges.extend(n2.edges)
		for edge in n1.edges:
			edge.nbr.edges.append(Edge(self, edge.wt))
			#print "removeEdge to n1:", n1, "\n", edge.nbr.edges
			removeEdge(edge.nbr.edges, n1)
			#print edge.nbr.edges, "\n"
		for edge in n2.edges:
			edge.nbr.edges.append(Edge(self, edge.wt))
			#print "removeEdge to n2:", n2, "\n", edge.nbr.edges
			removeEdge(edge.nbr.edges, n2)
			#print edge.nbr.edges, "\n"
		return self
	def remDups(self, nbr):
		"remove any extra edges between self and nbr"
		lastEdge = None
		removeList = []
		for edge in self.edges:
			if lastEdge is not None and edge.nbr is nbr:
				removeList.append(lastEdge) #self.edges.remove(lastEdge)
			if edge.nbr is nbr:
				lastEdge = edge
		for edge in removeList:
			self.edges.remove(edge)
		return self

			
class SPgraph:
	"self.graph is just a list of nodes"
	def __init__(self):
		self.graph = []
	def printMe(self):
		ctr = 0
		for node in self.graph:
			node.name = ctr
			ctr += 1
		print "s", self.s.name
		print "t", self.t.name
		for node in self.graph:
			print str(node.name) + ":",
			for edge in node.edges:
				print "(", edge.nbr.name, ",", edge.wt, ",", edge.poly, "),",
			print "\n",
		#deal with integral if graph is collapsed
		if len(self.graph) == 2:
			print "Evaluating integral:"
			intPoly = self.s.edges[0].poly.integral()
			print "value:", intPoly(0.0,1.0)
			print "thresh p:", findPolyThresh(self.s.edges[0].poly)
	def copyMe(self):
		"create a duplicate of the graph (usefull if want to keep a copy of graph after stcollapse"
		newG = SPgraph()
		for node in self.graph:
			newG.graph.append(Node())
		for i in range(len(self.graph)):
			for edge in self.graph[i].edges:
				toNode = newG.graph[self.graph.index(edge.nbr)]
				newEdge = Edge(toNode, edge.wt)  #avoid edge poly in this one
				#newEdge = Edge(toNode, edge.wt, edge.poly)  #node/wt/poly
				newG.graph[i].edges.append(newEdge)
		newG.s = newG.graph[self.graph.index(self.s)]
		newG.t = newG.graph[self.graph.index(self.t)]
		return newG
	
	def outputPoly(self, fileName="out.dat"):
		"assumes stcollapse has not be run, calcs values on the fly"
		if len(self.graph) == 2:
			print "ERROR: DON'T call stCollapse first"
			return
		fp = open(fileName, "w")
		xLst = [i/100.0 for i in range(101)]
		pVals = evalPolys(self, xLst)
		dVals = calDerivatives(xLst, pVals)
		ddVals = calDerivatives(xLst, dVals)
		integral = calIntegral(xLst,pVals)
		for i in xrange(len(xLst)-2):
			fp.write(str(xLst[i]))
			fp.write("\t")
			fp.write(str(pVals[i]))
			fp.write("\t")
			fp.write(str(dVals[i]))
			fp.write("\t")
		       	fp.write(str(ddVals[i]))
			fp.write("\n")
		fp.close()

		copyGraph = self.copyMe()
		copyGraph.stCollapse()

		fp2 = open("plot.gplot", "w")
		if GRAPH2DERIV:
			fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"polynomial\","
				  +"\"" + str(fileName) + "\" using 1:3 title \"derivative\","
				  +"\"" + str(fileName) + "\" using 1:4 title \"2nd derivative\","
				  + str(copyGraph.s.edges[0].wt) + " title \"ER\","
				  + str(integral) + " title \"integral\"")
		else:
       			fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"polynomial\","
				  +"\"" + str(fileName) + "\" using 1:3 title \"derivative\","
				  + str(1.0 - copyGraph.s.edges[0].wt) + " title \"1.0 - ER\","
				  + str(integral) + " title \"integral\"")
		fp2.close()
		print "ER:", copyGraph.s.edges[0].wt
		print "Approx ER:", 1.0 - integral
		os.system('gnuplot -persist plot.gplot')
# 	def outputPoly(self, fileName="out.dat"):
# 		"assumes that stCollapse has already been run, writes out values of st poly from 0.0 to 1.0 in .001 increments"
# 		if len(self.graph) != 2:
# 			print "ERROR: call stCollapse first"
# 			return
# 		fp = open(fileName, "w")
# 		step = .001
# 		for i in xrange(0, int(1.0/.001)):
# 			val = i * step
# 			fp.write(str(val))
# 			fp.write("\t")
# 			fp.write(str(self.graph[0].edges[0].poly(val)))
# 			fp.write("\n")
# 		fp.close()
# 		fp2 = open("plot.gplot", "w")
# 		fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"polynomial\"")
# 		fp2.close()
# 		os.system('gnuplot -persist plot.gplot')
	

	def seriesCollapse(self, node):
		"given a node in a graph of degree 2, eliminate node c = c + c"
		if len(node.edges) != 2:
			raise Exception("node is not degree 2!")
		node1 = node.edges[0].nbr
		node2 = node.edges[1].nbr
		newWt = node.edges[0].wt + node.edges[1].wt
		newPoly = None
		newPolyVal = None
		if not node.edges[0].poly is None:
			newPoly = node.edges[0].poly * node.edges[1].poly
		if not node.edges[0].polyVal is None:
			newPolyVal = node.edges[0].polyVal * node.edges[1].polyVal
		
		removeEdge(node1.edges, node)
		removeEdge(node2.edges, node)
		#print "newWt:", newWt
		node1.edges.append(Edge(node2, newWt, newPoly, newPolyVal))
		node2.edges.append(Edge(node1, newWt, newPoly, newPolyVal))
		self.graph.remove(node)
	def parallelCollapse(self, node1, node2):
		"given two nodes with two edges together, collapse edges c = (1/(1/c + 1/c))"
		oldWts = []
		oldPolys = []
		oldPolyVals = []
		for edge in node1.edges:
			if edge.nbr is node2:
				#if edge.wt <= 0.0:
				#	print "wt:", edge.wt
				oldWts.append(edge.wt)
				oldPolys.append(edge.poly)
				oldPolyVals.append(edge.polyVal)
		if len(oldWts) < 2:
			raise(Exception("there are not at least 2 edges between node1 and node2"))
		# FOR NOW ONLY join first two edges so we can deal with polynomial!
		#newWtBottom = 0.0
		#for wt in oldWts:
		#	newWtBottom += 1.0/wt
		#newWt = 1.0/newWtBottom
		
		newWt = 1.0/(1.0/oldWts[0] + 1.0/oldWts[1])
		if not oldPolys[0] is None:
			newPoly = (oldPolys[0] + oldPolys[1]) - (oldPolys[0] * oldPolys[1])
		else:
			newPoly = None
		if not oldPolyVals[0] is None:
			newPolyVal = (oldPolyVals[0] + oldPolyVals[1]) - (oldPolyVals[0] * oldPolyVals[1])
		else:
			newPolyVal = None
		#!!! if we are 
		#self.printMe()
		removeEdge(node1.edges, node2, 2)
		removeEdge(node2.edges, node1, 2)
		#self.printMe()
		node1.edges.append(Edge(node2, newWt, newPoly, newPolyVal))
		node2.edges.append(Edge(node1, newWt, newPoly, newPolyVal))
		#self.printMe()
		#raise(Exception())
		#print "did pcollapse\n\n\n\n\n"
	def stCollapse(self):
		"""essentially we calculate effective resistance of st edge in graph:
		we collapse all other edges using seriesCollapse & parallelCollapse
		until we are left with just the st edge with effective resistance as its wt"""
		while len(self.graph) > 2 or len(self.s.edges) > 1:
			for node in self.graph:
				if node is self.s or node is self.t:
					continue
				if len(node.edges) == 2:
					self.seriesCollapse(node)
			for node1 in self.graph:
				for edge in node1.edges:
					if countEdges(node1.edges, edge.nbr) > 1:
						self.parallelCollapse(node1, edge.nbr)
						break
					

		
def basicSPgraph():
	"create k2"
	graph = SPgraph()
	graph.graph = [Node(), Node()]
	graph.graph[0].adjoin([Edge(graph.graph[1], 1.0)])
	graph.s = graph.graph[0]
	graph.t = graph.graph[1]
	return graph

def randSPgraph(size, probParallel=.5, maxParallel=3):
	"""randomly  put together S-P graph randomly, call randGraph wrapper!
	size: recursive depth
	probParallel: probability to combine graphs parallel instead of in series
	"""
	if size == 0:
		graph = basicSPgraph()
		return graph
	else:
		graph1 = randSPgraph(size - 1, probParallel, maxParallel)
		randFloat = random.random()   #between 0 and 1.0
		merged = None
		if randFloat > probParallel:      # combine by series
			newGraph = randSPgraph(size - 1, probParallel, maxParallel)
			merged = seriesMerge(graph1, newGraph)
		else:
			merged = graph1
			for i in xrange(random.randrange(1, maxParallel - 1, 1)):
				newGraph = randSPgraph(size - 1, probParallel, maxParallel)
				merged = parallelMerge(merged, newGraph)
		return merged
		
def randGraph(size, probParallel=.5, maxParallel=3):
	"""wrapper for randSPgraph: ensures that there is st edge and that st edge has poly of 0
	see randSPgraph for parameters"""
	graph = randSPgraph(size, probParallel, maxParallel)
	basic = basicSPgraph()
	merged = parallelMerge(basic, graph)
	for edge in merged.s.edges:
		if edge.nbr is merged.t:
			edge.poly = Polynomial.Polynomial("0")
	for edge in merged.t.edges:
		if edge.nbr is merged.s:
			edge.poly = Polynomial.Polynomial("0")
	return merged

def parallelMerge(g, h):
	"fuse both s's and t's together, eliminate extra edges"
	mergedS = Node()
	mergedS.merge(g.s, h.s)
	mergedT = Node()
	mergedT.merge(g.t, h.t)
	graph = SPgraph()
	graph.graph.append(mergedS)
	graph.graph.append(mergedT)
	graph.graph.extend(g.graph)
	graph.graph.extend(h.graph)
	graph.graph.remove(g.s)
	graph.graph.remove(g.t)
	graph.graph.remove(h.s)
	graph.graph.remove(h.t)
	graph.s = mergedS
	graph.t = mergedT
	# eliminate double edge s-t
	graph.s.remDups(graph.t)
	graph.t.remDups(graph.s)
	return graph
		
def seriesMerge(g, h):
	"fuse t of g with s of h"
	mergedNode = Node()             # must merge first before using extend (so that side effects work correctly)
	mergedNode.merge(g.t, h.s)
	graph = SPgraph()
	graph.graph.extend(g.graph)
	graph.graph.extend(h.graph)
	graph.s = g.s
	graph.t = h.t
	#print "g.t", g.t, "h.s", h.s
	#print "mergedNode", mergedNode
	graph.graph.append(mergedNode)
	graph.graph.remove(g.t)
	graph.graph.remove(h.s)
	return graph

def evalPolys(g, valLst):
	"""g is a graph
	valLst is a list of float at which to evaluate polynomial"""
	outVals = []
	for v in valLst:
		newG = g.copyMe()
		for node in newG.graph:
			for edge in node.edges:
				# polyVal's on s and t must be 0
				if (edge.nbr is newG.s and node is newG.t) or (edge.nbr is newG.t and node is newG.s):
					edge.polyVal = 0.0
				else:
					edge.polyVal = v
		newG.stCollapse()
		outVals.append(newG.s.edges[0].polyVal)
	return outVals

def calIntegral(xLst, yLst):
	"given a list of x's and y's"
	out = 0.0
	for i in xrange(len(xLst)-1):
		out += (yLst[i+1] + yLst[i])/2.0 * (xLst[i+1]-xLst[i])
	return out

def calDerivatives(xLst, yLst):
	out = []
	for i in xrange(len(yLst)-1):
		out.append((yLst[i+1] - yLst[i])/(xLst[i+1] - xLst[i]))
	return out



def runTest(numGraphs, size, probParallel, maxParallel):
	"randomly pick bunch of graphs and figure out difference between two probs"
	graphs = []
	approxERs = []
	ERs = []
	for i in xrange(numGraphs):
		graphs.append(randGraph(size, probParallel, maxParallel))
	for graph in graphs:
		graph.stCollapse()
		ERs.append(graph.graph[0].edges[0].wt)
		intPoly = graph.graph[0].edges[0].poly.integral()
		approxERs.append(1.0 - intPoly(0.0,1.0))
	maxOvershoot = [0.0, 0.0, 0.0]    # 1st entry is overshoo, 2nd is ER, 3rd is approxER
	maxUndershoot = [100.0, 0.0, 0.0]
	for i in xrange(len(ERs)):
		overshoot = approxERs[i] / ERs[i]
		if overshoot > maxOvershoot[0]:
			maxOvershoot[0] = overshoot
			maxOvershoot[1] = ERs[i]
			maxOvershoot[2] = approxERs[i]
		if overshoot < maxUndershoot[0]:
			maxUndershoot[0] = overshoot
			maxUndershoot[1] = ERs[i]
			maxUndershoot[2] = approxERs[i]
	print "maxOvershoot:", maxOvershoot
	print "maxUndershoot:", maxUndershoot
	
def seriesGraph(n):
	"create a series of length n"
	basics = []
	for i in xrange(n):
		basics.append(basicSPgraph())
	out = basics[0]

	for i in xrange(1, n):
		out = seriesMerge(out, basics[i])
	return out

def parallelComp(lensLst):
	"get all the seriesGraph's according to the lensLst, then parallel them all together"
	series = []
	for length in lensLst:
		series.append(seriesGraph(length))
	series.append(seriesGraph(1))
	out = series[0]
	for i in xrange(1,len(series)):
		out = parallelMerge(out, series[i])
	return out

def calcSimpleParER(numParallels, lengthSeries):
	"""return the effective resistence of the graph formed by parallel combining 'numParallels' series graphs of length 'lengthSeries'
	numParallels does not include the length one edge that we are looking at"""
	return float(lengthSeries) / (numParallels + float(lengthSeries))

def binomial(n,k):
	"from http://neverland.ncssm.edu/~morrison/tcm05.html"
	p = 1.0
	for j in range(0,k):
		p = p*float(n - j)/float(j + 1)
	return p

def calcSimpleParApprox(numParallels, lengthSeries):
	"""return the approximation of the effective resistence of the graph formed by parallel combining 'numParallels' series graphs of length 'lengthSeries'
	numParallels does not include the length one edge that we are looking at """
	out = 0.0
	k = float(lengthSeries)
	for i in xrange(1, numParallels+1):
		out += ((-1)**(i+1)) * binomial(numParallels, i) / float(k * i + 1)
	return 1.0 - out

def parallelCompTest2(length, numGraphs, fileName="out.dat"):
	"same as parallelCompTest, except it calculates everything without actually dealing with graphs"
	approxs = []
	reals = []
	for test in xrange(1,numGraphs):
		approxs.append(calcSimpleParApprox(test, length))
		reals.append(calcSimpleParER(test, length))
	print "approxs", approxs
	print "reals", reals
	fp = open(fileName, "w")
	for i in xrange(len(approxs)):
		fp.write(str(i+1) + "\t" +  str(approxs[i]/reals[i]) + "\n")
	fp.close()
	fp2 = open(fileName + ".gplot", "w")
	fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"approx/real\"")
	fp2.close()
	os.system("gnuplot -persist " + fileName + ".gplot")


def parallelCompTest(length, numGraphs, fileName="out.dat"):
	"""length is how long the series graphs are
	numGraphs is how many graphs we are putting together"""
	approxs = []
	reals = []
	for test in xrange(1,numGraphs):
		lensLst = [length for i in range(test)]
		g = parallelComp(lensLst)
		xLst = [i/100.0 for i in range(101)]
		pVals = evalPolys(g,xLst)
		approx = 1.0 - calIntegral(xLst,pVals)
		g.stCollapse()
		real = g.s.edges[0].wt
		approxs.append(approx)
		reals.append(real)
	print "approxs", approxs
	print "reals", reals
	fp = open(fileName, "w")
	for i in xrange(len(approxs)):
		fp.write(str(i+1) + "\t" +  str(approxs[i]/reals[i]) + "\n")
	fp.close()
	fp2 = open("plot.gplot", "w")
	fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"approx/real\"")
	fp2.close()
	os.system('gnuplot -persist plot.gplot')


def parallelCompTestInterp(length1, length2, fileName="out.dat"):
	"""length1 and length2 are different series graph lengths
	NIXED:numGraphs is how many graphs we are putting together (divided between lenth1 and length2)"""
	approxs = []
	reals = []
	if length1 > length2:
		tmp = length1
		length1 = length2
		length2 = tmp

	for test in xrange(0, length2 - length1 + 1):
		lensLst1 = [length1 for i in range(test)]
		lensLst2 =  [length2 for i in range(length2 - length1 - test)]
		lensLst = lensLst1
		lensLst.extend(lensLst2)
		print "length lensLst", lensLst
		g = parallelComp(lensLst)
		xLst = [i/100.0 for i in range(101)]
		pVals = evalPolys(g,xLst)
		approx = 1.0 - calIntegral(xLst,pVals)
		g.stCollapse()
		real = g.s.edges[0].wt
		approxs.append(approx)
		reals.append(real)
	print "approxs", approxs
	print "reals", reals
	fp = open(fileName, "w")
	for i in xrange(len(approxs)):
		fp.write(str(i) + "\t" +  str(approxs[i]/reals[i]) + "\n")
	fp.close()
	fp2 = open("plot.gplot", "w")
	fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"approx/real l1:" + str(length1) + " l2:" +str(length2)+ "\"")
	fp2.close()
	os.system('gnuplot -persist plot.gplot')

def parallelCompTestRand(minlength, maxlength, numGraphs, fileName="out.dat"):
	"""length1 and length2 are different series graph lengths
	numGraphs is how many graphs we are putting together (divided between lenth1 and length2)"""
	approxs = []
	reals = []
	for test in xrange(0, numGraphs + 1):
		lensLst = [int(random.uniform(minlength, maxlength)) for i in range(test)]

		g = parallelComp(lensLst)
		xLst = [i/100.0 for i in range(101)]
		pVals = evalPolys(g,xLst)
		approx = 1.0 - calIntegral(xLst,pVals)
		g.stCollapse()
		real = g.s.edges[0].wt
		approxs.append(approx)
		reals.append(real)
	print "approxs", approxs
	print "reals", reals
	fp = open(fileName, "w")
	for i in xrange(len(approxs)):
		fp.write(str(i) + "\t" +  str(approxs[i]/reals[i]) + "\n")
	fp.close()
	fp2 = open("plot.gplot", "w")
	fp2.write("plot \"" + str(fileName) + "\" using 1:2 title \"approx/real minl:" + str(minlength) + " maxl:" +str(maxlength)+ "\"")
	fp2.close()
	os.system('gnuplot -persist plot.gplot')
