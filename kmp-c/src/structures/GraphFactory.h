#pragma once

#include <string>
#include "structures/Graph.h"
#include "util/aligned.h"

class GraphFactory
{
public:
	
	static Graph fromArrays(aligned_vector<int> u, aligned_vector<int> v, aligned_vector<double> w);
	static Graph fromFile(const std::string &filename);

};