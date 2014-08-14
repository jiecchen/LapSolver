#include "PartialCholeskyOrder.h"

PartialCholeskyOrder::PartialCholeskyOrder(const Graph &g) : graph(g) {
	n = g.nv;
}

PartialCholeskyOrder::~PartialCholeskyOrder() {
}
