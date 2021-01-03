#ifndef OUTPUT_FOR_GRAPH_GRAPH
#define OUTPUT_FOR_GRAPH_GRAPH

#include <vector>
#include <cinttypes>

#include "generator_A.lib/graph.hpp"
#include "generator_A.lib/random_number.hpp"
#include "generator_A.lib/constraints.hpp"

void output_for_graphgraph(const graph<edge> &G, std::vector<double> x, std::vector<double> y, FILE *fp) {
	
	size_t N = G.size();
	assert(N == x.size());
	assert(N == y.size());
	
	random_number rnd(283);
  std::vector<edge> edges;
	std::vector<size_t> nodes(N);
	std::vector<size_t> inv(N);
	for(size_t i = 0; i < N; i++) nodes[i] = i;
  for(size_t i=0; i<N; i++) {
    for(auto e : G[i]) {
      size_t u, v;
      u = e.src, v = e.dst;
      if(u > v) continue;
      edges.emplace_back(e);
    }
  }

  size_t M = edges.size();
  fprintf(fp, "%zu %zu\n", N, M);

	// print points
	for(size_t i = 0; i < N; i++) {
		inv[nodes[i]] = i;
	}
	for(size_t i = 0; i < N; i++) {
		fprintf(fp, "%f %f\n", x[inv[i]] * 100, y[inv[i]] * 100);
	}

	// print edges
  for(auto e : edges) {
    size_t u, v; edge::cost_type c;
    u = nodes[e.src] + INDEXED, v = nodes[e.dst] + INDEXED, c = e.cost;
		if((1 << 1) & e.bit or (1 << 2) & e.bit or (1 << 3) & e.bit or (1 << 4) & e.bit) {
			fprintf(fp, "%zu %zu %" PRId64 " 10\n", u, v, c);
		} else {
			fprintf(fp, "%zu %zu %" PRId64 " 1\n", u, v, c);
		}
  }

  fflush(fp);
}
#endif
