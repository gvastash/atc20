#ifndef RANDOM_ORDER_GENERATOR_HPP
#define RANDOM_ORDER_GENERATOR_HPP
#include <cstdio>
#include "generator_B.lib/random_number.hpp"
#include "generator_B.lib/constraints.hpp"

void generate_random_order(size_t N, random_number &rnd, FILE *fp) {

	size_t cur_order_id = 0;
	for (size_t t = 0; t <= T_LAST; t++) {
		size_t NOrder = 0;
		if (t <= T_LAST) {
			int x = rnd.next_int(0, 1000);
			if (x < P_TRANS_CONST*1000 ) {
				NOrder = 1;
			}
		}
		fprintf(fp, "%zu\n", NOrder);
		for (size_t i = 0; i < NOrder; i++) {
			size_t from = 0;
			size_t to = 0;
			while (from == to) {
				from = rnd.next_uint(0, N - 1);
				to = rnd.next_uint(0, N - 1);
			}
			fprintf(fp, "%zu %zu %zu %zu\n", cur_order_id+1, from+1, to+1, t);
			cur_order_id++;
		}
	}

	return;
	}
#endif
