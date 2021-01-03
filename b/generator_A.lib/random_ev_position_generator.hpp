#ifndef RANDOM_EV_POSITION_GENERATOR_HPP
#define RANDOM_EV_POSITION_GENERATOR_HPP

#include <cstdio>
#include <vector>
#include "generator_A.lib/constraints.hpp"
#include "generator_A.lib/random_number.hpp"

void generate_random_ev_position(size_t N, size_t N_EV, random_number &rnd, FILE *fp) {

  assert(N_EV <= N);
  std::vector<size_t> x(N_EV);

  std::vector<size_t> permutation(N);
  std::iota(permutation.begin(), permutation.end(), 1);
  rnd.shuffle(permutation.begin(), permutation.end());
  for (size_t i = 0; i < N_EV; i++) {
    x[i] = permutation[i];
  }
  
  std::sort(x.begin(), x.end());

  for (size_t i = 0; i < N_EV; i++) {
    fprintf(fp, "%zu\n", x[i]);
  }

  return;
}
#endif
