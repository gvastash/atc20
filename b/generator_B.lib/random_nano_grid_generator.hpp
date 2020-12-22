#ifndef RANDOM_NANO_GRID_GENERATOR_HPP
#define RANDOM_NANO_GRID_GENERATOR_HPP

#include <cstdio>
#include <vector>
#include "generator_B.lib/constraints.hpp"
#include "generator_B.lib/random_number.hpp"

std::pair<std::vector<size_t>, std::vector<size_t>> generate_random_nano_grid(size_t N, random_number &rnd, FILE *fp) {

  assert(N_GRID  <= N);
  std::vector<size_t> x(N_GRID);
  std::vector<size_t> pattern(N_GRID);
  std::vector<size_t> nanos(N_GRID);
  std::vector<size_t> patterns(N_GRID);

  std::vector<size_t> permutation(N);
  std::iota(permutation.begin(), permutation.end(), 1);
  rnd.shuffle(permutation.begin(), permutation.end());
  for (size_t i = 0; i < N_GRID; i++) {
    x[i] = permutation[i];
    pattern[i] = rnd.next_uint(1, N_PATTERN);
  }

  std::sort(x.begin(), x.end());

  for (size_t i = 0; i < N_GRID; i++) {
    nanos[i] = x[i];
    patterns[i] = pattern[i];
    fprintf(fp, "%zu %zu\n", x[i], pattern[i]);
  }

  return std::make_pair(nanos, patterns);
}
#endif
