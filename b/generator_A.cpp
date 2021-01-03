#include <cstdio>
#include <random>
#include "generator_A.lib/random_graph_generator.hpp"
#include "generator_A.lib/random_nano_grid_generator.hpp"
#include "generator_A.lib/random_ev_position_generator.hpp"
#include "generator_A.lib/random_expected_state_of_power_generator.hpp"
#include "generator_A.lib/random_actual_state_of_power_generator.hpp"
#include "generator_A.lib/output_undirected_graph.hpp"
#include "generator_A.lib/graph_relabeling.hpp"
#include "generator_A.lib/random_order_generator.hpp"
#include "generator_A.lib/cmdline.h"

uint_fast64_t seed_1 = 23587163955511535UL;

int main(int argc, const char** argv) {
  int mode = 1;

  cmdline::parser options;
  options.add<uint_fast64_t>("seed", 's', "Seed of pseudo-random sequence.", false, seed_1);
  options.add<size_t>("day-type", 'd', "Fix day type. The value must be 0, 1, 2 or 3.", false);
  options.add("random", 'r', "Obtain seed from random device.");
  options.add("deterministic", 'D', "Remove randomness from supply and demand for  electric power");
  if (mode == 0)
    options.add<int>("pw-type", 'p', "select pw type. The value must be 0 (A) or 1 (B).", true);
  options.parse_check(argc, argv);

  if (options.exist("random"))
    seed_1 = std::random_device()();
  else
    seed_1 = options.get<uint_least64_t>("seed");

  int fixed_day_type = options.exist("day-type") ? options.get<size_t>("day-type") : -1;
  assert(-1 <= fixed_day_type and fixed_day_type <= 3);
  const bool remove_randomness = options.exist("deterministic");
  size_t pw_type = mode == 0 ? options.get<int>("pw-type") : mode - 1;
  if (mode == 0) mode = pw_type + 1;
  assert(pw_type == 0 or pw_type == 1);

  fprintf(stdout, "%zu\n", mode == 1 ? 1 : N_REST);

  random_number rnd(seed_1);
  const size_t N = rnd.next_uint(MIN_N, MAX_N);
  const size_t M = rnd.next_uint(std::ceil(1.5 * N), 2 * N);
  assert((1.5 * N) <= M and M <= 2 * N);

  graph<edge> G = generate_random_graph(N, M, rnd);
  G = graph_relabeling(G, rnd);
  output_undirected_graph(G, rnd, stdout);

  const size_t DAY_TYPE = fixed_day_type == -1 ? rnd.next_uint(MIN_DAY_TYPE, MAX_DAY_TYPE) : fixed_day_type;
  fprintf(stdout, "%zu\n", DAY_TYPE);

  output_expected_state_of_power(DAY_TYPE, pw_type, stdout);

  fprintf(stdout, "%zu %zu %zu %zu\n", N_GRID, C_INIT_GRID, C_MAX_GRID, V_MAX_GRID);

  std::vector<size_t> nanos;
  std::vector<size_t> patterns;
  std::tie(nanos, patterns) = generate_random_nano_grid(N, rnd, stdout);

  const size_t N_EV = rnd.next_uint(MIN_EV, MAX_EV);
  fprintf(stdout, "%zu %zu %zu %zu %zu %zu\n", N_EV, C_INIT_EV, C_MAX_EV, V_MAX_EV, N_MAX_TRANS, DELTA_MOVE_EV);

  generate_random_ev_position(N, N_EV, rnd, stdout);

  fprintf(stdout, "%.02f %zu\n", P_TRANS_CONST, T_LAST);
  fprintf(stdout, "%zu %.02f %d %d\n", P_TRANS, GAMMA, S_ELE_REF, S_TRANS_REF);
  fprintf(stdout, "%zu\n", T_MAX);

  for (size_t i = 0; i < (mode == 1 ? 1 : N_REST); i++) {
    generate_random_order(N, rnd, stdout);

    output_actual_state_of_power(rnd, DAY_TYPE, pw_type, nanos, patterns, remove_randomness, stdout);
  }

  fflush(stdout);

  return 0;
}
