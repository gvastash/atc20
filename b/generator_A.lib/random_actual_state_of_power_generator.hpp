#ifndef OUTPUT_ACTUAL_STATE_OF_POWER
#define OUTPUT_ACTUAL_STATE_OF_POWER

#include "generator_A.lib/random_number.hpp"
#include "generator_A.lib/constraints.hpp"

#include <cinttypes>
#include <vector>
#include <random>

void output_actual_state_of_power(random_number &rnd, size_t DAY_TYPE, size_t pw_type, std::vector<size_t> &nanos, std::vector<size_t> &patterns, bool remove_randomness, FILE *fp) {


  assert(T_MAX % N_DIV == 0);
  size_t len = T_MAX / N_DIV; 

  std::normal_distribution<> dist(0, std::sqrt(DISPERSION));
  std::default_random_engine engine(rnd.next_uint());

  std::uniform_real_distribution<double> unit_random(0, 1);
  const double p_event = DAY_TYPE == 1 || DAY_TYPE == 3 ? 0.1 : 0;
  std::vector<std::vector<bool>> event_pos(N_DIV, std::vector<bool>(N_GRID, false)); 
  std::vector<size_t> selector(N_GRID); std::iota(selector.begin(), selector.end(), 0);
  for(size_t k = 0; k < N_DIV; ++k){
    if (unit_random(engine) < p_event){
      std::shuffle(selector.begin(), selector.end(), engine);
      for(size_t i = 0; i < N_EVENT; ++i)
        event_pos[k][selector[i]] = true;
    }
  }

  for (size_t i = 0; i < T_MAX; i++) {
    size_t kukaku = i / len;
    assert(kukaku < N_DIV);

    for (size_t j = 0; j < N_GRID; j++) {
      size_t weather = (DAY_TYPE <= 1 ? 0 : 1);
      size_t pattern = patterns[j] - 1;
      size_t id = nanos[j];

      int pw = PW[pw_type][weather][pattern][kukaku]; 
      if(!remove_randomness){
        if(event_pos[kukaku][j])
          pw += (DAY_TYPE == 1 ? -1 : DAY_TYPE == 3 ? 1 : 0)*DELTA_EVENT;

        int add = round(dist(engine));
        pw += add;
      }
      fprintf(fp, "%zu %d", id, pw); 
      if (j + 1 == N_GRID) {
        fprintf(fp, "\n");
      } else {
        fprintf(fp, " ");
      }
    }
  }
  
  fflush(fp);
}

#endif
