#ifndef OUTPUT_EXPECTED_STATE_OF_POWER
#define OUTPUT_EXPECTED_STATE_OF_POWER

#include "generator_A.lib/random_number.hpp"
#include "generator_A.lib/constraints.hpp"

#include <cinttypes>
#include <vector>
#include <string>


// pw_type, DAY_TYPE, PATTARN, pw_1, pd_2, ..., pw_{N_DIV}
int PW[2][2][3][20] = {
  {
    {
      {
        -247,-209,-163,-57,5,119,146,200,208,204,190,131,69,22,-95,-139,-173,-165,-97,-68
      },
      {
        8,43,88,129,166,93,82,49,56,53,81,60,32,-2,-42,-166,-209,-197,-174,-150
      },
      {
        -156,-116,-115,-89,-52,-26,-1,24,48,60,50,20,-43,-78,-115,-114,-157,-140,-123,-118
      },
    },
    {
      {
        27,36,47,40,73,82,51,52,56,54,10,7,-32,-52,-137,-146,-145,-89,-61,-32
      },
      {
        -38,-31,-22,-14,-7,16,-13,-58,-57,-57,-60,-64,-29,-36,-44,-53,-61,-41,-18,6
      },
      {
        -122,-110,-105,-112,-64,-63,-56,-44,-64,-50,30,16,-25,-32,-37,-41,-50,-24,-47,-42
      },
    }
  },
  {
    {
      {
        -207,-169,-123,-17,45,159,186,240,248,244,230,171,109,62,-55,-99,-133,-125,-57,-28
      },
      {
        48,83,128,169,206,133,122,89,96,93,121,100,72,38,-2,-126,-169,-157,-134,-110
      },
      {
        -116,-76,-75,-49,-12,14,39,64,88,100,90,60,-3,-38,-75,-74,-117,-100,-83,-78
      },
    },
    {
      {
        67,76,87,80,113,122,91,92,96,94,50,47,8,-12,-97,-106,-105,-49,-21,8
      },
      {
        2,9,18,26,33,56,27,-18,-17,-17,-20,-24,11,4,-4,-13,-21,-1,22,46
      },
      {
        -82,-70,-65,-72,-24,-23,-16,-4,-24,-10,70,56,15,8,3,-1,-10,16,-7,-2
      },
    }
  }

};


void output_expected_state_of_power(size_t DAY_TYPE, size_t pw_type, FILE* fp) {

  std::string P_EVENT = (DAY_TYPE == 1 or DAY_TYPE == 3 ? "0.1" : "0.0");
  fprintf(fp, "%zu %zu %zu %s %zu\n", N_DIV, N_PATTERN, DISPERSION, P_EVENT.c_str(), DELTA_EVENT);
  for (size_t i = 0; i < N_PATTERN; i++) {
    for (size_t j = 0; j < N_DIV; j++) {
      int wether = (DAY_TYPE <= 1 ? 0 : 1);
      int pw = PW[pw_type][wether][i][j];
      if (j + 1 == N_DIV) {
        fprintf(fp, "%d\n", pw);
      }
      else {
        fprintf(fp, "%d ", pw);
      }
    }
  }

  fflush(fp);
}

#endif
