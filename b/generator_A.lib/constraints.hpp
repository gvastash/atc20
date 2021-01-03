#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include <cstdint>
#include <cstddef>
#include <climits>

const size_t N_REST = 5;

const size_t MIN_N = 225;
const size_t MAX_N = 225;

const size_t MIN_DAY_TYPE = 0;
const size_t MAX_DAY_TYPE = 3;
const size_t N_DIV = 20;
const size_t N_PATTERN = 3;
const size_t DISPERSION = 100;
const size_t DELTA_EVENT = 1000;
const size_t N_EVENT = 3;

const int MIN_PW = -999;
const int MAX_PW = 999;

const size_t N_GRID = 20;
const size_t C_INIT_GRID = 25000;
const size_t C_MAX_GRID = 50000;
const size_t V_MAX_GRID = 800;

const size_t MIN_EV = 15;
const size_t MAX_EV = 25;
const size_t C_INIT_EV = 12500;
const size_t C_MAX_EV = 25000;
const size_t V_MAX_EV = 400;
const size_t N_MAX_TRANS = 4;
const size_t DELTA_MOVE_EV = 50;

const double P_TRANS_CONST = 0.7;
const size_t T_LAST = 900;
const size_t P_TRANS = 3000;
const double GAMMA = 2;
const int S_ELE_REF = -1500000;
const int S_TRANS_REF = -1900000;
const size_t T_MAX = 1000;

const size_t MIN_C = 1;
const size_t MAX_C = 100;

const size_t INDEXED = 1;

const size_t MAX_DEGREE = 5;

const size_t MIN_N_ORDER = 0;
const size_t MAX_N_ORDER = 1;

const size_t SCALE_ROAD = 2;

const double LIMIT_TIME_PER_ONECE = 100.0;

#endif
