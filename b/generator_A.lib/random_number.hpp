#ifndef RANDOM_NUMBER_HPP
#define RANDOM_NUMBER_HPP
#include <random>
#include <cassert>
#include <cstdint>
#include <limits>
#include <cmath>

class random_number : private std::mt19937_64 {
public:
  using uint64 = uint64_t;
  using uint32 = uint32_t;
  using int64 = int64_t;
  using int32 = int32_t;
  using std::mt19937_64::mt19937_64;
  using std::mt19937_64::operator();

  uint64 next_ulong();                   // [0, 2^{64})
  uint64 next_ulong(uint64 m);           // [0, m)
  uint64 next_ulong(uint64 a, uint64 b); // [a, b]

  uint32 next_uint();                    // [0, 2^{32})
  uint32 next_uint(uint32 m);            // [0, m)
  uint32 next_uint(uint32 a, uint32 b);  // [a, b]

  int64 next_long();                     // [0, 2^{63})
  int64 next_long(int64 m);              // [0, m)
  int64 next_long(int64 a, int64 b);     // [a, b]

  int32 next_int();                      // [0, 2^{31})
  int32 next_int(int32 m);               // [0, m)
  int32 next_int(int32 a, int32 b);      // [a, b]

  double next_double();                  // [0, 1)
  double next_double_by_exponential(double mean);
  double next_double_by_power_law(double eps, double median);
  
  template <typename random_access_iterator>
  void shuffle(random_access_iterator l, random_access_iterator r); // [l, r)
};

using uint64 = random_number::uint64;
using uint32 = random_number::uint32;
using int64 = random_number::int64;
using int32 = random_number::int32;

uint64 random_number::next_ulong() {
  return (*this)();
}
uint64 random_number::next_ulong(uint64 m) {
  assert(m > 0);
  return next_ulong() % m;
}
uint64 random_number::next_ulong(uint64 a, uint64 b) {
  return next_ulong(b - a + 1) + a;
}

uint32 random_number::next_uint() {
  return static_cast<uint32>(next_ulong() & 0xFFFFFFFFUL);
}
uint32 random_number::next_uint(uint32 m) {
  assert(m > 0);
  return next_uint() % m;
}
uint32 random_number::next_uint(uint32 a, uint32 b) {
  return next_uint(b - a + 1) + a;
}

int64 random_number::next_long() {
  return static_cast<int64>(next_ulong() & 0x7FFFFFFFFFFFFFFFUL);
}
int64 random_number::next_long(int64 m) {
  assert(m > 0);
  return next_long() % m;
}
int64 random_number::next_long(int64 a, int64 b) {
  return next_long(b - a + 1) + a;
}

int32 random_number::next_int() {
  return static_cast<int32>(next_ulong() & 0x7FFFFFFFUL);
}
int32 random_number::next_int(int32 m) {
  assert(m > 0);
  return next_int() % m;
}
int32 random_number::next_int(int32 a, int32 b) {
  return next_int(b - a + 1) + a;
}

double random_number::next_double() {
  int32 p = next_int();
  return 1.0 * p / std::numeric_limits<int32>::max();
}

double random_number::next_double_by_exponential(double mean) {
  while(true) {
    double r = next_double();
    if(r < 1.0) return -1.0 * mean * log(1.0 - r);
  }
  assert(false);
  return -1;
}

double random_number::next_double_by_power_law(double eps, double median) {
  while(true) {
    double r = next_double();
    if(r < 1.0) {
      double mu = median / (pow(2.0, 1.0 / eps) - 1);
      return mu * (pow(1.0 - r, -1.0 / eps) - 1);
    }
  }
  assert(false);
  return -1;
}

template <typename random_access_iterator>
void random_number::shuffle(random_access_iterator l, random_access_iterator r) {
  for(random_access_iterator itr=l; itr!=r; itr++) {
    std::iter_swap(l + next_uint(itr - l + 1), itr);
  }
}

#endif
