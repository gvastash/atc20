#ifndef NODE_HPP
#define NODE_HPP
#include <cstdio>
//#include "../utils/constraints.hpp"

class node {
public:
  size_t node_id;
  size_t max_cap;
  size_t cap;
  int type; // -1: normal grid, [0, MAX_TYPE]: nano_grid
  int pw_actual, pw_excess, pw_buy;


  node() = default;
  node(size_t max_cap, size_t cap, int type) :
    max_cap(max_cap), cap(cap), type(type) {}
  node(size_t max_cap, size_t cap, int type, size_t node_id) :
    node_id(node_id), max_cap(max_cap), cap(cap), type(type) {
    pw_actual = 0;
    pw_excess = 0;
    pw_buy = 0;
  }


  bool is_grid() {
    return type != 0;
  }

  bool charge(int amount) {
    if (!is_grid()) return false;
    if (amount < 0) return false;
    cap = std::min(max_cap, cap + amount);
    return true;
  }

  bool consume(int amount) {
    if (!is_grid()) return false;
    if (amount < 0) return false;
    if ((size_t)amount > cap) return false;
    cap -= amount;
    return true;
  }

  void init_pw() {
    pw_actual = 0;
    pw_excess = 0;
    pw_buy = 0;
  }

  void print(FILE* fp) {
    fprintf(fp, "%zu %zu %d %d %d\n", node_id + 1, cap, pw_actual, pw_excess, pw_buy);
  }
};
#endif
