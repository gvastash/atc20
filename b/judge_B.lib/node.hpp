#ifndef NODE_HPP
#define NODE_HPP
#include <cstdio>

class node {
public:
  size_t node_id;
  int max_cap;
  int cap;
  int type; // 0: normal grid, [1, MAX_TYPE]: nano_grid
  int pw_actual,pw_excess,pw_buy;


  node() = default;
  node(int max_cap, size_t cap, int type) : 
    max_cap(max_cap), cap(cap), type(type) {}
  node(int max_cap, size_t cap, int type, size_t node_id) : 
	  node_id(node_id),max_cap(max_cap),cap(cap),type(type) {
		  pw_actual = 0;
		  pw_excess = 0;
		  pw_buy = 0;
	  }


  bool is_grid(){
	  return type != 0;
  }

  bool charge(int amount) {
    if( !is_grid() ) return false;
    if (amount < 0) return false;
    cap = std::min(max_cap, cap + amount);
    return true;
  }

  bool consume(int amount) {
    if( !is_grid() ) return false;
    if (amount < 0) return false;
    if ( amount > cap) return false;
    cap -= amount;
    return true;
  }

  void init_pw(){
	  pw_actual = 0;
	  pw_excess = 0;
	  pw_buy = 0;
  }

  void print(FILE *fp) {
    fprintf(fp, "%zu %d %d %d %d\n",  node_id+1, cap,pw_actual, pw_excess, pw_buy);
  }
};
#endif
