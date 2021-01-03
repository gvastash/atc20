#ifndef VEHICLE_HPP
#define VEHICLE_HPP


#include <set>
#include "judge_A.lib/position.hpp"

class vehicle {
public:
  int max_cap;
  int cap;
  position pos;
  std::set<size_t> orders; 



  vehicle() = default;
  vehicle(int max_cap, int cap, position pos) : 
    max_cap(max_cap), cap(cap), pos(pos) {}
  
  bool charge(int amount) {
    if (amount < 0) return false;
    cap = std::min(max_cap, cap + amount);
    return true;
  }
  
  bool consume(int amount) {
    if (amount < 0 ) return false;
    if ( amount > cap) return false;
    cap -= amount;
    return true;
  }

  bool set_position(position p) {
    pos = p;
    return true;
  }


  int now_vertex() {
    return pos.now_vertex();
  }


  void add_order(size_t id){// this vehicle pick up about the order/
	  orders.insert(id);
  }
  void erase_order(size_t id ){
	  orders.erase(id);
  }


};
#endif
