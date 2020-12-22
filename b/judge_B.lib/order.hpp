#ifndef ORDER_HPP
#define ORDER_HPP

class order {
public:
  size_t order_id;//0-index
  size_t start_time;
  int  end_time;
  int from;//0-index
  int to;//0-index
  int vehicle_id; 

  order() = default;
  order(size_t order_id, size_t start_time, int from, int to) : 
    order_id(order_id), start_time(start_time), end_time(-1), from(from), to(to), vehicle_id(-1) {}


  void print(FILE *fp) {
    fprintf(fp, "%zu %d %d %d %zu\n", order_id+1, from+1, to+1, vehicle_id == -1 ? 0 : 1 , start_time) ;
  }

};
#endif
