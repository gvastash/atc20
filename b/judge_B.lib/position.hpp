#ifndef POSITION_HPP
#define POSITION_HPP

#include "judge_B.lib/graph.hpp" 

class position {
public:
	size_t from;
	size_t to;
	size_t distance;
	size_t edge_distance;

	position() = default;
	position(size_t from, size_t to, size_t distance, size_t edge_distance) :
		from(from), to(to), distance(distance), edge_distance(edge_distance) {}

	bool on_vertex() {
		if (distance == 0 || distance == edge_distance) {// 0 is from , edge_distance is to
			return true;
		}
		else {
			return false;
		}
	}

	int now_vertex() {
		if (on_vertex()) {
			if (distance == 0) {
				return from;
			}
			else {
				return to;
			}
		}
		return -1;
	}
	size_t rest_distance() {
		return edge_distance - distance;
	}

};
#endif
