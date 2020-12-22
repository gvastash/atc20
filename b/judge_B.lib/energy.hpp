#ifndef ENERGY_HPP
#define ENERGY_HPP

class energy{
public:
	size_t id;//nano_id
	int amount;

	energy() = default;
	energy(size_t id,int amount): id(id),amount(amount){}

};


#endif
