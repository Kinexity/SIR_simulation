#pragma once
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <vector>
#include <atomic>
#include <map>
#include <mutex>
#include <intrin.h>
#include <iterator>
#include "Population.h"

enum class Status {
	Susceptible,
	Infected,
	Recovered,
	Dead,
	Exposed
};

inline std::array<int_fast64_t, 4> to_array(Status stat) {
	switch (stat) {
	case Status::Susceptible:
		return { 1, 0, 0, 0 };
	case Status::Infected:
		return { 0, 1, 0, 0 };
	case Status::Recovered:
		return { 0, 0, 1, 0 };
	case Status::Dead:
		return { 0, 0, 0, 1 };
	default:
		return { 0, 0, 0, 0 };
	}
}

class Individual {
	friend class Population;
private:
	Population&
		population;
	const size_t
		index;
	std::vector<std::reference_wrapper<Individual>>
		bonded_members;
	Status
		infection_status = Status::Susceptible;
	bool
		positive = false;
	std::mutex
		mut_ex;
	size_t
		k,
		members_infected = 0,
		days_infected = 0;
public:
	Individual(Population& population_ref_arg, uint_fast64_t index_arg);
	~Individual() = default;
	void
		create_bond(uint_fast64_t member_index),
		create_bonds(),
		infect(),
		update_status(),
		clear_bonds();
	bool
		set_positive(),
		try_infect();
	Status
		get_status();
};

#endif // !INDIVIDUAL_H