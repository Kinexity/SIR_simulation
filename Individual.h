#pragma once
inline constexpr size_t status_count = 5;
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <vector>
#include <atomic>
#include <map>
#include <mutex>
#include <intrin.h>
#include <iterator>
#include "Population.h"
#include "xoshiro256pp.h"

enum class Status {
	Susceptible,
	Exposed,
	Infected,
	Recovered,
	Dead
};

inline std::array<int_fast64_t, status_count> to_array(Status stat) {
	switch (stat) {
	case Status::Susceptible:
		return { 1, 0, 0, 0 , 0 };
	case Status::Exposed:
		return { 0, 1, 0, 0 , 0 };
	case Status::Infected:
		return { 0, 0, 1, 0 , 0 };
	case Status::Recovered:
		return { 0, 0, 0, 1 , 0 };
	case Status::Dead:
		return { 0, 0, 0, 0 , 1 };
	default:
		return { 0, 0, 0, 0 , 0 };
	}
}

class Individual {
	friend class Population;
private:
	xoshiro256pp
		rnd;
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
	size_t
		k,
		days_to_state_change = 0;
	std::atomic<size_t>
		infected_neighbours = 0;
	size_t
		infected_neighbours_old = 0;
public:
	Individual(Population& population_ref_arg, uint_fast64_t index_arg);
	~Individual() = default;
	void
		create_bond(uint_fast64_t member_index),
		recreate_bond(uint_fast64_t member_index),
		create_bonds(),
		try_to_get_infected(),
		initial_infect(),
		update_status(),
		clear_bonds(),
		increment_infected_neighbours(),
		decrement_infected_neighbours();
	Status
		get_status();
};

#endif // !INDIVIDUAL_H