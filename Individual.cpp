#include "Individual.h"
#include "Population.h"
#include <set>
#include <string>

template <class T>
double uni_fast(T var) {
	return double(var) / std::numeric_limits<T>::max();
}

Individual::Individual(Population& population_ref_arg, uint_fast64_t index_arg) :
	population(population_ref_arg), index(index_arg) {
	k = std::max(population.rnd_norm(), size_t(1));
}

void Individual::create_bond(uint_fast64_t member_index) {
	if (member_index != index) {
		bonded_members.push_back(*(population.members[member_index]));
	}
	if (bonded_members.size() == k) {
		static std::mutex mtx;
		std::unique_lock<std::mutex> lck(mtx);
		std::erase(population.bondable_members, index);
	}
}

void Individual::decrement_infected_neighbours() {
	infected_neighbours--;
}

void Individual::increment_infected_neighbours() {
	infected_neighbours++;
}

void Individual::initial_infect() {
	infection_status = Status::Infected;
	auto drawn = (double)rnd() / std::numeric_limits<uint64_t>::max();
	std::for_each(std::execution::par_unseq, bonded_members.begin(), bonded_members.end(), [&](Individual& elem) { elem.increment_infected_neighbours();  });
	days_to_state_change = std::ceil(std::log(drawn) / population.disease_stats.recovery_probability_log);
}

void Individual::recreate_bond(uint_fast64_t member_index) {
	if (member_index != index) {
		bonded_members.push_back(*population.members[member_index]);
	}
}

void sample_fast(std::vector<size_t>& arr, std::set<size_t>& smpl_arr, xoshiro256pp& gen, size_t sample_size) {
	if (sample_size < arr.size()) {
		do {
			smpl_arr.insert(gen() % arr.size());
		} while (smpl_arr.size() < sample_size);
	}
	else {
		for (auto x : arr) {
			smpl_arr.insert(x);
		}
	}
}

bool smpl_std = true;

void Individual::create_bonds() {
	std::erase(population.bondable_members, index);
	std::set<size_t> chosen_members;
	//chosen_members.resize((k - bonded_members.size() % (population.bondable_members.size() + 1)));
	//if (smpl_std) {
	//	std::sample(
	//		population.bondable_members.begin(),
	//		population.bondable_members.end(),
	//		chosen_members.begin(),
	//		(k - bonded_members.size() % (population.bondable_members.size() + 1)),
	//		rnd_local);
	//}
	//else {
	sample_fast(population.bondable_members, chosen_members, rnd, (k - bonded_members.size() % (population.bondable_members.size() + 1)));
	//}
	//std::set<size_t> chosen_members;
	//while (chosen_members.size() < (k - bonded_members.size()) && (chosen_members.size() < population.population_size - index - 1)) {
	//	chosen_members.insert(index + 1 + rnd_local() % (population.population_size - index - 1));
	//}
	for (auto second_member_index : chosen_members) {
		population.create_bond(index, second_member_index);
	}
}

void Individual::clear_bonds() {
	bonded_members.clear();
}

void Individual::update_status() {
	if (infection_status == Status::Infected || infection_status == Status::Exposed) {
		if (--days_to_state_change == 0) {
			if (infection_status == Status::Infected) {
				if ((population.simulation_model_type == Model::SIRD || population.simulation_model_type == Model::SEIRD)) {
					infection_status = Status::Dead;
					population.simulation_stats.infected--;
					population.simulation_stats.dead++;
					std::for_each(std::execution::par_unseq, bonded_members.begin(), bonded_members.end(), [&](Individual& elem) { elem.decrement_infected_neighbours();  });
					bonded_members.clear();
				}
				else if (true) {
					if (population.simulation_model_type == Model::SIS) {
						infection_status = Status::Susceptible;
						population.simulation_stats.infected--;
						population.simulation_stats.suspectible++;
					}
					else {
						infection_status = Status::Recovered;
						population.simulation_stats.infected--;
						population.simulation_stats.recovered++;
						bonded_members.clear();
					}
					positive = false;
				}
				std::for_each(std::execution::par_unseq, bonded_members.begin(), bonded_members.end(), [&](Individual& elem) { elem.decrement_infected_neighbours();  });
			}
			else {
				if (true) {
					infection_status = Status::Infected;
					std::for_each(std::execution::par_unseq, bonded_members.begin(), bonded_members.end(), [&](Individual& elem) { elem.increment_infected_neighbours();  });
					population.simulation_stats.exposed--;
					population.simulation_stats.infected++;
				}
			}
		}
	}
	else if (infection_status == Status::Susceptible && positive) {
		population.simulation_stats.suspectible--;
		if (population.simulation_model_type != Model::SEIR && population.simulation_model_type != Model::SEIRD) {
			infection_status = Status::Infected;
			population.simulation_stats.infected++;
			std::for_each(std::execution::par_unseq, bonded_members.begin(), bonded_members.end(), [&](Individual& elem) { elem.increment_infected_neighbours();  });
		}
		else {
			infection_status = Status::Exposed;
			population.simulation_stats.exposed++;
		}
		auto drawn = uni_fast(rnd());
		days_to_state_change = std::ceil(std::log(drawn) / population.disease_stats.recovery_probability_log);
		//std::cout << (std::to_string(days_to_state_change) + "\n");
	}
}

Status Individual::get_status() {
	return infection_status;
}

void Individual::try_to_get_infected() {
	if (infection_status == Status::Susceptible && infected_neighbours > 0) {
		auto prob = uni_fast(rnd());
		positive = (prob < 1 - std::pow(1 - population.disease_stats.infection_or_exposition_probability, infected_neighbours));
		auto is_SS = (population.simulation_model_type == Model::SIS);
		std::erase_if(bonded_members, [&](Individual& elem) {
			Status st = elem.get_status();
			return (is_SS ? st == Status::Dead : st != Status::Susceptible);
			});
	}
}
