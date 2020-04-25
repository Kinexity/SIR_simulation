#include "Individual.h"
#include "Population.h"

Individual::Individual(Population& population_ref_arg, uint_fast64_t index_arg) :
	population(population_ref_arg), index(index_arg) {
	k = population.rnd_norm();
}

void Individual::create_bond(uint_fast64_t member_index) {
	if (member_index != index) {
		bonded_members.push_back(*population.members[member_index]);
	}
	if (bonded_members.size() == k) {
		static std::mutex mtx;
		std::unique_lock<std::mutex> lck(mtx);
		std::erase(population.bondable_members, index);
	}
}

void Individual::create_bonds() {
	static std::mt19937 rnd_local{ std::random_device{}() };
	std::erase(population.bondable_members, index);
	std::vector<size_t> chosen_members;
	chosen_members.resize((k - bonded_members.size() % (population.bondable_members.size() + 1)));
	std::sample(
		population.bondable_members.begin(),
		population.bondable_members.end(),
		chosen_members.begin(),
		(k - bonded_members.size() % (population.bondable_members.size() + 1)),
		rnd_local);
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
		auto drawn = population.rnd_uni();
		if (infection_status == Status::Infected) {
			if ((population.simulation_model_type == Model::SIRD || population.simulation_model_type == Model::SEIRD) && drawn < population.disease_stats.death_probability) {
				infection_status = Status::Dead;
				population.simulation_stats.infected--;
				population.simulation_stats.dead++;
				bonded_members.clear();
			}
			else if (drawn < population.disease_stats.death_probability + population.disease_stats.recovery_probability) {
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
		}
		else {
			if (drawn < population.disease_stats.exposition_to_infection_probability) {
				infection_status = Status::Infected;
				population.simulation_stats.exposed--;
				population.simulation_stats.infected++;
			}
		}
	}
	else if (infection_status == Status::Susceptible && positive) {
		population.simulation_stats.suspectible--;
		if (population.simulation_model_type != Model::SEIR && population.simulation_model_type != Model::SEIRD) {
			infection_status = Status::Infected;
			population.simulation_stats.infected++;
		}
		else {
			infection_status = Status::Exposed;
			population.simulation_stats.exposed++;
		}
	}
}

bool Individual::set_positive() {
	return !std::exchange((bool&)positive, true);
}

bool Individual::try_infect() {
	return get_status() == Status::Susceptible && population.rnd_uni() <  population.disease_stats.infection_or_exposition_probability && set_positive();
}

Status Individual::get_status() {
	return infection_status;
}

void Individual::infect() {
	if (infection_status == Status::Infected) {
		members_infected += std::transform_reduce(std::execution::par_unseq, bonded_members.begin(), bonded_members.end(), size_t(), std::plus<size_t>(), [&](Individual& elem)->size_t { return (elem.try_infect() ? 1 : 0);  });
		auto is_SS = (population.simulation_model_type == Model::SIS);
		std::erase_if(bonded_members, [&](Individual& elem) {
			Status st = elem.get_status();
			return (is_SS ? st == Status::Dead : st != Status::Susceptible);
		});
	}
}
