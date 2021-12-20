#pragma once
#ifndef POPULATION_H
#define POPULATION_H
#include <map>
#include <random>
#include <memory>
#include <algorithm>
#include <numeric>
#include <utility>
#include <execution>
#include <array>
#include <iostream>
#include <set>
#include <fstream>
#include <string>
#include <filesystem>
#include "Individual.h"
#include "Model.h"
#include "Disease.h"
#include "C_Time_Counter.h"
#include "xoshiro256pp.h"

class Population {
	friend class Individual;
	friend class Bond;
private:
	class {
	private:
		xoshiro256pp gen;
	public:
		double operator()() { return double(gen())/std::numeric_limits<uint64_t>::max(); }
	} rnd_uni;
	class internal_rnd_norm {
	private:
		xoshiro256pp gen;
		std::normal_distribution<double> norm;
	public:
		internal_rnd_norm(double k_mean, double k_stddev);
		size_t operator()() { return std::max(static_cast<int_fast64_t>(std::round(norm(gen))), int_fast64_t(1)); }
	} rnd_norm;
	std::vector<std::shared_ptr<Individual>>
		members;
	const Disease
		disease_stats;
	const size_t
		population_size;
	const size_t
		start_infected_number;
	const Model
		simulation_model_type;
	std::vector<size_t>
		bondable_members;
	const std::string
		grid_file_name;
	struct {
		std::atomic<int_fast64_t>
			suspectible = 0,
			exposed = 0,
			infected = 0,
			recovered = 0,
			dead = 0;
	} simulation_stats;
	void
		create_bond(uint_fast64_t member_1_index, uint_fast64_t member_2_index),
		infect();
	PCL::C_Time_Counter
		tc;
public:
	Population(Model simulation_model_type_arg, size_t population_size_arg, size_t patient_zero_number_arg, double k_mean, double k_stddev, Disease_ext disease_stats_arg);
	~Population() = default;
	void
		initialize_simulation(),
		build_grid(),
		save_grid(std::string filename),
		load_grid(std::string filename);
	std::array<std::vector<int_fast64_t>, status_count>
		simulate(std::stringstream& ss, std::stringstream& ss2);
};

#endif // !POPULATION_H