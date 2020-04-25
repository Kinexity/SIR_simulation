#pragma once
#ifndef POPULATION_H
#define POPULATION_H
#include <map>
#include <random>
#include <memory>
#include <algorithm>
#include <numeric>
#include <execution>
#include <array>
#include <iostream>
#include <set>
#include <fstream>
#include <filesystem>
#include "Individual.h"
#include "Model.h"
#include "Disease.h"
#include "C_Random.h"
#include "C_Time_Counter.h"

class Population {
	friend class Individual;
	friend class Bond;
private:
	class {
	private:
		std::random_device rd;
		std::mt19937 gen{ rd() };
		std::uniform_real_distribution<> uni{ 0.0, std::nextafter(1.0, 2.0) };
	public:
		double operator()() { return uni(gen); }
	} rnd_uni;
	class internal_rnd_norm {
	private:
		std::random_device rd;
		std::mt19937 gen{ rd() };
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
	void
		create_bond(uint_fast64_t member_1_index, uint_fast64_t member_2_index),
		infect();
public:
	Population(Model simulation_model_type_arg, size_t population_size_arg, size_t patient_zero_number_arg, double k_mean, double k_stddev, Disease disease_stats_arg);
	~Population() = default;
	void
		initialize_simulation(),
		build_grid(),
		save_grid(std::string filename),
		load_grid(std::string filename);
	bool
		test_grid();
	std::array<std::vector<int_fast64_t>, 4>
		simulate();
};

#endif // !POPULATION_H