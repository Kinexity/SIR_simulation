#include "Population.h"
#include <utility>
#include <format>

void Population::create_bond(uint_fast64_t member_1_index, uint_fast64_t member_2_index) {
	members[member_1_index]->create_bond(member_2_index);
	members[member_2_index]->create_bond(member_1_index);
}

void Population::infect() {
	std::vector<std::shared_ptr<Individual>> first_infected;
	std::sample(members.begin(), members.end(), std::back_inserter(first_infected), start_infected_number, std::mt19937{ std::random_device{}() });
	simulation_stats.infected = start_infected_number;
	simulation_stats.suspectible = members.size() - start_infected_number;
	std::for_each(std::execution::par_unseq, first_infected.begin(), first_infected.end(), [&](std::shared_ptr<Individual> member) { member->initial_infect(); });
}

Population::Population(Model simulation_model_type_arg, size_t population_size_arg, size_t start_infected_number_arg, double k_mean_arg, double k_stddev_arg, Disease_ext disease_stats_arg) :
	simulation_model_type(simulation_model_type_arg), population_size(population_size_arg), start_infected_number(start_infected_number_arg), disease_stats(dis_int_conv(disease_stats_arg)), rnd_norm(k_mean_arg, k_stddev_arg),
	grid_file_name(std::string("grid_") + std::format("{:x}", std::hash<std::string>{}(std::to_string(population_size) + "_" + std::to_string(k_mean_arg) + "_" + std::to_string(k_stddev_arg))) + ".bin") {}

void Population::initialize_simulation() {
	if (!std::filesystem::exists(grid_file_name)) {
		tc.start();
		members.clear();
		for (size_t index = 0; index < population_size; index++) {
			members.push_back(std::make_shared<Individual>(*this, index));
		}
		infect();
		tc.stop();
		std::cout << "Members initialization time: " << tc.measured_timespan().count() << '\n';
		tc.start();
		build_grid();
		tc.stop();
		std::cout << "Bonding time: " << tc.measured_timespan().count() << '\n';
		save_grid(grid_file_name);
	}
	else {
		tc.start();
		members.clear();
		load_grid(grid_file_name);
		infect();
		tc.stop();
		std::cout << "Loading time: " << tc.measured_timespan().count() << '\n';
	}
}

void Population::build_grid() {
	bondable_members.resize(population_size);
	std::iota(bondable_members.begin(), bondable_members.end(), 0);
	size_t sum = 0;
	PCL::C_Time_Counter tc_loc_bonds;
	tc_loc_bonds.start();
	for (auto& member : members) {
		member->create_bonds();
		sum += member->bonded_members.size();
		if ((std::distance(&*members.begin(), &member) + 1) % 10000 == 0) {
			tc_loc_bonds.stop();
			std::cout << "Bonds finished: " << std::distance(&*members.begin(), &member) + 1 << "	Time: " << tc_loc_bonds.measured_timespan().count() << '\n';
			tc_loc_bonds.start();
		}
	}
	tc_loc_bonds.stop();
	bondable_members.clear();
}

void Population::load_grid(std::string filename) {
	std::vector<size_t> grid_data;
	std::fstream file;
	grid_data.resize(std::filesystem::file_size(filename) / sizeof(size_t));
	file.open(filename, std::ios::binary | std::ios::in);
	file.read((char*)grid_data.data(), grid_data.size() * sizeof(size_t));
	file.close();
	members.resize(grid_data[0]);
	std::vector<size_t> members_indices;
	members_indices.resize(members.size());
	std::iota(members_indices.begin(), members_indices.end(), 0);
	std::vector<std::vector<size_t>> bonded_members;
	for (size_t index = 1; index < grid_data.size(); index += (grid_data[index] + 1)) {
		bonded_members.push_back(std::vector<size_t>{grid_data.begin() + index + 1, grid_data.begin() + index + grid_data[index] + 1});
	}
	std::for_each(std::execution::par_unseq, members_indices.begin(), members_indices.end(), [&](auto member_index) {
		members[member_index] = std::make_shared<Individual>(*this, bonded_members[member_index].size());
		});
	std::for_each(std::execution::par_unseq, members_indices.begin(), members_indices.end(), [&](auto member_index) {
		for (auto bond_member_index : bonded_members[member_index]) {
			members[member_index]->recreate_bond(bond_member_index);
		}
		});
	std::cout << "Grid loaded from " << filename << '\n';
}

void Population::save_grid(std::string filename) {
	std::vector<size_t> grid_data;
	grid_data.push_back(members.size());
	for (auto member : members) {
		grid_data.push_back(member->bonded_members.size());
		for (auto& bond : member->bonded_members) {
			grid_data.push_back(bond.get().index);
		}
	}
	std::fstream file;
	file.open(filename, std::ios::binary | std::ios::trunc | std::ios::out);
	file.write((char*)grid_data.data(), grid_data.size() * sizeof(size_t));
	file.close();
	std::cout << "Grid saved to " << filename << '\n';
}

std::array<std::vector<int_fast64_t>, status_count> Population::simulate(std::stringstream& ss) {
	tc.start();
	ss.clear();
	std::array<std::vector<int_fast64_t>, status_count> result_stats;
	size_t sim_time = 0;
	auto update_arr = [&] {
		ss << sim_time;
		std::cout << sim_time++;
		result_stats[0].push_back(simulation_stats.suspectible);
		result_stats[1].push_back(simulation_stats.exposed);
		result_stats[2].push_back(simulation_stats.infected);
		result_stats[3].push_back(simulation_stats.recovered);
		result_stats[4].push_back(simulation_stats.dead);
		for (auto& vec : result_stats) {
			ss << "	" << vec.back();
			std::cout << "	" << vec.back();
		}
		ss << '\n';
		std::cout << '\n';
	};
	update_arr();
	auto last = members.end();
	do {
		if (result_stats[2].back() != 0) {
			std::for_each(std::execution::par_unseq, members.begin(), last, [&](std::shared_ptr<Individual> elem) { elem->try_to_get_infected(); });
		}
		std::for_each(std::execution::par_unseq, members.begin(), last, [&](std::shared_ptr<Individual> elem) { elem->update_status(); });
		update_arr();
		//last = std::remove_if(std::execution::par_unseq, members.begin(), last, [](std::shared_ptr<Individual>& i) {return i->get_status() == Status::Recovered; });
	} while (!(simulation_stats.exposed == 0 && simulation_stats.infected == 0 || simulation_stats.suspectible == 0));
	tc.stop();
	std::cout << "Simulation time: " << tc.measured_timespan().count() << '\n';
	return result_stats;
}

Population::internal_rnd_norm::internal_rnd_norm(double k_mean, double k_stddev) : norm(k_mean, k_stddev) {}
