#include "Population.h"

void Population::create_bond(uint_fast64_t member_1_index, uint_fast64_t member_2_index) {
	members[member_1_index]->create_bond(member_2_index);
	members[member_2_index]->create_bond(member_1_index);
}

void Population::infect() {
	std::vector<std::shared_ptr<Individual>> first_infected;
	std::sample(members.begin(), members.end(), std::back_inserter(first_infected), start_infected_number, std::mt19937{ std::random_device{}() });
	std::for_each(std::execution::par_unseq, first_infected.begin(), first_infected.end(), [&](std::shared_ptr<Individual> member) { member->infection_status = Status::Infected; });
}

Population::Population(Model simulation_model_type_arg, size_t population_size_arg, size_t start_infected_number_arg, double k_mean_arg, double k_stddev_arg, Disease disease_stats_arg) :
	simulation_model_type(simulation_model_type_arg), population_size(population_size_arg), start_infected_number(start_infected_number_arg), disease_stats(disease_stats_arg), rnd_norm(k_mean_arg, k_stddev_arg) {}

void Population::initialize_simulation() {
	if (false) {
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
		save_grid("grid.bin");
	}
	else {
		tc.start();
		members.clear();
		load_grid("grid.bin");
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

bool Population::test_grid()
{
	return false;
}

void Population::load_grid(std::string filename) {
	std::vector<size_t> grid_data;
	std::fstream file;
	grid_data.resize(std::filesystem::file_size(filename) / sizeof(size_t));
	file.open(filename, std::ios::binary | std::ios::in);
	file.read((char*)grid_data.data(), grid_data.size() * sizeof(size_t));
	file.close();
	members.resize(grid_data[0]);
	std::vector<std::vector<size_t>> bonded_members;
	for (size_t index = 1; index < grid_data.size(); index += (grid_data[index] + 1)) {
		bonded_members.push_back(std::vector<size_t>{grid_data.begin() + index + 1, grid_data.begin() + index + grid_data[index] + 1});
	}
	for (size_t member_index = 0; member_index < members.size(); member_index++) {
		members[member_index] = std::make_shared<Individual>(*this, bonded_members[member_index].size());
	}
	for (size_t member_index = 0; member_index < members.size(); member_index++) {
		for (auto bond_member_index : bonded_members[member_index]) {
			members[member_index]->create_bond(bond_member_index);
		}
	}
	std::cout << "Grid loaded from " << filename << '\n';
}

void Population::save_grid(std::string filename) {
	std::vector<size_t> grid_data;
	grid_data.push_back(members.size());
	for (auto member : members) {
		grid_data.push_back(member->bonded_members.size());
		for (auto& bond: member->bonded_members) {
			grid_data.push_back(bond.get().index);
		}
	}
	std::fstream file;
	file.open(filename, std::ios::binary | std::ios::trunc | std::ios::out);
	file.write((char*)grid_data.data(), grid_data.size() * sizeof(size_t));
	file.close();
	std::cout << "Grid saved to " << filename << '\n';
}

std::array<std::vector<int_fast64_t>, 4> Population::simulate() {
	tc.start();
	std::array<std::vector<int_fast64_t>, 4> result_stats;
	auto trn_arr = [&](std::shared_ptr<Individual> elem) { return to_array(elem->get_status()); };
	auto red_arr = [&](std::array<int_fast64_t, 4> a, std::array<int_fast64_t, 4> b) {
		std::transform(std::execution::par_unseq, a.begin(), a.end(), b.begin(), a.begin(), std::plus<>());
		return a;
	};
	size_t sim_time = 0;
	auto update_arr = [&] {
		auto res = std::transform_reduce(std::execution::par_unseq, members.begin(), members.end(), std::array<int_fast64_t, 4>{0}, red_arr, trn_arr);
		std::cout << sim_time++;
		for (int i = 0; i < 4; i++) {
			std::cout << "	" << res[i];
			result_stats[i].push_back(res[i]);
		}
		std::cout << '\n';
	};
	update_arr();
	do {
		if (result_stats[1].back() != 0) {
			std::for_each(std::execution::par_unseq, members.begin(), members.end(), [&](std::shared_ptr<Individual> elem) { elem->infect(); });
		}
		std::for_each(std::execution::par_unseq, members.begin(), members.end(), [&](std::shared_ptr<Individual> elem) { elem->update_status(); });
		update_arr();
	} while (result_stats[1].back() != 0);
	tc.stop();
	std::cout << "Simulation time: " << tc.measured_timespan().count() << '\n';
	return result_stats;
}

Population::internal_rnd_norm::internal_rnd_norm(double k_mean, double k_stddev) : norm(k_mean, k_stddev) {}
