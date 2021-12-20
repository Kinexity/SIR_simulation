// SIR_simulation.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <fstream>
#include "Population.h"
#include <matplot/matplot.h>
#include <sstream>
//namespace pt = matplotlibcpp;

int main() {
	std::array<std::vector<int_fast64_t>, status_count> res;
	std::stringstream ss, ss2;
	size_t N = 10000000;
	do {
		ss.clear();
		Population pop{ Model::SIR, N, 10, 100, 5, {0.0004, 0.004, 0.003, 0.1} };
		pop.initialize_simulation();
		res = pop.simulate(ss, ss2);
	} while (res[3].back() < std::sqrt(N) && false);
	std::vector<int_fast64_t> time_stamps;
	std::vector<int_fast64_t> sum;
	sum.resize(res[0].size());
	for (int i = 1; i < 4; i++) {
		std::transform(std::execution::par_unseq, sum.begin(), sum.end(), res[i].begin(), sum.begin(), std::plus<>());
	}
	time_stamps.resize(res[0].size());
	std::iota(time_stamps.begin(), time_stamps.end(), 0);
	std::string str[status_count] = { "y-", "C1-","r-","g-","k-" };
	std::string legend[status_count] = { "Zagorożeni", "Narażeni","Zarażeni","Wyzdrowiali","Zmarli" };
	std::string str_s = "b-";
	std::cout << std::filesystem::current_path();
	std::fstream("sim_data.txt", std::ios::trunc | std::ios::out) << ss.str();
	std::fstream("sim_histo_data.txt", std::ios::trunc | std::ios::out) << ss2.str();
}