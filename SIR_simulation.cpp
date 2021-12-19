﻿// SIR_simulation.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <fstream>
#include "Population.h"
#include <matplot/matplot.h>
#include <sstream>
//namespace pt = matplotlibcpp;

int main() {
	std::array<std::vector<int_fast64_t>, status_count> res;
	std::stringstream ss;
	size_t N = 100000;
	do {
		ss.clear();
		Population pop{ Model::SIR, N, 1, 100, 5, {0.0004, 0.004, 0.003, 0.1} };
		pop.initialize_simulation();
		res = pop.simulate(ss);
		std::cout << (double)pop.get()/(N - res[0].back()) << '\n';
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
	std::fstream f("sim_data.txt", std::ios::trunc | std::ios::out);
	f << ss.str();
	//for (int type = 0; type < status_count; type++) {
	//	matplot::plot(/*legend[type],*/ time_stamps, res[type], str[type]);
	//}
	//matplot::show();
	//pt::legend();
	//pt::show();
	//for (int type = 0; type < status_count; type++) {
	//	pt::/*named_*/semilogy(/*legend[type],*/ time_stamps, res[type], str[type]);
	//}
	//pt::semilogy(time_stamps, sum, str_s);
	//pt::legend();
	//pt::show();
}