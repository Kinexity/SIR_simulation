// SIR_simulation.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <fstream>
#include <matplotlibcpp.h>
#include "Population.h"
namespace pt = matplotlibcpp;

int main() {
	std::array<std::vector<int_fast64_t>, status_count> res;
	do {
		Population pop{ Model::SEIRD, 100000, 1, 100, 5, {0.006, 0.06, 0.003, 0.1} };
		pop.initialize_simulation();
		res = pop.simulate();
	} while (res[3].back() < 400);
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
	if (true) {
		for (int type = 0; type < status_count; type++) {
			pt::/*named_*/plot(/*legend[type],*/ time_stamps, res[type], str[type]);
		}
	}
	else {
		for (int type = 0; type < status_count; type++) {
			pt::named_semilogy(legend[type], time_stamps, res[type], str[type]);
		}
		pt::semilogy(time_stamps, sum, str_s);
	}
	pt::legend();
	pt::show();
}