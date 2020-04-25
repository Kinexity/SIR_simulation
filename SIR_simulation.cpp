// SIR_simulation.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <fstream>
#include <matplotlibcpp.h>
#include "Population.h"
namespace pt = matplotlibcpp;

int main() {
	std::array<std::vector<int_fast64_t>, 4> res;
	do {
		Population pop{ Model::SIR, 1000000, 10, 100, 5, {0.0005, 0.06, 0.003} };
		pop.initialize_simulation();
		res = pop.simulate();
	} while (res[2].back() < 400);
	std::vector<size_t> time_stamps;
	std::vector<size_t> sum;
	sum.resize(res[0].size());
	for (int i = 1; i < 4; i++) {
		std::transform(sum.begin(), sum.end(), res[i].begin(), sum.begin(), std::plus<>());
	}
	time_stamps.resize(res[0].size());
	std::iota(time_stamps.begin(), time_stamps.end(), 0);
	std::string str[4] = { "y-","r-","g-","k-" };
	std::string legend[4] = { "Zagorożeni","Zarażeni","Wyzdrowiali","Zmarli" };
	for (int type = 0; type < 4; type++) {
		pt::semilogy(time_stamps, res[type], str[type]);
	}
	std::string str_s = "b-";
	pt::semilogy(time_stamps, sum, str_s);
	pt::legend();
	pt::show();
}