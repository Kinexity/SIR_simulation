#pragma once
#ifndef DISEASE_H
#define DISEASE_H
#include <cmath>

struct Disease_ext {
	double
		infection_or_exposition_probability = 0.0,
		recovery_probability = 0.0,
		death_probability = 0.0,
		exposition_to_infection_probability = 0.0;
};


struct Disease {
	double
		infection_or_exposition_probability = 0.0,
		infection_or_exposition_probability_log = 0.0,
		recovery_probability = 0.0,
		recovery_probability_log = 0.0,
		death_probability = 0.0,
		death_probability_log = 0.0,
		exposition_to_infection_probability = 0.0,
		exposition_to_infection_probability_log = 0.0;
};

inline Disease dis_int_conv(Disease_ext dis) {
	return Disease{ dis.infection_or_exposition_probability,
		std::log(1. - dis.infection_or_exposition_probability),
		dis.recovery_probability,
		std::log(1. - dis.recovery_probability),
		dis.death_probability,
		std::log(1. - dis.death_probability),
		dis.exposition_to_infection_probability,
		std::log(1. - dis.exposition_to_infection_probability) };
};

#endif // !DISEASE_H


