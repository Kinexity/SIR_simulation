#pragma once
#ifndef DISEASE_H
#define DISEASE_H

struct Disease {
	double
		infection_or_exposition_probability = 0.0,
		recovery_probability = 0.0,
		death_probability = 0.0,
		exposition_to_infection_probability = 0.0;
};

#endif // !DISEASE_H


