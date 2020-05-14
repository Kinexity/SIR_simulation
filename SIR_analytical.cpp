#include "SIR_analytical.h"

double I_Simpson(std::function<double(double)> fn, double x_p, double x_k, size_t n) {
	auto delta = (x_k - x_p) / (n - 1);
	static std::vector<double> increments;
	increments.resize(n - 1);
	std::transform(std::execution::par_unseq, increments.begin(), increments.end(), increments.begin(), [&](auto& elem) {return delta * (std::distance(&*increments.begin(), &elem) + 1); });
	return delta / 3. * (fn(x_p) + fn(x_k) + 2 * std::transform_reduce(std::execution::par_unseq, increments.begin(), increments.end(), 0.0, std::plus<double>(), [&](auto& increment) { return (2 - (std::distance(&*increments.begin(), &increment) + 1) % 2) * fn(x_p + increment); }));
}

double integrate(std::function<double(double)> fn, double x_p, double x_k) {
	size_t n_0 = 10, n_1 = 0;
	constexpr auto eps = 1e-6;
	auto delta = std::numeric_limits<double>::max();
	auto fn_b = std::bind(I_Simpson, fn, x_p, x_k, _1);
	auto I0 = fn_b(n_0), I1 = 0.0;
	do {
		n_1 = 2 * (n_0 - 1) + 1;
		I1 = fn_b(n_1);
		delta = std::abs(1 - I0 / I1);
		if (delta < eps) {
			break;
		}
		n_0 = n_1;
		I0 = I1;
	} while (true);
	return I1;
}

double t(double u, double inv_gamma, double x_0, double beta) {
	return inv_gamma * integrate([&](double x) {return 1 / (x * (std::log(x) - beta * inv_gamma * x_0 * x - beta * inv_gamma)); }, u, 1);
}

std::array<double, 3> fit() {
	DFTaskPtr task;

	return {};
}