#include "SIR_analytical.h"
#include <gsl/gsl_multifit_nlinear.h>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <iostream>
//#include <curve_fit.h>
#include <fstream>

using boost::math::quadrature::trapezoidal;

double I_Simpson(std::function<double(double)> fn, double x_p, double x_k, size_t n) {
	auto delta = (x_k - x_p) / (n - 1);
	static std::vector<double> increments;
	increments.resize(n - 1);
	std::transform(std::execution::par_unseq, increments.begin(), increments.end(), increments.begin(), [&](auto& elem) {return delta * (std::distance(&*increments.begin(), &elem) + 1); });
	return delta / 3. * (fn(x_p) + fn(x_k) + 2 * std::transform_reduce(std::execution::par_unseq, increments.begin(), increments.end(), 0.0, std::plus<double>(), [&](auto& increment) { return (2 - (std::distance(&*increments.begin(), &increment) + 1) % 2) * fn(x_p + increment); }));
}

double integrate(std::function<double(double)> fn, double x_p, double x_k) {
	size_t n_0 = 10, n_1 = 0;
	constexpr auto eps = 1e-4;
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

double t(double z, double inv_gamma, double x_0, double beta, double offset) {
	auto u = std::exp(-beta * inv_gamma * z);
	if (z < 1.) {
		return inv_gamma * trapezoidal([&](double x) {return 1. / (x * (std::log(x) - beta * inv_gamma * x_0 * x + beta * inv_gamma)); }, u, 1., 1e-6, 20) + offset;
	}
	else {
		return -inv_gamma * trapezoidal([&](double x) {return 1. / (x * (std::log(x) - beta * inv_gamma * x_0 * x + beta * inv_gamma)); }, 1., u, 1e-6, 20) + offset;
	}
}

std::array<double, 3> fit() {
	//DFTaskPtr task;
	//auto status = dfdNewTask1D(&task, nx, x, xhint, ny, y, yhint);
	//Modify the task parameters.
	//status = dfdEditPPSpline1D(task, s_order, c_type, bc_type, bc, ic_type, ic, scoeff, scoeffhint);
	//Perform Data Fitting spline - based computations.You may reiterate steps 2 - 3 as needed.
	//status = dfdInterpolate1D(task, estimate, method, nsite, site, sitehint, ndorder,dorder, datahint, r, rhint, cell);
	//Destroy the task or tasks.
	//status = dfDeleteTask(&task);
	std::fstream file;
	std::vector<double> t_data, z_data;
	file.open("C:\\Users\\26kuba05\\source\\repos\\SpacePhysicsSimulator\\conv.txt", std::ios::in);
	std::string line;
	while (std::getline(file, line)) { 
		int day;
		double recovered_prob;
		std::stringstream(line) >> day >> recovered_prob;
		std::cout << day << "	" << recovered_prob << '\n';
		t_data.push_back((double)day);
		z_data.push_back(recovered_prob);
	}
	file.close();
	std::function fn = t;
	//curve_fit(fn, z_data, t_data, {});
	return {};
}