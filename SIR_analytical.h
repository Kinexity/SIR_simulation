#pragma once
#ifndef SIR_ANALYTICAL_H
#define SIR_ANALYTICAL_H
#include <functional>
#include <algorithm>
#include <execution>
#include <utility>
#include <array>
#include <fstream>
#include <memory>
using namespace std::placeholders;

std::array<double, 3> fit();

#endif // !SIR_ANALYTICAL_H