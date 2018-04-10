#pragma once
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

int64_t Shenks(uint64_t a, uint64_t b, uint64_t m);
int64_t Pollard(uint64_t alpha, uint64_t beta, uint64_t N);
int64_t Polig_Helman(uint64_t a, uint64_t b, uint64_t p);
int64_t index(uint64_t a, uint64_t b, uint64_t p, int n);

