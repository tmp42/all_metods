#pragma once
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

mp::cpp_int powmod(mp::cpp_int a, mp::cpp_int b, mp::cpp_int m);
mp::cpp_int gcdex(mp::cpp_int a, mp::cpp_int b, mp::cpp_int &x, mp::cpp_int &y);
mp::cpp_int ReverseElement(mp::cpp_int a, mp::cpp_int N);
mp::cpp_int NOD(mp::cpp_int a, mp::cpp_int b);
mp::cpp_int phi(mp::cpp_int n);