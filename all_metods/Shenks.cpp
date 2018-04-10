#include <iostream>
#include <math.h>
#include "methods.h"
#include <map>
#include "functions.h"
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;
using int128_t = boost::multiprecision::int128_t;

using namespace std;



int64_t Shenks(uint64_t a, uint64_t b, uint64_t m) {
	int128_t n = (int128_t)sqrt(m) + 1;
	try
	{
		int128_t *r;
		r = new int128_t[n.convert_to<size_t>()];
		delete[] r;
	}
	catch (...)
	{
		return -2;
	}
	map<int128_t,int128_t> vals;
	for (int128_t i = n; i >= 1; --i)
		vals[int128_t(powmod(a, i * n, m))] = i;
	for (int128_t i = 0; i <= n; ++i) {
		int128_t cur = (int128_t(powmod(a, i, m) * b)) % m;
		if (vals.count(cur)) {
			int128_t ans = vals[cur] * n - i;
			if (ans < m)
				return uint64_t(ans);
		}
	}
	return -1;
}