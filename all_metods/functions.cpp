#include "functions.h"
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

mp::cpp_int powmod(mp::cpp_int a, mp::cpp_int b, mp::cpp_int m) {
	mp::cpp_int res = 1;

	while (b > 0)
		if (b & 1) {
			res = (res * a) % m;
			--b;
		}
		else {
			a = (a * a) % m;
			b >>= 1;
		}
		return res % m;
}
mp::cpp_int gcdex(mp::cpp_int a, mp::cpp_int b, mp::cpp_int &x, mp::cpp_int &y) {
	if (b == 0) {
		x = 1;
		y = 0;
		return a;
	}
	mp::cpp_int x1, y1;
	mp::cpp_int d1 = gcdex(b, a % b, x1, y1);
	x = y1;
	y = x1 - (a / b) * y1;
	return d1;
}

mp::cpp_int ReverseElement(mp::cpp_int a, mp::cpp_int N) {
	mp::cpp_int x, y, d;
	d = gcdex(a, N, x, y);
	if (d != 1) {
		return -1;
	}
	else
	{
		if (x < 0)
			x += N;
		return x;
	}
}

mp::cpp_int NOD(mp::cpp_int a, mp::cpp_int b)
{
	while (a > 0 && b > 0)

		if (a > b)
			a %= b;

		else
			b %= a;

	return a + b;
}

mp::cpp_int phi(mp::cpp_int n) {
	mp::cpp_int result = n;
	for (mp::cpp_int i = 2; i*i <= n; ++i)
		if (n % i == 0) {
			while (n % i == 0)
				n /= i;
			result -= result / i;
		}
	if (n > 1)
		result -= result / n;
	return result;
}