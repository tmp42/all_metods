#include <iostream>
#include "methods.h"
#include <vector>
#include "functions.h"

using namespace std;



void new_xab( mp::cpp_int& x,mp::cpp_int& a,mp::cpp_int& b, mp::cpp_int alpha, mp::cpp_int  beta,  mp::cpp_int N, mp::cpp_int n) {
	switch (x.convert_to<size_t>() % 3) {
	case 0: x = x*x     % N;  a = a * 2 % n;  b = b * 2 % n;  break;
	case 1: x = x*alpha % N;  a = (a + 1) % n;                  break;
	case 2: x = x*beta  % N;                  b = (b + 1) % n;  break;
	}
}

int64_t Pollard(uint64_t alpha, uint64_t beta, uint64_t N) {
	mp::cpp_int x = 1, a = 0, b = 0;
	mp::cpp_int X = x, A = a, B = b;
	mp::cpp_int ev, eu, nod, chk, n = N-1, h = 0;

	for (int i = 1; i < n; ++i) {
		new_xab(x, a, b,alpha,beta,N,n);
		new_xab(X, A, B, alpha, beta, N, n);
		new_xab(X, A, B, alpha, beta, N, n);
		if (x == X)
		{
			ev = B - b;
			if (ev < 0)
				ev += n;
			eu = a - A;
			if (eu < 0)
				eu += n;
			nod = NOD(ev, n);
			ev /= nod;
			eu /= nod;
			n /= nod;
			chk = ReverseElement(ev, n);
			if (chk < 0)
				chk += n;
			x = chk*eu%n;
			while (h < nod)
			{
				if (powmod(alpha, x, N) != beta)
					x += n;
				h++;
			}
			return int64_t(x);
			break;
		}
	}
	return 0;
}