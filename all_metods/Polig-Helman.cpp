#include <iostream>
#include <map>
#include "methods.h"
#include "functions.h"

using namespace std;

map<mp::cpp_int, mp::cpp_int> factorize(mp::cpp_int x) {
	map<mp::cpp_int, mp::cpp_int> factors;

	for (mp::cpp_int i = 2; i <= sqrt(x); i++) {
		while (x % i == 0) {
			factors[i]++;
			x /= i;
		}
	}

	if (x != 1) {
		factors[x]++;
	}

	return factors;
}

uint64_t china_alg(map<mp::cpp_int, mp::cpp_int> factor)
{
	mp::cpp_int k = 0;
	mp::cpp_int o = 1;
	auto it = factor.begin();
	mp::cpp_int u,e,h,x=0,w=0;

	while (w<factor.size())
	{
		o *= it->first;
		it++;
		w++;
	}
	it = factor.begin();
	for (mp::cpp_int i = 0; i < factor.size(); i++)
	{
		e = o / it->first;
		h=ReverseElement(e , it->first);
		if (h < 0)
			h = h+ it->first;
		x += it->second * e*h;
		x %= o;
		it++;
	}
	return uint64_t(x);
}

int64_t Polig_Helman(uint64_t a, uint64_t b, uint64_t p)
{
	mp::cpp_int l = phi(p);
	auto ps = factorize(l);
	auto it = ps.begin();
	auto it2 = ps.begin();
	auto it3 = ps.begin();
	auto it4 = ps.end();
	mp::cpp_int **r;
	r = new mp::cpp_int*[ps.size()];
	try
	{
		for (mp::cpp_int i = 0; i < ps.size(); i++)
		{
			r[i.convert_to<size_t>()] = new mp::cpp_int[(it->first).convert_to<size_t>()];
			it++;
		}
	}
	catch(...)
	{
		return -2;
	}
	mp::cpp_int t = 0;
	mp::cpp_int bk = 0;
	mp::cpp_int x0 = 0;
	mp::cpp_int ui = 0;
	mp::cpp_int step = 1;
	it = ps.begin();
	map<mp::cpp_int, mp::cpp_int> factorik;

	for (mp::cpp_int i = 0; i < ps.size(); i++)
	{
		for (mp::cpp_int j = 0; j < it->first; j++)
			r[i.convert_to<size_t>()][j.convert_to<size_t>()] = powmod(a, j*(p - 1) / it->first, p);
		it++;
	}

	it = ps.begin();
	//тут считает x
	do {
		for (mp::cpp_int h = 0; h < it->second; h++)
		{
			bk = powmod(mp::cpp_int(b)*mp::pow(mp::cpp_int(a), int(-x0)), (p - 1) / (pow(it->first, int((h + 1)))), p);
			for (mp::cpp_int i = 0; i < ps.size(); i++)
			{
				for (mp::cpp_int j = 0; j < it2->first; j++)
				{
					if (r[i.convert_to<size_t>()][j.convert_to<size_t>()] == bk)
					{
						x0 = j;
						ui += x0*pow(it->first, int(h));
					}
				}
				it2++;
			}
			it2 = ps.begin();
		}

		x0 = 0;
		factorik[pow(it->first, int(it->second))] = ui;
		it++;
		ui = 0;
		it2 = ps.begin();
		t++;
	} while (t<ps.size());
	delete[] r;
	return china_alg(factorik);
}