#include <iostream>
#include <vector>
#include <random>
#include "methods.h"
#include "functions.h"

using namespace std;
using int128_t = boost::multiprecision::int128_t;

vector<vector<int128_t>> gauss(int128_t **a, int128_t *y, int128_t n, int128_t p)
{
	//int128_t **mas = new int128_t *[n];
	vector<vector<int128_t>> mas(n.convert_to<size_t>());
	int128_t max = -1;
	int k, index;
	int128_t s = 0; //предыдущще гцд
	int128_t rez = 0; //сумарное гцд
	int128_t w = 0;//для подсчета гцд
	int128_t ost = 0; //ответы для дерева
	bool check = true;
	bool new_check = true;
	int128_t gcd;
	//x = new int128_t[n];
	k = 0;
	while (k < n - 1)
	{
		// Поиск строки с обратным элементом
		index = k;
		for (int i = k; i < n&& check; i++)
		{
			if (ReverseElement(abs(a[i][k]), p) != -1)
			{
				max = abs(a[i][k]);
				index = i;
				check = false;
			}
		}

		//если нет обратного элемента
		if (max == -1)
		{
			// Поиск строки с максимальным a[i][k]
			max = abs(a[k][k]);
			index = k;
			for (int i = k + 1; i < n; i++)
			{
				if (abs(a[i][k]) > max)
				{
					max = abs(a[i][k]);
					index = i;
				}
			}
			// Перестановка строк

			for (int j = 0; j < n; j++)
			{
				int128_t temp = a[k][j];
				a[k][j] = a[index][j];
				a[index][j] = temp;
			}
			int128_t temp = y[k];
			y[k] = y[index];
			y[index] = temp;

			for (int i = k + 1; i < n; i++)
			{
				int128_t temp = a[i][k];
				if (abs(temp) < 0) continue; // для нулевого коэффициента пропустить
				for (int j = k; j < n; j++)
				{
					a[i][j] = int128_t(max*a[i][j] - a[k][j] * temp) % p;
					if (a[i][j] < 0)
						do {
							a[i][j] += p;
						} while (a[i][j] < 0);
				}
				y[i] = int128_t(max*y[i] - y[k] * temp) % p;
				if (y[i] < 0)
					do {
						y[i] += p;
					} while (y[i] < 0);
			}
		}
		//если существует обратный элемент
		else
		{
			for (int j = 0; j < n; j++)
			{
				int128_t temp = a[k][j];
				a[k][j] = a[index][j];
				a[index][j] = temp;
			}
			int128_t temp = y[k];
			y[k] = y[index];
			y[index] = temp;


			for (int i = 0; i < n; i++)
				a[k][i] = int128_t(a[k][i] * ReverseElement(max, p)) % p;
			y[k] = int128_t(y[k] * ReverseElement(max, p)) % p;

			// Нормализация уравнений
			for (int i = k + 1; i < n; i++)
			{
				int128_t temp = a[i][k];
				if (abs(a[i][k]) < 0) continue; // для нулевого коэффициента пропустить
				for (int j = k; j < n; j++)
				{
					a[i][j] = int128_t(a[i][j] - temp * a[k][j]) % p;
					if (a[i][j] < 0)
						do {
							a[i][j] += p;
						} while (a[i][j] < 0);
				}
				y[i] = int128_t(y[i] - temp*y[k]) % p;
				if (y[i] < 0)
					do {
						y[i] += p;
					} while (y[i] < 0);
			}
		}
		k++;
		max = -1;
		check = true;
	}

	// обратная подстановка

	for (k = int(n - 1); k >= 0; k--)
	{
		gcd = int128_t(NOD(a[k][k], p));
		rez = gcd*s;

		if (k == n - 1)
		{
			mas[k] = vector<int128_t>(gcd.convert_to<size_t>());
			for (int i = 0; i < gcd; i++)
				mas[k][i] = int128_t(y[k] / gcd*ReverseElement(a[k][k] / gcd, p / gcd) % (p / gcd) + i*(p / gcd));
		}
		else
		{
			mas[k] = vector<int128_t>(rez.convert_to<size_t>());
			for (int i = 0; i < rez; i++)
			{
				ost = y[k];
				for (int r = int(n - 1); r > k; r--)
				{
					ost -= a[k][r] * mas[r][int(i / (rez / mas[r].size()))] % p;
				}
				if (ost < 0)
					ost += p;
				mas[k][i] = int128_t( ost / gcd*ReverseElement(a[k][k] / gcd, p / gcd) % (p / gcd) + w*(p / gcd));
				if (mas[k][i] < 0)
					mas[k][i] += p;
				w++;
				if (w >= gcd)
					w = 0;
			}
		}
		s = mas[k].size();
	}
	return mas;
}

int64_t index(uint64_t a, uint64_t b, uint64_t p, int n)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution <unsigned long long> dis(1, uint64_t(p));
	int128_t x = 0;
	vector<int64_t> vector1{ 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019};
	vector<int128_t> vector3(n);
	int128_t **mas = new int128_t *[n];
	int128_t *result = new int128_t[n];
	int128_t r = 0;
	int e = n - 1;
	int128_t s = 1;
	int128_t t = 1;
	int128_t l = 0;
	bool zax = false;
	int128_t check = 0;
	int128_t g = 0;
	bool new_check = true;
	bool ew = true;

	for (int i = 0; i<n; i++)
		mas[i] = new int128_t[n];


	for (int i = 0; i <n; i++)
		for (int j = 0; j <n; j++)
			mas[i][j] = 0;
	//выбор элемента, который еще не встречался
	for (int i = 0; i < n; i++)
	{
		do
		{
			do
			{
				s = dis(gen);

				for (int i = 0; i < vector3.size(); i++)
				{
					if (NOD(vector3[i], s) != 1 && vector3[i] != 1 && vector3[i] != 0)
						check++;
				}
				if (check != 0)
				{
					check = 0;
					zax = true;
				}
				else
				{
					r = int128_t( powmod(a, s, p));
					if (r != 1)
					{
						result[i] = s;
						zax = false;
						vector3[i] = s;
					}
				}
			} while (zax);
			//возведение в степень и разложение числа
			r = int128_t(powmod(a, s, p));
			result[i] = s;
			while (e != -1)
			{
				if (r%vector1[e] == 0)
				{
					r /= vector1[e];
					mas[i][e]++;
				}
				else
					e--;
			}

			//если элемент не ддо конца разложился
			if (r != 1)
			{
				for (int j = 0; j < n; j++)
					mas[i][j] = 0;
				result[i] = 0;
				vector3[i] = 0;
				//	i--;
			}
			//проверка получившийся системы
			else if (i == n - 1)
			{
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						if (mas[j][i] != 0)
							g++;
					}
					//если получился столбец из 0
					if (g == 0)
						new_check = false;

					//если получился столбец из одного элемента
					if (g == 1)
					{
						for (int f = 0; f < n; f++)
						{
							if (mas[f][i] != 0)
							{
								if (ReverseElement(mas[f][i], p - 1) == -1)
									new_check = false;
							}
						}
					}
					g = 0;
				}
				//обнуляем последнюю строку
				if (!new_check)
				{
					new_check = true;
					for (int j = 0; j < n; j++)
						mas[i][j] = 0;
					result[i] = 0;
					i--;
				}
			}
			e = n - 1;
			zax = true;
		} while (r != 1);
	}
	//решение системы
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << mas[i][j] << " ";
		cout << result[i];
		cout << endl;
	}
	cout << endl << endl << endl;
	auto c = gauss(mas, result, n, p - 1);
	s = 0;
	int128_t koef = 0;
	int128_t sum = 0;
	//разложение числа
	do
	{
		for (int i = 0; i < n; i++)
			mas[0][i] = 0;
		do {
			s = dis(gen);
			r = int128_t(b*powmod(a, s, p) % p);
		} while (r == 1);

		while (e != -1)
		{
			if (r%vector1[e] == 0)
			{
				r /= vector1[e];
				mas[0][e]++;
			}
			else
				e--;
		}
		e = n - 1;
	} while (r != 1);
	cout << endl << s << endl << endl;

	for (int h = 0; h < c[0].size() && ew == true; h++)
	{
		for (int j = 0; j < n; j++)
			if (mas[0][j] != 0)
				x += (mas[0][j] * c[j][h / (c[0].size() / c[j].size())]) % (p - 1);
		x = (x - s) % (p - 1);
		if (powmod(a, x, p) == b)
			ew = false;
		else
			x = 0;
	}
	delete[] mas,result;
	return uint64_t(x);
}