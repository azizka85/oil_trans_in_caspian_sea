#include <iostream>

#include <slae/direct/tridiagonal.h>

using namespace std;

using namespace SLAE::Direct;

double maxAbsVal(vector<double>& v) {
	double maxVal = 0;

	for (auto& val : v) {
		if (abs(val) > maxVal) {
			maxVal = abs(val);
		}
	}

	return maxVal;
}

int main() {
	int n = 100;

	double l = 1;
	double l1 = 0;

	double c = -2;
	double c0 = 1;
	double c1 = 1;

	double r = 1;
	double r0 = 0;

	vector<double> u(n, 0);
	vector<double> d(n, 0);

	d[0] = 0;
	d[n - 1] = 1;

	auto d1 = d;

	Tridiagonal::solve(l, l1, c, c0, c1, r, r0, d, u);

	span<double> su(u);
	span<double> sd(d1);

	Tridiagonal::residual(n, l, l1, c, c0, c1, r, r0, su, sd);

	cout << "Max absolute value of residual: " << maxAbsVal(d1) << endl;

	vector<double> al(n, 1);

	al[0] = 0;
	al[n - 1] = 0;

	vector<double> ac(n, -2);

	ac[0] = 1;
	ac[n - 1] = 1;

	vector<double> ar(n, 1);

	ar[0] = 0;
	ar[n - 1] = 0;

	d = vector<double>(n, 0);

	d[0] = 0;
	d[n - 1] = 1;

	d1 = d;	

	Tridiagonal::solve(al, ac, ar, d, u);
	Tridiagonal::residual(n, al, ac, ar, su, sd);

	cout << "Max absolute value of residual: " << maxAbsVal(d1) << endl;

	l = 1;
	l1 = 2;

	c = -2;
	c0 = 1;
	c1 = -2;

	r = 1;
	r0 = 0;

	d = vector<double>(n, 0);

	d[0] = 1;
	d[n - 1] = 0;

	d1 = d;

	Tridiagonal::solve(l, l1, c, c0, c1, r, r0, d, u);
	Tridiagonal::residual(n, l, l1, c, c0, c1, r, r0, su, sd);

	cout << "Max absolute value of residual: " << maxAbsVal(d1) << endl;

	al[0] = 0;
	al[n - 1] = 2;

	ac[0] = 1;
	ac[n - 1] = -2;

	ar[0] = 0;
	ar[n - 1] = 0;

	d = vector<double>(n, 0);

	d[0] = 1;
	d[n - 1] = 0;

	d1 = d;

	Tridiagonal::solve(al, ac, ar, d, u);
	Tridiagonal::residual(n, al, ac, ar, su, sd);

	cout << "Max absolute value of residual: " << maxAbsVal(d1) << endl;

	return 0;
}