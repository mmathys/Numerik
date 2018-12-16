#include <cmath>
#include <iostream>
#include <functional>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

double GaussLegendre5(const std::function<double(double)> &f, double a, double b)
{
	const std::vector<double> c = {-0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593};
	const std::vector<double> w = {0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505};

	VectorXd cTilde = VectorXd::Map(c.data(), c.size());
	VectorXd wTilde = VectorXd::Map(w.data(), w.size());
	VectorXd _w = (0.5 * (b - a)) * wTilde;
	VectorXd _c = 0.5 * (VectorXd::Ones(c.size()) - cTilde) * a + 0.5 * (VectorXd::Ones(c.size()) + cTilde) * b;
	VectorXd _f = _c.unaryExpr(f);
	VectorXd elements = _w.cwiseProduct(_f);
	return elements.sum();
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x)
{
	// convert into eigen vector
	VectorXd _x = VectorXd::Map(x.data(), x.size());
	double m = _x.size() - 1;
	double a = _x(0);
	double b = _x(m); // = last

	double part1 = (1.0 / 6) * (_x(1) - _x(0)) * f(a);

	VectorXd sum1 = (_x.tail(m - 1) - _x.head(m - 1)).cwiseProduct(_x.segment(1, m - 1).unaryExpr(f));
	double part2 = (1.0 / 6) * sum1.sum();

	VectorXd sum2_first = _x.tail(m) - _x.head(m);
	VectorXd sum2_last = (0.5 * (_x.tail(m) + _x.head(m))).unaryExpr(f);
	VectorXd sum2 = sum2_first.cwiseProduct(sum2_last);
	double part3 = (2.0 / 3) * sum2.sum();

	double part4 = (1.0 / 6) * (_x(m) - _x(m - 1)) * f(b);

	return part1 + part2 + part3 + part4;
}

std::vector<double> LinSpace(int n, double a, double b)
{
	std::vector<double> x(n);
	double d = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i)
	{
		x[i] = i * d + a;
	}
	return x;
}

double f(double x)
{
	return std::sqrt(x);
}

double F(double x)
{
	double y = std::sqrt(x);
	return 2.0 / 3.0 * y * y * y;
}

int main()
{
	int n = 5;
	double a = 0.0;
	double b = 1.0;

	std::cout << "Gauss-Legendre: " << GaussLegendre5(f, a, b) << std::endl;
	std::cout << "Simpson: " << CompositeSimpson(f, LinSpace(n, a, b)) << std::endl;
	std::cout << "Exact value: " << F(b) - F(a) << std::endl;

	return 0;
}
