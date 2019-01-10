#include <cmath>
#include <iostream>
#include <functional>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

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
	int size = x.size();
	VectorXd X = VectorXd::Map(x.data(), x.size());
	VectorXd sum1 = (1. / 6) * (X.tail(size - 2) - X.head(size - 2)).cwiseProduct(X.segment(1, size - 2).unaryExpr(f));
	VectorXd sum2 = (2. / 3) * (X.tail(size - 1) - X.head(size - 1)).cwiseProduct((0.5 * (X.tail(size - 1) + X.head(size - 1))).unaryExpr(f));
	double part1 = (1. / 6) * (X(1) - X(0)) * f(X(0));
	double part2 = sum1.sum();
	double part3 = sum2.sum();
	double part4 = (1. / 6) * (X(size - 1) - X(size - 2)) * f(X(size - 1));

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
