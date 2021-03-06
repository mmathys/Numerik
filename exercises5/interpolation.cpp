#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct Newton
{
	Newton(const Eigen::VectorXd &x) : _x(x), _a(x.size()) {}
	void Interpolate(const Eigen::VectorXd &y);
	double operator()(double x) const;

  private:
	Eigen::VectorXd _x; // nodes
	Eigen::VectorXd _a; // coefficients
};

// Compute the coefficients in the Newton basis.
void Newton::Interpolate(const Eigen::VectorXd &y)
{
	int n = _x.size() - 1;

	Eigen::MatrixXd A(n + 1, n + 1);
	for (int i = 0; i < n + 1; i++)
	{
		A(i, 0) = 1;
		int t = _x(i);
		for (int j = 1; j <= i; j++)
		{
			A(i, j) = A(i, j - 1) * (_x(i) - _x(j - 1));
		}
	}

	//cout << A << endl;
	//cout << _a.size() << endl;
	//cout << y.size() << endl;

	_a = A.partialPivLu().solve(y);
	//cout << _a << endl;

	// TODO: Task (a)
	// ...
	// compute _a
	// ...
}

// Evaluate the interpolant at x.
double Newton::operator()(double x) const
{
	int n = _a.size() - 1;
	double v = _a(n);
	for (int i = n - 1; i >= 0; i--)
	{
		v *= (x - _x(i));
		v += _a(i);
	}
	return v; // dummy
}

struct Lagrange
{
	Lagrange(const Eigen::VectorXd &x);
	void Interpolate(const Eigen::VectorXd &y) { _y = y; }
	double operator()(double x) const;

  private:
	Eigen::VectorXd _x; // nodes
	Eigen::VectorXd _l; // weights
	Eigen::VectorXd _y; // coefficients
};

// Compute the weights l for given nodes x.
Lagrange::Lagrange(const Eigen::VectorXd &x) : _x(x), _l(x.size()), _y(x.size())
{
	int n = _l.size();
	for (int i = 0; i < n; i++)
	{
		_l(i) = 1;
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				continue;
			_l(i) *= 1 / (_x(i) - _x(j));
		}
	}
}

// Evaluate the interpolant at x.
double Lagrange::operator()(double x) const
{
	double q = 0;
	int n = _y.size() - 1;

	// calculate w
	double w = 1;
	for (int i = 0; i <= n; i++)
	{
		w *= (x - _x(i));
	}

	for (int i = 0; i <= n; i++)
	{
		double L = w * _l(i) / (x - _x(i));
		q += _y(i) * L;
	}
	return q; // dummy
}

// Runge function
Eigen::VectorXd r(const Eigen::VectorXd &x)
{
	return (1.0 / (1.0 + 25.0 * x.array() * x.array())).matrix();
}

int main()
{
	int n = 5;
	Eigen::VectorXd x;
	x.setLinSpaced(5, -1.0, 1.0);
	Eigen::VectorXd y = r(x);

	Newton p(x);
	p.Interpolate(y); // correct result: p._a = [0.0384615, 0.198939, 1.5252, -3.31565, 3.31565]

	Lagrange q(x); // correct result: p._l = [0.666667, -2.66667, 4, -2.66667, 0.666667]
	q.Interpolate(y);

	// Compute difference of p and q.
	int m = 22;
	double offset = 0.08333333333;
	x.setLinSpaced(m, -1.0 + offset, 1.0 - offset);
	double norm2 = .0;
	for (int i = 0; i < m; ++i)
	{
		double d = p(x(i)) - q(x(i));
		norm2 += d * d;
	}

	// By uniquenss of the interpolation polynomial, we expect p = q.
	std::cout << "This number should be close to zero: " << norm2 << std::endl;

	return 0;
}
