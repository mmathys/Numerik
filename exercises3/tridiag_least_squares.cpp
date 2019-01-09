#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

/* @brief 
 * @param[in] z An $n$-dimensional vector containing one side of input data
 * @param[in] c An $n$-dimensional vector containing the other side of input data
 * @param[out] x The vector of parameters $(\alpha,\beta)$, intercept and slope of the line fitted
 */
Eigen::Vector2d lsqEst(const Eigen::VectorXd &z, const Eigen::VectorXd &c)
{
	Eigen::Vector2d x;

	int n = z.size();
	cout << "size = " << n << endl;

	cout << "Vector z = " << endl;
	cout << z << endl;

	cout << "Vector c = " << endl;
	cout << c << endl;

	// normal equation
	MatrixXd A(n, n);
	for (int i = 0; i < n; i++)
	{
		int diagSize = A.diagonal(i).size();
		VectorXd diag = z(i) * VectorXd::Ones(diagSize);
		A.diagonal(i) = diag;
		A.diagonal(-i) = diag;
	}

	cout << "Matrix A = " << endl;
	cout << A << endl;

	// bad normal equation
	//VectorXd c2 = A.transpose() * c;
	//VectorXd A2 = A.transpose() * A;
	//x = A2.fullPivLu().solve(c2);

	return x;
}

int main()
{
	int n = 10;
	Eigen::VectorXd z(n), c(n);
	for (int i = 0; i < n; ++i)
	{
		z(i) = i + 1;
		c(i) = n - i;
	}

	Eigen::Vector2d x = lsqEst(z, c);

	std::cout << "alpha = " << x(0) << std::endl;
	std::cout << "beta = " << x(1) << std::endl;

	return 0;
}
