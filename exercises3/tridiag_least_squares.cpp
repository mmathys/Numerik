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
	assert(z.size() == c.size());

	// construct A
	MatrixXd A(n, 2);
	A.col(0) = z;
	A.block(0, 1, n - 1, 1) = z.tail(n - 1);
	A.block(1, 1, n - 1, 1) += z.head(n - 1);

	cout << "A =\n"
		 << A << endl;

	return (A.transpose() * A).fullPivLu().solve(A.transpose() * c);
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
