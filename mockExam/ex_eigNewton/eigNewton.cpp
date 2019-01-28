#include <iostream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void eigNewton(const Eigen::MatrixXd &A, double tol, int maxItr, Eigen::VectorXd &z)
{
	// initial guess: z
	int m = z.size();
	int n = m - 1;
	cout << A.size() << endl;
	cout << n << endl;
	assert(A.size() == n);

	cout << MatrixXd::Identity(n, n) << endl;

	for (int i = 0; i < maxItr; i++)
	{
		VectorXd x = z.head(z.size() - 1);
		VectorXd lambda = z.tail(1);

		// Build DF
		MatrixXd DF(m, m);
		DF.block(0, 0, n, n) = A - lambda * MatrixXd::Identity(n, n);
	}
}

int main()
{
	int n = 2;
	int m = n + 1;

	Eigen::VectorXd x = Eigen::VectorXd::Ones(n);
	x << 1.0, 0.0;
	double lambda = 0.0;

	Eigen::VectorXd z(m);
	z.head(n) = x;
	z(n) = lambda;

	Eigen::MatrixXd A = Eigen::MatrixXd::Ones(n, n);
	A << 2, 1, 1, 2;

	eigNewton(A, .0001, 10, z);

	std::cout << "eigenvector:" << std::endl;
	std::cout << z.head(n) << std::endl;
	std::cout << "eigenvalue: " << z(n) << std::endl;

	return 0;
}
