#include <iostream>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace Eigen;

// QR decomposition using Householder reflections
void householderQR(const MatrixXd &A, MatrixXd &R, MatrixXd &Q)
{

	size_t m = A.rows();
	size_t n = A.cols();

	HouseholderQR<MatrixXd> QR = A.householderQr();
	Q = QR.householderQ();
	R = QR.matrixQR().triangularView<Upper>();
	// If A: m x n, then Q: m x m and R: m x n.
}

// solves the LLSQ using normal equation method
void llsq_normal(const MatrixXd &A, const VectorXd &b, VectorXd &x)
{

	size_t n = A.cols();

	VectorXd b_ = A.transpose() * b;
	x = (A.transpose() * A).fullPivLu().solve(b_);
	std::cout << "\n\n"
			  << x.transpose() << std::endl;
}

// solves the LLSQ using QR decomposition
void llsq_qr(const MatrixXd &A, const VectorXd &b, VectorXd &x)
{

	size_t m = A.rows();
	size_t n = A.cols();

	MatrixXd R, Q;

	householderQR(A, R, Q);
	std::cout << "QR decomposition: R =" << std::endl
			  << R << "\n"
			  << std::endl;
	std::cout << "QR decomposition: Q =" << std::endl
			  << Q << std::endl;

	VectorXd b_ = (Q.transpose() * b).head(n);
	MatrixXd R_ = R.block(0, 0, n, n);

	x = R_.triangularView<Upper>().solve(b_);
	std::cout << x.transpose() << std::endl;
}

// solves the LLSQ using SVD decomposition
double llsq_gsol(const MatrixXd &A, const VectorXd &b, VectorXd &x)
{

	double llsq_err = 0;

	size_t m = A.rows();
	size_t n = A.cols();
	JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV); // SVD decomposition

	int r = 0;
	double tol = 1e-12;
	VectorXd sv = svd.singularValues();
	for (int i = 0; i < sv.size(); i++)
	{ // number of non-zero singular values
		if (fabs(sv(i)) > tol)
			r++;
	}

	// Compute the matrices in Moore-Penrose pseudoinverse
	MatrixXd invD_(r, r);
	for (int i = 0; i < r; i++)
		invD_(i, i) = 1.0 / sv(i);

	MatrixXd U = svd.matrixU();
	MatrixXd U_ = U.block(0, 0, m, r);

	MatrixXd V = svd.matrixV();
	MatrixXd V_ = V.block(0, 0, n, r);

	//MatrixXd D = sv.asDiagonal();
	//std::cout << "D =" << std::endl << D << "\n" << std::endl;
	//std::cout << "U =" << std::endl << U << "\n" << std::endl;
	//std::cout << "V =" << std::endl << V << std::endl;

	x = V_ * invD_ * U_.transpose() * b; // Generalized solution

	std::cout << "Sigma^-1 =" << std::endl
			  << invD_ << std::endl;

	std::cout << "x =" << std::endl
			  << x << std::endl;

	llsq_err = (U.block(0, r, m, m - r).transpose() * b).norm(); // least-squares error

	return llsq_err;
}

void run_llsq_gsol()
{

	size_t m = 5;
	size_t n = 4;

	VectorXd t(m), b(m);
	t << 11, 32, 77, 121, 152;
	b << 9, 10, 12, 14, 15;

	MatrixXd A(m, n);
	for (int i = 0; i < m; i++)
	{
		A(i, 0) = 1;
		A(i, 1) = sin(2 * M_PI * t(i) / 366);
		A(i, 2) = cos(2 * M_PI * t(i) / 366);
		A(i, 3) = sin(M_PI * t(i) / 366) * cos(M_PI * t(i) / 366);
	}
	std::cout << "A =" << std::endl
			  << A << "\n"
			  << std::endl;

	VectorXd x(n);
	double llsq_err;
	llsq_err = llsq_gsol(A, b, x);

	std::cout << "x: " << x.transpose() << std::endl;
	std::cout << "llsq error: " << llsq_err << std::endl;
}

/***** TESTING ******/
// "" Do NOT CHANGE the routines below ""

void runTests()
{

	bool success = true;
	double TOL = 1e-6;

	size_t m = 5;
	size_t n = 4;

	VectorXd t_test(m), b_test(m);
	t_test << 11, 32, 77, 121, 152;
	b_test << 10, 11, 12, 15, 16;

	MatrixXd A_test(m, n);
	A_test << 1, 0.187718565437199, 0.982222856682841, 0.0938592827185993,
		1, 0.522132581076825, 0.85286433140216, 0.261066290538412,
		1, 0.969178065439254, 0.246361274293315, 0.484589032719627,
		1, 0.87448095783104, -0.485059846195196, 0.43724047891552,
		1, 0.507415093293846, -0.861701759948068, 0.253707546646923;

	std::cout << "A_test =" << std::endl
			  << A_test << "\n"
			  << std::endl;

	VectorXd x_test(n), x_ans(n);
	double llsq_err_test, llsq_err_ans = 0.67829;
	x_ans << 13.455, -0.238115, -3.21762, -0.119057;
	llsq_err_test = llsq_gsol(A_test, b_test, x_test);
	if ((x_ans - x_test).norm() > TOL * x_ans.norm())
	{
		std::cout << "\nTest llsq_gsol FAILED.\n";
		std::cerr << "Your answer:\n"
				  << x_test << "\n\n"
				  << "Correct answer:\n"
				  << x_ans << "\n\n";
		success = false;
	}

	if (success)
	{
		std::cout << "\n All tests PASSED.\n"
				  << std::endl;
	}
}

int main(int argc, char **argv)
{

	int runTests_flag = 0;

	if (argc == 2)
	{
		runTests_flag = std::stoi(argv[1]);
	}

	if (runTests_flag)
	{
		std::cout << "\nRun tests ...." << std::endl;
		runTests();
	}

	std::cout << "\nRun generalized linear least squares solver ...." << std::endl;
	run_llsq_gsol();
}
