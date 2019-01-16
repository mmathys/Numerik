#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>

using namespace Eigen;
using namespace std;

// routine: llsq_gsol
// " Computes the Generalized linear least squares solution "
// (in)  A: model matrix
// (in)  b: rhs vector
// (out)  x: solution vector
double llsq_gsol(const MatrixXd &A, const VectorXd &b, VectorXd &x)
{
    double llsq_err = 0;

    cout << A << endl;
    cout << b << endl;
    cout << x << endl;

    cout << "SVD..." << endl;
    // SVD
    Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();
    VectorXd sv = svd.singularValues();
    // tol
    double tol = 1e-12;
    sv = sv.unaryExpr([&tol](double e) {
        if (e > tol)
            return e;
        else
            return 0.;
    });
    MatrixXd Sigma = sv.asDiagonal();

    cout << "moore penrose..." << endl;
    // calculate x with moore-penrose pseudo inverse
    x = V * Sigma.inverse() * U.transpose() * b;

    cout << Sigma.inverse() << endl;

    cout << "x = " << endl
         << x << endl;

    cout << "error calculation..." << endl;
    MatrixXd err = (A * x) - b;
    llsq_err = err.norm();

    return llsq_err;
}

// routine: run_llsq_gsol
// " Executes the Generalized LLSQ solver "
void run_llsq_gsol()
{
    double PI = 3.14159265359;
    VectorXd k(5);
    k << 11, 32, 77, 121, 152;
    VectorXd b(5);
    b << 9, 10, 12, 14, 15;
    VectorXd _a = VectorXd::Ones(5);
    VectorXd _b = k.unaryExpr([&PI](int i) {
        return std::sin(2 * PI / 366. * i);
    });
    VectorXd _c = k.unaryExpr([&PI](int i) {
        return std::cos(2 * PI / 366. * i);
    });
    VectorXd _d = k.unaryExpr([&PI](int i) {
        return std::cos(PI / 366. * i) * std::sin(PI / 366. * i);
    });

    MatrixXd A(5, 4);
    A.col(0) = _a;
    A.col(1) = _b;
    A.col(2) = _c;
    A.col(3) = _d;

    VectorXd x;
    double err = llsq_gsol(A, b, x);
    cout << err << endl;
    // TODO: Solve exercise 2.d
}

/***** TESTING ******/
// "" Do NOT CHANGE the routines below ""

void runTests()
{

    bool success = true;
    double TOL = 1e-5;

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
