#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;

//!
//! \brief buildDistanceLSQMatrix Efficiently build the system matrix
//! mapping positions to distances.
//! \param n Number of points (including $x_1$).
//! \return The system matrix $A$.
//!

SparseMatrix<double> buildDistanceLSQMatrix(int n)
{
    SparseMatrix<double> A(n * (n - 1) / 2, n - 1);
    //cout << "n=" << n << endl;
    //cout << "m=" << n * (n - 1) / 2 << endl;

    int c = 0;
    for (int j = 1; j <= n - 1; j++)
    {
        for (int i = j + 1; i <= n; i++)
        {
            int row = c++;
            cout << i << ", " << j << " at row " << row << endl;

            A.insert(row, i - 1 - 1) = 1.;

            if (j - 1 - 1 >= 0)
            {
                A.insert(row, j - 1 - 1) = -1.;
            }
        }
    }

    //cout << "debug A = " << endl;
    //cout << A.toDense() << endl;

    // exercise c.
    //cout << "c)" << endl;
    //cout << A.transpose() * A << endl;

    // TODO efficiently build and return the system matrix A
    A.makeCompressed();
    return A;
}

//!
//! \brief estimatePointsPositions Return positions (without $x_1$).
//! The coordinate $x_1$ is assumed $x_1 = 0$.
//! \param D An $n \times n$ anti-symmetric matrix of distances.
//! \return Vector of positions $x_2, \dots, x_n$.
//!

VectorXd estimatePointsPositions(const MatrixXd &D)
{
    VectorXd x;
    int n = D.cols();

    // TODO compute and return the vector of positions $[x_2, \dots, x_n]$ given
    // the matrix of distances and assuming $x_1 = 0$.
    // Make the implementation as efficient as possible.

    // calculate b.
    VectorXd b(n * (n - 1) / 2);

    int c = 0;
    for (int j = 1; j <= n - 1; j++)
    {
        for (int i = j + 1; i <= n; i++)
        {
            assert(c < b.size());
            cout << "d " << i << "," << j << " = " << D(i - 1, j - 1) << endl;
            b(c) = D(i - 1, j - 1);
            c++;
        }
    }

    //cout << "debug b=" << endl;
    //cout << b << endl;

    // calculate A
    MatrixXd A = buildDistanceLSQMatrix(n);
    VectorXd b2 = A.transpose() * b;
    //cout << "b2 size: " << b2.size() << ", n is " << n << endl;

    VectorXd u = VectorXd::Ones(n - 1);
    VectorXd v = -u;

    SparseMatrix<double> Diag(n - 1, n - 1);
    for (int i = 0; i < n - 1; i++)
    {
        Diag.insert(i, i) = 5;
    }

    SparseLU<SparseMatrix<double>> solver(Diag);

    assert(b2.size() == n - 1);
    VectorXd z = solver.solve(b2);
    assert(u.size() == n - 1);
    VectorXd w = solver.solve(u);
    double alpha = 1.0 + v.dot(w);

    x = (z - w * v.dot(z) / alpha);

    return x;
}

int main()
{

    int n = 5;

    // PART 1: build and print system matrix A
    std::cout << "**** PART 1 ****" << std::endl;
    std::cout << "The Matrix A is:"
              << std::endl
              << buildDistanceLSQMatrix(n)
              << std::endl;

    // PART 2: solve the LSQ system and find positions
    std::cout << "**** PART 2 ****" << std::endl;
    // Vector of positions
    n = 5;

    // Build D
    MatrixXd D(n, n);
    D << 0, -2.1, -3, -4.2, -5,
        2.1, 0, -0.9, -2.2, -3.3,
        3, 0.9, 0, -1.3, -1.1,
        4.2, 2.2, 1.3, 0, -1.1,
        5, 3.3, 1.1, 1.1, 0;
    std::cout << "The matrix D is:"
              << std::endl
              << D
              << std::endl;

    // Find out positions
    VectorXd x_recovered = estimatePointsPositions(D);
    std::cout << "The positions [x_2, ..., x_n] obtained from the LSQ system are:"
              << std::endl
              << x_recovered
              << std::endl;
}
