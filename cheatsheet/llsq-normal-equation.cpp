// solves the LLSQ using normal equation method
void llsq_normal(const MatrixXd &A, const VectorXd &b, VectorXd &x)
{

    size_t n = A.cols();

    VectorXd b_ = A.transpose() * b;
    x = (A.transpose() * A).fullPivLu().solve(b_);
    std::cout << "\n\n"
              << x.transpose() << std::endl;
}