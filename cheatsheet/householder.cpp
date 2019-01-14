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