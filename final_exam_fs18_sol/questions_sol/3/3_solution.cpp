#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;


// routine: evalF
// " Evaluates the Gaussian distribution function "
// (in)  w: weight
// (in)  mu: mean
// (in)  sigma: variance
// (out) $gaussian(x, mu, sigma)$
double gaussian(const double x, const double mu, const double sigma) {
    return 1./(sqrt(2*M_PI)*sigma)*exp( -(x - mu)*(x - mu)/(2*sigma*sigma) );
}


// routine: evalF
// " Computes $F(x;S)$ "
// (in)  S: randomized samples of weights
// (in)  x: model parameters
// (out) $F(x; S)$
double evalF(const VectorXd& S, const VectorXd& x) {

    int n = x.size();
    int M = S.size();
    double F;
    
    F = 0;
    for (int j=0; j<M; j++) {
        double w = S(j);
        F += log(gaussian(w, x(0), x(1)));
    }
    F *= -1;

  return F;
}


// routine: evalGradF
// " Computes the Gradient of $F(x;S)$ "
// (in)  S: randomized samples of weights
// (in)  x: model parameters
// (out) gradF: Gradient of $F(x;S)$ wrt $x$.
VectorXd evalGradF(const VectorXd& S, const VectorXd& x) {
    
    int n = x.size();
    int M = S.size();
    VectorXd gradF(n);

    gradF.setZero();
    for (int j=0; j<M; j++) {
        double w = S(j);
        gradF(0) += (w - x(0))/(x(1)*x(1));
        gradF(1) += -1./x(1) + (w - x(0))*(w - x(0))/pow(x(1),3);
    }
    gradF *= -1;
    
    return gradF;
}


// routine: evalHessF
// " Computes the Hessian of $F(x;S)$ "
// (in)  S: randomized samples of weights
// (in)  x: model parameters
// (out) hessF: Hessian of $F(x;S)$ wrt $x$.
MatrixXd evalHessF(const VectorXd& S, const VectorXd& x) {
    
    int n = x.size();
    int M = S.size();
    MatrixXd hessF(n,n);
    
    hessF.setZero();
    for (int j=0; j<M; j++) {
        double w = S(j);
        hessF(0,0) += -1./(x(1)*x(1));
        hessF(1,0) += -2*(w - x(0))/pow(x(1),3);
        hessF(0,1) += -2*(w - x(0))/pow(x(1),3);
        hessF(1,1) += 1./(x(1)*x(1)) - 3*(w - x(0))*(w - x(0))/pow(x(1),4);
    }
    hessF *= -1;
        
    return hessF;
}


// routine: newtonOpt
// " Solves the unconstrained minimization problem using Newton method "
// (in) S: randomized samples of weights
// (in)  tol: tolerance
// (in)  maxItr: maximum iterations
// (in/out) x: model parameters
void newtonOpt(const VectorXd& S, const double tol, const int maxItr, VectorXd& x) {

    const int n = x.size();
    
    VectorXd dx;
    VectorXd gradF(n);
    MatrixXd hessF(n,n);
  
    for (int itr=0; itr<maxItr; itr++) {
        
        gradF = evalGradF(S, x);
        hessF = evalHessF(S, x);
    			
        dx = hessF.fullPivLu().solve(gradF);
        std::cout << itr << "\t" << dx.norm() << std::endl;
        if (dx.norm() <= tol)
            break;
    
        x -= dx;
    }
	
}


// routine: runNewtonOpt
// " Executes the Newton Optimization "
// (in) S: randomized samples of weights
void runNewtonOpt(const VectorXd& S) {

    double tol = 1e-8; // tolerance
    int maxItr = 100; // maximum iterations
  
    VectorXd x(2);
    x(0) = 20; x(1) = 3; // initial guess
    newtonOpt(S, tol, maxItr, x);
    std::cout << "Mean and variance:\n" << x.transpose() << std::endl;
}



/***** TESTING ******/
// "" Do NOT CHANGE the routines below ""

void runTests() {
    
    bool success = true;
    double TOL = 1e-8;
    
    int M1 = 8;
    VectorXd S1(M1);
    S1 << -3, -2, -1, 0, 1, 2, 3, 4;
    VectorXd x1(2);
    x1 << 0, 2;
    
    // Test evalGradF
    VectorXd gradF_test(2), gradF_ans(2);
    gradF_ans << -1, -1.5;
    gradF_test = evalGradF(S1, x1);
    if((gradF_ans-gradF_test).norm() > TOL * gradF_ans.norm()){
		std::cout << "\nTest evalGradF FAILED.\n";
		std::cerr << "Your answer:\n" << gradF_test << "\n\n" << "Correct answer:\n" << gradF_ans << "\n\n";
		success = false;
	}
	
	// Test evalHessF
	MatrixXd hessF_test(2,2), hessF_ans(2,2);
    hessF_ans << 2, 1, 1, 6.25;
    hessF_test = evalHessF(S1, x1);
    if((hessF_ans-hessF_test).norm() > TOL * hessF_ans.norm()){
		std::cout << "\nTest evalHessF FAILED.\n";
		std::cerr << "Your answer:\n" << hessF_test << "\n\n" << "Correct answer:\n" << hessF_ans << "\n\n";
		success = false;
	}
	
    // Test newtonOpt
    int M2 = 7;
    VectorXd S2(M2);
    S2 << -3, -2, -1, 0, 1, 2, 3;
    VectorXd x2_test(2), x2_ans(2);
    x2_ans << 0, 2;
    
    int maxItr = 100;    
    x2_test << 0, 0.5;
    newtonOpt(S2, TOL, maxItr, x2_test);
    if((x2_ans-x2_test).norm() > TOL * x2_ans.norm()){
		std::cout << "\nTest newtonOpt FAILED.\n";
		std::cerr << "Your answer:\n" << x2_test << "\n\n" << "Correct answer:\n" << x2_ans << "\n\n";
		success = false;
	}
	
	if (success) {
	    std::cout << "\n All tests PASSED.\n" << std::endl;
	}
    
}


int main(int argc, char** argv) {
    
    int runTests_flag=0;
    
    if (argc==2) {
        runTests_flag = std::stoi(argv[1]);
    }
    
    if (runTests_flag) {
        std::cout << "\nRun tests ...." << std::endl;
        runTests();
    }
    
    std::cout << "\nRun Newton optimization for given data set ...." << std::endl;
    int M = 20;
    VectorXd S(M);
    S << 12.0, 15.2, 16.1, 17.0, 17.8, 18.1, 18.9, 19.0, 19.2, 20.0,\
         20.2, 21.0, 21.4, 22.1, 22.5, 23.5, 23.8, 24.0, 26.0, 28.0;
    
    runNewtonOpt(S);
    
}

// END OF FILE
