// Likelihood lower bound for the variational approximation
// the random effects have a diagonal covariance matrix
#include <TMB.hpp>

template<class Type>
    Type objective_function<Type>::operator() ()
{
    using namespace Eigen;
    
    // define the data objects
    DATA_VECTOR(y);
    DATA_MATRIX(X);
    DATA_SPARSE_MATRIX(S);
    DATA_SPARSE_MATRIX(Q);
    DATA_SCALAR (logdetQ); // log determinant of the neighbor matrix.
    DATA_VECTOR(wt);
    DATA_IVECTOR(s_indx);
    DATA_SPARSE_MATRIX(means_mat);
    
    // define the parameter objects
    PARAMETER(ltau);
    //PARAMETER(ltau_A);
    PARAMETER_VECTOR(M);
    PARAMETER_VECTOR(logV);
    PARAMETER_VECTOR(beta);
    
    // define variables
    int r = S.cols();
    int n = X.rows();
    int p = X.cols();
    int q = means_mat.rows();
    vector<Type> eta(n);
    Type floatr = r;
    Type floatq = q;
    
    // The new version (diagonal V)
    Type prod;
    vector<Type> v(n);
    for (int i=0; i<n; i++) {
        if (s_indx(i) >= 0)
            v(i) = exp(logV(s_indx(i)));
    }
    
    // calculate the linear predictor
    eta = X * beta + S * M;
    
    Type quad = M.matrix().transpose() * Q * M.matrix();
    
    vector<Type> A = means_mat * M.matrix();
    
    
    // negative log-lik of the not-Poisson part:
    Type nll = (exp(ltau) * (quad + (Q.diagonal().array() * exp(logV)).sum()) - floatr * ltau - logdetQ) / 2.0; // Add expectation of the log prior on the random effects
    nll -= sum(dnorm(A, 0.0, exp(-ltau/2.0), true)); 
    
    nll -= r/2.0 * (1.0 + log(2.0*PI)) + logV.sum() / 2.0; // Subtract the entropy term
    
    // Poisson part of the likelihood
    for (int i=0; i<n; i++)
        nll -= wt(i) * (y(i)*eta(i) - exp(eta(i) + v(i)/2.0));
    
    return nll;
    }
