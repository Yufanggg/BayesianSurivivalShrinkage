//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  // Recursive B-spline basis function
  real bspline_basis(real x, vector knots, int degree, int i) {
    if (degree == 0) {
      return (knots[i] <= x && x < knots[i + 1]) ? 1.0 : 0.0;
    } else {
      real left = (x - knots[i]) / (knots[i + degree] - knots[i]);
      real right = (knots[i + degree + 1] - x) / (knots[i + degree + 1] - knots[i + 1]);
      return left * bspline_basis(x, knots, degree - 1, i) +
             right * bspline_basis(x, knots, degree - 1, i + 1);
    }
  }

  // Function to generate B-spline basis matrix
  matrix generate_bspline_basis(vector x, vector knots, int degree, int intercept) {
    int N = num_elements(x);
    int num_basis = num_elements(knots) + degree - 1;
    matrix[N, num_basis] B;

    for (j in 1:num_basis) {
      for (i in 1:N) {
        B[i, j] = bspline_basis(x[i], knots, degree, j);
      }
    }

    if (intercept) {
      matrix[N, num_basis + 1] B_intercept;
      B_intercept[, 1] = rep_vector(1, N);
      B_intercept[, 2:(num_basis + 1)] = B;
      return B_intercept;
    } else {
      return B;
    }
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;          // Number of data points
  vector[N] x;             // Predictor variable
  int<lower=1> num_knots;  // Number of knots
  vector[num_knots] knots; // Knots for B-splines
  int<lower=0> degree;     // Degree of B-splines
  int<lower=0, upper=1> intercept; // Include intercept
}

transformed data {
  matrix[N, num_knots + degree - 1 + intercept] B;
  B = generate_bspline_basis(x, knots, degree, intercept);
}

parameters {
  vector[num_knots + degree - 1 + intercept] a; // Coefficients for B-splines
  real<lower=0> sigma;     // Standard deviation of the error term
}

model {
  vector[N] y_hat;
  
  // Linear predictor
  y_hat = B * a;
  
  // Likelihood
  y ~ normal(y_hat, sigma);
  
  // Priors
  a ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
}



data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}

