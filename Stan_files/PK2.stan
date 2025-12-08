functions {
  array[] real ode(real t, array[] real y, array[] real theta, array[] real x_r, array[] int x_i) {
    // Define the state variables
    real Plasma = y[1];
    real KidneyPt = y[2];
    
    real k_PlasKid = theta[1];
    real k_ePlas = theta[2];
    real k_KidPlas = theta[3];
    real k_eKid = theta[4];
    
    // Differential equations of the model
    real dPlasma = KidneyPt * k_KidPlas - Plasma * (k_PlasKid + k_ePlas);
    real dKidneyPt = Plasma * k_PlasKid - KidneyPt * (k_KidPlas + k_eKid);

    // Return the derivatives
    return {dPlasma, dKidneyPt};
  }
}

data {
  int<lower=1> N1;                // Number of time points for first dataset
  int<lower=1> N2;                // Number of time points for second dataset
  real<lower=0> t0;               // Initial time
  array[N1] real ts1;             // Time points for first dataset
  array[N2] real ts2;             // Time points for second dataset
  array[N1,1] real y1;              // Observed state for first dataset
  array[N2,1] real y2;              // Observed state for second dataset
  array[N1,1] real sigma1;          // Known standard deviations of the observations for first dataset
  array[N2,1] real sigma2;          // Known standard deviations of the observations for second dataset
  real<lower=0> KidneyPt0;        // Initial state value for KidneyPt
}

transformed data {
  array[0] real x_r;
  array[0] int x_i;
}

parameters {
  real<lower=0> k_PlasKid;
  real<lower=0> k_ePlas;
  real<lower=0> k_KidPlas;
  real<lower=0> k_eKid;
  real<lower=0> Plasma0;
  real<lower=0> scale;
}

transformed parameters {
  array[N1, 2] real z_hat1;
  array[N2, 2] real z_hat2;
  array[2] real y0 = {Plasma0, KidneyPt0};
  array[4] real theta = {k_PlasKid, k_ePlas, k_KidPlas, k_eKid};
  z_hat1 = integrate_ode_rk45(ode, y0, t0, ts1, theta, x_r, x_i);
  z_hat2 = integrate_ode_rk45(ode, y0, t0, ts2, theta, x_r, x_i);
}

model {
  // Priors
  k_PlasKid ~ exponential(1);
  k_ePlas ~ exponential(1);
  k_KidPlas ~ exponential(1);
  k_eKid ~ exponential(1);
  Plasma0 ~ normal(50,1);
  scale ~ normal(4,0.25);

  // Likelihood
  for (n in 1:N1) {
    y1[n] ~ normal(z_hat1[n, 1], sigma1[n]);
  }

  for (n in 1:N2) {
    real z_hat2_scaled = z_hat2[n, 1] * scale;
    y2[n] ~ normal(z_hat2_scaled, sigma2[n]);
  }
}




