functions {
  array[] real ode(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    // Define the state variables
    real Plasma = y[1];
    real KidneyPt = y[2];
    real AccuPt = y[3];
    real DD = y[4];
    real CD = y[5];
    real INF = y[6];
    real KF = y[7];
    real p = y[8];

    real k_kidDD = theta[1];
    real d_DD = theta[2];
    real maxdeath = theta[3];
    real hillDD = theta[4];
    real k_hillnec = theta[5];
    real d_CD = theta[6];
    real h1 = theta[7];
    real k_CDINF = theta[8];
    real d_INF = theta[9];
    real k_INFKF = theta[10];
    real hillKF = theta[11];
    real d_p = theta[12];
    real k_PlasKid = theta[13];
    real k_ePlas = theta[14];
    real k_KidPlas = theta[15];
    real k_eKid = theta[16];
    real k_KidAccu = theta[17];
    real k_AccuKid = theta[18];
    real scale = theta[19];

    /** Below you can fill in the differential equations of the model */
    real dPlasma = KidneyPt * k_KidPlas - Plasma * (k_PlasKid + k_ePlas);
    real dKidneyPt = AccuPt * k_AccuKid + Plasma * k_PlasKid - KidneyPt * (k_KidPlas + k_KidAccu + k_eKid);
    real dAccuPt = KidneyPt * k_KidAccu - AccuPt * k_AccuKid;
    real dDD = (k_kidDD * scale*(KidneyPt + AccuPt))/(hillDD+scale*(KidneyPt + AccuPt)) - DD * d_DD;
    real dCD = (maxdeath*DD^h1/(k_hillnec^h1+DD^h1)) * (1-CD) - d_CD*CD*INF; 
    real dINF = k_CDINF*CD - d_INF*INF;
    real dKF = (k_INFKF*INF^(p+1))/(hillKF^(p+1) + INF^(p+1)) * (1-KF);
    real dp = -d_p * p;

    return {dPlasma, dKidneyPt, dAccuPt, dDD, dCD, dINF, dKF, dp};
  }
}

data {
  real<lower=0> t0;
  int<lower=0> NDD; // Number of time points DD 
  int<lower=0> NINF; // Number of time points INF
  int<lower=0> NCD; // Number of time points CD
  int<lower=0> NKF; // Number of time points KF
  int<lower=0> NTOT; // Total number of timepoints
  array[NDD] real tsdd; // Time points DD
  array[NINF] real tsinf; // Time points INF
  array[NCD] real tscd; // Time points CD 
  array[NKF] real tskf; // Time points KF 
  array[NTOT] real tot; // Total time points
  array[NDD] int idd;
  array[NINF] int iinf;
  array[NCD] int icd;
  array[NKF] int ikf;
  array[NDD, 1] real ydd; // Observed state DD
  array[NINF, 1] real yinf; // Observed state INF
  array[NCD, 1] real ycd; // Observed state CD 
  array[NKF,1] real ykf; // Observed state KF
  array[NDD, 1] real sigmadd; // standard deviations DD -0.5 AT LAST TIMEPOINT
  array[NINF, 1] real sigmainf; // standard deviations INF
  array[NCD, 1] real sigmacd; // standard deviations CD
  array[NKF,1] real sigmakf; // standard deviations KF
  real Plasma0;
  real Kid0;
  real Accu0;
  real DD0;
  real CD0;
  real KF0;
  real<lower=5, upper = 10> p0;
  real<lower=0> d_p;
  real<lower=0> k_PlasKid;
  real<lower=0> k_ePlas;
  real<lower=0> k_KidPlas;
  real<lower=0> k_eKid;
  real<lower=0> k_KidAccu;
  real<lower=0> k_AccuKid;
  real<lower=0> scale;
}

transformed data {
  array[0] real x_r;
  array[0] int x_i;
}

parameters {
  real<lower=0> k_kidDD; 
  real<lower=0> d_DD;
  real<lower=0> maxdeath;
  real<lower=0> hillDD;
  real<lower=0> k_hillnec;
  real<lower=0> d_CD;
  real<lower=1,upper=10> h1;
  real<lower=0> k_CDINF;
  real<lower=0> d_INF;
  real<lower=0> k_INFKF;
  real<lower=0> hillKF;
  real<lower=0, upper = 0.005> INF0;
  //real<lower=0, upper = 15> p0;  //comment if you want to fix the initial hill exponent with decay
  // Degrees of freedom for Student's t-distribution
  real<lower=2> nu_DD;
  real<lower=2> nu_CD;
  real<lower=2> nu_INF;
  real<lower=2> nu_KF;
}

transformed parameters {
  array[NDD] real y_DD;
  array[NINF] real y_INF;
  array[NCD] real y_CD;
  array[NKF] real y_KF;
  array[NTOT,8] real y_hat;
  array[8] real y0 = {Plasma0, Kid0, Accu0, DD0, CD0, INF0, KF0, p0};
  array[19] real theta = {k_kidDD, d_DD, maxdeath, hillDD, k_hillnec, d_CD, h1,
  k_CDINF, d_INF, k_INFKF, hillKF, d_p, k_PlasKid, k_ePlas, k_KidPlas,
  k_eKid, k_KidAccu, k_AccuKid, scale};
  
  // Integrate the ODE system
  y_hat = integrate_ode_rk45(ode, y0, t0, tot, theta, x_r, x_i);
  y_DD = y_hat[idd, 4];
  y_CD = y_hat[icd, 5];
  y_INF = y_hat[iinf, 6];
  y_KF = y_hat[ikf, 7];
}

model {
  // Priors for parameters
  k_kidDD ~ normal(1.644928e-01, 0.5);  
  d_DD ~ normal(5.475403e-02, 0.05);  
  maxdeath ~ normal(1.292922e+01, 5.116602e-01);
  hillDD ~ normal(40.138282925, 2.919770e+00);
  k_hillnec ~ normal(4.078113e+00, 1);
  d_CD ~ normal(6.063535e-03, 0.005);
  h1 ~ normal(8, 1);   
  //k_CDINF ~ normal(2.603806e-01, 0.5);
  k_CDINF ~ normal(0.5, 0.25);
  //d_INF ~ normal(1.255164e-02, 0.05);
  d_INF ~ normal(0.5, 0.25);
  k_INFKF ~ exponential(0.5);
  hillKF ~ exponential(0.5);
  INF0 ~ normal(0.0005, 0.00025);
  //p0 ~ normal(7.5, 1.75);
  // Priors for degrees of freedom
  nu_DD ~ gamma(2, 0.1);
  nu_CD ~ gamma(2, 0.1);
  nu_INF ~ gamma(2, 0.1);
  nu_KF ~ gamma(2, 0.1);

  // Likelihood using Student's t-distribution
  for (n in 1 : NDD) {
    ydd[n, 1] ~ student_t(nu_DD, y_DD[n], sigmadd[n, 1]);
  }
  
  for (n in 1 : NCD) {
    ycd[n, 1] ~ student_t(nu_CD, y_CD[n], sigmacd[n, 1]);
  }
  
  for (n in 1 : NINF) {
    yinf[n, 1] ~ student_t(nu_INF, y_INF[n], sigmainf[n, 1]);
  }
  
  for (n in 1 : NKF) {
    ykf[n, 1] ~ student_t(nu_KF, y_KF[n], sigmakf[n, 1]);
  }
}

generated quantities {
  vector[NDD] logLikelihood_dd;
  vector[NCD] logLikelihood_cd;
  vector[NINF] logLikelihood_inf;
  vector[NKF] logLikelihood_kf;
  vector[NTOT] logLikelihood_total;

  // Compute log-likelihoods for each observation using Student's t-distribution
  for (n in 1 : NDD) {
    logLikelihood_dd[n] = student_t_lpdf(ydd[n, 1] | nu_DD, y_DD[n], sigmadd[n, 1]);
  }
  for (n in 1 : NCD) {
    logLikelihood_cd[n] = student_t_lpdf(ycd[n, 1] | nu_CD, y_CD[n], sigmacd[n, 1]);
  }
  for (n in 1 : NINF) {
    logLikelihood_inf[n] = student_t_lpdf(yinf[n, 1] | nu_INF, y_INF[n], sigmainf[n, 1]);
  }
  for (n in 1 : NKF) {
    logLikelihood_kf[n] = student_t_lpdf(ykf[n, 1] | nu_KF, y_KF[n], sigmakf[n, 1]);
  }

  // Initialize logLikelihood_total with zeros
  for (n in 1 : NTOT) {
    logLikelihood_total[n] = 0;
  }

  // Accumulate log-likelihoods into total likelihood
  for (n in 1 : NDD) {
    logLikelihood_total[idd[n]] += logLikelihood_dd[n];
  }
  for (n in 1 : NCD) {
    logLikelihood_total[icd[n]] += logLikelihood_cd[n];
  }
  for (n in 1 : NINF) {
    logLikelihood_total[iinf[n]] += logLikelihood_inf[n];
  }
  for (n in 1 : NKF) {
    logLikelihood_total[ikf[n]] += logLikelihood_kf[n];
  }
}
