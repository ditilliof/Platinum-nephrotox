functions {
  array[] real ode(real t, array[] real y, array[] real theta,
                   array[] real  x_r, array[] int  x_i) {
    real Qapi = y[1];
    real Qbas = y[2];
    real Qcell = y[3];
    real Qinter = y[4];
    real DD = y[5];
    real NEC = y[6];

    real Vapi = theta[1];
    real Vbas = theta[2];
    real Vcell = theta[3];
    real Ncell = theta[4];
    real Fouta = theta[5];
    real Fina = theta[6];
    real Kmet = theta[7];
    real Foutb = theta[8];
    real Finb = theta[9];
    real degDD = theta[10];
    real hillDD = theta[11];
    real k_cDD = theta[12];
    real k_inter = theta[13];
    real maxdeath = theta[14];
    real h = theta[15];
    real k_hillnec = theta[16];
    real p = theta[17];


    real dQapi = Ncell * (Fouta * Qcell/(Ncell*Vcell) - Fina * Qapi/Vapi - Kmet * Qcell/(Ncell*Vcell)); // change in time in amount cisplatin in apical medium in ug/h
    real dQbas = Ncell * (Foutb * Qcell/(Ncell*Vcell) - Finb * Qbas/Vbas); // change in time in amount cisplatin in basolateral medium in ug/h
    real dQcell = Ncell * (Fina * Qapi/Vapi + Finb * Qbas/Vbas - (Fouta + Foutb) * Qcell/(Ncell*Vcell)); // change in amount of cisplatin in cells in ug/h
    real dQinter = k_inter * (Qcell-Qinter);
    real dDD = -degDD * DD + (k_cDD * Qinter)/(hillDD + Qinter);
    real dNEC = ((maxdeath * DD^h)/(k_hillnec^h + DD^h))*(p - NEC); 

    return {dQapi, dQbas, dQcell, dQinter, dDD, dNEC};
  }
}

data {
  int<lower=1> N_DD; // Number of time points DD
  int<lower=1> N_NEC; // Number of time points NEC
  int<lower=1> M;
  real t0;
  array[N_DD] real ts; // time points DD
  array[N_NEC] real ts1; // time points NEC
  array[N_DD + N_NEC] real tot;
  array[N_DD] int<lower=1, upper=100> i_ts;
  array[N_NEC] int<lower=1, upper=100> i_ts1;
  array[N_DD, M] real y_DD; 
  array[N_NEC, M] real y_NEC;
  array[N_DD, M] real<lower=0> sigma_DD;
  array[N_NEC, M] real<lower=0> sigma_NEC;
  real NEC0;
  real<lower=0> Vapi;
  real<lower=0> Vbas;
  real<lower=0> Vcell;
  real<lower=0> Ncell;
  real<lower=0> Fouta;
  real<lower=0> Fina;
  real<lower=0> Kmet;
  real<lower=0> Foutb;
  real<lower=0> Finb;
  real<lower=0> Qapi3;
  real<lower=0> Qbas3;
  real<lower=0> Qapi4;
  real<lower=0> Qbas4;
  real<lower=0> Qapi5;
  real<lower=0> Qbas5;
  real<lower=0> Qapi6;
  real<lower=0> Qbas6;
  real<lower=0> Qapi7;
  real<lower=0> Qbas7;
  real<lower=0> Qapi8;
  real<lower=0> Qbas8;
  real<lower=0> Qapi9;
  real<lower=0> Qbas9;
  real<lower=0> Qcell0;
  real<lower=0, upper=0.5> DD0;
  real<lower=0> Qinter0;
}
transformed data {
  array[0] int x_i;
  array[0] real x_r;
}
parameters {
  real<lower=0> k_cDD;
  real<lower=0> degDD;
  real<lower = 0> hillDD;
  real<lower=0> k_inter;
  real<lower=0> maxdeath;
  real<lower=0> k_hillnec;
  real<lower=0.7, upper=1> p;
  real<lower=1, upper=100> h;
  real<lower=2> nu_DD;  // Degrees of freedom for Student's t-distribution for DD
  real<lower=2> nu_NEC; // Degrees of freedom for Student's t-distribution for NEC
}

transformed parameters {
  array[N_DD + N_NEC, 6] real sol3;
  array[N_DD + N_NEC, 6] real sol4;
  array[N_DD + N_NEC, 6] real sol5;
  array[N_DD + N_NEC, 6] real sol6;
  array[N_DD + N_NEC, 6] real sol7;
  array[N_DD + N_NEC, 6] real sol8;
  array[N_DD + N_NEC, 6] real sol9;
  
  array[N_DD] real z_DD3;
  array[N_NEC] real z_NEC3;
  array[N_DD] real z_DD4;
  array[N_NEC] real z_NEC4;
  array[N_DD] real z_DD5;
  array[N_NEC] real z_NEC5;
  array[N_DD] real z_DD6;
  array[N_NEC] real z_NEC6;
  array[N_DD] real z_DD7;
  array[N_NEC] real z_NEC7;
  array[N_DD] real z_DD8;
  array[N_NEC] real z_NEC8;
  array[N_DD] real z_DD9;
  array[N_NEC] real z_NEC9;
  
  {
     array[17] real theta = {Vapi, Vbas, Vcell, Ncell, Fouta, Fina, Kmet, Foutb,
     Finb, degDD, hillDD, k_cDD, k_inter, maxdeath, h, k_hillnec, p};
     
     array[6] real y3 = {Qapi3, Qbas3, Qcell0, Qinter0, DD0, NEC0}; 
     array[6] real y4 = {Qapi4, Qbas4, Qcell0, Qinter0, DD0, NEC0};
     array[6] real y5 = {Qapi5, Qbas5, Qcell0, Qinter0, DD0, NEC0};
     array[6] real y6 = {Qapi6, Qbas6, Qcell0, Qinter0, DD0, NEC0};
     array[6] real y7 = {Qapi7, Qbas7, Qcell0, Qinter0, DD0, NEC0};
     array[6] real y8 = {Qapi8, Qbas8, Qcell0, Qinter0, DD0, NEC0};
     array[6] real y9 = {Qapi9, Qbas9, Qcell0, Qinter0, DD0, NEC0};
    
    // Solve ODEs
    sol3 = integrate_ode_rk45(ode, y3, t0, tot, theta, x_r, x_i);
    sol4 = integrate_ode_rk45(ode, y4, t0, tot, theta, x_r, x_i);
    sol5 = integrate_ode_rk45(ode, y5, t0, tot, theta, x_r, x_i);
    sol6 = integrate_ode_rk45(ode, y6, t0, tot, theta, x_r, x_i);
    sol7 = integrate_ode_rk45(ode, y7, t0, tot, theta, x_r, x_i);
    sol8 = integrate_ode_rk45(ode, y8, t0, tot, theta, x_r, x_i);
    sol9 = integrate_ode_rk45(ode, y9, t0, tot, theta, x_r, x_i);
    
    // Extract solutions for DD and NEC
    z_DD3 = sol3[i_ts, 5];
    z_DD4 = sol4[i_ts, 5];
    z_DD5 = sol5[i_ts, 5];
    z_DD6 = sol6[i_ts, 5];
    z_DD7 = sol7[i_ts, 5];
    z_DD8 = sol8[i_ts, 5];
    z_DD9 = sol9[i_ts, 5];
    
    z_NEC3 = sol3[i_ts1, 6];
    z_NEC4 = sol4[i_ts1, 6];
    z_NEC5 = sol5[i_ts1, 6];
    z_NEC6 = sol6[i_ts1, 6];
    z_NEC7 = sol7[i_ts1, 6];
    z_NEC8 = sol8[i_ts1, 6];
    z_NEC9 = sol9[i_ts1, 6];
  }
}


model {
  // Priors
 
  k_cDD ~ lognormal(-1.25, 0.5);   
  degDD ~ lognormal(-3.0, 0.5);    
  hillDD ~ lognormal(-6.0, 0.5);   
  k_inter ~ lognormal(-2.5, 0.5);  
  maxdeath ~ normal(1.0, 0.3);     
  k_hillnec ~ normal(5.0, 0.5);    
  p ~ beta(2, 2);                  
  h ~ normal(20, 3);               

  // Priors for degrees of freedom (Student's t-distribution)
  nu_DD ~ exponential(1);          
  nu_NEC ~ exponential(1);        

  // Likelihood using Student's t-distribution
  for (n in 1 : N_DD) {
    y_DD[n, 1] ~ student_t(nu_DD, z_DD3[n], sigma_DD[n, 1]);
    y_DD[n, 2] ~ student_t(nu_DD, z_DD4[n], sigma_DD[n, 2]);
    y_DD[n, 3] ~ student_t(nu_DD, z_DD5[n], sigma_DD[n, 3]);
    y_DD[n, 4] ~ student_t(nu_DD, z_DD6[n], sigma_DD[n, 4]);
    y_DD[n, 5] ~ student_t(nu_DD, z_DD7[n], sigma_DD[n, 5]);
    y_DD[n, 6] ~ student_t(nu_DD, z_DD8[n], sigma_DD[n, 6]);
    y_DD[n, 7] ~ student_t(nu_DD, z_DD9[n], sigma_DD[n, 7]);
  }
  
  for (n in 1 : N_NEC) {
    y_NEC[n, 1] ~ student_t(nu_NEC, z_NEC3[n], sigma_NEC[n, 1]);
    y_NEC[n, 2] ~ student_t(nu_NEC, z_NEC4[n], sigma_NEC[n, 2]);
    y_NEC[n, 3] ~ student_t(nu_NEC, z_NEC5[n], sigma_NEC[n, 3]);
    y_NEC[n, 4] ~ student_t(nu_NEC, z_NEC6[n], sigma_NEC[n, 4]);
    y_NEC[n, 5] ~ student_t(nu_NEC, z_NEC7[n], sigma_NEC[n, 5]);
    y_NEC[n, 6] ~ student_t(nu_NEC, z_NEC8[n], sigma_NEC[n, 6]);
    y_NEC[n, 7] ~ student_t(nu_NEC, z_NEC9[n], sigma_NEC[n, 7]);
  }
}

generated quantities {
  vector[N_DD] logLikelihood_dd3;
  vector[N_NEC] logLikelihood_cd3;
  vector[N_DD] logLikelihood_dd4;
  vector[N_NEC] logLikelihood_cd4;
  vector[N_DD] logLikelihood_dd5;
  vector[N_NEC] logLikelihood_cd5;
  vector[N_DD] logLikelihood_dd6;
  vector[N_NEC] logLikelihood_cd6;
  vector[N_DD] logLikelihood_dd7;
  vector[N_NEC] logLikelihood_cd7;
  vector[N_DD] logLikelihood_dd8;
  vector[N_NEC] logLikelihood_cd8;
  vector[N_DD] logLikelihood_dd9;
  vector[N_NEC] logLikelihood_cd9;
  
  vector[N_DD + N_NEC] logLikelihood_total3;
  vector[N_DD + N_NEC] logLikelihood_total4;
  vector[N_DD + N_NEC] logLikelihood_total5;
  vector[N_DD + N_NEC] logLikelihood_total6;
  vector[N_DD + N_NEC] logLikelihood_total7;
  vector[N_DD + N_NEC] logLikelihood_total8;
  vector[N_DD + N_NEC] logLikelihood_total9;

  // Compute log-likelihoods for each observation
  for (n in 1 : N_DD) {
    logLikelihood_dd3[n] = student_t_lpdf(y_DD[n, 1] | nu_DD, z_DD3[n], sigma_DD[n, 1]);
    logLikelihood_dd4[n] = student_t_lpdf(y_DD[n, 2] | nu_DD, z_DD4[n], sigma_DD[n, 2]);
    logLikelihood_dd5[n] = student_t_lpdf(y_DD[n, 3] | nu_DD, z_DD5[n], sigma_DD[n, 3]);
    logLikelihood_dd6[n] = student_t_lpdf(y_DD[n, 4] | nu_DD, z_DD6[n], sigma_DD[n, 4]);
    logLikelihood_dd7[n] = student_t_lpdf(y_DD[n, 5] | nu_DD, z_DD7[n], sigma_DD[n, 5]);
    logLikelihood_dd8[n] = student_t_lpdf(y_DD[n, 6] | nu_DD, z_DD8[n], sigma_DD[n, 6]);
    logLikelihood_dd9[n] = student_t_lpdf(y_DD[n, 7] | nu_DD, z_DD9[n], sigma_DD[n, 7]);
  }

  for (n in 1 : N_NEC) {
    logLikelihood_cd3[n] = student_t_lpdf(y_NEC[n, 1] | nu_NEC, z_NEC3[n], sigma_NEC[n, 1]);
    logLikelihood_cd4[n] = student_t_lpdf(y_NEC[n, 2] | nu_NEC, z_NEC4[n], sigma_NEC[n, 2]);
    logLikelihood_cd5[n] = student_t_lpdf(y_NEC[n, 3] | nu_NEC, z_NEC5[n], sigma_NEC[n, 3]);
    logLikelihood_cd6[n] = student_t_lpdf(y_NEC[n, 4] | nu_NEC, z_NEC6[n], sigma_NEC[n, 4]);
    logLikelihood_cd7[n] = student_t_lpdf(y_NEC[n, 5] | nu_NEC, z_NEC7[n], sigma_NEC[n, 5]);
    logLikelihood_cd8[n] = student_t_lpdf(y_NEC[n, 6] | nu_NEC, z_NEC8[n], sigma_NEC[n, 6]);
    logLikelihood_cd9[n] = student_t_lpdf(y_NEC[n, 7] | nu_NEC, z_NEC9[n], sigma_NEC[n, 7]);
  }
  
  // Initialize logLikelihood_total with zeros
  for (n in 1 : N_DD + N_NEC) {
    logLikelihood_total3[n] = 0;
    logLikelihood_total4[n] = 0;
    logLikelihood_total5[n] = 0;
    logLikelihood_total6[n] = 0;
    logLikelihood_total7[n] = 0;
    logLikelihood_total8[n] = 0;
    logLikelihood_total9[n] = 0;
  }

  // Sum log-likelihood components
  for (n in 1 : N_DD) {
    logLikelihood_total3[i_ts[n]] += logLikelihood_dd3[n];
    logLikelihood_total4[i_ts[n]] += logLikelihood_dd4[n];
    logLikelihood_total5[i_ts[n]] += logLikelihood_dd5[n];
    logLikelihood_total6[i_ts[n]] += logLikelihood_dd6[n];
    logLikelihood_total7[i_ts[n]] += logLikelihood_dd7[n];
    logLikelihood_total8[i_ts[n]] += logLikelihood_dd8[n];
    logLikelihood_total9[i_ts[n]] += logLikelihood_dd9[n];
  }
  for (n in 1 : N_NEC) {
    logLikelihood_total3[i_ts1[n]] += logLikelihood_cd3[n];
    logLikelihood_total4[i_ts1[n]] += logLikelihood_cd4[n];
    logLikelihood_total5[i_ts1[n]] += logLikelihood_cd5[n];
    logLikelihood_total6[i_ts1[n]] += logLikelihood_cd6[n];
    logLikelihood_total7[i_ts1[n]] += logLikelihood_cd7[n];
    logLikelihood_total8[i_ts1[n]] += logLikelihood_cd8[n];
    logLikelihood_total9[i_ts1[n]] += logLikelihood_cd9[n];
  }
}
