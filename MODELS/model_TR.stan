functions {
  real switch_lock(real t, real t_lock, real r_lock, real shift_lock, real m_lock) {
    return(r_lock+(1-r_lock)/(1+exp(m_lock*(t-t_lock-shift_lock))));
  }
  real switch_relax(real t, real t_relax, real r_relax, real r_end, real m_relax) {
    return(r_relax+1./(1/(r_end-r_relax)+exp(-m_relax*(t-t_relax-(6*(1/m_relax)*(1-m_relax)+1/m_relax)))));
  }
  real switch_asymp(real t, real t_test, real r_test, real shift_test, real m_test) {
    return(r_test+1./(1/(1-r_test)+exp(-m_test*(t-t_test-shift_test))));
  }
  // THIS IS FOR POLYNOMIAL FIT
  // real switch_asymp(real t, real p1, real p2, real p3, real p4) {
  //   real out; 
  //   out = p4*t^3 + p3*t^2 + p2*t + p1;
  //   if(out<0){
  //     out=0;
  //   }
  //   return(out);
  // }
  
  real[] SEIR(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[]  x_i 
  ) {
    real tswitch_1 = x_r[1]; // lockdown time 1
    real trelax    = x_r[2];
    real r_end     = x_r[3];
    real m_relax   = x_r[4]; //
    real tasymp    = x_r[5];
    // real p1        = x_r[5];
    // real p2        = x_r[6];
    // real p3        = x_r[7];
    // real p4        = x_r[8];
    int  r_c       = x_i[1]; //the only fixed integer is the reduction in inf due to being chronic
    real p_tswitch_1;
    real p_relax;
    real p_asymp;
    real dydt[11];
    real init[2]; // inital intercept for S and E
    
    real pi; // number of cases at t0
    real R0; //
    real tau;
    real gamma_s;
    real gamma_H;
    real gamma_ICU;
    real eps_H_ICU;
    real eps_H_x;
    real eps_ICU_x;
    real eps_H;
    real r_d_s; // detection rate of symp
    real r_d_a; // detection rate of asymp
    real r_r_s; // recovery rate of asymp
    real r_r_a; // recovery rate of symp
    real r_lock_1; // reduction in transmission rate after lockdown
    real r_test;
    real shift_lock_1;
    real shift_test;
    real m_lock_1;
    real m_test;
    real coeffR;
    
    // Free parameters
    pi          = theta[1];
    R0          = theta[2];
    tau         = theta[3];
    gamma_s     = theta[4];
    gamma_H     = theta[5];
    gamma_ICU   = theta[6];
    eps_H_ICU   = theta[7];
    eps_H_x     = theta[8];
    eps_ICU_x   = theta[9];
    eps_H       = theta[10];
    r_d_s       = theta[11];
    r_d_a       = theta[12];
    r_r_s       = theta[13];
    r_r_a       = theta[14];
    r_lock_1    = theta[15];
    r_test      = theta[16];
    shift_lock_1= theta[17];
    shift_test  = theta[18];
    m_lock_1    = theta[19];
    m_test      = theta[20];

    p_tswitch_1 = switch_lock(t,tswitch_1,r_lock_1,shift_lock_1,m_lock_1);
    p_asymp     = switch_asymp(t,tasymp,r_test,shift_test,m_test);
    
    if(m_relax>0){ //m_relax=0 no relaxation
    p_relax   = switch_relax(t,trelax,r_lock_1,r_end,m_relax);
    if(p_tswitch_1>p_relax){
      coeffR = p_tswitch_1;
    }
    else{
      coeffR = p_relax;
    }
    }else{
      p_relax= 0;
      coeffR = p_tswitch_1;
    }
    
    // Initial conditions
    init[1] = 1-pi; // -> sensitives, the actual initial contidition (this is for speed, check below)
    init[2] = pi; // -> exposed
    // SEIR SYSTEM //
    
    // S
    dydt[1] = - coeffR*(y[1]+init[1])*(R0*gamma_s*y[3]); 
    // E
    dydt[2] = + coeffR*(y[1]+init[1])*(R0*gamma_s*y[3]) - tau*(y[2]+init[2]); 
    // I
    dydt[3] = + tau*(y[2]+init[2]) - gamma_s*y[3];
    // H
    dydt[4] = + eps_H*r_d_s*gamma_s*y[3] - gamma_H*y[4]; // to ICU, recovery, and death
    // ICU
    dydt[5] = + gamma_H*eps_H_ICU*y[4] - gamma_ICU*y[5]; // death and recovery
    // R -> all recovereds in the population
    dydt[6] = + gamma_H*(1-eps_H_ICU-eps_H_x)*y[4] + gamma_ICU*(1-eps_ICU_x)*y[5] + (1-eps_H*r_d_s)*gamma_s*y[3]; 
    // X -> acc deaths 
    dydt[7] = + gamma_H*eps_H_x*y[4] + gamma_ICU*eps_ICU_x*y[5]; //assumption : deaths only from ICU
    // DUMMY COMPARTMENTS (FOR COUNTING) //
    // C -> recorded cases (total) - CUMULATIVE
    dydt[8]  = + (r_d_s + p_asymp*r_d_a)*gamma_s*y[3];
    // C1 -> recorded cases (out of hospital, symptomatic) - CUMULATIVE
    dydt[9]  = + r_d_s*(1-eps_H)*gamma_s*y[3];
    // C2 -> recorded cases (out of hospital, asymptomatic) - CUMULATIVE
    dydt[10] = + p_asymp*r_d_a*gamma_s*y[3];
    // RC -> recorded recovereds : from hospital + from icu + not hosp (therefore not icu) but detected - CUMULATIVE
    dydt[11] = + gamma_H*(1-eps_H_ICU-eps_H_x)*y[4] + gamma_ICU*(1-eps_ICU_x)*y[5] + r_r_s*y[9] + r_r_a*y[10]; //
    
    return(dydt);
  }
}

data {
  // data
  int pop_t; // total population
  real tswitch_1; 
  real trelax;
  real r_end;
  real m_relax; // real gamma_c;
  real tasymp;
  int D; // number of days with reported incidence
  int k_daily_cases[D];
  int k_icu[D];
  int k_daily_deaths[D];
  int k_daily_recov[D];

  // priors
  real p_pi[2];
  real p_R0[2];
  real p_tau[2];
  real p_gamma_s[2];
  real p_gamma_H[2];
  real p_gamma_ICU[2];
  real p_eps_H_ICU[2];
  real p_eps_H_x[2];
  real p_eps_ICU_x[2];
  real p_eps_H[2];
  real p_r_d_s[2];
  real p_r_d_a[2];
  real p_r_r_s[2];
  real p_r_r_a[2];
  real p_r_lock_1[2];
  real p_r_test[2];
  real p_phi;
  
  // Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S; // total simulation time
  int E; // extra simulation time (S=D+E)
  real ts[D]; // time bins
  real ts_pred[S];
  
  // controls
  int r_c;
}

transformed data {
  real x_r[5]; 
  int x_i[1]; // this is for r_c
  real init[11] = rep_array(1e-9,11); // initial values -> this is for speed, keep it like that
  x_r[1] = tswitch_1;
  x_r[2] = trelax;
  x_r[3] = r_end;
  x_r[4] = m_relax;// x_r[6] = gamma_c;
  x_r[5] = tasymp;
  // x_r[5] = p1;
  // x_r[6] = p2;
  // x_r[7] = p3;
  // x_r[8] = p4;
  x_i[1] = r_c;
}

parameters{
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=1, upper=5> R0;
  real<lower=0, upper=1> tau; 
  real<lower=0, upper=1> gamma_s;  
  real<lower=0, upper=1> gamma_H; 
  real<lower=0, upper=1> gamma_ICU; 
  real<lower=0, upper=1> eps_H_ICU; 
  real<lower=0, upper=1> eps_H_x; 
  real<lower=0, upper=1> eps_ICU_x;
  real<lower=0, upper=1> eps_H;
  real<lower=0, upper=0.5> r_d_s;
  real<lower=0, upper=0.3> r_d_a;
  real<lower=0, upper=1> r_r_s;
  real<lower=0, upper=1> r_r_a;
  real<lower=0, upper=1> r_lock_1;
  real<lower=0, upper=1> r_test;
  real<lower=0, upper=1> m_lock_raw_1; // slope of quarantine implementation
  real<lower=0, upper=1> m_test_raw; // slope of test implementation
  real<lower=0> shift_lock_1; // shift of quarantine implementation
  real<lower=0> shift_test; // shift of test implementation
  real<lower=1> phi[4]; // dispersion parameters
}
transformed parameters {
  real m_lock_1 = m_lock_raw_1+0.5;
  real m_test   = m_test_raw+0.05; // not sure if necessary
  real theta[20]; // vector of parameters
  real y[D,11]; // raw ODE output
  vector[D] output_ICU;
  vector[D] output_cumC;
  vector[D] output_cumX;
  vector[D] output_cumRH;
  // outcomes
  vector[D] output_k_icu; 
  vector[D] output_k_daily_deaths; 
  vector[D] output_k_daily_cases;
  vector[D] output_k_daily_recov;
  
  theta = {pi,R0,tau,gamma_s,gamma_H,gamma_ICU,eps_H_ICU,eps_H_x,
  eps_ICU_x,eps_H,r_d_s,r_d_a,r_r_s,r_r_a,
  r_lock_1,r_test,shift_lock_1,shift_test,
  m_lock_1,m_test};
  // run ODE solver
  y = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances 
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(i in 1:D) {
      output_ICU[i]  = (y[i,5]+1.0E-9)*pop_t;
      output_cumC[i] = (y[i,8]+1.0E-9)*pop_t;
      output_cumX[i] = (y[i,7]+1.0E-9)*pop_t;
      output_cumRH[i]= (y[i,11]+1.0E-9)*pop_t;
      // fitted outputs
      output_k_icu[i] = output_ICU[i]; 
      output_k_daily_cases[i] = i==1 ? output_cumC[i] : 1.0E-9*pop_t + output_cumC[i] - output_cumC[i-1]; 
      output_k_daily_deaths[i] = i==1 ? output_cumX[i] : 1.0E-9*pop_t + output_cumX[i] - output_cumX[i-1]; 
      output_k_daily_recov[i] = i==1 ? output_cumRH[i] : 1.0E-9*pop_t + output_cumRH[i] - output_cumRH[i-1]; 
    }
}

model {
  // priors
  pi ~ beta(p_pi[1],p_pi[2]);// R0 ~ uniform(p_R0[1],p_R0[2]);
  R0 ~ normal(p_R0[1],p_R0[2]);
  tau~ beta(p_tau[1],p_tau[2]);
  gamma_s ~ beta(p_gamma_s[1],p_gamma_s[2]);
  gamma_H ~ beta(p_gamma_H[1],p_gamma_H[2]);
  gamma_ICU ~ beta(p_gamma_ICU[1],p_gamma_ICU[2]);
  eps_H_ICU ~ beta(p_eps_H_ICU[1],p_eps_H_ICU[2]);
  eps_H_x ~ beta(p_eps_H_x[1],p_eps_H_x[2]);
  eps_ICU_x ~ beta(p_eps_ICU_x[1],p_eps_ICU_x[2]);
  eps_H ~ beta(p_eps_H[1],p_eps_H[2]);
  r_d_s ~ beta(p_r_d_s[1],p_r_d_s[2]);
  r_d_a ~ beta(p_r_d_a[1],p_r_d_a[2]);
  r_r_s ~ beta(p_r_d_s[1],p_r_d_s[2]);
  r_r_a ~ beta(p_r_r_a[1],p_r_r_a[2]);
  r_lock_1 ~ beta(p_r_lock_1[1],p_r_lock_1[2]);
  phi ~ exponential(p_phi);
  m_lock_raw_1 ~ beta(1,1); 
  m_test_raw ~ beta(1,1); 
  shift_lock_1 ~ exponential(1/15.0);
  shift_test ~ exponential(1/15.0);
  // likelihood
  for(i in 1:D) {
    target += neg_binomial_2_lpmf( k_icu[i] | output_k_icu[i], phi[1]);
    target += neg_binomial_2_lpmf( k_daily_cases[i] | output_k_daily_cases[i], phi[2]);
    target += neg_binomial_2_lpmf( k_daily_deaths[i] | output_k_daily_deaths[i], phi[3]);
    target += neg_binomial_2_lpmf( k_daily_recov[i] | output_k_daily_recov[i], phi[4]);
  }
}

generated quantities{
  
  real y_pred[S,11]; // raw ODE output
  
  vector[S] comp_S;
  vector[S] comp_E;
  vector[S] comp_Is;
  vector[S] comp_H;
  vector[S] comp_ICU;
  vector[S] comp_R;
  vector[S] comp_X;
  vector[S] comp_C;
  vector[S] comp_RH;
  vector[S] comp_diffX;
  vector[S] comp_diffC;
  vector[S] comp_diffRH;
  vector[S] RE_raw;
  vector[S] asymp_detect;
  
  real p_tswitch_1;
  real p_asymp;
  real p_relax;
  real coeffR;
  
  int fitted_k_daily_cases[D];
  int fitted_k_icu[D];
  int fitted_k_daily_deaths[D];
  int fitted_k_daily_recov[D];
  
  int predicted_k_daily_cases[E];
  int predicted_k_icu[E];
  int predicted_k_daily_deaths[E];
  int predicted_k_daily_recov[E];
  
  real predicted_daily_cases[S];
  real predicted_cum_cases[S];
  real predicted_current_hospit[S];
  real predicted_current_icu[S];
  real predicted_daily_deaths[S];
  real predicted_cum_deaths[S];
  real predicted_cum_recov[S];
  real predicted_daily_recov[S];
  
  
  y_pred = integrate_ode_bdf(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts_pred, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3);
    
    for(i in 1:S) {
      
      comp_S[i]  = (y_pred[i,1] + 1 - pi)*pop_t;
      comp_E[i]  = (y_pred[i,2] + pi)*pop_t;
      comp_Is[i] = (y_pred[i,3]+1.0E-9)*pop_t;
      comp_H[i]  = (y_pred[i,4]+1.0E-9)*pop_t;
      comp_ICU[i]= (y_pred[i,5]+1.0E-9)*pop_t;
      comp_R[i]  = (y_pred[i,6]+1.0E-9)*pop_t;
      comp_X[i]  = (y_pred[i,7]+1.0E-9)*pop_t;
      comp_C[i]  = (y_pred[i,8]+1.0E-9)*pop_t;
      comp_RH[i] =  (y_pred[i,11]+1.0E-9)*pop_t;
      // lagged differences of cumulative compartments (daily deaths and daily cases)
      comp_diffX[i] = i==1 ? comp_X[i] : 1.0E-9*pop_t + comp_X[i] - comp_X[i-1];
      comp_diffC[i] = i==1 ? comp_C[i] : 1.0E-9*pop_t + comp_C[i] - comp_C[i-1];
      comp_diffRH[i] = i==1 ? comp_RH[i] : 1.0E-9*pop_t + comp_RH[i] - comp_RH[i-1];
    }
    
    for(i in 1:D) {
      fitted_k_icu[i]            = neg_binomial_2_rng( output_k_icu[i], phi[1]);
      fitted_k_daily_cases[i]    = neg_binomial_2_rng( output_k_daily_cases[i], phi[2]);
      fitted_k_daily_deaths[i]   = neg_binomial_2_rng( output_k_daily_deaths[i], phi[3]);
      fitted_k_daily_recov[i]    = neg_binomial_2_rng( output_k_daily_recov[i], phi[4]);
    }
    for(i in 1:E) {
      predicted_k_icu[i]           = neg_binomial_2_rng( comp_ICU[D+i], phi[1]);
      predicted_k_daily_cases[i]   = neg_binomial_2_rng( comp_diffC[D+i], phi[2]);
      predicted_k_daily_deaths[i]  = neg_binomial_2_rng( comp_diffX[D+i], phi[3]);
      predicted_k_daily_recov[i]   = neg_binomial_2_rng( comp_diffRH[D+i], phi[4]);
    }
    for(i in 1:S) {
      predicted_current_hospit[i] = comp_H[i];
      predicted_current_icu[i]    = comp_ICU[i];
      predicted_daily_deaths[i]   = comp_diffX[i];
      predicted_cum_deaths[i]     = comp_X[i];
      predicted_daily_cases[i]    = comp_diffC[i];
      predicted_cum_cases[i]      = comp_C[i];
      predicted_cum_recov[i]      = comp_RH[i];
      predicted_daily_recov[i]    = comp_diffRH[i];
      
      p_tswitch_1 = switch_lock(i,tswitch_1,r_lock_1,shift_lock_1,m_lock_1);
      p_asymp     = switch_asymp(i,tasymp,r_test,shift_test,m_test);
      
      asymp_detect[i] = p_asymp*r_d_a;
      
      if(m_relax>0){ //m_relax=0 no relaxation
      p_relax   = switch_relax(i,trelax,r_lock_1,r_end,m_relax);
      if(p_tswitch_1>p_relax){
        coeffR = p_tswitch_1;
      }
      else{
        coeffR = p_relax;
      }
      }else{
        p_relax= 0;
        coeffR = p_tswitch_1;
      }
      RE_raw[i] = coeffR*R0*(comp_S[i]/comp_S[1]);
    }
    
}
