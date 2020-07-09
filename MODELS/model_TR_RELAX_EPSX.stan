functions {
  real switch_lock(real t, real t_lock, real r_lock, real shift_lock, real m_lock) {
    return(r_lock+(1-r_lock)/(1+exp(m_lock*(t-t_lock-shift_lock))));
  }
  real switch_relax(real t, real t_relax, real r_relax, real m_relax, real mult, real r_end) {
    return(r_relax+1./(1/(r_end-r_relax)+exp(-m_relax*(t-t_relax-(mult*(1/m_relax)*(1-m_relax)+1/m_relax)))));
  }
  real switch_epsx(real t, real c_epsx, real k_epsx) {
    return(c_epsx*exp(-t/k_epsx));
  }
  // real switch_relax(real t, real t_relax, real r_relax, real s_relax, real m_relax, real r_end) {
    // return(r_relax+1./(1/(r_end-r_relax)+exp(-m_relax*(t-t_relax-s_relax))));
  //}
    real switch_asymp(real t, real K, real mu, real sig) {
      real out;
      real piVal = 3.141593;
      out = K*(1/(t*sig*sqrt(2*piVal)))*exp(-((log(t)-mu)/(sqrt(2)*sig))^2);
      if(is_nan(out)){
        out = 0;
      }
      if(out<0){
        out = 0;
      }
      return(out);
    }
    
    real[] SEIR(real t,
    real[] y,
    real[] theta,
    real[] x_r,
    int[]  x_i 
    ) {
      real tlock_1   = x_r[1]; // lockdown time 1
      real trelax_1  = x_r[2];
      real r_end     = x_r[3];
      real K         = x_r[4];
      real mu        = x_r[5];
      real sig       = x_r[6];
      int  mult      = x_i[1]; //the only fixed integer is the reduction in inf due to being chronic
      real p_lock_1;
      real p_relax_1;
      real p_asymp;
      real dydt[8];
      real init[2]; // inital intercept for S and E
      
      real pi; // number of cases at t0
      real R0; //
      real tau;
      real gamma_s;
      real gamma_H;
      real gamma_ICU;
      real eps_H_ICU;
      // real eps_H_x;
      real eps_ICU_x;
      real r_d_s; // detection rate of symp
      real r_d_a; // detection rate of asymp
      real r_lock_1; // reduction in transmission rate after lockdown
      real shift_lock_1;
      real m_lock_1;
      real r_relax_1; // reduction in transmission rate after lockdown
      real m_relax_1;   
      real k_epsx;
      real coeffR;
      real r_test;
      real eps_ICU_x_mod;
      
      
      // Free parameters
      pi          = theta[1];
      R0          = theta[2];
      tau         = theta[3];
      gamma_s     = theta[4];
      gamma_H     = theta[5];
      gamma_ICU   = theta[6];
      eps_H_ICU   = theta[7];
      // eps_H_x     = theta[8];
      eps_ICU_x   = theta[8];
      r_d_s       = theta[9];
      r_d_a       = theta[10];
      r_lock_1    = theta[11];
      shift_lock_1= theta[12];
      m_lock_1    = theta[13];
      m_relax_1   = theta[14];
      k_epsx      = theta[15];
      
      p_lock_1    = switch_lock(t,tlock_1,r_lock_1,shift_lock_1,m_lock_1);
      p_relax_1   = switch_relax(t,trelax_1,r_lock_1,m_relax_1,mult,r_end);
      p_asymp     = switch_asymp(t,K,mu,sig);
      
      if(p_relax_1>p_lock_1){
        coeffR = p_relax_1;
      }
      else{
        coeffR = p_lock_1;
      }
      
      // time dependent death rate
      eps_ICU_x_mod = switch_epsx(t, eps_ICU_x, k_epsx);
      
      // Initial conditions
      init[1] = 1-pi; // -> sensitives, the actual initial contidition (this is for speed, check below)
      init[2] = pi; // -> exposed
      
      // GENERALIZED SEIR SYSTEM //
      // S
      dydt[1] = - coeffR*(y[1]+init[1])*(R0*gamma_s*y[3]); 
      // E
      dydt[2] = + coeffR*(y[1]+init[1])*(R0*gamma_s*y[3]) - tau*(y[2]+init[2]); 
      // I
      dydt[3] = + tau*(y[2]+init[2]) - gamma_s*y[3];
      // H
      dydt[4] = + r_d_s*gamma_s*y[3] - gamma_H*y[4]; // to ICU, recovery, and death
      // ICU
      dydt[5] = + gamma_H*eps_H_ICU*y[4] - gamma_ICU*y[5]; // death and recovery
      // R -> all recovereds in the population
      // dydt[6] = + gamma_H*(1-eps_H_ICU-eps_H_x)*y[4] + gamma_ICU*(1-eps_ICU_x)*y[5] + (1-r_d_s)*gamma_s*y[3]; 
      dydt[6] = + gamma_H*(1-eps_H_ICU)*y[4] + gamma_ICU*(1-eps_ICU_x_mod)*y[5] + (1-r_d_s)*gamma_s*y[3]; 
      // X -> acc deaths 
      // dydt[7] = + gamma_H*eps_H_x*y[4] + gamma_ICU*eps_ICU_x*y[5]; 
      dydt[7] = + gamma_ICU*eps_ICU_x_mod*y[5]; //assumption : deaths only from ICU
      // DUMMY COMPARTMENTS (FOR COUNTING) //
      // C -> recorded cases (total) - CUMULATIVE
      dydt[8]  = + (r_d_s + p_asymp*r_d_a)*gamma_s*y[3];
      return(dydt);
    }
}

data {
  // data
  int pop_t; // total population
  real tlock_1; 
  real trelax_1;
  real r_end;
  real K;
  real mu;
  real sig;
  int mult;
  int D; // number of days with reported incidence
  int k_daily_cases[D];
  int k_icu[D];
  int k_daily_deaths[D];
  
  // priors
  real p_pi[2];
  real p_R0[2];
  real p_tau[2];
  real p_gamma_s[2];
  real p_gamma_H[2];
  real p_gamma_ICU[2];
  real p_eps_H_ICU[2];
  // real p_eps_H_x[2];
  real p_eps_ICU_x[2];
  real p_r_d_s[2];
  real p_r_d_a[2];
  real p_r_lock_1[2];
  real p_r_mr[2];
  real p_k_epsx[2];
  real p_phi;
  
  // Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S; // total simulation time
  int E; // extra simulation time (S=D+E)
  real ts[D]; // time bins
  real ts_pred[S];
  
}

transformed data {
  real x_r[6]; 
  int x_i[1]; // this is for mult
  real init[8] = rep_array(1e-9,8); // initial values -> this is for speed, keep it like that
  x_r[1] = tlock_1;
  x_r[2] = trelax_1;
  x_r[3] = r_end;
  x_r[4] = K;
  x_r[5] = mu;
  x_r[6] = sig;
  x_i[1] = mult;
}

parameters{
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=1, upper=5> R0;
  real<lower=0, upper=1> tau; 
  real<lower=0, upper=1> gamma_s;  
  real<lower=0, upper=1> gamma_H; 
  real<lower=0, upper=1> gamma_ICU; 
  real<lower=0, upper=1> eps_H_ICU; 
  // real<lower=0, upper=1> eps_H_x; 
  real<lower=0, upper=1> eps_ICU_x;
  real<lower=0, upper=1> r_d_s;
  real<lower=0, upper=1> r_d_a;
  real<lower=0, upper=r_end> r_lock_1;
  real<lower=0, upper=1> m_lock_raw_1; // slope of quarantine implementation
  real<lower=0, upper=1> m_relax_raw_1; // slope of quarantine implementation
  real<lower=0> shift_lock_1; // shift of quarantine implementation
  real<lower=1> phi[3]; // dispersion parameters
  real<lower=0> k_epsx; 
}
transformed parameters {
  real m_lock_1  = m_lock_raw_1+0.1;
  real m_relax_1 = m_relax_raw_1+0.1;
  real theta[15]; // vector of parameters
  real y[D,8]; // raw ODE output
  vector[D] output_ICU;
  vector[D] output_cumC;
  vector[D] output_cumX;
  // outcomes
  vector[D] output_k_icu; 
  vector[D] output_k_daily_deaths; 
  vector[D] output_k_daily_cases;
  
  //   theta = {pi,R0,tau,gamma_s,gamma_H,gamma_ICU,eps_H_ICU,eps_H_x
  // eps_ICU_x,r_d_s,r_d_a,
  // r_lock_1,shift_lock_1,m_lock_1,m_relax_1};
  
  theta = {pi,R0,tau,gamma_s,gamma_H,gamma_ICU,eps_H_ICU,
  eps_ICU_x,r_d_s,r_d_a,
  r_lock_1,shift_lock_1,m_lock_1,m_relax_1,k_epsx};
  
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
      // fitted outputs
      output_k_icu[i] = output_ICU[i]; 
      output_k_daily_cases[i] = i==1 ? output_cumC[i] : 1.0E-9*pop_t + output_cumC[i] - output_cumC[i-1]; 
      output_k_daily_deaths[i] = i==1 ? output_cumX[i] : 1.0E-9*pop_t + output_cumX[i] - output_cumX[i-1]; 
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
  // eps_H_x ~ beta(p_eps_H_x[1],p_eps_H_x[2]);
  k_epsx ~ normal(p_k_epsx[1],p_k_epsx[2]);
  eps_ICU_x ~ beta(p_eps_ICU_x[1],p_eps_ICU_x[2]);
  r_d_s ~ beta(p_r_d_s[1],p_r_d_s[2]);
  r_d_a ~ beta(p_r_d_a[1],p_r_d_a[2]);
  r_lock_1 ~ beta(p_r_lock_1[1],p_r_lock_1[2]);
  phi ~ exponential(p_phi);
  m_lock_raw_1 ~ beta(1,1); 
  shift_lock_1 ~ exponential(1/15.0);
  m_relax_raw_1 ~ beta(p_r_mr[1],p_r_mr[2]); 

  // likelihood
  for(i in 1:D) {
    target += neg_binomial_2_lpmf( k_icu[i] | output_k_icu[i], phi[1]);
    target += neg_binomial_2_lpmf( k_daily_cases[i] | output_k_daily_cases[i], phi[2]);
    target += neg_binomial_2_lpmf( k_daily_deaths[i] | output_k_daily_deaths[i], phi[3]);
  }
}

generated quantities{
  
  real y_pred[S,8]; // raw ODE output
  
  vector[S] comp_S;
  vector[S] comp_E;
  vector[S] comp_Is;
  vector[S] comp_H;
  vector[S] comp_ICU;
  vector[S] comp_R;
  vector[S] comp_X;
  vector[S] comp_C;
  vector[S] comp_diffX;
  vector[S] comp_diffC;
  vector[S] RE_raw;
  vector[S] asymp_detect;
  
  real p_lock_1;
  real p_relax_1;
  real p_asymp;
  real coeffR;
  
  int fitted_k_daily_cases[D];
  int fitted_k_icu[D];
  int fitted_k_daily_deaths[D];
  
  int predicted_k_daily_cases[E];
  int predicted_k_icu[E];
  int predicted_k_daily_deaths[E];
  
  real predicted_daily_cases[S];
  real predicted_cum_cases[S];
  real predicted_current_hospit[S];
  real predicted_current_icu[S];
  real predicted_daily_deaths[S];
  real predicted_cum_deaths[S];
  
  
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
      // lagged differences of cumulative compartments (daily deaths and daily cases)
      comp_diffX[i] = i==1 ? comp_X[i] : 1.0E-9*pop_t + comp_X[i] - comp_X[i-1];
      comp_diffC[i] = i==1 ? comp_C[i] : 1.0E-9*pop_t + comp_C[i] - comp_C[i-1];
    }
    
    for(i in 1:D) {
      fitted_k_icu[i]            = neg_binomial_2_rng( output_k_icu[i], phi[1]);
      fitted_k_daily_cases[i]    = neg_binomial_2_rng( output_k_daily_cases[i], phi[2]);
      fitted_k_daily_deaths[i]   = neg_binomial_2_rng( output_k_daily_deaths[i], phi[3]);
    }
    for(i in 1:E) {
      predicted_k_icu[i]           = neg_binomial_2_rng( comp_ICU[D+i], phi[1]);
      predicted_k_daily_cases[i]   = neg_binomial_2_rng( comp_diffC[D+i], phi[2]);
      predicted_k_daily_deaths[i]  = neg_binomial_2_rng( comp_diffX[D+i], phi[3]);
    }
    for(i in 1:S) {
      predicted_current_hospit[i] = comp_H[i];
      predicted_current_icu[i]    = comp_ICU[i];
      predicted_daily_deaths[i]   = comp_diffX[i];
      predicted_cum_deaths[i]     = comp_X[i];
      predicted_daily_cases[i]    = comp_diffC[i];
      predicted_cum_cases[i]      = comp_C[i];
      
      p_lock_1  = switch_lock(i,tlock_1,r_lock_1,shift_lock_1,m_lock_1);
      p_relax_1 = switch_relax(i,trelax_1,r_lock_1,m_relax_1,mult,r_end);
      p_asymp   = switch_asymp(i,K,mu,sig);
      
      asymp_detect[i] = p_asymp*r_d_a;
      
      if(p_relax_1>p_lock_1){
        coeffR = p_relax_1;
      }
      else{
        coeffR = p_lock_1;
      }
      RE_raw[i] = coeffR*R0*(comp_S[i]/comp_S[1]);
    }
    
}
