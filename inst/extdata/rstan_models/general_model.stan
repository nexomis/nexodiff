functions {
  real sign(real x) {
    real y;
    y = fabs(x)/x;
    return y;
  }
  
  //Using Skew Cauchy's pdf by retro-engineering the sn package
  real log_skew_cauchy(real x, real mu, real sigma, real alpha) {
    real y;
    real z;
    z = (x-mu)/sigma;
    y = log(1/((pi()*sigma*(1+z^2))))+log(1+alpha*z/sqrt(1+z^2*(1+alpha^2)));
    return y;
  }
}


data {
  int<lower=0> N;
  vector[N] values;  
  
  // ordinal vector of design (group) : 1 for NT VS NT, 2 for NT VS T
  int<lower=0> group_design[N];
    
  // distribution options, 0 for Gaussian, 1 for Cauchy, 2 for Logistic
  int<lower = 0, upper = 2> distribution; 
  
  // options for the model, 0 for False, 1 for True
  int<lower = 0, upper = 1> base_alpha_equals_0; 
  int<lower = 0, upper = 1> alpha_shift_equals_0; 
  int<lower = 0, upper = 1> base_mu_equals_0;
  int<lower = 0, upper = 1> mu_shift_equals_0;
  int<lower = 0, upper = 1> sigma_ratio_equals_1;   
}

parameters {
  real mu_base[1-base_mu_equals_0];
  real mu_shift[1-mu_shift_equals_0];
  
  real alpha_base[(base_alpha_equals_0 == 1 || distribution == 2)?0:1];
  real alpha_shift[(alpha_shift_equals_0 == 1 || distribution == 2)?0:1];
  
  real<lower=0> sigma_base;
  real<lower=0> sigma_ratio[1-sigma_ratio_equals_1];
}

transformed parameters {
  real mu_1[1-base_mu_equals_0];
  real mu_2[1-mu_shift_equals_0]; 
  
  real alpha_1[(base_alpha_equals_0 == 1 || distribution == 2)?0:1];
  real alpha_2[(alpha_shift_equals_0 == 1 || distribution == 2)?0:1];
  
  real<lower=0> sigma_1;
  real<lower=0> sigma_2[1-sigma_ratio_equals_1];
  
  
  if(base_mu_equals_0 == 0 && mu_shift_equals_0 == 0){
    mu_1[1] = mu_base[1];
    mu_2[1] = mu_base[1] + mu_shift[1];
  }
  if(base_mu_equals_0 == 1 && mu_shift_equals_0 == 0)
  {
      mu_2[1] = 0 + mu_shift[1];
  }
  if(base_mu_equals_0 == 0 && mu_shift_equals_0 == 1)
  {
      mu_1[1] = mu_base[1];
  }
  
  if(distribution != 2){
    if(base_alpha_equals_0 == 0 && alpha_shift_equals_0 == 0){
      alpha_1[1] = alpha_base[1];
      alpha_2[1] = alpha_base[1] + alpha_shift[1];
    }
    if(base_alpha_equals_0 == 1 && alpha_shift_equals_0 == 0)
    {
        alpha_2[1] = 0 + alpha_shift[1];
    }
    if(base_alpha_equals_0 == 0 && alpha_shift_equals_0 == 1)
    {
        alpha_1[1] = alpha_base[1];
    }
  }
  
  sigma_1 = sigma_base;
  if(sigma_ratio_equals_1 == 0){
    sigma_2[1] = sigma_base * sigma_ratio[1];
  }
}

model {
  //Priors
  if(base_mu_equals_0 == 0)
  {
    target += normal_lpdf(mu_1[1] | 0, 1);
  }
  if(mu_shift_equals_0 == 0)
  {
    target += normal_lpdf(mu_2[1] | 0, 10);
  }

  target += cauchy_lpdf(sigma_1 | 1,2);
  if(sigma_ratio_equals_1 == 0){
    target += cauchy_lpdf(sigma_2[1] | 1,2);
  }

  if(distribution != 2){
    if(base_alpha_equals_0 == 0)
    {
      target += normal_lpdf(alpha_1[1] | 0, 1);
    }
    if(alpha_shift_equals_0 == 0)
    {
      target += normal_lpdf(alpha_2[1] | 0, 10);
    }
  }


  
  for (i in 1:N)
  {
    if(distribution == 0){
      if(group_design[i] == 1){
        target += skew_normal_lpdf(values[i] | (base_mu_equals_0 == 1)?0:mu_1[1], 
        sigma_1,
        (base_alpha_equals_0 == 1)?0:alpha_1[1]);
      }else{
        target += skew_normal_lpdf(values[i] | (mu_shift_equals_0 == 1)?0:mu_2[1], 
        (sigma_ratio_equals_1 == 1)?1:sigma_2[1],
        (alpha_shift_equals_0 == 1)?0:alpha_2[1]);
      }
    }
    
    if(distribution == 1){

        if(group_design[i] == 1){
          target += log_skew_cauchy(values[i], (base_mu_equals_0 == 1)?0:mu_1[1],
          sigma_1,
          (base_alpha_equals_0 == 1)?0:alpha_1[1]);
        }else{
          target += log_skew_cauchy(values[i], (mu_shift_equals_0 == 1)?0:mu_2[1],
          (sigma_ratio_equals_1 == 1)?1:sigma_2[1],
          (alpha_shift_equals_0 == 1)?0:alpha_2[1]);
        }
    }
    if(distribution == 2){
      if(group_design[i] == 1){
        target += logistic_lpdf(values[i] | (base_mu_equals_0 == 1)?0:mu_1[1], 
        sigma_1);
      }else{
        target += logistic_lpdf(values[i] | (mu_shift_equals_0 == 1)?0:mu_2[1], 
        (sigma_ratio_equals_1 == 1)?1:sigma_2[1]);
      }
    }
  }
}

