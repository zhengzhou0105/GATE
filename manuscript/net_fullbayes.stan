data {
  // Overall
  int<lower=1> N;                // Number of chemicals
  
  // Assay counts
  int<lower=0> A;                // Number of in vitro GM assays
  int<lower=0> B;                // Number of in vivo GM assays
  int<lower=0> C;                // Number of in vitro CA assays
  int<lower=0> D;                // Number of in vivo CA assays

  // Individual Data (-1 for missing, 0 for negative, 1 for positive)
  array[N, A] int<lower=-1, upper=1> assay_GM_vt;
  array[N, B] int<lower=-1, upper=1> assay_GM_vv;
  array[N, C] int<lower=-1, upper=1> assay_CA_vt;
  array[N, D] int<lower=-1, upper=1> assay_CA_vv;
  array[N] int<lower=-1, upper=1> RC;   // Individual rodent carcinogenicity
  
  // External Validation Data
  int<lower=0> N_1;       // Number of known RC
  int<lower=0> N_0;       // Number of known non-RC
  array[A + B] int<lower=0, upper=N_1> k_plus_1_GM; 
  array[A + B] int<lower=0, upper=N_0> k_minus_0_GM;
  array[C + D] int<lower=0, upper=N_1> k_plus_1_CA; 
  array[C + D] int<lower=0, upper=N_0> k_minus_0_CA;
  
  // Priors
  vector[2] po_GM;       // Potency of GM
  vector[2] po_CA;       // Potency of CA
  vector[2] prior_RC;    // Prevalence of RC
  vector[2] prior_vt;    // Prior for in vitro assays
  vector[2] prior_vv;    // Prior for in vivo assays
}
transformed data {
  // Aggregate counts for convenience
  int<lower=0> n_GM = A + B;
  int<lower=0> n_CA = C + D;
  
  // Reconstruct unified data matrices
  array[N, n_GM] int<lower=-1, upper=1> assay_GM;
  array[N, n_CA] int<lower=-1, upper=1> assay_CA;
  
  for (i in 1:N) {
    // Merge in vitro (vt) and in vivo (vv) GM assays
    if (A > 0) assay_GM[i, 1:A] = assay_GM_vt[i];
    if (B > 0) assay_GM[i, (A+1):n_GM] = assay_GM_vv[i];
    
    // Merge in vitro (vt) and in vivo (vv) CA assays
    if (C > 0) assay_CA[i, 1:C] = assay_CA_vt[i];
    if (D > 0) assay_CA[i, (C+1):n_CA] = assay_CA_vv[i];
  }
}
parameters {
  // Endpoint potency P(Beta | RC)
  real<lower=0, upper=1> p_GM1_1; // P(GM=1 | RC=1)
  real<lower=0, upper=1> p_GM0_0; // P(GM=0 | RC=0)
  real<lower=0, upper=1> p_CA1_1; // P(CA=1 | RC=1)
  real<lower=0, upper=1> p_CA0_0; // P(CA=0 | RC=0)
  
  // Prior probability of a chemical being an RC (needed for missing RC data)
  real<lower=0, upper=1> pi_RC;  // P(RC=1)
  
  // Assay accuracy (Split by in vitro and in vivo)
  vector<lower=0, upper=1>[A] se_GM_vt; 
  vector<lower=0, upper=1>[B] se_GM_vv; 
  vector<lower=0, upper=1>[A] sp_GM_vt; 
  vector<lower=0, upper=1>[B] sp_GM_vv; 
  
  vector<lower=0, upper=1>[C] se_CA_vt; 
  vector<lower=0, upper=1>[D] se_CA_vv; 
  vector<lower=0, upper=1>[C] sp_CA_vt; 
  vector<lower=0, upper=1>[D] sp_CA_vv;
}
transformed parameters {
  // Aggregate parameter vectors back to original shapes
  vector<lower=0, upper=1>[n_GM] se_GM = append_row(se_GM_vt, se_GM_vv);
  vector<lower=0, upper=1>[n_GM] sp_GM = append_row(sp_GM_vt, sp_GM_vv);
  vector<lower=0, upper=1>[n_CA] se_CA = append_row(se_CA_vt, se_CA_vv);
  vector<lower=0, upper=1>[n_CA] sp_CA = append_row(sp_CA_vt, sp_CA_vv);
  
  // L_ext: empirical anchoring
  // Computed true probabilities
  vector<lower=0,upper=1>[n_GM] p_plus_1_GM;        
  vector<lower=0,upper=1>[n_GM] p_minus_0_GM;
  vector<lower=0,upper=1>[n_CA] p_plus_1_CA;
  vector<lower=0,upper=1>[n_CA] p_minus_0_CA;
  
  // eq:short, Deterministic mapping of accuracy & potency to summary statistics from external validation data
  p_plus_1_GM  = se_GM * p_GM1_1 + (1-sp_GM) * (1-p_GM1_1);
  p_minus_0_GM = (1-se_GM) * (1-p_GM0_0) + sp_GM * p_GM0_0;
  p_plus_1_CA  = se_CA * p_CA1_1 + (1-sp_CA) * (1-p_CA1_1);
  p_minus_0_CA = (1-se_CA) * (1-p_CA0_0) + sp_CA * p_CA0_0;
}
model {
  // Priors
  p_GM1_1 ~ beta(po_GM[1],po_GM[2]);
  p_GM0_0 ~ beta(po_GM[1],po_GM[2]);
  p_CA1_1 ~ beta(po_CA[1],po_CA[2]);
  p_CA0_0 ~ beta(po_CA[1],po_CA[2]);
  pi_RC ~ beta(prior_RC[1],prior_RC[2]);
  se_GM_vt ~ beta(prior_vt[1],prior_vt[2]);
  sp_GM_vt ~ beta(prior_vt[1],prior_vt[2]);
  se_CA_vt ~ beta(prior_vt[1],prior_vt[2]);
  sp_CA_vt ~ beta(prior_vt[1],prior_vt[2]);
  se_GM_vv ~ beta(prior_vv[1],prior_vv[2]);
  sp_GM_vv ~ beta(prior_vv[1],prior_vv[2]);
  se_CA_vv ~ beta(prior_vv[1],prior_vv[2]);
  sp_CA_vv ~ beta(prior_vv[1],prior_vv[2]);

  // Model
  // Part 1: Empirical anchoring on external data
  k_plus_1_GM  ~ binomial(N_1, p_plus_1_GM);
  k_minus_0_GM ~ binomial(N_0, p_minus_0_GM);
  k_plus_1_CA  ~ binomial(N_1, p_plus_1_CA);
  k_minus_0_CA ~ binomial(N_0, p_minus_0_CA);

  // Part 2: Convergence on endpoint by individual data
  for (i in 1:N) {
    // Because HMC does not allow int parameters
    // The latent endpoint values are computed through marginal probabilities
    // Log probability of assay conditional on latent states
    real log_p_GM_1 = 0; // log P(GM Assays | GM=1)
    real log_p_GM_0 = 0; // log P(GM Assays | GM=0)
    real log_p_CA_1 = 0; // log P(CA Assays | CA=1)
    real log_p_CA_0 = 0; // log P(CA Assays | CA=0)
    
    // Bottom up approach: P(A|E,se_A,sp_A)
    // Accumulate likelihood for observed GM assays
    for (j in 1:n_GM) {
      if (assay_GM[i, j] == 1) { // Positive assay result 
        log_p_GM_1 += log(se_GM[j]);  // P(X=+ | GM=1)
        log_p_GM_0 += log(1 - sp_GM[j]); // P(X=+|GM=0)
      } else if (assay_GM[i, j] == 0) { // Negative result
        log_p_GM_1 += log(1 - se_GM[j]); // P(X=-|GM=1)
        log_p_GM_0 += log(sp_GM[j]);     // P(X=-|GM=0)
      } // Missing values (-1) ignored, not contributing
    }
    
    // same for observed CA assays
    for (j in 1:n_CA) {
      if (assay_CA[i, j] == 1) { 
        log_p_CA_1 += log(se_CA[j]);
        log_p_CA_0 += log(1 - sp_CA[j]);
      } else if (assay_CA[i, j] == 0) { 
        log_p_CA_1 += log(1 - se_CA[j]);
        log_p_CA_0 += log(sp_CA[j]);
      }
    }
    
    // Top-down approach: P(RC,E,po_E)
    if (RC[i] == 1) {         // Positive RC
      vector[4] log_probs; // joint probability includes top down and bottom up
      log_probs[1] = log(1 - p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_0 + log_p_CA_0; // State 1: GM=0, CA=0
      log_probs[2] = log(p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_1 + log_p_CA_0; // State 2: GM=1, CA=0
      log_probs[3] = log(1 - p_GM1_1) + log(p_CA1_1) + log_p_GM_0 + log_p_CA_1; // State 3: GM=0, CA=1
      log_probs[4] = log(p_GM1_1) + log(p_CA1_1) + log_p_GM_1 + log_p_CA_1; // State 4: GM=1, CA=1
      
      target += log_sum_exp(log_probs);
      
    } else if (RC[i] == 0) {          // Negative RC
      vector[4] log_probs;
      log_probs[1] = log(p_GM0_0) + log(p_CA0_0) + log_p_GM_0 + log_p_CA_0; // same order as above
      log_probs[2] = log(1 - p_GM0_0) + log(p_CA0_0) + log_p_GM_1 + log_p_CA_0;
      log_probs[3] = log(p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_0 + log_p_CA_1;
      log_probs[4] = log(1 - p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_1 + log_p_CA_1;
      
      target += log_sum_exp(log_probs);  
      
    } else if (RC[i] == -1) {
      // If RC is unobserved, marginalize over RC states
      vector[4] log_probs_rc1;
      vector[4] log_probs_rc0;
      
      // Sample from priors instead
      log_probs_rc1[1] = log(pi_RC) + log(1 - p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_0  + log_p_CA_0;
      log_probs_rc1[2] = log(pi_RC) + log(p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_1 + log_p_CA_0;
      log_probs_rc1[3] = log(pi_RC) + log(1 - p_GM1_1) + log(p_CA1_1) + log_p_GM_0 + log_p_CA_1;
      log_probs_rc1[4] = log(pi_RC) + log(p_GM1_1) + log(p_CA1_1) + log_p_GM_1 + log_p_CA_1;
      
      log_probs_rc0[1] = log(1 - pi_RC) + log(p_GM0_0) + log(p_CA0_0) + log_p_GM_0 + log_p_CA_0;
      log_probs_rc0[2] = log(1 - pi_RC) + log(1 - p_GM0_0) + log(p_CA0_0) + log_p_GM_1 + log_p_CA_0;
      log_probs_rc0[3] = log(1 - pi_RC) + log(p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_0 + log_p_CA_1;
      log_probs_rc0[4] = log(1 - pi_RC) + log(1 - p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_1 + log_p_CA_1;
      
      // Sum probabilities across all 8 possible hidden realities (4 for RC=1, 4 for RC=0)
      target += log_sum_exp( log_sum_exp(log_probs_rc1), log_sum_exp(log_probs_rc0) );
    }
  }
}

generated quantities {
  // Log-likelihood for model evaluation
  vector[N] log_lik;
  // Posterior predictive values of RC for performance check
  vector<lower=0,upper=1>[N] prob_RC_1;           // P(RC=1 | Assays)
  array[N] int<lower=0,upper=1> RC_pred;
  // Posterior predictive values of GM and CA
  vector<lower=0,upper=1>[N] prob_GM;
  vector<lower=0,upper=1>[N] prob_CA;
  array[N] int<lower=0,upper=1> GM_pred;
  array[N] int<lower=0,upper=1> CA_pred;
  
  // Calculate likelihood for model evaluation
  // Use only the individual data
  // Repeat the same individual likelihood calculation in the model block
  for (i in 1:N) {
    real log_p_GM_1 = 0; 
    real log_p_GM_0 = 0; 
    real log_p_CA_1 = 0; 
    real log_p_CA_0 = 0; 
    
    for (j in 1:n_GM) {
      if (assay_GM[i, j] == 1) { 
        log_p_GM_1 += log(se_GM[j]);
        log_p_GM_0 += log(1 - sp_GM[j]);
      } else if (assay_GM[i, j] == 0) { 
        log_p_GM_1 += log(1 - se_GM[j]);
        log_p_GM_0 += log(sp_GM[j]);
      }
    }
    
    for (j in 1:n_CA) {
      if (assay_CA[i, j] == 1) { 
        log_p_CA_1 += log(se_CA[j]);
        log_p_CA_0 += log(1 - sp_CA[j]);
      } else if (assay_CA[i, j] == 0) { 
        log_p_CA_1 += log(1 - se_CA[j]);
        log_p_CA_0 += log(sp_CA[j]);
      }
    }
    
    if (RC[i] == 1) {
      vector[4] log_probs;
      log_probs[1] = log(1 - p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_0 + log_p_CA_0; // State 1: GM=0, CA=0
      log_probs[2] = log(p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_1 + log_p_CA_0; // State 2: GM=1, CA=0
      log_probs[3] = log(1 - p_GM1_1) + log(p_CA1_1) + log_p_GM_0 + log_p_CA_1; // State 3: GM=0, CA=1
      log_probs[4] = log(p_GM1_1) + log(p_CA1_1) + log_p_GM_1 + log_p_CA_1; // State 4: GM=1, CA=1
      
      log_lik[i] = log_sum_exp(log_probs);
      
    } else if (RC[i] == 0) {
      vector[4] log_probs;
      log_probs[1] = log(p_GM0_0) + log(p_CA0_0) + log_p_GM_0 + log_p_CA_0; // same order as above
      log_probs[2] = log(1 - p_GM0_0) + log(p_CA0_0) + log_p_GM_1 + log_p_CA_0;
      log_probs[3] = log(p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_0 + log_p_CA_1;
      log_probs[4] = log(1 - p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_1 + log_p_CA_1;
      
      log_lik[i] = log_sum_exp(log_probs);

    } else if (RC[i] == -1) {
      vector[4] log_probs_rc1;
      vector[4] log_probs_rc0;
      
      log_probs_rc1[1] = log(pi_RC) + log(1 - p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_0  + log_p_CA_0;
      log_probs_rc1[2] = log(pi_RC) + log(p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_1 + log_p_CA_0;
      log_probs_rc1[3] = log(pi_RC) + log(1 - p_GM1_1) + log(p_CA1_1) + log_p_GM_0 + log_p_CA_1;
      log_probs_rc1[4] = log(pi_RC) + log(p_GM1_1) + log(p_CA1_1) + log_p_GM_1 + log_p_CA_1;
      
      log_probs_rc0[1] = log(1 - pi_RC) + log(p_GM0_0) + log(p_CA0_0) + log_p_GM_0 + log_p_CA_0;
      log_probs_rc0[2] = log(1 - pi_RC) + log(1 - p_GM0_0) + log(p_CA0_0) + log_p_GM_1 + log_p_CA_0;
      log_probs_rc0[3] = log(1 - pi_RC) + log(p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_0 + log_p_CA_1;
      log_probs_rc0[4] = log(1 - pi_RC) + log(1 - p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_1 + log_p_CA_1;
      
      // sum over marginalized RC states
      vector[4] log_probs_marginal; // To hold the sum over RC states
      for (k in 1:4) {
        log_probs_marginal[k] = log_sum_exp(log_probs_rc1[k], log_probs_rc0[k]);
      }
      log_lik[i] = log_sum_exp(log_probs_marginal);
    }
    
    // Calculate the joint log-probabilities of (RC, Assays, Latent States)
    // For generating predictive values of RC and endpoints
    vector[4] log_joint_rc1;
    vector[4] log_joint_rc0;
    
    // Holding out RC states, use pi_RC prior instead.
    // RC = 1 | prior, GM, CA
    log_joint_rc1[1] = log(pi_RC) + log(1 - p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_0  + log_p_CA_0; // GM=0 CA=0
    log_joint_rc1[2] = log(pi_RC) + log(p_GM1_1) + log(1 - p_CA1_1) + log_p_GM_1 + log_p_CA_0; // GM=1 CA=0
    log_joint_rc1[3] = log(pi_RC) + log(1 - p_GM1_1) + log(p_CA1_1) + log_p_GM_0 + log_p_CA_1; // GM=0 CA=1
    log_joint_rc1[4] = log(pi_RC) + log(p_GM1_1) + log(p_CA1_1) + log_p_GM_1 + log_p_CA_1; // GM=1 CA=1
    // RC = 0 | prior GM CA
    log_joint_rc0[1] = log(1 - pi_RC) + log(p_GM0_0) + log(p_CA0_0) + log_p_GM_0  + log_p_CA_0; // GM=0 CA=0
    log_joint_rc0[2] = log(1 - pi_RC) + log(1 - p_GM0_0) + log(p_CA0_0) + log_p_GM_1  + log_p_CA_0; // GM=1 CA=0
    log_joint_rc0[3] = log(1 - pi_RC) + log(p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_0 + log_p_CA_1; // GM=0 CA=1
    log_joint_rc0[4] = log(1 - pi_RC) + log(1 - p_GM0_0) + log(1 - p_CA0_0) + log_p_GM_1 + log_p_CA_1; // GM=1 CA=1
    
    // Marginalize over the four states to get the unnormalized log-prob of RC
    real log_marginal_rc1 = log_sum_exp(log_joint_rc1);
    real log_marginal_rc0 = log_sum_exp(log_joint_rc0);
    
    // Normalize to get P(RC=1 | Assays) using bayes theorom
    prob_RC_1[i] = exp(log_marginal_rc1 - log_sum_exp(log_marginal_rc1, log_marginal_rc0));
    // Simulate using bernoulli
    RC_pred[i] = bernoulli_rng(prob_RC_1[i]);
    
    // Marginalize over RC to get the probability of the four states
    vector[4] log_marginal_states;
    for(k in 1:4){
      log_marginal_states[k] = log_sum_exp(log_joint_rc1[k],log_joint_rc0[k]);
    }
    vector[4] prob_states = softmax(log_marginal_states);
    prob_GM[i] = prob_states[2] + prob_states[4];
    prob_CA[i] = prob_states[3] + prob_states[4];
    GM_pred[i] = bernoulli_rng(prob_GM[i]);
    CA_pred[i] = bernoulli_rng(prob_CA[i]);
  }
}

