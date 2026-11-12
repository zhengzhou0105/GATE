# Settings--------

# Load packages----
library(cmdstanr)
library(readr)
library(tidyverse)
library(readxl)
library(bayesplot)
library(loo)
library(pROC)
library(tidyverse)
library(readr)
library(pROC)
library(BayesLCA)
library(extraDistr)

# Functions----

# Fit a network using bnlearn
dtm_bn <- function(K,p_prior){
  # configure
  p_X <- rep(p_prior, K)          # Prior probability P(Xk = 1) for each Xk
  
  # node names X and Y
  X_names <- paste0("X", 1:K)
  nodes <- c(X_names, "Y")
  
  # network
  parents_str <- paste(X_names, collapse = ":")
  model_str <- paste0("[", paste(X_names, collapse = "]["),
                      "][Y|", parents_str, "]")
  dag <- model2network(model_str)
  
  # define CPT
  
  
  # Priors for X1 to XK
  cpds <- list()
  for (i in seq_len(K)) {
    xname <- X_names[i]
    cpds[[xname]] <- matrix(
      c(1 - p_X[i], p_X[i]), ncol = 2,
      dimnames = list(NULL, c("0", "1"))
    )
  }
  
  # CPT for Y (deterministic)
  Y_dimnames <- c(list(Y = c("0", "1")), lapply(X_names, function(n) c("0", "1")))
  names(Y_dimnames)[-1] <- X_names
  Y_table <- expand.grid(rev(Y_dimnames[-1])) # all X combinations
  Y_vals <- apply(Y_table, 1, function(row) as.numeric(any(row == "1")))
  p_Y <- rbind(1 - Y_vals, Y_vals)  # rows: Y=0, Y=1
  # Now reshape and name
  dim(p_Y) <- c(2, rep(2, K))  # 2 for Y, 2^K for parents
  dimnames(p_Y) <- Y_dimnames
  cpds[["Y"]] <- p_Y
  
  # Fit network
  fitted <- custom.fit(dag, cpds)
  
  return(fitted)
}

# Sample from a fitted bayesian network
bn_sample <-function(fitted,n_samples){
  # Simulate data
  samples <- rbn(fitted, n = n_samples)
  
  # Summarize Y
  freq_Y <- as.data.frame(table(samples$Y))
  colnames(freq_Y) <- c("Y", "Count")
  freq_Y$Proportion <- freq_Y$Count / sum(freq_Y$Count)
  return(freq_Y[2,3])  # Return proportion of Y=1
}

# Compute posterior probability based on combination

point_prob_GM_CA <- function(combo_GM, combo_CA,
                             p_GM, p_CA,
                             se_GM, sp_GM,
                             se_CA, sp_CA) {
  # combo_GM : vector of length n_GM with 0,1,NA
  # combo_CA : vector of length n_CA with 0,1,NA
  n_GM <- length(se_GM)
  n_CA <- length(ss_CA)
  # All other arguments are scalars/vectors of posterior means
  
  # Likelihoods for the four latent states
  L11 <- 1; L10 <- 1; L01 <- 1; L00 <- 1
  
  # GM assays under GM=1 (for L11 and L10)
  for (j in 1:n_GM) {
    if (!is.na(combo_GM[j])) {
      if (combo_GM[j] == 1) {
        L11 <- L11 * se_GM[j]
        L10 <- L10 * se_GM[j]
      } else {
        L11 <- L11 * (1 - se_GM[j])
        L10 <- L10 * (1 - se_GM[j])
      }
    }
  }
  # GM assays under GM=0 (for L01 and L00)
  for (j in 1:n_GM) {
    if (!is.na(combo_GM[j])) {
      if (combo_GM[j] == 1) {
        L01 <- L01 * (1 - sp_GM[j])
        L00 <- L00 * (1 - sp_GM[j])
      } else {
        L01 <- L01 * sp_GM[j]
        L00 <- L00 * sp_GM[j]
      }
    }
  }
  # CA assays under CA=1 (for L11 and L01)
  for (j in 1:n_CA) {
    if (!is.na(combo_CA[j])) {
      if (combo_CA[j] == 1) {
        L11 <- L11 * se_CA[j]
        L01 <- L01 * se_CA[j]
      } else {
        L11 <- L11 * (1 - se_CA[j])
        L01 <- L01 * (1 - se_CA[j])
      }
    }
  }
  # CA assays under CA=0 (for L10 and L00)
  for (j in 1:n_CA) {
    if (!is.na(combo_CA[j])) {
      if (combo_CA[j] == 1) {
        L10 <- L10 * (1 - sp_CA[j])
        L00 <- L00 * (1 - sp_CA[j])
      } else {
        L10 <- L10 * sp_CA[j]
        L00 <- L00 * sp_CA[j]
      }
    }
  }
  
  # Prior probabilities of the four latent states
  prior11 <- p_GM * p_CA
  prior10 <- p_GM * (1 - p_CA)
  prior01 <- (1 - p_GM) * p_CA
  prior00 <- (1 - p_GM) * (1 - p_CA)
  
  # Normalizing constant
  norm <- prior00 * L00 + prior01 * L01 + prior10 * L10 + prior11 * L11
  
  # Posterior probabilities
  prob_GM   <- (prior10 * L10 + prior11 * L11) / norm
  prob_CA   <- (prior01 * L01 + prior11 * L11) / norm
  prob_Geno <- 1 - (prior00 * L00) / norm
  
  return(c(GM = prob_GM, CA = prob_CA, Geno = prob_Geno))
}


## extract a parameter from cmdstanr fit object
posterior_extract_byname <- function(fit,name){
  # extract posterior draws from cmdstanr fit object
  fit$draws(format="df") %>% as_tibble %>% select(contains(name))
}

## Plot posterior predictive density plots
cmdstanr_posterior_predictive <- function(stanfit,y,ypredname="fd",ndraws=50){
  # Identify valid observations
  id_valid <- which(!is.na(y) & y != -1)
  
  # Filter valid obs
  y_valid <- y[id_valid]
  nobs <- length(y_valid)
  
  # Extract fitted values
  df_ypred <- posterior_extract_byname(stanfit,ypredname) 
  
  # Subset the observations with valid values
  df_pred_valid <- df_ypred[, id_valid, drop = FALSE]
  
  # Take some random samples
  df_sample <- slice_sample(df_pred_valid,n = ndraws)
  
  # Prepare data for plotting
  df_sample_plot <- pivot_longer(df_sample,cols = everything(),
                                 names_to = "obs",values_to = "y") 
  df_sample_plot$iter <- rep(1:ndraws,each=nobs)
  df_raw <- data.frame(
    y = y_valid,iter = rep(0,nobs)
  )
  df_plot <- bind_rows(df_sample_plot,df_raw) %>%
    mutate(type =ifelse(iter == 0,"y","ypred"))
  
  # Plotting
  g0 <- ggplot(data=df_plot,
               aes(x=y,group=iter))+
    geom_density(aes(color = type))+
    geom_density(data = subset(df_plot,iter ==0),
                 aes(x=y),linewidth = 1)+
    scale_color_manual(values = c("y"="black","ypred"="lightgrey"))+
    theme_classic()
  print(g0)
}

cmdstanr_ROC <- function(
    y,stanfit,probpredname,probabilistic=T,ndraws = 50
){
  # Use only observations with valid Y values
  id_valid <- which(!is.na(y) & y != -1)
  nvalid <- length(id_valid)
  y_valid <- y[id_valid]  # extract valid Y
  
  # Extract predicted probability
  prob_pred <- posterior_extract_byname(stanfit,probpredname)
  prob_pred_valid <- prob_pred[, id_valid, drop = FALSE]  # Use only valid obs
  prob_pred_valid <- as.data.frame(prob_pred_valid)
  
  if(probabilistic){
    # Draw random samples
    if(nvalid < ndraws){
      warning("Specified ndraws=",ndraws," is larger than number of valid observations = ",nvalid,". Use all observations.")
      draw_indices <- seq_len(nvalid)
    } else {
      draw_indices <- sample(seq_len(nvalid),size = ndraws,replace = F)
    }
    
    # ROC for each drawn observation
    roc_list <- list()
    auc_vec <- vector()
    
    for(i in seq_along(draw_indices)){
      idx <- draw_indices[i]
      suppressMessages(     # Suppress setting level messages
        roc_id <- roc(
          response=y_valid,
          predictor = prob_pred_valid[idx,] %>% as.numeric(),
          auc=T
        )
      )
      
      # Create a plot base using the 1st iter, then add lines
      if(i == 1){
        g0 <- ggroc(roc_id,col="blue")
        print(g0)
        # plot(roc_id,
        #      # smooth=T,
        #      # main= paste0("Probabilistic ROC curves (",length(draw_indices)," iterations)"),
        #      col="blue",
        #      xlab = "False Positive Rate",
        #      ylab = "True Positive Rate")
      } else {
        # lines(roc_id, col="blue")
        g1 <- g0 + ggroc(roc_id)
        print(g1)
      }
      
      roc_list[[i]] <- roc_id
      auc_vec[i] <- str_remove(roc_id$auc,"Area under the curve: ") %>% as.numeric()
    }
    
    # Add captions
    g1 <- g1 + 
      labs(
        caption = paste0("AUC (mean=",round(mean(auc_vec),2),
                         ", 95% CI= [",round(quantile(auc_vec,0.05),2),", ",
                         round(quantile(auc_vec,0.95),2),"])")
      )
    # mtext(
    #   paste0("Observed N_pos=", sum(y_valid == 1),
    #              ", Observed N_neg=", sum(y_valid == 0)),
    #       side = 4, line = 0,adj=0.99,cex=0.9,
    # )
    # mtext(
    #   "Blue lines are iterations from MCMC",
    #       side = 4,line = 1,adj=0.99,cex=0.9)
    # mtext(paste0(            # below x axis
    #   "AUC (mean=",round(mean(auc_vec),2),
    #   ", 95% CI= [",round(quantile(auc_vec,0.05),2),", ",
    #   round(quantile(auc_vec,0.95),2),"])"),
    #   side = 1,line = 3.5,cex=0.9,adj=0.99
    # )
    
  } else {
    prob_pred_mean <- apply(prob_pred_valid,2,mean)
    suppressMessages(     # Suppress setting level messages
      roc_mean <- roc(
        response=y_valid,
        predictor = prob_pred_mean,
        auc=T
      )
    )
    auc_mean <- str_remove(roc_mean$auc,"Area under the curve: ")
    g0 <- ggroc(roc_mean,col="blue")+
      labs(caption = paste0("AUC=",auc_mean))
    print(g0)
    # plot(
    #   roc_mean,
    #   col="blue",
    #   xlab = "False Positive Rate",
    #   ylab = "True Positive Rate")
    # mtext(
    
    #   paste("AUC=",auc_mean),
    #   side = 1,line = 3.5,cex=0.9,adj=0.99
    # )
    roc_list <- list(roc_mean)
  }
  
  # Return the list of ROC
  return(roc_list)
}

# softmax function
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

# Numerically stable version of log_sum_exp from stan
log_sum_exp <- function(x) {
  max.x <- max(x)
  return(max.x + log(sum(exp(x - max.x))))
}

# Function to predict endpoint probability given performance parameter and assay data
pred_prob_E <- function(
    data_GM,data_CA, # all vectors
    se_GM,sp_GM,se_CA,sp_CA,pars_po){
  if(length(pars_po) != 5) warnings("Potency vector must have 5 values as pi_RC,p_GM1_1,p_GM0_0,p_CA1_1,p_CA0_0")
  n_GM <- length(se_GM)  # number of GM assays
  n_CA <- length(se_CA)  # number of CA assays
  
  # Compute endpoint likelihood
  log_p_GM_1 <- 0         # temporary storage
  log_p_GM_0 <- 0
  for(i in seq_len(n_GM)){   # for each assay
    # Assay results
    if(data_GM[i] == 1){     # update based on the assay result
      log_p_GM_1 <- log_p_GM_1 + log(se_GM[i])
      log_p_GM_0 <- log_p_GM_0 + log(1-sp_GM[i])
    } else if(data_GM[i] == 0){
      log_p_GM_1 <- log_p_GM_1 + log(1-se_GM[i])
      log_p_GM_0 <- log_p_GM_0 + log(sp_GM[i])
    }
  }
  p_GM_L <- softmax(c(log_p_GM_1,log_p_GM_0))
  
  log_p_CA_1 <- 0         # temporary storage
  log_p_CA_0 <- 0
  for(j in seq_len(n_CA)){
    # Assay results
    if(data_CA[j] == 1){     # update based on the assay result
      log_p_CA_1 <- log_p_CA_1 + log(se_CA[j])
      log_p_CA_0 <- log_p_CA_0 + log(1-sp_CA[j])
    } else if(data_CA[j] == 0){
      log_p_CA_1 <- log_p_CA_1 + log(1-se_CA[j])
      log_p_CA_0 <- log_p_CA_0 + log(sp_CA[j])
    }
  }
  p_CA_L <- softmax(c(log_p_CA_1,log_p_CA_0))
  
  message(
    paste0(
      "If not consider the endpoint prevalence among known chemicals, the probabilities of GM and CA are calculated from the likelihood directly as ",
      p_GM_L," and ",p_CA_L
    )
  )
  
  # Update likelihood by priors
  # Prior prevalence of endpoints from the network
  log_joint_rc1 <- vector(length = 4)
  log_joint_rc0 <- vector(length = 4)
  
  ## RC not observed, use prior pi_RC
  # pars_po: pi_RC,p_GM1_1,p_GM0_0,p_CA1_1,p_CA0_0
  ## RC= 1 | prior GM CA
  log_joint_rc1[1] <- log(pars_po[1]) + log(1 - pars_po[2]) + log(1 - pars_po[4]) + log_p_GM_0  + log_p_CA_0 ## GM=0 CA=0
  log_joint_rc1[2] <- log(pars_po[1]) + log(pars_po[2]) + log(1 - pars_po[4]) + log_p_GM_1 + log_p_CA_0 ## GM=1 CA=0
  log_joint_rc1[3] <- log(pars_po[1]) + log(1 - pars_po[2]) + log(pars_po[4]) + log_p_GM_0 + log_p_CA_1 ## GM=0 CA=1
  log_joint_rc1[4] <- log(pars_po[1]) + log(pars_po[2]) + log(pars_po[4]) + log_p_GM_1 + log_p_CA_1 ## GM=1 CA=1
  ## RC = 0 | prior GM CA
  log_joint_rc0[1] <- log(1 - pars_po[1]) + log(pars_po[3]) + log(pars_po[5]) + log_p_GM_0  + log_p_CA_0 ## GM=0 CA=0
  log_joint_rc0[2] <- log(1 - pars_po[1]) + log(1 - pars_po[3]) + log(pars_po[5]) + log_p_GM_1  + log_p_CA_0 ## GM=1 CA=0
  log_joint_rc0[3] <- log(1 - pars_po[1]) + log(pars_po[3]) + log(1 - pars_po[5]) + log_p_GM_0 + log_p_CA_1 ## GM=0 CA=1
  log_joint_rc0[4] <- log(1 - pars_po[1]) + log(1 - pars_po[3]) + log(1 - pars_po[5]) + log_p_GM_1 + log_p_CA_1 ## GM=1 CA=1
  
  ## Marginalize over the four states to get the unnormalized log-prob of RC
  log_marginal_rc1 <- log_sum_exp(log_joint_rc1)
  log_marginal_rc0 <- log_sum_exp(log_joint_rc0)
  ## Normalize to get P(RC=1 | Assays) using bayes theorom
  p_RC <- exp(log_marginal_rc1 - log_sum_exp(c(log_marginal_rc1, log_marginal_rc0)))
  
  log_marginal_states <- vector(length = 4)
  for(k in 1:4){
    log_marginal_states[k] = log_sum_exp(c(log_joint_rc1[k],log_joint_rc0[k]))
  }
  prob_states <- softmax(log_marginal_states)
  p_GM <- prob_states[2] + prob_states[4]
  p_CA <- prob_states[3] + prob_states[4]
  return(c(p_GM,p_CA,p_RC))
}

pred_iter_E_prob <- function(
    data_GM,  # data are vectors 
    data_CA,
    se_GM,sp_GM,    # parameter are dataframes
    se_CA,sp_CA,
    pars_po){
  niter <- nrow(se_GM)
  mtx_pred <- matrix(0,niter,3)
  for(i in seq_len(niter)){
    se_GM_i <- as.numeric(se_GM[i,])
    sp_GM_i <- as.numeric(sp_GM[i,])
    se_CA_i <- as.numeric(se_CA[i,])
    sp_CA_i <- as.numeric(sp_CA[i,])
    pars_po_i <- as.numeric(pars_po[i,])
    suppressWarnings(suppressMessages(
      mtx_pred[i,] <- pred_prob_E(
        se_GM = se_GM_i,
        sp_GM = sp_GM_i,
        data_GM = data_GM,
        se_CA = se_CA_i,
        sp_CA = sp_CA_i,
        data_CA = data_CA,
        pars_po = pars_po_i
      )
    ))
  }
  return(mtx_pred)
}


# find_combo_matches <- function(row, raw, assays) {
#   results <- list()
#   for (a in assays) {
#     if (is.na(row[[a]])) {
#       for (v in c(T, F)) {
#         candidate <- row
#         candidate[[a]] <- v
#         # Search for matching rows (by assays only)
#         for (i in 1:nrow(raw)) {
#           match_assays <- all(
#             sapply(assays, function(ac) {
#               (is.na(raw[i, ac]) && is.na(candidate[[ac]])) ||
#                 (!is.na(raw[i, ac]) && !is.na(candidate[[ac]]) && raw[i, ac] == candidate[[ac]])
#             })
#           )
#           if (match_assays) {
#             results[[length(results) + 1]] <- list(
#               N = row[["N_vt"]],
#               N_next = raw[i, "N_vt"],
#               Geno = row[["Geno_TRUE"]],
#               Geno_next = raw[i, "Geno_TRUE"],
#               Assay = a,
#               Value = v,
#               combo_next = raw[i, "combo"]
#             )
#           }
#         }
#       }
#     }
#   }
#   return(results)
# }

match_assay_combo <- function(assay, y) {
  # Validate inputs and setup
  y <- as_vector(y)
  if (nrow(assay) != length(y)) {
    stop("The number of rows in 'assay' must match the length of 'y'.")
  }
  
  assay <- as.data.frame(assay)
  assay_names <- colnames(assay)
  
  # Count the number of non-NA assays (-1 represents NA)
  num_non_na <- rowSums(assay != -1)
  
  # Create a combined data frame with original row index and y
  df <- tibble(
    assay,
    y = y,
    orig_idx = 1:nrow(assay),
    num_non_na = num_non_na
  ) %>%
    arrange(num_non_na)   # sort by number of non-NA assays ascending
  
  # Prepare a list to collect result rows
  result_list <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    curr_count <- df$num_non_na[i]
    curr_vals <- df[i, assay_names, drop = FALSE]
    curr_y <- df$y[i]
    curr_combo <- df$orig_idx[i]
    
    # Current row data (assay values, y, combo, N, N_prev)
    base_row <- curr_vals %>%
      mutate(
        y = curr_y,
        combo = curr_combo,
        N = curr_count,
        N_prev = ifelse(curr_count == 0, 0, curr_count - 1)
      )
    
    # If no non-NA assays, it has no previous combo: one row with y0 = y
    if (curr_count == 0) {
      result_list[[i]] <- base_row %>%
        mutate(
          combo_prev = NA_real_,
          new_Assay = NA_character_,
          new_Value = NA_real_,
          y0 = curr_y
        )
      next
    }
    
    # Candidates: rows above current with exactly curr_count - 1 non-NA assays
    candidate_indices <- if (i == 1) integer(0) else which(df$num_non_na[1:(i-1)] == curr_count - 1)
    
    # If no candidates, still output one row but with NA for previous info? 
    # Original didn't handle this, but we'll produce a row with NA previous combo and y0 = NA (or maybe keep as NA, but that might be confusing). 
    # The validation of previous combo is optional; better to output the current combo with NAs.
    if (length(candidate_indices) == 0) {
      result_list[[i]] <- base_row %>%
        mutate(
          combo_prev = NA_real_,
          new_Assay = NA_character_,
          new_Value = NA_real_,
          y0 = NA_real_
        )
      next
    }
    
    matches <- list()  # to collect rows for all valid previous combos
    for (j in candidate_indices) {
      prev_vals <- df[j, assay_names]
      
      # Which assays were performed in the candidate
      non_na_mask <- prev_vals != -1
      
      # Overlapping assays must match
      prev_sub <- unlist(prev_vals[non_na_mask])
      curr_sub <- unlist(curr_vals[non_na_mask])
      
      if (length(prev_sub) == 0 || all(prev_sub == curr_sub)) {
        # Found a matching previous combo
        new_mask <- (prev_vals == -1) & (curr_vals != -1)
        new_col_idx <- which(new_mask)
        
        # If exactly one new assay (should be the case when counts differ by 1)
        if (length(new_col_idx) == 1) {
          new_assay <- assay_names[new_col_idx]
          new_val <- as.numeric(curr_vals[new_col_idx])
        } else {
          new_assay <- NA_character_
          new_val <- NA_real_
        }
        
        match_row <- base_row %>%
          mutate(
            combo_prev = df$orig_idx[j],
            new_Assay = new_assay,
            new_Value = new_val,
            y0 = df$y[j]
          )
        matches <- c(matches, list(match_row))
      }
    }
    
    # If no valid matches found, output a row with NAs
    if (length(matches) == 0) {
      result_list[[i]] <- base_row %>%
        mutate(
          combo_prev = NA_real_,
          new_Assay = NA_character_,
          new_Value = NA_real_,
          y0 = NA_real_
        )
    } else {
      result_list[[i]] <- bind_rows(matches)
    }
  }
  
  final_df <- bind_rows(result_list)
  return(final_df)
}

# function to update positive negative and equivacal in the data
update_values <- function(df, columns, Pvalue, Nvalue, Evalue) {
  for (col in columns) {
    # Get the current column
    vals <- as.character(df[[col]])
    # Use grepl for each rule
    is_plus <- grepl("\\+", vals) | grepl("\\+",vals) | grepl("P",vals)
    is_minus <- grepl("\\-", vals) |grepl("\\−",vals) | grepl("N",vals)    # Use standard minus!
    is_E <- grepl("E", vals) | grepl("equiv",vals)
    # Assign NA by default
    new_vals <- rep(NA, length(vals))
    # Apply rules
    new_vals[is_plus] <- Pvalue
    new_vals[is_minus] <- Nvalue
    new_vals[is_E] <- Evalue
    df[[col]] <- new_vals
  }
  return(df)
}  

# Normalize enligsh greek letters back to greek letters
normalize_greek <- function(chemical_name) {
  # Greek letter replacements (only when standalone parts)
  greek_replacements <- c(
    # Pattern matches:
    # - beginning of string or non-letter character before
    # - Greek letter name
    # - non-letter character or end of string after
    "(^|[^a-zA-Z])alpha([^a-zA-Z]|$)" = "\\1α\\2",
    "(^|[^a-zA-Z])beta([^a-zA-Z]|$)" = "\\1β\\2",
    "(^|[^a-zA-Z])gamma([^a-zA-Z]|$)" = "\\1γ\\2",
    "(^|[^a-zA-Z])delta([^a-zA-Z]|$)" = "\\1δ\\2",
    "(^|[^a-zA-Z])epsilon([^a-zA-Z]|$)" = "\\1ε\\2",
    "(^|[^a-zA-Z])zeta([^a-zA-Z]|$)" = "\\1ζ\\2",
    "(^|[^a-zA-Z])eta([^a-zA-Z]|$)" = "\\1η\\2",
    "(^|[^a-zA-Z])theta([^a-zA-Z]|$)" = "\\1θ\\2",
    "(^|[^a-zA-Z])iota([^a-zA-Z]|$)" = "\\1ι\\2",
    "(^|[^a-zA-Z])kappa([^a-zA-Z]|$)" = "\\1κ\\2",
    "(^|[^a-zA-Z])lambda([^a-zA-Z]|$)" = "\\1λ\\2",
    "(^|[^a-zA-Z])mu([^a-zA-Z]|$)" = "\\1μ\\2",
    "(^|[^a-zA-Z])nu([^a-zA-Z]|$)" = "\\1ν\\2",
    "(^|[^a-zA-Z])xi([^a-zA-Z]|$)" = "\\1ξ\\2",
    "(^|[^a-zA-Z])omicron([^a-zA-Z]|$)" = "\\1ο\\2",
    "(^|[^a-zA-Z])pi([^a-zA-Z]|$)" = "\\1π\\2",
    "(^|[^a-zA-Z])rho([^a-zA-Z]|$)" = "\\1ρ\\2",
    "(^|[^a-zA-Z])sigma([^a-zA-Z]|$)" = "\\1σ\\2",
    "(^|[^a-zA-Z])tau([^a-zA-Z]|$)" = "\\1τ\\2",
    "(^|[^a-zA-Z])upsilon([^a-zA-Z]|$)" = "\\1υ\\2",
    "(^|[^a-zA-Z])phi([^a-zA-Z]|$)" = "\\1φ\\2",
    "(^|[^a-zA-Z])chi([^a-zA-Z]|$)" = "\\1χ\\2",
    "(^|[^a-zA-Z])psi([^a-zA-Z]|$)" = "\\1ψ\\2",
    "(^|[^a-zA-Z])omega([^a-zA-Z]|$)" = "\\1ω\\2"
  )
  
  # Perform replacements for each Greek letter pattern
  for (pattern in names(greek_replacements)) {
    chemical_name <- gsub(pattern, greek_replacements[pattern], chemical_name, ignore.case = FALSE)
  }
  
  return(chemical_name)
}

# Manual computation of posterior prob
compute_posterior <- function(fit, data) {
  # Convert data to matrix and ensure numeric values
  data <- as.matrix(data)
  storage.mode(data) <- "numeric"
  
  # Extract model parameters
  classprob <- fit$classprob
  itemprob <- fit$itemprob
  n_classes <- length(classprob)
  n_obs <- nrow(data)
  n_items <- ncol(data)
  
  # Initialize log-posterior matrix
  log_posterior <- matrix(0, nrow = n_obs, ncol = n_classes)
  
  for (k in 1:n_classes) {
    for (j in 1:n_items) {
      # Identify non-missing observations
      non_na <- !is.na(data[, j])
      
      # Only process if there are valid observations
      if (sum(non_na) > 0) {
        # Extract relevant values
        x <- data[non_na, j]
        p <- itemprob[k, j]
        
        # Calculate log-probabilities manually
        log_probs <- ifelse(x == 1, log(p), log(1 - p))
        
        # Add to log-posterior
        log_posterior[non_na, k] <- log_posterior[non_na, k] + log_probs
      }
    }
    # Add class prior
    log_posterior[, k] <- log_posterior[, k] + log(classprob[k])
  }
  
  # Normalize to probabilities
  max_log <- apply(log_posterior, 1, max)
  posterior <- exp(log_posterior - max_log)
  posterior / rowSums(posterior)
}