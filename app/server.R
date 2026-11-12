library(shiny)
library(shinyjs)
library(readr)

## Utility functions -------
# Softmax function
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par - 1)) {
    Lk <- max(par1[k + 1], Lk) + log1p(exp(-abs(par1[k + 1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

# Numerically stable version of log_sum_exp from Stan
log_sum_exp <- function(x) {
  max.x <- max(x)
  return(max.x + log(sum(exp(x - max.x))))
}

# Function to predict endpoint probability given performance parameter and assay data
pred_prob_E <- function(data_GM, data_CA, se_GM, sp_GM, se_CA, sp_CA, pars_po){
  if(length(pars_po) != 5) warning("Potency vector must have 5 values as pi_RC, p_GM1_1, p_GM0_0, p_CA1_1, p_CA0_0")
  n_GM <- length(se_GM)
  n_CA <- length(se_CA)
  
  # Compute GM endpoint likelihood
  log_p_GM_1 <- 0
  log_p_GM_0 <- 0
  for(i in seq_len(n_GM)){
    if(data_GM[i] == 1){
      log_p_GM_1 <- log_p_GM_1 + log(se_GM[i])
      log_p_GM_0 <- log_p_GM_0 + log(1 - sp_GM[i])
    } else if(data_GM[i] == 0){
      log_p_GM_1 <- log_p_GM_1 + log(1 - se_GM[i])
      log_p_GM_0 <- log_p_GM_0 + log(sp_GM[i])
    }
  }
  p_GM_L <- softmax(c(log_p_GM_1, log_p_GM_0))
  
  # Compute CA endpoint likelihood
  log_p_CA_1 <- 0
  log_p_CA_0 <- 0
  for(j in seq_len(n_CA)){
    if(data_CA[j] == 1){
      log_p_CA_1 <- log_p_CA_1 + log(se_CA[j])
      log_p_CA_0 <- log_p_CA_0 + log(1 - sp_CA[j])
    } else if(data_CA[j] == 0){
      log_p_CA_1 <- log_p_CA_1 + log(1 - se_CA[j])
      log_p_CA_0 <- log_p_CA_0 + log(sp_CA[j])
    }
  }
  p_CA_L <- softmax(c(log_p_CA_1, log_p_CA_0))
  
  # Update likelihood by priors
  log_joint_rc1 <- vector(length = 4)
  log_joint_rc0 <- vector(length = 4)
  
  # pars_po: pi_RC, p_GM1_1, p_GM0_0, p_CA1_1, p_CA0_0
  # RC = 1 | prior GM CA
  log_joint_rc1[1] <- log(pars_po[1]) + log(1 - pars_po[2]) + log(1 - pars_po[4]) + log_p_GM_0 + log_p_CA_0 # GM=0 CA=0
  log_joint_rc1[2] <- log(pars_po[1]) + log(pars_po[2])     + log(1 - pars_po[4]) + log_p_GM_1 + log_p_CA_0 # GM=1 CA=0
  log_joint_rc1[3] <- log(pars_po[1]) + log(1 - pars_po[2]) + log(pars_po[4])     + log_p_GM_0 + log_p_CA_1 # GM=0 CA=1
  log_joint_rc1[4] <- log(pars_po[1]) + log(pars_po[2])     + log(pars_po[4])     + log_p_GM_1 + log_p_CA_1 # GM=1 CA=1
  
  # RC = 0 | prior GM CA
  log_joint_rc0[1] <- log(1 - pars_po[1]) + log(pars_po[3])     + log(pars_po[5])     + log_p_GM_0 + log_p_CA_0 # GM=0 CA=0
  log_joint_rc0[2] <- log(1 - pars_po[1]) + log(1 - pars_po[3]) + log(pars_po[5])     + log_p_GM_1 + log_p_CA_0 # GM=1 CA=0
  log_joint_rc0[3] <- log(1 - pars_po[1]) + log(pars_po[3])     + log(1 - pars_po[5]) + log_p_GM_0 + log_p_CA_1 # GM=0 CA=1
  log_joint_rc0[4] <- log(1 - pars_po[1]) + log(1 - pars_po[3]) + log(1 - pars_po[5]) + log_p_GM_1 + log_p_CA_1 # GM=1 CA=1
  
  # Marginalize over the four states to get unnormalized log-prob of RC
  log_marginal_rc1 <- log_sum_exp(log_joint_rc1)
  log_marginal_rc0 <- log_sum_exp(log_joint_rc0)
  p_RC <- exp(log_marginal_rc1 - log_sum_exp(c(log_marginal_rc1, log_marginal_rc0)))
  
  log_marginal_states <- vector(length = 4)
  for(k in 1:4){
    log_marginal_states[k] <- log_sum_exp(c(log_joint_rc1[k], log_joint_rc0[k]))
  }
  prob_states <- softmax(log_marginal_states)
  p_GM <- prob_states[2] + prob_states[4]
  p_CA <- prob_states[3] + prob_states[4]
  return(c(p_GM, p_CA, p_RC))
}

# Helper to evaluate a probability against conclusive thresholds
eval_endpoint_status <- function(prob, ep_label, upper_t, lower_t) {
  if (prob > upper_t) {
    list(conclusive = TRUE, 
         msg = sprintf("%s: Probability = %s (> %s) — Endpoint is CONFIRMED (Positive).", 
                       ep_label, fmt_num(prob), fmt_num(upper_t)))
  } else if (prob < lower_t) {
    list(conclusive = TRUE, 
         msg = sprintf("%s: Probability = %s (< %s) — Endpoint is REJECTED (Negative).", 
                       ep_label, fmt_num(prob), fmt_num(lower_t)))
  } else {
    list(conclusive = FALSE, 
         msg = sprintf("%s: Probability = %s (between %s and %s) — Inconclusive; more testing is needed.", 
                       ep_label, fmt_num(prob), fmt_num(lower_t), fmt_num(upper_t)))
  }
}

# Helper to parse select inputs
parse_result <- function(res) {
  if (res == "Positive") return(1)
  if (res == "Negative") return(0)
  return(-1)
}

# Global display precision setting
DISPLAY_DIGITS <- 3

# Helper function to format any numeric vector to DISPLAY_DIGITS
fmt_num <- function(x, digits = DISPLAY_DIGITS) {
  sprintf(paste0("%.", digits, "f"), x)
}

## Backend Server -------
server <- function(input, output, session) {
  
  # Global reactive state
  rv <- reactiveValues(
    mtx_pred_expand = NULL,
    current_combo_idx = NULL,
    assay_test = NULL,
    Expand = NULL,
    par_mean = NULL,
    assay_names = c("Ames", "MLA", "HPRT", "MNvt", "CAvt"),
    thresh_upper = 0.975,
    thresh_lower = 0.025
  )
  
  # Load data files safely (checks both 'data/' subfolder and current directory)
  load_app_data <- function() {
    if (is.null(rv$Expand)) {
      expand_file <- if (file.exists("data/vt_expand.csv")) {
        "data/vt_expand.csv"
      } else if (file.exists("vt_expand.csv")) {
        "vt_expand.csv"
      } else {
        NULL
      }
      if (is.null(expand_file)) {
        showNotification("Error: Could not find 'vt_expand.csv' in 'data/' or app folder.", type = "error", duration = 10)
        return(FALSE)
      }
      rv$Expand <- read_csv(expand_file, show_col_types = FALSE)
    }
    
    if (is.null(rv$par_mean)) {
      par_candidates <- c(
        "data/par_mean_raw_1.1.csv",
        "par_mean_raw_1.1.csv",
        "data/par_mean_raw.csv",
        "par_mean_raw.csv"
      )
      existing_par <- par_candidates[file.exists(par_candidates)]
      if (length(existing_par) == 0) {
        showNotification("Error: Could not find 'par_mean_raw_1.1.csv' in 'data/' or app folder.", type = "error", duration = 10)
        return(FALSE)
      }
      rv$par_mean <- read_csv(existing_par[1], show_col_types = FALSE)
    }
    return(TRUE)
  }
  
  # Updated Start Button Handler
  observeEvent(input$start_btn, {
    shinyjs::hide("page1")
    shinyjs::show("page2")
    load_app_data()
  })
  
  # Helper to compute the prediction matrix for all 243 combinations given pi_RC
  compute_prediction_matrix <- function(pi_RC_val) {
    load_app_data()
    par_mean <- rv$par_mean
    se_GM_mean <- par_mean[[2]][1:3]
    sp_GM_mean <- par_mean[[2]][7:9]
    se_CA_mean <- par_mean[[2]][13:14]
    sp_CA_mean <- par_mean[[2]][17:18]
    pars_po_mean <- par_mean[[2]][21:24]
    
    pars_po <- c(pi_RC_val, pars_po_mean)
    nobs <- nrow(rv$Expand)
    mtx <- matrix(0, nobs, 3)
    colnames(mtx) <- c("p_GM", "p_CA", "p_RC")
    
    for(n in seq_len(nobs)){
      suppressMessages(
        pred_vt <- pred_prob_E(
          data_GM = as.numeric(rv$Expand[n, 1:3]),
          data_CA = as.numeric(rv$Expand[n, 4:5]),
          se_GM = se_GM_mean, sp_GM = sp_GM_mean,
          se_CA = se_CA_mean, sp_CA = sp_CA_mean,
          pars_po = pars_po
        )
      )
      mtx[n, ] <- pred_vt
    }
    return(mtx)
  }
  
  # Threshold UI toggles
  observeEvent(input$keep_thresh_p2, {
    if (input$keep_thresh_p2 == "no") show("custom_thresh_p2_div") else hide("custom_thresh_p2_div")
  })
  
  observeEvent(input$keep_thresh_3b, {
    if (input$keep_thresh_3b == "no") show("custom_thresh_3b_div") else hide("custom_thresh_3b_div")
  })
  
  # Prior help modal & toggle
  observeEvent(input$help_prior, {
    showModal(modalDialog(
      title = "Prior Probability",
      "Prior belief affects posterior probability significantly when assays show conflicting results. See our paper for details.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$provide_prior, {
    if(input$provide_prior == "yes") show("prior_input_div") else hide("prior_input_div")
  })
  
  # ==========================================
  # Navigation Handlers
  # ==========================================
  observeEvent(input$start_btn, {
    load_app_data()
    hide("page1")
    show("page2")
  })
  
  observeEvent(input$next_btn, {
    # Save threshold chosen on Page 2 and sync with Page 3b
    if (input$keep_thresh_p2 == "yes") {
      rv$thresh_upper <- 0.975
      rv$thresh_lower <- 0.025
      updateRadioButtons(session, "keep_thresh_3b", selected = "yes")
    } else {
      rv$thresh_upper <- input$thresh_upper_p2
      rv$thresh_lower <- input$thresh_lower_p2
      updateRadioButtons(session, "keep_thresh_3b", selected = "no")
      updateNumericInput(session, "thresh_upper_3b", value = input$thresh_upper_p2)
      updateNumericInput(session, "thresh_lower_3b", value = input$thresh_lower_p2)
    }
    
    hide("page2")
    if (input$app_intent == "a") {
      show("page3a")
    } else {
      show("page3b")
    }
  })
  
  observeEvent(input$back_btn_3a, {
    hide("page3a")
    show("page2")
  })
  
  observeEvent(input$back_btn_3b, {
    hide("page3b")
    show("page2")
  })
  
  # Transition from Path A -> Path B when user wants follow-up
  observeEvent(input$yes_to_3b_btn, {
    updateSelectInput(session, "ames", selected = "NA")
    updateSelectInput(session, "mla", selected = "NA")
    updateSelectInput(session, "hprt", selected = input$hprt_3a)
    updateSelectInput(session, "mnvt", selected = "NA")
    updateSelectInput(session, "cavt", selected = input$cavt_3a)
    
    hide("page3a")
    hide("followup_div")
    hide("ask_followup_prompt_div")
    output$results_ui <- renderUI({ NULL })
    show("page3b")
  })
  
  # Transitions to Ending Page (Page 4)
  observeEvent(input$no_to_end_3a_btn, {
    hide("page3a")
    show("page4")
  })
  
  observeEvent(input$no_to_end_3b_btn, {
    hide("page3b")
    show("page4")
  })
  
  observeEvent(input$to_end_from_followup_btn, {
    hide("page3b")
    show("page4")
  })
  
  observeEvent(input$restart_btn, {
    hide("page4")
    hide("more_testing_3a_div")
    hide("ask_followup_prompt_div")
    hide("followup_div")
    output$results_3a_ui <- renderUI({ NULL })
    output$results_ui <- renderUI({ NULL })
    show("page1")
  })
  
  # ==========================================
  # Path A: Starting Assay Calculation
  # ==========================================
  observeEvent(input$calc_3a_btn, {
    mtx <- compute_prediction_matrix(0.5)
    
    test_arr_3a <- c(
      -1, # Ames
      -1, # MLA
      parse_result(input$hprt_3a),
      -1, # MNvt
      parse_result(input$cavt_3a)
    )
    
    match_idx <- which(apply(rv$Expand[, 1:5], 1, function(r) all(r == test_arr_3a)))
    req(length(match_idx) > 0)
    idx <- match_idx[1]
    
    p_GM <- mtx[idx, "p_GM"]
    p_CA <- mtx[idx, "p_CA"]
    
    gm_eval <- eval_endpoint_status(p_GM, "Gene Mutation P(GM)", rv$thresh_upper, rv$thresh_lower)
    ca_eval <- eval_endpoint_status(p_CA, "Chromosome Aberration P(CA)", rv$thresh_upper, rv$thresh_lower)
    
    output$results_3a_ui <- renderUI({
      tagList(
        wellPanel(
          h4("Calculated Probabilities & Threshold Evaluation"),
          p(sprintf("Active Conclusive Thresholds: Positive > %s | Negative < %s", 
                    fmt_num(rv$thresh_upper), fmt_num(rv$thresh_lower))),
          tags$ul(
            tags$li(strong(gm_eval$msg)),
            tags$li(strong(ca_eval$msg))
          ),
          if (gm_eval$conclusive && ca_eval$conclusive) {
            p(strong("Both Gene Mutation and Chromosome Aberration endpoints have reached conclusive decisions!"), style = "color: green;")
          } else {
            p(strong("More testing is needed to reach a conclusive decision."), style = "color: #d9534f;")
          }
        )
      )
    })
    
    if (!gm_eval$conclusive || !ca_eval$conclusive) {
      show("more_testing_3a_div")
    } else {
      hide("more_testing_3a_div")
    }
  })
  
  # ==========================================
  # Path B: Existing Data Calculation
  # ==========================================
  observeEvent(input$calc_btn, {
    # Update active thresholds from Path B settings
    if (input$keep_thresh_3b == "yes") {
      rv$thresh_upper <- 0.975
      rv$thresh_lower <- 0.025
    } else {
      rv$thresh_upper <- input$thresh_upper_3b
      rv$thresh_lower <- input$thresh_lower_3b
    }
    
    pi_RC <- if(input$provide_prior == "yes") input$pi_RC_input else 0.5
    rv$mtx_pred_expand <- compute_prediction_matrix(pi_RC)
    
    test_arr <- c(
      parse_result(input$ames),
      parse_result(input$mla),
      parse_result(input$hprt),
      parse_result(input$mnvt),
      parse_result(input$cavt)
    )
    rv$assay_test <- test_arr
    
    match_idx <- which(apply(rv$Expand[, 1:5], 1, function(r) all(r == rv$assay_test)))
    
    if(length(match_idx) > 0) {
      rv$current_combo_idx <- match_idx[1]
      
      p_GM <- rv$mtx_pred_expand[rv$current_combo_idx, "p_GM"]
      p_CA <- rv$mtx_pred_expand[rv$current_combo_idx, "p_CA"]
      p_RC <- rv$mtx_pred_expand[rv$current_combo_idx, "p_RC"]
      
      gm_eval <- eval_endpoint_status(p_GM, "Gene Mutation P(GM)", rv$thresh_upper, rv$thresh_lower)
      ca_eval <- eval_endpoint_status(p_CA, "Chromosome Aberration P(CA)", rv$thresh_upper, rv$thresh_lower)
      rc_eval <- eval_endpoint_status(p_RC, "Rodent Carcinogenicity P(RC)", rv$thresh_upper, rv$thresh_lower)
      
      show_rc_flag <- (input$show_rc == "yes")
      
      all_conclusive <- gm_eval$conclusive && ca_eval$conclusive && (!show_rc_flag || rc_eval$conclusive)
      
      output$results_ui <- renderUI({
        tagList(
          wellPanel(
            h4("Current Probabilities & Conclusive Evaluation"),
            p(sprintf("Active Conclusive Thresholds: Positive > %s | Negative < %s", 
                      fmt_num(rv$thresh_upper), fmt_num(rv$thresh_lower))),
            tags$ul(
              tags$li(strong(gm_eval$msg)),
              tags$li(strong(ca_eval$msg)),
              if (show_rc_flag) tags$li(strong(rc_eval$msg))
            ),
            if (all_conclusive) {
              tagList(
                p(strong("All displayed endpoints have reached a conclusive decision!"), style = "color: green;"),
                actionButton("conclusive_finish_btn", "Finish", class = "btn-primary")
              )
            } else {
              p(strong("More testing is needed for at least one displayed endpoint."), style = "color: #d9534f;")
            }
          )
        )
      })
      
      hide("followup_div")
      if (!all_conclusive) {
        show("ask_followup_prompt_div")
      } else {
        hide("ask_followup_prompt_div")
      }
      
    } else {
      output$results_ui <- renderUI({
        p("Error: Could not match the selected pattern in the Expand combinations.")
      })
    }
  })
  
  observeEvent(input$conclusive_finish_btn, {
    hide("page3b")
    show("page4")
  })
  
  # ==========================================
  # Follow-up Recommendation Logic (GM, CA, RC)
  # ==========================================
  observeEvent(input$ask_followup, {
    req(rv$Expand, rv$mtx_pred_expand, rv$current_combo_idx)
    
    # Helper to find the best next assay for a specific endpoint & allowed assay subset
    get_best_followup <- function(ep_col, allowed_indices, ep_title) {
      untested_allowed <- intersect(which(rv$assay_test == -1), allowed_indices)
      
      if (length(untested_allowed) == 0) {
        return(tagList(
          h5(strong(ep_title)),
          p("All relevant assays for this endpoint have already been tested. No follow-up assay can be suggested.")
        ))
      }
      
      curr_N <- sum(rv$assay_test != -1)
      tested_idx <- which(rv$assay_test != -1)
      
      # Find candidate rows where exactly 1 new assay from allowed_indices is tested
      valid_next_idx <- which(apply(rv$Expand[, 1:5], 1, function(r) {
        if (sum(r != -1) != (curr_N + 1)) return(FALSE)
        if (length(tested_idx) > 0 && !all(r[tested_idx] == rv$assay_test[tested_idx])) return(FALSE)
        new_assay_idx <- which(r != -1 & rv$assay_test == -1)
        return(length(new_assay_idx) == 1 && (new_assay_idx %in% allowed_indices))
      }))
      
      if (length(valid_next_idx) == 0) {
        return(tagList(
          h5(strong(ep_title)),
          p("No valid follow-up combinations found.")
        ))
      }
      
      candidate_mtx <- rv$mtx_pred_expand[valid_next_idx, , drop = FALSE]
      max_idx_local <- which.max(candidate_mtx[, ep_col])
      min_idx_local <- which.min(candidate_mtx[, ep_col])
      
      row_max <- valid_next_idx[max_idx_local]
      row_min <- valid_next_idx[min_idx_local]
      
      added_max <- rv$assay_names[which(rv$Expand[row_max, 1:5] != rv$assay_test & rv$assay_test == -1)][1]
      added_min <- rv$assay_names[which(rv$Expand[row_min, 1:5] != rv$assay_test & rv$assay_test == -1)][1]
      
      tagList(
        h5(strong(ep_title)),
        tags$ul(
          tags$li(sprintf("Best follow-up towards CONFIRMING %s: Add %s (reaches %.3f if positive)", 
                          ep_title, added_max, candidate_mtx[max_idx_local, ep_col])),
          tags$li(sprintf("Best follow-up towards REJECTING %s: Add %s (reaches %.3f if negative)", 
                          ep_title, added_min, candidate_mtx[min_idx_local, ep_col]))
        )
      )
    }
    
    # 1. GM recommendation: Ames (1), MLA (2), HPRT (3)
    gm_rec_ui <- get_best_followup("p_GM", c(1, 2, 3), "Gene Mutation")
    
    # 2. CA recommendation: MNvt (4), CAvt (5)
    ca_rec_ui <- get_best_followup("p_CA", c(4, 5), "Chromosome Aberration")
    
    # 3. RC recommendation: All 5 assays (1:5), shown only if user asked to see P(RC)
    rc_rec_ui <- if (input$show_rc == "yes") {
      tagList(
        hr(),
        get_best_followup("p_RC", 1:5, "Rodent Carcinogenicity")
      )
    } else {
      NULL
    }
    
    output$followup_ui <- renderUI({
      wellPanel(
        gm_rec_ui,
        hr(),
        ca_rec_ui,
        rc_rec_ui
      )
    })
    
    show("followup_div")
  })
}