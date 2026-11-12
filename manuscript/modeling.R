## Model--------
mod_fullnet <- cmdstan_model(
  "/manuscript/net_fullbayes.stan"
)

## priors-------
prior_uninform <- list(
  po_GM = c(1,1),
  po_CA = c(1,1),
  prior_RC = c(1,1),
  prior_vt = c(1,1),
  prior_vv = c(1,1)
)

prior_inform <- list(
  po_GM = c(10,1),
  po_CA = c(10,1),
  prior_RC = c(10,1),
  prior_vt = c(3,1),
  prior_vv = c(10,1)
)

## Sampling------

init_fun <- function() {
  list(
    p_GM1_1  = 0.8, p_GM0_0  = 0.8,
    p_CA1_1 = 0.8, p_CA0_0 = 0.8,
    pi_RC = 0.7,
    se_GM_vt = rep(0.7, length(name_vt_GM)), sp_GM_vt = rep(0.7, length(name_vt_GM)),
    se_GM_vv = rep(0.7, length(name_vv_GM)), sp_GM_vv = rep(0.7, length(name_vv_GM)),
    se_CA_vt = rep(0.7, length(name_vt_CA)), sp_CA_vt = rep(0.7, length(name_vt_CA)),
    se_CA_vv = rep(0.7, length(name_vv_CA)), sp_CA_vv = rep(0.7, length(name_vv_CA))
  )
}

fit_raw_uninform <- mod_fullnet$sample(
  data = append(prior_uninform,data_list),
  init = init_fun,
  timestamp = F,random = F
)

fit_raw_uninform$save_output_files(
  dir = "fits",basename = "fit_raw_uninformative",
  init = init_fun,
  timestamp = F,random = F
)

## Predictive check-----
loo_raw_uninform <- fit_raw_uninform$loo()

loo_raw_inform <- fit_raw_inform$loo()

loo_compare(loo_raw_uninform,loo_raw_inform)   # raw data
