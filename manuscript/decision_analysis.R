## Extract parameter posterior------

# List parameters of interest

pars_print <- c(
  "se_GM","sp_GM","se_CA","sp_CA",
  "p_GM1_1","p_GM0_0","p_CA1_1","p_CA0_0","pi_RC"
)

# Parameter names for display

name_GM <- c(name_vt_GM,name_vv_GM)
name_CA <- c("MNvt","CAvt","MNvv","CAvv")

pars_raw <- fit_raw_uninform$summary(pars_print)

df_mean_raw <- tibble(
"Assay"= c(name_GM,name_CA),
"se_mean" = rbind(
  pars_raw[1:6,2],
  pars_raw[13:16,2]
) %>% unlist(),
"sp_mean" = rbind(
  pars_raw[7:12,2],
  pars_raw[17:20,2]
)   %>% unlist()
) %>% 
  mutate_if(is.numeric,round,digits=3)

df_low_raw <- tibble(
  "Assay"= c(name_GM,name_CA),
  "se_low" = rbind(
    pars_raw[1:6,6],
    pars_raw[13:16,6]
  ) %>% unlist(),
  "sp_low" = rbind(
    pars_raw[7:12,6],
    pars_raw[17:20,6]
  )   %>% unlist()
) %>% 
  mutate_if(is.numeric,round,digits=3)

df_up_raw <- tibble(
  "Assay"= c(name_GM,name_CA),
  "se_up" = rbind(
    pars_raw[1:6,7],
    pars_raw[13:16,7]
  ) %>% unlist(),
  "sp_up" = rbind(
    pars_raw[7:12,7],
    pars_raw[17:20,7]
  )   %>% unlist()
) %>% 
  mutate_if(is.numeric,round,digits=3)

df_sd_raw <- tibble(
  "Assay"= c(name_GM,name_CA),
  "se_sd" = rbind(
    pars_raw[1:6,4],
    pars_raw[13:16,4]
  ) %>% unlist(),
  "sp_sd" = rbind(
    pars_raw[7:12,4],
    pars_raw[17:20,4]
  )   %>% unlist()
) %>% 
  mutate_if(is.numeric,round,digits=3)

## Predict by all assays------
fit <- fit_raw_uninform
draw_se_GM <- posterior_extract_byname(fit,c("se_GM_vt","se_GM_vv"))
draw_sp_GM <- posterior_extract_byname(fit,c("sp_GM_vt","sp_GM_vv"))
draw_se_CA <- posterior_extract_byname(fit,c("se_CA_vt","se_CA_vv"))
draw_sp_CA <- posterior_extract_byname(fit,c("sp_CA_vt","sp_CA_vv"))
pars_po <- posterior_extract_byname(fit,
                                    c("pi_RC","p_GM1_1","p_GM0_0","p_CA1_1","p_CA0_0")                   
)
assay_GM <- cbind(assay_GM_vt,assay_GM_vv)
assay_CA <- cbind(assay_CA_vt,assay_CA_vv)

nobs <- nrow(assay_GM)

# Only use the average
se_GM_mean <- apply(draw_se_GM,2,mean)
sp_GM_mean <- apply(draw_sp_GM,2,mean)
se_CA_mean <- apply(draw_se_CA,2,mean)
sp_CA_mean <- apply(draw_sp_CA,2,mean)
pars_po_mean <- apply(pars_po,2,mean)
mtx_pred_all <- matrix(0,nobs,3)
colnames(mtx_pred_all) <- c("p_GM","p_CA","p_RC")
for(n in seq_len(nobs)){
  suppressMessages(  
    pred_all <- pred_prob_E(
      data_GM = assay_GM[n,],
      data_CA = assay_CA[n,],
      se_GM = se_GM_mean,sp_GM = sp_GM_mean,
      se_CA = se_CA_mean,sp_CA = sp_CA_mean,
      pars_po = pars_po_mean
    )
  )
  mtx_pred_all[n,1] <- pred_all[1]
  mtx_pred_all[n,2] <- pred_all[2]
  mtx_pred_all[n,3] <- pred_all[3]
}

## Predict by in vitro only------
fit <- fit_raw_uninform
draw_se_GM <- posterior_extract_byname(fit,"se_GM_vt")
draw_sp_GM <- posterior_extract_byname(fit,"sp_GM_vt")
draw_se_CA <- posterior_extract_byname(fit,"se_CA_vt")
draw_sp_CA <- posterior_extract_byname(fit,"sp_CA_vt")
pars_po <- posterior_extract_byname(fit,
                                    c("pi_RC","p_GM1_1","p_GM0_0","p_CA1_1","p_CA0_0")                   
)
nobs <- nrow(assay_GM_vt)

# Only use the average value

se_GM_mean <- apply(draw_se_GM,2,mean)
sp_GM_mean <- apply(draw_sp_GM,2,mean)
se_CA_mean <- apply(draw_se_CA,2,mean)
sp_CA_mean <- apply(draw_sp_CA,2,mean)
pars_po_mean <- apply(pars_po,2,mean)
mtx_pred_vt <- matrix(0,nobs,3)
colnames(mtx_pred_vt) <- c("p_GM","p_CA","p_RC")
for(n in seq_len(nobs)){
  suppressMessages(  
    pred_vt <- pred_prob_E(
      data_GM = assay_GM_vt[n,],
      data_CA = assay_CA_vt[n,],
      se_GM = se_GM_mean,sp_GM = sp_GM_mean,
      se_CA = se_CA_mean,sp_CA = sp_CA_mean,
      pars_po = pars_po_mean
    )
  )
  mtx_pred_vt[n,1] <- pred_vt[1]
  mtx_pred_vt[n,2] <- pred_vt[2]
  mtx_pred_vt[n,3] <- pred_vt[3]
}

## Predict by in vivo only------
fit <- fit_raw_uninform
draw_se_GM <- posterior_extract_byname(fit,"se_GM_vv")
draw_sp_GM <- posterior_extract_byname(fit,"sp_GM_vv")
draw_se_CA <- posterior_extract_byname(fit,"se_CA_vv")
draw_sp_CA <- posterior_extract_byname(fit,"sp_CA_vv")
pars_po <- posterior_extract_byname(fit,
                                    c("pi_RC","p_GM1_1","p_GM0_0","p_CA1_1","p_CA0_0")                   
)
nobs <- nrow(assay_GM_vv)

# Only use the average value

se_GM_mean <- apply(draw_se_GM,2,mean)
sp_GM_mean <- apply(draw_sp_GM,2,mean)
se_CA_mean <- apply(draw_se_CA,2,mean)
sp_CA_mean <- apply(draw_sp_CA,2,mean)
pars_po_mean <- apply(pars_po,2,mean)
mtx_pred_vv <- matrix(0,nobs,3)
colnames(mtx_pred_vv) <- c("p_GM","p_CA","p_RC")
for(n in seq_len(nobs)){
  suppressMessages(  
    pred_vv <- pred_prob_E(
      data_GM = assay_GM_vv[n,],
      data_CA = assay_CA_vv[n,],
      se_GM = se_GM_mean,sp_GM = sp_GM_mean,
      se_CA = se_CA_mean,sp_CA = sp_CA_mean,
      pars_po = pars_po_mean
    )
  )
  mtx_pred_vv[n,1] <- pred_vv[1]
  mtx_pred_vv[n,2] <- pred_vv[2]
  mtx_pred_vv[n,3] <- pred_vv[3]
}

## Stepwise adaptation-----
nobs_expand <- nrow(Expand)
fit <- fit_raw_uninform
draw_se_GM <- posterior_extract_byname(fit,"se_GM_vt")
draw_sp_GM <- posterior_extract_byname(fit,"sp_GM_vt")
draw_se_CA <- posterior_extract_byname(fit,"se_CA_vt")
draw_sp_CA <- posterior_extract_byname(fit,"sp_CA_vt")
pars_po <- posterior_extract_byname(fit,
                                    c("pi_RC","p_GM1_1","p_GM0_0","p_CA1_1","p_CA0_0")                   
)
# Use only average values
se_GM_mean <- apply(draw_se_GM,2,mean)
sp_GM_mean <- apply(draw_sp_GM,2,mean)
se_CA_mean <- apply(draw_se_CA,2,mean)
sp_CA_mean <- apply(draw_sp_CA,2,mean)
pars_po_mean <- apply(pars_po,2,mean)
mtx_pred_expand <- matrix(0,nobs_expand,3)
colnames(mtx_pred_expand) <- c("p_GM","p_CA","p_RC")
for(n in seq_len(nobs_expand)){
  suppressMessages(  
    pred_vt <- pred_prob_E(
      data_GM = Expand[n,1:3],
      data_CA = Expand[n,4:5],
      se_GM = se_GM_mean,sp_GM = sp_GM_mean,
      se_CA = se_CA_mean,sp_CA = sp_CA_mean,
      pars_po = pars_po_mean
    )
  )
  mtx_pred_expand[n,1] <- pred_vt[1]
  mtx_pred_expand[n,2] <- pred_vt[2]
  mtx_pred_expand[n,3] <- pred_vt[3]
}
vt_expand_pred <- cbind(Expand,mtx_pred_expand)    # Expand and mtx_pred match?

df_GM <- df_expand %>%
  filter(MNvt == -1 & CAvt == -1)
results_GM <- match_assay_combo(
  assay = df_GM %>% select(1:3),
  y= df_GM %>% select(p_GM)
)

df_CA <- df_expand %>%
  filter(Ames == -1 & MLA == -1 & HPRT == -1)


results_CA <- match_assay_combo(
  assay = df_CA %>% select(4:5),
  y= df_CA %>% select(p_CA)
)

results_RC <- match_assay_combo(
  assay = df_expand %>% select(1:5),
  y= df_expand %>% select(p_RC)
)
