## RC by all assays----

id_valid <- which(!is.na(RC) & RC != -1) %>% unname()
nvalid <- length(id_valid)
y_valid <- RC[id_valid] %>% unname()  # extract valid Y
prob_pred_valid <- mtx_pred_all[id_valid,3] %>% as_vector()  # Use only valid obs
roc_raw <- roc(response=y_valid,predictor=prob_pred_valid)
roc_raw$auc
ggroc(roc_raw)

## RC by in vitro AUC-------

prob_pred_valid <- mtx_pred_vt[id_valid,3] %>% as_vector()
ROC_vt <- roc(response=y_valid,predictor=prob_pred_valid)
ROC_vt$auc

ggroc(ROC_vt)

## RC by in vivo AUC-------

prob_pred_valid <- mtx_pred_vv[id_valid,3] %>% as_vector()
ROC_vv <- roc(response=y_valid,predictor=prob_pred_valid)
ROC_vv$auc

ggroc(ROC_vv)

## RC in vitro vs in vivo vs all------
# get auc values
roc_raw$auc   
ROC_vt$auc    
ROC_vv$auc    

# plotting
g0 <- ggroc(
  list(
    "In Vitro"=ROC_vt,"In Vivo"=ROC_vv,
    "In Vitro & In Vivo"=roc_raw
  )
)+
  geom_segment(aes(x=1,xend=0,y=0,yend=1),col="grey",linetype="dashed")+
  labs(caption = "AUC:                                            \n In Vitro=0.573 \n In Vivo = 0.759 \n In Vivo & In Vitro=0.643")+
  theme_classic(base_size = 15)+
  theme(legend.position = "bottom")

## Assay accuracy contour plot----

# Use computed mean from network.R
df_mean <- df_mean_raw
df_mean_plot <- df_mean %>%
  mutate(
    "Endpoint" = c(rep_len(c("Gene Mutation"),c(6)),rep_len("Chromosome Aberration",4)),
    "Type" = c(rep("in vitro",3),rep("in vivo",3),rep("in vitro",2),rep("in vivo",2))
  )
g1 <- ggplot(df_mean_plot,aes(x=se_mean,y=sp_mean,
                              label=Assay,col=Type))+
  scale_color_manual(values = 2:11)+
  facet_wrap(~ Endpoint)+
  geom_point(aes(col=Type,shape = Type))+
  geom_text(size = 5,nudge_y=0.02,show.legend = F)+
  annotate(geom="path",
           x=1-0.1*cos(seq(0,pi,length.out=100)),
           y=1-0.1*sin(seq(0,pi,length.out=100)),
           color="gray",linetype="dashed",linewidth=0.5)+
  annotate(geom="path",
           x=1-0.2*cos(seq(0,pi,length.out=100)),
           y=1-0.2*sin(seq(0,pi,length.out=100)),
           color="gray",linetype="dashed",linewidth=0.5)+
  annotate(geom="path",
           x=1-0.3*cos(seq(0,pi,length.out=100)),
           y=1-0.3*sin(seq(0,pi,length.out=100)),
           color="gray",linetype="dashed",linewidth=0.5)+
  annotate(geom="path",
           x=1-0.4*cos(seq(0,pi,length.out=100)),
           y=1-0.4*sin(seq(0,pi,length.out=100)),
           color="gray",linetype="dashed",linewidth=0.5)+
  annotate(geom="path",
           x=1-0.5*cos(seq(0,pi,length.out=100)),
           y=1-0.5*sin(seq(0,pi,length.out=100)),
           color="gray",linetype="dashed",linewidth=0.5)+
  annotate(geom="text",x=0.51,y=0.5,label="Poor")+
  annotate(geom="text",x=0.97,y=1,label="Excellent")+
  coord_cartesian(xlim = c(0.5,1),ylim = c(0.5,1))+
  labs(x="Sensitivity",y="Specificity",
       caption = "Poor: no better than a random guess \n Excellent: perfect accuracy")+
  theme_classic(base_size = 20)+
  theme(legend.position = "bottom")

## In Vivo Endpoint by all AUC check-------

# GM
ROC_GM <- roc(
  response = df_Luijten25_DRP09$GM,
  predictor = mtx_pred_all[,1] %>% as_vector()
)

# CA
ROC_CA <- roc(
  response = df_Luijten25_DRP09$GM,
  predictor = mtx_pred_all[,1] %>% as_vector()
)

## In Vivo  Endpoint by in vitro AUC--------
ROC_GM_vt <- roc(
  response = df_Luijten25_DRP09$GM,
  predictor = mtx_pred_vt[,1] %>% as_vector()
)
ggroc(ROC_GM_vt)


ROC_CA_vt <- roc(
  response = df_Luijten25_DRP09$CA,
  predictor = mtx_pred_vt[,2] %>% as_vector()
)
ggroc(ROC_CA_vt)

## In Vivo  Endpoint by in vitro AUC--------
ROC_GM_vv <- roc(
  response = df_Luijten25_DRP09$GM,
  predictor = mtx_pred_vv[,1] %>% as_vector()
)
ggroc(ROC_GM_vv)

ROC_CA_vv <- roc(
  response = df_Luijten25_DRP09$CA,
  predictor = mtx_pred_vv[,2] %>% as_vector()
)
ggroc(ROC_CA_vv)

## In Vivo  Endpoint in vitro vs in vivo vs all-------

# plotting
g2 <- ggroc(list(         # Gene Mutation
  "In Vitro"=ROC_GM_vt,"In Vivo"=ROC_GM_vv,
  "In Vitro & In Vivo"=ROC_GM
))+
  geom_segment(aes(x=1,xend=0,y=0,yend=1),col="grey",linetype="dashed")+
  labs(caption = "AUC:                                     \n In Vitro=0.6884 \n In Vivo=0.9832 \n In Vitro & In Vivo=0.9213",
       subtitle = "Gene Mutation")+
  theme_classic(base_size = 20)+
  theme(legend.position = "bottom")

g3 <- ggroc(list(         # Chromosome abberation
  "In Vitro"=ROC_CA_vt,"In Vivo"=ROC_CA_vv,
  "In Vitro & In Vivo"=ROC_CA
))+
  geom_segment(aes(x=1,xend=0,y=0,yend=1),col="grey",linetype="dashed")+
  labs(caption = "AUC:                                        \n In Vitro=0.7295 \n In Vivo=0.9931 \n In Vivo & In Vivo=0.9365",
       subtitle = "Chromosome Aberration")+
  theme_classic(base_size = 20)+
  theme(legend.position = "bottom")

## Parallel coordinates GM------
head(results_GM)
colnames(results_GM)

results_plotting <- results_GM %>%
  mutate(
    across(c(1:3, 10), ~ case_when(
      . == 1   ~ "+",
      . == 0   ~ "-",
      . == -1  ~ NA_character_,
      TRUE     ~ as.character(.)
    )),
    y0 = ifelse(N_prev == 0, 0.5,y0)
  ) %>%
  mutate(new_Value = as.character(new_Value)) %>%
  filter(!is.na(new_Assay))

# For color mapping: assign unique color per assay
assay_colors <- c(
  Ames = "#E41A1C",
  MLA = "#377EB8",
  HPRT = "#4DAF4A"
)

# Best combo for TRUE: 0.5, HPRT=0.916, MLA=0.981, Ames=0.996
# Best combo for FALSE: 0.5, MLA = 0.109, HPRT = 0.0196, Ames=0.0084

# Plot
fontsize <- 6
g4 <- ggplot(results_plotting) +
  # facet_wrap(~Prior,ncol=1)+
  geom_segment(aes(x = N_prev, y = y0, xend = N, yend = y,
                   color = new_Assay, linetype = new_Value),
               linewidth  = 0.8,alpha = 0.5) +
  # text labels for best combo
  annotate(geom = "text",label="None,P=0.5",
           x=0.15,y=0.5,size=fontsize)+
  annotate(geom = "text",label="Add HPRT=+ \n P=0.916",
           x=1,y=0.97,size=fontsize)+
  annotate(geom = "text",label="Add MLA=+ \n P=0.981",
           x=2,y=1.03,size=fontsize)+
  annotate(geom = "text",label="Add Ames=+ \n P=0.996",
           x=2.8,y=1.05,size=fontsize)+
  annotate(geom = "text",label="Add MLA=- \n P=0.109",
           x=1,y=0.04,size=fontsize)+
  annotate(geom = "text",label="Add HPRT=- \n P=0.020",
           x=2,y=-0.02,size=fontsize)+
  annotate(geom = "text",label="Add Ames=- \n P=0.008",
           x=2.8,y=-0.04,size=fontsize)+
  scale_color_manual(
    values = assay_colors
  )+
  scale_linetype_manual(
    values=c("+"="solid","-"="dotted")
  )+
  labs(
    x = "Number of assays tested",
    y = "P(GM=TRUE)",
    color="Assay",
    linetype = "Assay Result",
    caption = "At each number of assays tested\n Text labels are upper and lower bounds of P(GM=TRUE)\n And the newly added assay"
  ) +
  theme_minimal(base_size = 20)+
  theme(legend.position = "right")

## Parallel coordinates CA------
head(results_CA)
colnames(results_CA)

results_plotting <- results_CA %>%
  mutate(across(c(1:3, 9), ~ case_when(
    . == 1   ~ "+",
    . == 0   ~ "-",
    . == -1  ~ NA_character_,
    TRUE     ~ as.character(.)
  )),
  y0 = ifelse(N_prev == 0, 0.5,y0)) %>%
  mutate(new_Value = as.character(new_Value),
         y = as.numeric(y)) %>%
  filter(!is.na(new_Assay))

# For color mapping: assign unique color per assay
assay_colors <- c(
  MNvt = "#984EA3",
  CAvt = "#FF7F00"
)

# Best combo for TRUE: 0.5, MNvt=0.832, CAvt=0.941
# Best combo for FALSE: 0.5, CAvt=0.132, MNvt=0.027

# Plot
fontsize <- 6
g5 <- ggplot(results_plotting) +
  # facet_wrap(~Prior,ncol=1)+
  geom_segment(aes(x = N_prev, y = y0, xend = N, yend = y,
                   color = new_Assay, linetype = new_Value),
               linewidth  = 0.8,alpha = 0.5) +
  # text labels for best combo
  annotate(geom = "text",label="None,P=0.5",
           x=0.1,y=0.5,size=fontsize)+
  annotate(geom = "text",label="Add MNvt=+ \n P=0.832",
           x=1,y=0.92,size=fontsize)+
  annotate(geom = "text",label="Add CAvt=+ \n P=0.941",
           x=1.9,y=1.01,size=fontsize)+
  annotate(geom = "text",label="Add CAvt=- \n P=0.132",
           x=1,y=0.08,size=fontsize)+
  annotate(geom = "text",label="Add MNvt=- \n P=0.027",
           x=1.9,y=0.00,size=fontsize)+ 
  scale_color_manual(
    values = assay_colors
  )+
  scale_linetype_manual(
    values=c("+"="solid","-"="dotted")
  )+
  labs(
    x = "Number of assays tested",
    y = "P(CA=TRUE)",
    color="Assay",
    linetype = "Assay Result",
    caption = "At each number of assays tested\n Text labels are upper and lower bounds of P(CA=TRUE)\n And the newly added assay"
  ) +
  theme_minimal(base_size = 20)+
  theme(legend.position = "right")

## Parallel coordinates RC------
head(results_RC)
colnames(results_RC)

results_plotting <- results_RC %>%
  mutate(across(c(1:5, 12), ~ case_when(
    . == 1   ~ "+",
    . == 0   ~ "-",
    . == -1  ~ NA_character_,
    TRUE     ~ as.character(.)
  ))) %>%
  mutate(new_Value = as.character(new_Value)) %>%
  filter(!is.na(new_Assay))

# For color mapping: assign unique color per assay
assay_colors <- c(
  Ames = "#E41A1C",
  MLA = "#377EB8",
  HPRT = "#4DAF4A",
  MNvt = "#984EA3",
  CAvt = "#FF7F00"
)

# Best combo for TRUE: 535, HPRT=0.910,MNvt=0.977,MLA=0.994,CAvt=0.998,Ames=0.999
# Best combo for FALSE: 0.535, MLA = 0.128, CAvt = 0.028, HPRT=0.008,MNvt=0.003,Ames=0.001

# Plot
fontsize <- 6
g6 <- ggplot(results_plotting) +
  # facet_wrap(~Prior,ncol=1)+
  geom_segment(aes(x = N_prev, y = y0, xend = N, yend = y,
                   color = new_Assay, linetype = new_Value),
               linewidth  = 0.8,alpha = 0.5) +
  # text labels for best combo
  annotate(geom = "text",label="None,P=0.535",
           x=0.4,y=0.55,size=fontsize)+
  annotate(geom = "text",label="Add HPRT=+ \n P=0.910",
           x=1,y=0.93,size=fontsize)+
  annotate(geom = "text",label="Add MNvt=+ \n P=0.977",
           x=2,y=1.01,size=fontsize)+
  annotate(geom = "text",label="Add MLA=+ \n P=0.994",
           x=3,y=1.03,size=fontsize)+
  annotate(geom = "text",label="Add CAvt=+ \n P=0.998",
           x=4,y=1.03,size=fontsize)+
  annotate(geom = "text",label="Add Ames=+ \n P=0.999",
           x=5,y=1.03,size=fontsize)+
  annotate(geom = "text",label="Add MLA=- \n P=0.128",
           x=1,y=0.09,size=fontsize)+
  annotate(geom = "text",label="Add CAvt=- \n P=0.028",
           x=2,y=0.00,size=fontsize)+
  annotate(geom = "text",label="Add HPRT=- \n P=0.008",
           x=3,y=-0.03,size=fontsize)+
  annotate(geom = "text",label="Add MNvt=- \n P=0.003",
           x=4,y=-0.03,size=fontsize)+
  annotate(geom = "text",label="Add Ames=- \n P=0.001",
           x=5,y=-0.03,size=fontsize)+
  
  scale_color_manual(
    values = assay_colors
  )+
  scale_linetype_manual(
    values=c("+"="solid","-"="dotted")
  )+
  labs(
    x = "Number of assays tested",
    y = "P(RC=TRUE)",
    color="Assay",
    linetype = "Assay Result",
    caption = "At each number of assays tested\n Text labels are upper and lower bounds of P(RC=TRUE)\n And the newly added assay"
  ) +
  theme_minimal(base_size = 20)+
  theme(legend.position = "right")
