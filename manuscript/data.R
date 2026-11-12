# Settings-----
dir_data <- file.path("/shared/data")
dir_fig <- file.path("/shared/figs")

source(file.path("/manuscript","utilities.R"))

# dir_fit <- file.path(getwd(),"fit")
Pvalue <- "TRUE"
Nvalue <- "FALSE"
Evalue <-  Pvalue


# Read data----
df_individual <- read_excel(
  file.path(dir_data,"Luijte25_TRAID09.xlsx"),
  sheet = "Luijte25_TRAID09",
  col_types = c("text","text",rep("logical",13))) 

df_ext <- read_excel(
  file.path(dir_data,"Luijte25_TRAID09.xlsx"),
  sheet = "external") %>%
  rename(
    vt_MN = MNvt,
    vv_MN = MNvv,
    vt_CA = CAvt,
    vv_CA = CAvv,
    comet = Comet)

Expand <- read_csv(
  file.path(dir_data,"vt_Expand.csv")
)

# Processing----

## Extract individual chemical data
df_Luijten25_DRP09 <- df_individual %>%
  # change logical to 0,1, NA to -1
  mutate(
    across(3:13,
           ~ replace(as.integer(.),is.na(.),-1)
    )
  )

name_vt_GM <- c("Ames","MLA","HPRT") 
name_vv_GM <- c("TGR","comet","PigA")
name_vt_CA <- c("vt_MN","vt_CA")
name_vv_CA <- c("vv_MN","vv_CA")
name_RC <- "CarcR"

assay_GM_vt = df_Luijten25_DRP09 %>% select(any_of(name_vt_GM)) %>%
  mutate(across(everything(),~as.integer(.))) %>%
  as.matrix()
assay_GM_vv = df_Luijten25_DRP09 %>% select(any_of(name_vv_GM)) %>%
  mutate(across(everything(),~as.integer(.))) %>%
  as.matrix()
assay_CA_vt = df_Luijten25_DRP09 %>% select(any_of(name_vt_CA)) %>%
  mutate(across(everything(),~as.integer(.))) %>%
  as.matrix()
assay_CA_vv = df_Luijten25_DRP09 %>% select(any_of(name_vv_CA)) %>%
  mutate(across(everything(),~as.integer(.))) %>%
  as.matrix()
RC = df_Luijten25_DRP09 %>% select(all_of(name_RC)) %>%
  mutate(across(everything(),~as.integer(.))) %>%
  unlist()
N_1 <- 200                      # external validation data
N_0 <- 200

## Compute in vivo endpoint using alternative "precautionary majority" rule
df_Luijten25_DRP09_alt <- df_individual %>%
  mutate(
    GM = case_when(
      is.na(TGR) & is.na(PigA) & is.na(comet) ~ NA,
      rowSums(cbind(TGR == FALSE, PigA == FALSE, comet == FALSE), na.rm = TRUE) >
        rowSums(cbind(TGR == TRUE, PigA == TRUE, comet == TRUE), na.rm = TRUE) ~ FALSE,
      (TGR | PigA | comet) ~ TRUE,
      .default = FALSE
    ),
    CA = case_when(
      is.na(MNvv) & is.na(CAvv) ~ NA,
      (!MNvv & !CAvv) ~ FALSE,
      (MNvv | CAvv) ~ TRUE,
      .default = FALSE
    )
  )

## Prepare data as a list for stan
k_plus_1_GM <- N_1 * df_ext[1,c(name_vt_GM,name_vv_GM)] %>% as.numeric() / 100 %>% as.integer()
k_minus_0_GM <- N_0 * df_ext[2,c(name_vt_GM,name_vv_GM)] %>% as.numeric() / 100 %>% as.integer()
k_plus_1_CA <- N_1 * df_ext[1,c(name_vt_CA,name_vv_CA)] %>% as.numeric() / 100 %>% as.integer()
k_minus_0_CA <- N_0 * df_ext[2,c(name_vt_CA,name_vv_CA)] %>% as.numeric() / 100 %>% as.integer()

data_list <- list(
  N = nrow(df_Luijten25_DRP09),
  A = length(name_vt_GM),
  B = length(name_vv_GM),
  C = length(name_vt_CA),
  D = length(name_vv_CA),
  assay_GM_vt = assay_GM_vt,
  assay_GM_vv = assay_GM_vv,
  assay_CA_vt = assay_CA_vt,
  assay_CA_vv = assay_CA_vv,
  # RC = RC,
  # RC = RC_NGC,
  N_1 = N_1,
  N_0 = N_0,
  k_plus_1_GM = k_plus_1_GM,
  k_minus_0_GM = k_minus_0_GM,
  k_plus_1_CA = k_plus_1_CA,
  k_minus_0_CA = k_minus_0_CA
)
