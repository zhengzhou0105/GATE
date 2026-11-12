message("Please place the files according to the layout described in the README")

source("data.R")
message("The data is ready!")

message("Begin sampling for Bayesian inference. This may take a while...")
source("modeling.R")
message("Sampling complete!")

source("decision_analysis.R")

message("Figures are saved as objects. Save them on your local drive if desired.")
source("plotting.R")