library(dplyr)

# processor  cross validation   workers   runtime
# 9950X      10x10              10        20 min          
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
acmaOptions = tibble(fitHeight = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbh = TRUE,
                     fitDbhMixed = TRUE,
                     recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 10)

source("trees/height-diameter/ACMA3.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("ACMA job ran for ", format(Sys.time() - jobStartTime), "."))