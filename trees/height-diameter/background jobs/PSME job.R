library(dplyr)

# processor  cross validation   workers   height runtime   DBH runtime
# 9900X      10x10              10        9.2 min          2.5 + 2.9 min
jobStartTime = Sys.time()
psmeOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbhPrimary = TRUE,
                     fitDbhMixed = TRUE,
                     fitPhysioGams = TRUE,
                     recalcPreferredModels = FALSE)

source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 10)

source("trees/height-diameter/PSME.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("PSME job ran for ", format(Sys.time() - jobStartTime), "."))