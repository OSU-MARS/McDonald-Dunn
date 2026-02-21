library(dplyr)

# processor  cross validation   workers   runtime
# 9950X      10x10              10        12 min
jobStartTime = Sys.time()
otherOptions = tibble(fitHeight = TRUE, 
                      fitHeightMixed = TRUE,
                      fitDbh = TRUE,
                      fitDbhMixed = TRUE,
                      recalcPreferredModels = FALSE)

source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential) # mixed effects may need to run sequential at lower numbers of cross validation folds due to internal errors in future
future::plan(future::multisession, workers = 10)

source("trees/height-diameter/other.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("other species job ran for ", format(Sys.time() - jobStartTime), "."))