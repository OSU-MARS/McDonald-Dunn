# processor  cross validation   workers   runtime
# 9950X      10x10              10        12 min
#            5x200              8         23 + 46 + 38 + 59 = 2h46m
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
#htDiaOptions$crossValidation = "randomByCruiseRecord"
message(paste0(htDiaOptions$folds, "x", htDiaOptions$repetitions, " cross validation ", htDiaOptions$crossValidation, "..."))

otherOptions = tibble(fitHeight = TRUE,
                      fitHeightMixed = TRUE,
                      fitDbh = TRUE,
                      fitDbhMixed = TRUE,
                      recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential) # mixed effects may need to run sequential at lower numbers of cross validation folds due to internal errors in future
future::plan(future::multisession, workers = 8)

source("trees/height-diameter/other.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("other species job ran for ", format(Sys.time() - jobStartTime), "."))