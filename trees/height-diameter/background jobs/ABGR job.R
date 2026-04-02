# processor  cross validation   workers   runtime
# 9950X      10x10              10        19 m
#            5x200              8         1: + 1:19 + 1:13 + 48
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")

htDiaOptions$folds = 5
htDiaOptions$repetitions = 100
htDiaOptions$crossValidation = "randomByCruiseRecord" # blockedByStand, randomByCruiseRecord

message(paste0(htDiaOptions$folds, "x", htDiaOptions$repetitions, " cross validation ", htDiaOptions$crossValidation, "..."))

abgrOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbhPrimary = TRUE,
                     fitDbhMixed = TRUE, 
                     recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential) # mixed effects may need to run sequential at lower numbers of cross validation folds due to internal errors in future
future::plan(future::multisession, workers = 8)

source("trees/height-diameter/ABGR.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("ABGR job ran for ", format(Sys.time() - jobStartTime), "."))