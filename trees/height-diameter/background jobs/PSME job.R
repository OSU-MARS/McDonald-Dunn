# processor  cross validation   workers   height runtime   DBH runtime
# 9900X      10x10              10        9.2 min          2.5 + 2.9 min
#            5x200               8        3:43 + 6:52      3:49 + 5:30
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")

htDiaOptions$folds = 5
htDiaOptions$repetitions = 100
htDiaOptions$crossValidation = "randomByCruiseRecord" # blockedByStand, randomByCruiseRecord

message(paste0(htDiaOptions$folds, "x", htDiaOptions$repetitions, " cross validation ", htDiaOptions$crossValidation, "..."))

psmeOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbhPrimary = TRUE,
                     fitDbhMixed = TRUE,
                     fitPhysioGams = TRUE,
                     recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 8)

source("trees/height-diameter/PSME.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("PSME job ran for ", format(Sys.time() - jobStartTime), "."))