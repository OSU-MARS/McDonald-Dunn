# processor  cross validation   workers   runtime, hours
# 9950X      10x50              8         10.8
#            5x100              8         10.8
#            2x250              8         10.7
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")

htDiaOptions$folds = 2
htDiaOptions$repetitions = 250
htDiaOptions$crossValidation = "blockedByStand" # blockedByStand, randomByCruiseRecord

message(paste0("Bigleaf maple ", htDiaOptions$folds, "x", htDiaOptions$repetitions, " cross validation ", htDiaOptions$crossValidation, "..."))

acmaOptions = tibble(fitHeight = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbh = TRUE,
                     fitDbhMixed = TRUE,
                     recalcResultCollections = TRUE,
                     recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 8)

source("trees/height-diameter/ACMA3.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("ACMA job ran for ", format(Sys.time() - jobStartTime), "."))