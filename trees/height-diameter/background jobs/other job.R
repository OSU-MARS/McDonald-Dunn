# processor  cross validation   workers   runtime, hours
# 9950X      10x50              8         ~4
#            5x100              8         ~4
#            2x250              8         ~4
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")

htDiaOptions$folds = 10
htDiaOptions$repetitions = 50
htDiaOptions$crossValidation = "blockedByStand" # blockedByStand, randomByCruiseRecord

message(paste0("Other species ", htDiaOptions$folds, "x", htDiaOptions$repetitions, " cross validation ", htDiaOptions$crossValidation, "..."))

otherOptions = tibble(fitHeight = TRUE,
                      fitHeightMixed = TRUE,
                      fitDbh = TRUE,
                      fitDbhMixed = TRUE,
                      recalcResultCollections = TRUE,
                      recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential) # mixed effects may need to run sequential at lower numbers of cross validation folds due to internal errors in future
future::plan(future::multisession, workers = 8)

source("trees/height-diameter/other.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("other species job ran for ", format(Sys.time() - jobStartTime), "."))