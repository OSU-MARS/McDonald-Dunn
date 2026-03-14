# processor  cross validation   workers   runtime
# 9950X      10x10              10        20 min
#            5x200              8         1:41 + 1:43 + 1:15 + 27 = 5h11m
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
message(paste0(htDiaOptions$folds, "x", htDiaOptions$repetitions, " cross validation ", htDiaOptions$crossValidation, "..."))

acmaOptions = tibble(fitHeight = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbh = TRUE,
                     fitDbhMixed = TRUE,
                     recalcPreferredModels = FALSE)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 8)

source("trees/height-diameter/ACMA3.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("ACMA job ran for ", format(Sys.time() - jobStartTime), "."))