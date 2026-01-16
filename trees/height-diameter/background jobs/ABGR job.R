# processor  cross validation   workers   height runtime   DBH runtime
# 9900X      10x10              10        
jobStartTime = Sys.time()
abgrOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightMixed = FALSE,
                     fitDbhPrimary = FALSE,
                     fitDbhMixed = FALSE, 
                     recalcPreferredModels = FALSE)

source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential) # mixed effects may need to run sequential at lower numbers of cross validation folds due to internal errors in future
future::plan(future::multisession, workers = 10)

source("trees/height-diameter/ABGR.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("ABGR job ran for ", format(Sys.time() - jobStartTime), "."))