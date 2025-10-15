# processor  cross validation   workers   height runtime   DBH runtime
# 9900X      10x10              10        4.6 + 3.7 min          
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 10) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/ACMA3.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("ACMA job ran for ", format(Sys.time() - jobStartTime), "."))