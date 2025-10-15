# processor  cross validation   workers   height runtime   DBH runtime
# 9900X      10x10              10        3.4 + 1.3 min    2.3 min       
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential)
future::plan(future::multisession, workers = 10) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/other.R")
warnings()
dplyr::last_dplyr_warnings() 
print(paste0("other species job ran for ", format(Sys.time() - jobStartTime), "."))