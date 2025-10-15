# 9950X runtimes with 10x10 cross validation and 2 workers
# height: 
# DBH: 
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
#future::plan(future::sequential)
future::plan(future::multisession, workers = 10) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/other.R")
warnings()
print(paste0("other species job ran for ", format(Sys.time() - jobStartTime), "."))