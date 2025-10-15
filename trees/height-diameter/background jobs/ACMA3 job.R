# 9950X runtimes with 10x10 cross validation and 2 workers
# height: 
# DBH: 
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 2) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/ACMA3.R")
warnings()
print(paste0("ACMA job ran for ", format(Sys.time() - jobStartTime), "."))