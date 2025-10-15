# 9950X runtimes with 10x10 cross validation and eight workers
# height: 
# DBH: 
#
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 8) # increase worker count for mixed effects DBH GAMs

source("trees/height-diameter/PSME.R")
warnings()
print(paste0("PSME job ran for ", format(Sys.time() - jobStartTime), "."))