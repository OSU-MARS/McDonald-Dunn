### Overview
This repo contains research code pertaining to the McDonald-Dunn Research Forest. Since it doesn't contain the data (it's not something this research group has permission to release) it's unlikely to be of interest unless you're working on the forest, peer reviewing journal publication code, or are interested in replicating similar work elsewhere.

### Dependencies
[R](https://www.r-project.org/) is the primary tool used here, mainly via [RStudio](https://www.rstudio.com/) Desktop, and thus most code is in .R files.

### Height and diameter modeling
R scripts are in trees/height-diameter and require forest inventory data from out of repo (it's not ours to redistribute). Basic workflow is to run the scripts in the background jobs directory to fit models and then analyze the results via results.R. See notes at the top of results.R.