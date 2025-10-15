library(dplyr)
library(ggplot2)
library(ggtrendline)
library(gslnls)
library(mgcv)
library(patchwork)
library(readxl)
library(sf)
library(stringr)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             axis.title = element_text(size = 9),
                             legend.title = element_text(size = 9),
                             panel.border = element_blank(),
                             plot.title = element_text(size = 9)))
plotLetters = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)")

htDiaOptions = tibble(folds = 2,
                      repetitions = 2,
                      includeInvestigatory = FALSE, # default to excluding plotting and other add ons in species scripts
                      includeSetup = FALSE,
                      rangerThreads = 0.5 * future::availableCores(), 
                      retainModelThreshold = 5) # cross validation retains model objects if folds * repetitions is less than or equal to this threshold, e.g. 25 = retaining models up to and including 5x5 cross validation

create_fit_statistics = function(name, fittingMethod, fitSet)
{
  # https://github.com/tidyverse/tibble/issues/1569
  return(tibble::tibble_row(fitSet = fitSet, fitting = fittingMethod, name = name, 
                            effectiveDegreesOfFreedom = NA_real_, 
                            coefficients = vector("list", length = 1),
                            training = tibble(.rows = 1), validation = tibble(.rows = 1),
                            isConverged = NA_real_, significant = NA_real_,
                            fitTimeInS = NA_real_))
}

# currently supports fixed effect GAMs and mixed GAMs with a single random smooth
# mgcv::gam()'s default is family = gaussian() but family = Gamma() can be more accurate, partly due to preventing prediction of negative
# heights. However, mgcv's fitting of Gaussian random effects seems to compose poorly, if at all, with out of dataset prediction using 
# the gamma family.
#
# number of predictors   minimum k
# 1                      4
# 2                      4
# 3                      11
# 4                      16
# 5                      57
# 6                      85
# 7                      331
fit_gam = function(name, formula, data, family = gaussian(), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, randomTermIndex = NA_integer_, threads = 1, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = formula[2] # displays as height or dbh but compares as height() or dbh()
  
  if (responseVariable == "height()")
  {
    responseVariable = "height"
    allFitWeights = data$dbhWeight
  }
  else
  {
    if (responseVariable != "dbh()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    responseVariable = "DBH"
    allFitWeights = data$heightWeight
  }
  if (significant == FALSE)
  {
    # don't fit GAMs which are marked as not significant
    # This bypass is, if not required, desirable for smooths whose degrees of freedom exceed the amount of data, meaning they can't be
    # fit.
    nonsignificantFitStats = create_fit_statistics(name = name, fittingMethod = "gam", fitSet = "primary")
    nonsignificantFitStats$significant = significant
    return(nonsignificantFitStats)
  }
  
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using gam()..."))
  progressBar = progressr::progressor(steps = folds * repetitions)
  if ((folds == 1) & (repetitions == 1))
  {
    startFit = Sys.time()
    if (responseVariable == "height") 
    { 
      gamModel = gam(formula = formula, data = data, family = family, method = "REML", select = TRUE, weights = dbhWeight, nthreads = threads)
    } else { 
      gamModel = gam(formula = formula, data = data, family = family, method = "REML", select = TRUE, weights = heightWeight, nthreads = threads)
    }
    effectiveDegreesOfFreedom = sum(gamModel$edf)
    if (is.na(randomTermIndex))
    {
      validationPrediction = fitted(gamModel)
    } else {
      validationPrediction = predict(gamModel, validationData, exclude = paste0("s(", attr(gamModel$terms, "term.labels")[randomTermIndex], ")"), newdata.guaranteed = TRUE) # remove random smooth for consistency with cross validated cases
    }
    allFitStats = get_fit_statistics(name = name, fittingMethod = "gam", responseVariable = responseVariable, 
                                     trainingData = data, trainingPrediction = fitted(gamModel), effectiveDegreesOfFreedom = effectiveDegreesOfFreedom, 
                                     validationData = data, validationPrediction = validationPrediction,
                                     significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
    allFitStats$fitTimeInS = get_elapsed_time(startFit)
    progressBar()
    return(get_fit_return_value(allFit, allFitStats, returnModel))
  }

  fitFunction = function(dataFold)
  {
    startFit = Sys.time()
    trainingData = rsample::analysis(dataFold)
    if (responseVariable == "height") 
    { 
      gamModel = gam(formula = formula, data = trainingData, family = family, method = "REML", select = TRUE, weights = dbhWeight, nthreads = threads)
    } else { 
      gamModel = gam(formula = formula, data = trainingData, family = family, method = "REML", select = TRUE, weights = heightWeight, nthreads = threads)
    }
    effectiveDegreesOfFreedom = sum(gamModel$edf)
    
    validationData = rsample::assessment(dataFold)
    if (is.na(randomTermIndex))
    {
      validationPrediction = predict(gamModel, validationData)
    } else {
      validationPrediction = predict(gamModel, validationData, exclude = paste0("s(", attr(gamModel$terms, "term.labels")[randomTermIndex], ")"), newdata.guaranteed = TRUE) # random effect must be removed to make predict outside of the training dataset
    }
    
    fitStatistics = get_fit_statistics(name = name, fittingMethod = "gam", responseVariable = responseVariable, 
                                       trainingData = trainingData, trainingPrediction = fitted(gamModel), effectiveDegreesOfFreedom = effectiveDegreesOfFreedom, 
                                       validationData = validationData, validationPrediction = validationPrediction,
                                       significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
    # fitStatistics$coefficients not useful to set as GAM object is required
    fitStatistics$fitTimeInS = get_elapsed_time(startFit)
    progressBar()
    return(get_fit_return_value(gamModel, fitStatistics, returnModel))
  }
  
  splitsAndFits = rsample::group_vfold_cv(data, v = folds, repeats = repetitions, group = stand) %>% mutate(fit = furrr::future_map(splits, fitFunction, .options = furrr::furrr_options(seed = TRUE)))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

fit_gsl_nls = function(name, formula, data, start, control = gsl_nls_control(maxiter = 100), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = formula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using gsl_nls()..."))
  progressBar = progressr::progressor(steps = folds * repetitions)
  
  # gsl_nls() 1.4.1 documentation states weights should be a vector but this is incorrect, the name of column in data is required
  if (responseVariable == "height()")
  {
    responseVariable = "height"
    allFitNonlinear = gsl_nls(fn = formula, data = data, start = start, weights = dbhWeight, control = control)
  } else {
    if (responseVariable != "dbh()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    responseVariable = "DBH"
    allFitNonlinear = gsl_nls(fn = formula, data = data, start = start, weights = heightWeight, control = control)
  }
    
  startFit = Sys.time()
  if ((folds == 1) & (repetitions == 1))
  {
    allFitStats = get_fit_statistics(name = name, fittingMethod = "gsl_nls", responseVariable = responseVariable, 
                                     trainingData = data, trainingPrediction = fitted(allFitNonlinear), effectiveDegreesOfFreedom = length(coef(allFitNonlinear)) + 1, 
                                     validationData = data, validationPrediction = fitted(allFitNonlinear),
                                     significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
    allFitStats$coefficients[[1]] = bind_rows(coef(allFitNonlinear))
    allFitStats$fitTimeInS = get_elapsed_time(startFit)
    progressBar()
    return(get_fit_return_value(allFitNonlinear, allFitStats, returnModel))
  }
  
  allFitParameters = coef(allFitNonlinear)
  fitFunction = function(dataFold)
  {
    startFit = Sys.time()
    trainingData = rsample::analysis(dataFold)
    if (responseVariable == "height")
    {
      nonlinearModel = gsl_nls(fn = formula, data = trainingData, start = allFitParameters, weights = dbhWeight, control = control)
    } else {
      nonlinearModel = gsl_nls(fn = formula, data = trainingData, start = allFitParameters, weights = heightWeight, control = control)
    }
    
    validationData = rsample::assessment(dataFold)
    fitStatistics = get_fit_statistics(name = name, fittingMethod = "gsl_nls", responseVariable = responseVariable, 
                                       trainingData = trainingData, trainingPrediction = fitted(nonlinearModel), effectiveDegreesOfFreedom = length(coef(nonlinearModel)) + 1, 
                                       validationData = validationData, validationPrediction = predict(nonlinearModel, validationData),
                                       significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
    fitStatistics$coefficients[[1]] = bind_rows(coef(nonlinearModel))
    fitStatistics$fitTimeInS = get_elapsed_time(startFit)
    progressBar()
    return(get_fit_return_value(nonlinearModel, fitStatistics, returnModel))
  }

  splitsAndFits = rsample::group_vfold_cv(data, v = folds, repeats = repetitions, group = stand) %>% mutate(fit = furrr::future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

fit_lm = function(name, formula, data, folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = formula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using lm()..."))
  progressBar = progressr::progressor(steps = folds * repetitions)
  
  # separate height and diameter all fit cases and cross validation setups since height prediction uses offset
  if (responseVariable == "height()")
  {
    if ((folds == 1) & (repetitions == 1))
    {
      startFit = Sys.time()
      allFit = lm(formula = formula, data = data, offset = breastHeight, weights = dbhWeight)
      allFitStats = get_fit_statistics(name = name, fittingMethod = "lm", responseVariable = "height", 
                                       trainingData = data, trainingPrediction = fitted(allFit), effectiveDegreesOfFreedom = length(coef(allFit)) + 1, 
                                       validationData = data, validationPrediction = fitted(allFit),
                                       significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$coefficients[[1]] = bind_rows(coef(allFit))
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = rsample::analysis(dataFold)
      linearModel = lm(formula = formula, data = trainingData, offset = breastHeight, weights = dbhWeight) # specifying breastHeight as a constant fails with 'variable lengths differ'
      validationData = rsample::assessment(dataFold)
      fitStatistics = get_fit_statistics(name = name, fittingMethod = "lm", responseVariable = "height", 
                                         trainingData = trainingData, trainingPrediction = fitted(linearModel), effectiveDegreesOfFreedom = length(coef(linearModel)) + 1,
                                         validationData = validationData, validationPrediction = predict(linearModel, validationData), 
                                         significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      fitStatistics$coefficients[[1]] = bind_rows(coef(linearModel))
      fitStatistics$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(linearModel, fitStatistics, returnModel))
    }
  }
  else
  {
    if (responseVariable != "dbh()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    if ((folds == 1) & (repetitions == 1))
    {
      startFit = Sys.time()
      allFit = lm(formula = formula, data = data, weights = heightWeight)
      allFitStats = get_fit_statistics(name = name, fittingMethod = "lm", responseVariable = "DBH", 
                                       trainingData = data, trainingPrediction = fitted(allFit), effectiveDegreesOfFreedom = length(coef(allFit)) + 1, 
                                       validationData = data, validationPrediction = fitted(allFit), 
                                       significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$coefficients[[1]] = bind_rows(coef(allFit))
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = rsample::analysis(dataFold)
      linearModel = lm(formula = formula, data = trainingData, weights = heightWeight)
      validationData = rsample::assessment(dataFold)
      fitStatistics = get_fit_statistics(name = name, fittingMethod = "lm", responseVariable = "DBH", 
                                         trainingData = trainingData, trainingPrediction = fitted(linearModel), effectiveDegreesOfFreedom = length(coef(linearModel)) + 1, 
                                         validationData = validationData, validationPrediction = predict(linearModel, validationData), 
                                         significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      fitStatistics$coefficients[[1]] = bind_rows(coef(linearModel))
      fitStatistics$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(linearModel, fitStatistics, returnModel))
    }
  }
  
  splitsAndFits = rsample::group_vfold_cv(data, v = folds, repeats = repetitions, group = stand) %>% mutate(fit = furrr::future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

fit_nlme = function(name, modelFormula, data, fixedFormula, randomFormula, start, control = nlmeControl(maxIter = 100), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = modelFormula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using nlme()..."))
  progressBar = progressr::progressor(steps = folds * repetitions)
  
  if (responseVariable == "height()")
  {
    responseVariable = "height"
    # https://stackoverflow.com/questions/11778773/using-predict-in-a-function-call-with-nlme-objects-and-a-formula
    # Debatable if varFixed() should include treeCount. While variance is independent of the count nlme()'s offers no separate
    # mechanism capturing the number of observations. This does not appear to be an issue with varFixed() but might affect cases
    # such as varPower() where variance fitting is done within nlme() rather than externally as performed here.
    allFitNonlinear = do.call(nlme, list(model = modelFormula, data = data, fixed = fixedFormula, random = randomFormula, start = start, weights = varFixed(~1/dbhWeight), control = control))
  } else {
    if (responseVariable != "dbh()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    responseVariable = "DBH"
    allFitNonlinear = do.call(nlme, list(model = modelFormula, data = data, fixed = fixedFormula, random = randomFormula, start = start, weights = varFixed(~1/heightWeight), control = control))
  }
  
  startFit = Sys.time()
  if ((folds == 1) & (repetitions == 1))
  {
    allFitStats = get_fit_statistics(name = name, fittingMethod = "nlme", responseVariable = responseVariable, 
                                     trainingData = data, trainingPrediction = fitted(allFitNonlinear), effectiveDegreesOfFreedom = length(coef(allFitNonlinear)) + 1, 
                                     validationData = data, validationPrediction = predict(allFitNonlinear, newdata = data, level = 0),
                                     significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
    allFitStats$coefficients[[1]] = bind_rows(coef(allFitNonlinear))
    allFitStats$fitTimeInS = get_elapsed_time(startFit)
    progressBar()
    return(get_fit_return_value(allFitNonlinear, allFitStats, returnModel))
  }
  
  allFitParameters = allFitNonlinear$coefficients$fixed
  fitFunction = function(dataFold)
  {
    startFit = Sys.time()
    trainingData = rsample::analysis(dataFold)
    if (responseVariable == "height")
    {
      nonlinearModel = do.call(nlme, list(model = modelFormula, data = trainingData, fixed = fixedFormula, random = randomFormula, start = allFitParameters, weights = varFixed(~1/dbhWeight), control = control))
    } else {
      nonlinearModel = do.call(nlme, list(model = modelFormula, data = trainingData, fixed = fixedFormula, random = randomFormula, start = allFitParameters, weights = varFixed(~1/heightWeight), control = control))
    }
    
    validationData = rsample::assessment(dataFold)
    fitStatistics = get_fit_statistics(name = name, fittingMethod = "nlme", responseVariable = responseVariable, 
                                       trainingData = trainingData, trainingPrediction = fitted(nonlinearModel), effectiveDegreesOfFreedom = length(coef(nonlinearModel)) + 1, 
                                       validationData = validationData, validationPrediction = predict(nonlinearModel, validationData, level = 0),
                                       significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
    fitStatistics$coefficients[[1]] = bind_rows(coef(nonlinearModel))
    fitStatistics$fitTimeInS = get_elapsed_time(startFit)
    progressBar()
    return(get_fit_return_value(nonlinearModel, fitStatistics, returnModel))
  }
  
  splitsAndFits = rsample::group_vfold_cv(data, v = folds, repeats = repetitions, group = stand) %>% mutate(fit = furrr::future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

get_cross_validation_return_value = function(splitsAndFits, returnModel)
{
  # vfold_cv() uses "Fold" as prefix in id2, group_vfold_cv() uses "Resample"
  if (returnModel)
  {
    # if full models are retained retain also the training-validation data splits for each model fit
    return(splitsAndFits %>% 
             mutate(repetition = as.numeric(str_replace(id, "Repeat", "")), 
                    fold = as.numeric(str_replace(id2, "Resample", ""))) %>%
             select(-id, -id2) %>%
             relocate(repetition, fold, fit, splits))
  }
  # if only model statistics are kept then drop the data splits and bind together each fit's single 
  # row statistics tibble into one tibble with kr rows (k-folds x r-repetitions = kr fits)
  return(bind_rows(splitsAndFits$fit) %>% 
           mutate(repetition = as.numeric(str_replace(splitsAndFits$id, "Repeat", "")), 
                  fold = as.numeric(str_replace(splitsAndFits$id2, "Resample", ""))) %>%
           relocate(repetition, fold))
}

get_elapsed_time = function(start)
{
  return(as.numeric(difftime(Sys.time(), start, units = "secs")))
}

get_fit_return_value = function(model, fitStatistics, returnModel)
{
  if (returnModel)
  {
    model$stats = fitStatistics
    return(model)
  } else {
    return(fitStatistics)
  }
}

get_fit_statistics = function(name, fittingMethod, responseVariable, trainingData, trainingPrediction, effectiveDegreesOfFreedom, validationData, validationPrediction, fitSet = "primary", isConverged = TRUE, significant = TRUE, tDegreesOfFreedom = 8)
{
  fitStatistics = create_fit_statistics(name = name, fittingMethod = fittingMethod, fitSet = fitSet)
  fitStatistics$effectiveDegreesOfFreedom = effectiveDegreesOfFreedom
  fitStatistics$isConverged = isConverged
  fitStatistics$significant = significant

  if (responseVariable == "DBH")
  {
    trainingHeightDiameterRatio = trainingData$height / (0.01 * trainingPrediction)
    fitStatistics$training = get_goodness_of_fit_statistics(trainingData, trainingPrediction, trainingData$dbh, trainingData$heightWeight, effectiveDegreesOfFreedom, trainingHeightDiameterRatio, 0.04, trainingData$dbhMax, tDegreesOfFreedom = tDegreesOfFreedom)
    
    validationHeightDiameterRatio = validationData$height / (0.01 * validationPrediction)
    fitStatistics$validation = get_goodness_of_fit_statistics(validationData, validationPrediction, validationData$dbh, validationData$heightWeight, effectiveDegreesOfFreedom, validationHeightDiameterRatio, 0.04, validationData$dbhMax, tDegreesOfFreedom = tDegreesOfFreedom)
    
    validationResiduals = validationPrediction - validationData$dbh
    dbhByHeightClass = validationData %>% mutate(residuals = validationResiduals) %>%
      group_by(heightClass) %>%
      summarize(n = sum(treeCount),
                nNaturalRegen = sum(treeCount * (isPlantation == FALSE)),
                nPlantation = sum(treeCount * isPlantation),
                meanBiasPerTree = sum(treeCount * residuals) / n,
                meanBiasPerTreePct = 100 * sum(treeCount * residuals / dbh) / n,
                meanDbh = sum(treeCount * dbh) / n,
                meanNaturalRegenDbh = if_else(nNaturalRegen > 0, sum(treeCount * dbh * (isPlantation == FALSE)) / nNaturalRegen, NA_real_),
                meanPlantationDbh = if_else(nPlantation > 0, sum(treeCount * dbh * isPlantation) / nPlantation, NA_real_),
                minPlantationNaturalRegenN = min(nPlantation, nNaturalRegen),
                plantationEffect = meanPlantationDbh - meanNaturalRegenDbh,
                plantationEffectPct = 100 * plantationEffect / meanDbh,
                .groups = "drop") %>%
      filter(n > 0)
    
    fitStatistics$validation$mab = sum(dbhByHeightClass$n * abs(dbhByHeightClass$meanBiasPerTree)) / sum(dbhByHeightClass$n)
    fitStatistics$validation$mapb = sum(dbhByHeightClass$n * abs(dbhByHeightClass$meanBiasPerTreePct)) / sum(dbhByHeightClass$n)
    fitStatistics$validation$meanAbsolutePlantationEffect = sum(dbhByHeightClass$minPlantationNaturalRegenN * abs(dbhByHeightClass$plantationEffect), na.rm = TRUE) / sum(dbhByHeightClass$minPlantationNaturalRegenN * (is.na(dbhByHeightClass$plantationEffect) == FALSE), na.rm = TRUE)
    fitStatistics$validation$meanAbsolutePercentPlantationEffect = sum(dbhByHeightClass$minPlantationNaturalRegenN * abs(dbhByHeightClass$plantationEffectPct), na.rm = TRUE) / sum(dbhByHeightClass$minPlantationNaturalRegenN * (is.na(dbhByHeightClass$plantationEffect) == FALSE), na.rm = TRUE)
  } else if (responseVariable == "height") {
    trainingHeightDiameterRatio = trainingPrediction / (0.01 * trainingData$dbh)
    fitStatistics$training = get_goodness_of_fit_statistics(trainingData, trainingPrediction, trainingData$height, trainingData$dbhWeight, effectiveDegreesOfFreedom, trainingHeightDiameterRatio, 1.37, trainingData$heightMax, tDegreesOfFreedom = tDegreesOfFreedom)
    
    validationHeightDiameterRatio = validationPrediction / (0.01 * validationData$dbh)
    fitStatistics$validation = get_goodness_of_fit_statistics(validationData, validationPrediction, validationData$height, validationData$dbhWeight, effectiveDegreesOfFreedom, validationHeightDiameterRatio, 1.37, validationData$heightMax, tDegreesOfFreedom = tDegreesOfFreedom)
    
    validationResiduals = validationPrediction - validationData$height
    heightByDbhClass = validationData %>% mutate(residuals = validationResiduals) %>%
      group_by(dbhClass) %>%
      summarize(n = sum(treeCount),
                nNaturalRegen = sum(treeCount * (isPlantation == FALSE)),
                nPlantation = sum(treeCount * isPlantation),
                meanBiasPerTree = sum(treeCount * residuals) / n,
                meanBiasPerTreePct = 100 * sum(treeCount * residuals / height) / n,
                meanHeight = sum(treeCount * height) / n,
                meanNaturalRegenHeight = if_else(nNaturalRegen > 0, sum(treeCount * height * (isPlantation == FALSE)) / nNaturalRegen, NA_real_),
                meanPlantationHeight = if_else(nPlantation > 0, sum(treeCount * height * isPlantation) / nPlantation, NA_real_),
                minPlantationNaturalRegenN = min(nPlantation, nNaturalRegen),
                plantationEffect = meanPlantationHeight - meanNaturalRegenHeight,
                plantationEffectPct = 100 * plantationEffect / meanHeight,
                .groups = "drop") %>%
      filter(n > 0)
    
    fitStatistics$validation$mab = sum(heightByDbhClass$n * abs(heightByDbhClass$meanBiasPerTree)) / sum(heightByDbhClass$n)
    fitStatistics$validation$mapb = sum(heightByDbhClass$n * abs(heightByDbhClass$meanBiasPerTreePct)) / sum(heightByDbhClass$n)
    fitStatistics$validation$meanAbsolutePlantationEffect = sum(heightByDbhClass$minPlantationNaturalRegenN * abs(heightByDbhClass$plantationEffect), na.rm = TRUE) / sum(heightByDbhClass$minPlantationNaturalRegenN * (is.na(heightByDbhClass$plantationEffect) == FALSE), na.rm = TRUE)
    fitStatistics$validation$meanAbsolutePercentPlantationEffect = sum(heightByDbhClass$minPlantationNaturalRegenN * abs(heightByDbhClass$plantationEffectPct), na.rm = TRUE) / sum(heightByDbhClass$minPlantationNaturalRegenN * (is.na(heightByDbhClass$plantationEffect) == FALSE), na.rm = TRUE)
  } else {
    stop(paste0("Unhandled response variable ", responseVariable, "."))
  }

  return(fitStatistics)
}

get_goodness_of_fit_statistics = function(data, predicted, measured, weights, effectiveDegreesOfFreedom, heightDiameterRatio, minPrediction, maxPrediction, tDegreesOfFreedom)
{
  nObservations = sum(data$treeCount)
  residualDegreesOfFreedom = nObservations - effectiveDegreesOfFreedom
  residuals = predicted - measured
  
  # logLik.lm() (https://github.com/wch/r-source/blob/trunk/src/library/stats/R/logLik.R)
  #  log likelihood = 1/2 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w*res^2))))
  # logLik.nls() (https://github.com/wch/r-source/blob/trunk/src/library/stats/R/nls.R)
  #  log likelihood = -N/2 * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(regression$m$resid()^2)))
  #                 = 1/2 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(regression$m$resid()^2)))) if no weights are zero
  standardDeviation = sqrt(1/residualDegreesOfFreedom * sum(weights * residuals^2)) / sqrt(weights)
  logLikelihoodGaussian = sum(data$treeCount * dnorm(residuals, sd = standardDeviation, log = TRUE))
  logLikelihoodT = sum(data$treeCount * dt(residuals / standardDeviation, df = tDegreesOfFreedom, log = TRUE) - log(standardDeviation))
  
  naturalRegenIndices = which(data$isPlantation == FALSE)
  measuredNaturalRegen = measured[naturalRegenIndices]
  predictedNaturalRegen = predicted[naturalRegenIndices]
  residualsNaturalRegen = predictedNaturalRegen - measuredNaturalRegen
  naturalRegenTreeCount = data$treeCount[naturalRegenIndices]
  naturalRegenTreeCountTotal = sum(naturalRegenTreeCount)
  
  plantationIndices = which(data$isPlantation)
  measuredPlantation = measured[plantationIndices]
  predictedPlantation = predicted[plantationIndices]
  residualsPlantation = predictedPlantation - measuredPlantation
  plantationTreeCount = data$treeCount[plantationIndices]
  plantationTreeCountTotal = sum(plantationTreeCount)
  
  fitStatistics = tibble(n = nObservations,
                         aic = -2*logLikelihoodGaussian + 2 * effectiveDegreesOfFreedom, # calculate AIC and BIC manually
                         bic = -2*logLikelihoodGaussian + effectiveDegreesOfFreedom * log(nObservations),
                         aict = -2*logLikelihoodT + 2 * effectiveDegreesOfFreedom,
                         bict = -2*logLikelihoodT + effectiveDegreesOfFreedom * log(nObservations),
                         bias = sum(data$treeCount * residuals) / nObservations, # overall bias
                         mab = NA_real_,
                         mapb = NA_real_,
                         mae = sum(data$treeCount * abs(residuals)) / nObservations, # mean absolute error
                         mape = 100 * sum(data$treeCount * abs(residuals / measured)) / nObservations, # mean absolute percent error
                         nse = 1 - sum(data$treeCount * residuals^2) / sum(data$treeCount * (measured - sum(data$treeCount * measured) / nObservations)^2), # Nash-Sutcliffe model efficiency
                         rmse = sqrt(sum(data$treeCount * residuals^2) / nObservations), # root mean squared error
                         rmspe = 100 * sqrt(sum(data$treeCount * (residuals / measured)^2) / nObservations), # root mean squared percent error
                         nPredictionOutOfRange = sum(data$treeCount * (is.na(predicted) | (predicted < minPrediction) | (predicted > maxPrediction))),
                         nHtDiaRatioImplausible = sum(data$treeCount * ((heightDiameterRatio < data$heightDiameterRatioMin) | (heightDiameterRatio > data$heightDiameterRatioMax)), na.rm = TRUE),
                         biasNaturalRegen = sum(naturalRegenTreeCount * residualsNaturalRegen) / naturalRegenTreeCountTotal,
                         maeNaturalRegen = sum(naturalRegenTreeCount * abs(residualsNaturalRegen)) / naturalRegenTreeCountTotal,
                         mapeNaturalRegen = 100 * sum(naturalRegenTreeCount * abs(residualsNaturalRegen / measuredNaturalRegen)) / naturalRegenTreeCountTotal,
                         nseNaturalRegen = 1 - sum(naturalRegenTreeCount * residualsNaturalRegen^2) / sum(naturalRegenTreeCount * (measuredNaturalRegen - sum(naturalRegenTreeCount * measuredNaturalRegen) / naturalRegenTreeCountTotal)^2),
                         paeNaturalRegen = 100 * sum(naturalRegenTreeCount * abs(residualsNaturalRegen / measuredNaturalRegen)) / naturalRegenTreeCountTotal,
                         rmseNaturalRegen = sqrt(sum(naturalRegenTreeCount * residualsNaturalRegen^2) / naturalRegenTreeCountTotal),
                         rmspeNaturalRegen = 100 * sqrt(sum(naturalRegenTreeCount * (residualsNaturalRegen / measuredNaturalRegen)^2) / naturalRegenTreeCountTotal),
                         biasPlantation = sum(plantationTreeCount * residualsPlantation) / plantationTreeCountTotal,
                         maePlantation = sum(plantationTreeCount * abs(residualsPlantation)) / plantationTreeCountTotal,
                         mapePlantation = 100 * sum(plantationTreeCount * abs(residualsPlantation / measuredPlantation)) / plantationTreeCountTotal,
                         nsePlantation = 1 - sum(plantationTreeCount * residualsPlantation^2) / sum(plantationTreeCount * (measuredPlantation - sum(plantationTreeCount * measuredPlantation) / plantationTreeCountTotal)^2),
                         paePlantation = 100 * sum(plantationTreeCount * abs(residualsPlantation / measuredPlantation)) / plantationTreeCountTotal,
                         rmsePlantation = sqrt(sum(plantationTreeCount * residualsPlantation^2) / plantationTreeCountTotal),
                         rmspePlantation = 100 * sqrt(sum(plantationTreeCount * (residualsPlantation / measuredPlantation)^2) / plantationTreeCountTotal),
                         meanAbsolutePlantationEffect = NA_real_, 
                         meanAbsolutePercentPlantationEffect = NA_real_)
  return(fitStatistics)
}

impute_height = function(height, dbh)
{
  return(if_else(is.na(height), dbh, NA_real_))
}

predict_bootstrap_dbh = function(height, predictedHeight)
{
  return(if_else(is.na(height), predictedHeight, height))
}

to_fixed_coeffficients = function(crossValidationWithModels)
{
  return(bind_rows(lapply(crossValidationWithModels$fit, function(fit)
  {
    if (class(fit)[1] == "nlme")
    {
      return(fit$coefficients$fixed)
    } else {
      return(coef(fit))
    }
  })))
}

to_parameter_confidence_intervals = function(crossValidationWithModels)
{
  return(crossValidationWithModels %>% mutate(confidenceInterval = lapply(fit, function(fit) 
          { 
            if (class(fit)[1] == "nlme")
            {
              return(as_tibble(intervals(fit, level = 0.99, which = "fixed")$fixed, rownames = NA) %>% tibble::rownames_to_column(var = "parameter"))
            } else {
              return(as_tibble(confint(fit, level = 0.99), rownames = NA) %>% tibble::rownames_to_column(var = "parameter")) 
            }
          })) %>%
          select(-fit, -splits) %>% unnest(confidenceInterval) %>%
          arrange(parameter, repetition, fold))
}

stands2023gis = left_join(st_drop_geometry(st_read("GIS/McDunn inventory 2022.gdb", layer = "McDunn_Inv_Stands_Update_Net_2023", as_tibble = TRUE, quiet = TRUE)) %>% select(-starts_with("Shape")),
                          st_drop_geometry(st_read("GIS/McDunn_Data_2025.gdb", layer = "McDunn_Stands_2024_External", as_tibble = TRUE, quiet = TRUE)) %>% select(-ACRES, -INV_YR, -ORIGIN_YR, -starts_with("Shape")),
                          by = join_by(STANDID)) %>%
  filter(is.na(ALLOCATION) | (((ALLOCATION %in% c("Arboretum", "Lake", "Meadow", "Office", "Powerline", "Rental", "Reservoir", "Rock Pit", "Towers")) == FALSE) & (str_starts(ALLOCATION, "SLTP") == FALSE))) %>%
  select(STANDID, ACRES, OriginYr, HARVEST_YR, HARVEST_TY, INV_YR, BHAGE_2024, RESEARCH_A) %>%
  rename(stand = STANDID, area = ACRES, harvestYear = HARVEST_YR, harvestType = HARVEST_TY, allocation = RESEARCH_A, originYear = OriginYr) %>%
  mutate(area = area / 2.47105381, # convert to metric: acres -> ha
         harvestType = forcats::fct_recode(as.factor(harvestType), none = " ", clearcut = "Clearcut", other = "Other", thin = "Thin"),
         inventoryAge = if_else(is.na(BHAGE_2024), INV_YR - originYear, BHAGE_2024 + INV_YR - 2024), # 
         isPlantation = ((harvestType == "clearcut") | (harvestType == "thin") | (originYear > 1950)) & (is.na(allocation) | str_starts(allocation, "Uneven") == FALSE)) # TODO: refine transition to effective replanting if possible (introduced in 1941 Oregon Forest Conservation Act but not required until 1971 passage of Oregon Forest Practices Act)
stands2022fvs = read_xlsx("inventory/FVS.xlsx", sheet = "FVS_StandInit") %>% 
  select(STAND_ID, NUM_PLOTS, SAM_WT, INV_YEAR, BASAL_AREA_FACTOR, INV_PLOT_SIZE, BRK_DBH, SITE_SPECIES, SITE_INDEX) %>%
  rename(stand = STAND_ID, plotsInStand = NUM_PLOTS, areaNet = SAM_WT, inventoryYear = INV_YEAR, baf = BASAL_AREA_FACTOR, fixedRadiusExpansionFactor = INV_PLOT_SIZE, maxFixedDbh = BRK_DBH, siteSpecies = SITE_SPECIES, siteIndex = SITE_INDEX) %>%
  mutate(areaNet = areaNet / 2.47105381, # convert to metric: acres -> ha, ft²/ac to m²/ha, in -> cm, feet -> m, TPA -> TPH
         baf = 2.47105381 * 0.3048^2 * baf,
         maxFixedDbh = 2.54 * as.numeric(maxFixedDbh),
         siteIndex = 0.3048 * siteIndex,
         fixedRadiusExpansionFactor = 2.47105381 * as.numeric(fixedRadiusExpansionFactor),
         siteSpecies = forcats::fct_recode(as.factor(siteSpecies), ABGR = "GF", ACMA = "BM", FRLA = "FL", PSME = "DF"))
stands2022 = left_join(stands2022fvs, stands2023gis, by = join_by(stand)) %>%
  mutate(isPlantation = replace_na(isPlantation, FALSE)) # from inspection and aerial imagery review
plots2022 = st_drop_geometry(st_read("GIS/McDonald-Dunn.gpkg", layer = "inventory plots 2019 2020", as_tibble = TRUE, quiet = TRUE))

#plots2022 = st_read("GIS/McDonald-Dunn.gpkg", layer = "inventory plots 2019 2020", quiet = TRUE)
#plots2022 %<>% group_by(STAND) %>% arrange(PlotID) %>% mutate(fvsPlotID = row_number())
#st_write(plots2022, dsn = "GIS/McDonald-Dunn.gpkg", layer = "inventory plots 2019 2020", append = FALSE)

trees2022 = read_xlsx("inventory/FVS.xlsx", sheet = "FVS_TreeInit") %>% 
  select(-New_ID, -STANDPLOT_ID, -TAG_ID, -DIAMETER_HT, -DG, -HTG, -HT_TO_LIVE_CROWN, -DAMAGE2, -SEVERITY2, -DAMAGE3, -SEVERITY3, -DEFECT_CUBIC, -DEFECT_BOARD, -TREEVALUE, -PRESCRIPTION, -PV_CODE, -PV_REF_CODE, -TOPOCODE, -SLOPE, -ASPECT, -SITEPREP) %>%
  rename(stand = STAND_ID, plot = PLOT_ID, tree = TREE_ID, species = SPECIES, treeCount = TREE_COUNT, dbh = DIAMETER, height = HT, heightTopKill = HTTOPK, crownRatio = CRRATIO, age = AGE, history = HISTORY, damage = DAMAGE1, severity = SEVERITY1) %>%
  mutate(tree = as.numeric(tree),
         treeCount = pmax(treeCount, 1), # fix records with tree count of zero (n = 1)
         dbh = 2.54 * dbh, # convert to metric: inch -> cm, feet -> m
         height = 0.3048 * as.numeric(height),
         heightTopKill = 0.3048 * heightTopKill,
         species = forcats::fct_recode(as.factor(species), ABGR = "GF", ABGR = "Gf", ACMA = "BM", ACMA = "Bm", ALRU = "RA", ARME = "MA", CADE = "IC", CHLA = "PC", conifer = "OT", CONU = "DG", FRLA = "FL", FRLA = "OA", PIPO = "PP", POTR = "CW", PSME = "DF", PSME = "Df", PSME = "df", Prunus = "CH", QUGA = "WO", Salix = "WI", TABR = "PY", THPL = "RC", TSHE = "WH"),
         speciesModel = factor(if_else(species %in% c("PSME", "ABGR", "ACMA"), species, "other"), levels = c("PSME", "ABGR", "ACMA", "other")),
         history = forcats::fct_recode(as.factor(history), default = "1", ingrowth = "2", planted = "3", stumpSprout = "4", dead = "6", snag = "8"),
         damage = forcats::fct_recode(as.factor(replace_na(damage, 0)), none = "0", defect = "27", brokenTop = "96", deadTop = "97"),
         isConifer = species %in% c("ABGR", "CADE", "CHLA", "conifer", "PIPO", "PSME", "TABR", "THPL", "TSHE"),
         isLive = (history != "dead") & (history != "snag"), 
         isLiveUnbroken = isLive & (damage != "brokenTop"),
         uniquePlotID = factor(paste(stand, plot)),
         breastHeight = 1.37) %>% # m, used for offset in lm() height regressions
  filter(height >= 1.37) %>% # about 1% of records are for trees shorter than breast height with a 0.25 cm or 2.5 cm DBH dubbed in rather than having a NA DBH
  relocate(stand, plot, tree, species, dbh, height, treeCount, isLive, isLiveUnbroken)
trees2022 = left_join(left_join(trees2022, stands2022 %>% rename(standArea = area, standAreaNet = areaNet), by = join_by(stand)),
                      plots2022 %>% select(-PlotID, -Contractor, -Year, -cosAspect, -sinAspect),
                      by = join_by(stand == STAND, plot == fvsPlotID)) %>%
  mutate(basalAreaPerHectare = treeCount * if_else(is.na(dbh), if_else(height >= 1.37, fixedRadiusExpansionFactor * pi/4 * (0.01 * height / 150)^2, 0), # fallback: assume a height-diameter ratio of 150 for small stems (n = 1) missing DBH measurements
                                                               if_else(dbh >= maxFixedDbh, if_else(is.na(baf), fixedRadiusExpansionFactor * pi/4 * (0.01 * dbh)^2, baf), # fallback: use fixed radius if no stand BAF (n = 1)
                                                                                           fixedRadiusExpansionFactor * pi/4 * (0.01 * dbh)^2)), # m²/ha
         treesPerHectare = treeCount * if_else(is.na(dbh), fixedRadiusExpansionFactor,
                                               if_else(dbh >= maxFixedDbh, if_else(is.na(baf), fixedRadiusExpansionFactor, baf / (pi/4 * (0.01 * dbh)^2)), # fallback: use fixed radius if stand cruised without BAF (n = 1 tree) due to nearly all trees being fixed radius
                                                       fixedRadiusExpansionFactor))) %>%
  group_by(stand) %>%
  arrange(desc(isLiveUnbroken), desc(dbh), .by_group = TRUE) %>% # put largest diameter live trees first in each stand for calculating BAL (numbers sort before NA)
  mutate(#plotsInStand = length(unique(plot)), # if necessary the number of plots can be inferred from measurements but this will miss plots without any trees, potentially leading to density overestimation (unrelated note: nested fixed radius and BAF plots share same plot ID)
         basalAreaLarger = (cumsum(isLive * basalAreaPerHectare) - basalAreaPerHectare[1]) / plotsInStand, # m²/ha
         standBasalAreaPerHectare = sum(isLive * basalAreaPerHectare) / plotsInStand, # m²/ha
         standTreesPerHectare = sum(isLive * treesPerHectare) / plotsInStand, # stand's total trees per hectare
         standQmd = sqrt(standBasalAreaPerHectare / (pi/4 * 0.01^2 * standTreesPerHectare)), # quadratic mean diameter, cm
         relativeDiameter = dbh / standQmd,
         isMapped = any(is.na(elevation)) == FALSE) %>%
  filter(isMapped) %>% # exclude trees on plots not matched to GIS data
  # top height by tallest trees in stand, regardless of plot
  group_by(stand, plot, isLiveUnbroken) %>%
  arrange(desc(height), .by_group = TRUE) %>% 
  mutate(topHeightTph = if_else(isLiveUnbroken, pmin(cumsum(if_else(is.na(height), 0, treesPerHectare)), 100), NA_real_),
         topHeightWeight = if_else(isLiveUnbroken, pmax((topHeightTph - lag(topHeightTph, default = 0)) / treesPerHectare, 0), NA_real_),
         topHeight = if_else(isLiveUnbroken, sum(topHeightWeight * height, na.rm = TRUE) / sum(topHeightWeight, na.rm = TRUE), NA_real_),
         topHeightMask = if_else(isLiveUnbroken & (is.na(topHeightWeight) == FALSE), row_number() == n(), NA_real_)) %>%
  group_by(stand) %>%
  arrange(desc(isLiveUnbroken), desc(height), .by_group = TRUE) %>%
  mutate(topHeight = sum(topHeightMask * topHeightTph * topHeight, na.rm = TRUE) / sum(topHeightMask * topHeightTph, na.rm = TRUE),
         predictedHeight = impute_height(height, dbh),
         relativeHeight = if_else(is.na(height), predictedHeight, height) / topHeight,
         bootstrapDbh = predict_bootstrap_dbh(height, dbh), # TODO cm, predicted from height where available
         predictedBasalAreaPerHectare = treesPerHectare * pi/4 * (0.01 * if_else(is.na(bootstrapDbh), dbh, bootstrapDbh))^2, # m²/ha, based on diameter predicted from height whenever available (72.5% of trees)
         bootstrapStandBasalAreaPerHectare = sum(isLive * predictedBasalAreaPerHectare) / plotsInStand, # m²/ha
         basalAreaTaller = (cumsum(isLive * predictedBasalAreaPerHectare) - predictedBasalAreaPerHectare[1]) / plotsInStand,
         treesPerHectareTaller = cumsum(isLiveUnbroken * treesPerHectare) / plotsInStand) %>%
  ungroup()

heightClassBreaks = trees2022 %>% filter(isLiveUnbroken, is.na(height) == FALSE) %>%
  group_by(speciesModel) %>%
  group_modify(~{
    quantileBreaks = seq(0, 1, length.out = min(50, sum(.$treeCount) / (5 * 10))) # constrain maximum number of classes based on data availability: setting the max to n / (meanClassN*k) classes averages meanClassN samples per class in validation folds => primarily affects low n models
    return(tibble(heightBreaks = unique(ceiling(c(0, quantile(.$height, probs = quantileBreaks, na.rm = TRUE))))))
  }) %>%
  unstack(heightBreaks ~ speciesModel) # list of height class breaks in meters, named by species model group
dbhClassBreaks = trees2022 %>% filter(isLiveUnbroken, is.na(dbh) == FALSE) %>%
  group_by(speciesModel) %>%
  group_modify(~{
    quantileBreaks = seq(0, 1, length.out = min(50, sum(.$treeCount) / (5 * 10) - 3))
    return(tibble(dbhBreaks = unique(c(2.5 * c(0, 1.5, 2.5, 3.5), 2.5 * ceiling(quantile(.$dbh, probs = quantileBreaks, na.rm = TRUE) / 2.5) + 0.5 * 2.5))))
  }) %>%
  unstack(dbhBreaks ~ speciesModel) # list of DBH class breaks in cm, named by species model group

trees2022 %<>% group_by(speciesModel) %>%
  mutate(heightClass = cut(height, breaks = heightClassBreaks[[cur_group()$speciesModel]], labels = 0.5 * (head(heightClassBreaks[[cur_group()$speciesModel]], -1) + tail(heightClassBreaks[[cur_group()$speciesModel]], -1))),
         dbhClass = cut(dbh, breaks = dbhClassBreaks[[cur_group()$speciesModel]], labels = 0.5 * (head(dbhClassBreaks[[cur_group()$speciesModel]], -1) + tail(dbhClassBreaks[[cur_group()$speciesModel]], -1)))) %>%
  ungroup()
#trees2022 %>% filter(isLiveUnbroken) %>% group_by(is.na(dbh), is.na(height)) %>% summarize(records = n(), trees = sum(treeCount))

# use only unbroken and fully measured trees for model development
# Species limits and height-diameter ratio bounds from fits in data investigatory block below.
# Weights for heteroskedasticity are found from gsl_nls() regressions in the first setup block below.
psme2022 = trees2022 %>% filter(species == "PSME", is.na(dbh) == FALSE, is.na(height) == FALSE, height >= 1.37, isLiveUnbroken) %>%
  #mutate(dbhWeight = treeCount / (0.699 * dbh^0.863), heightWeight = treeCount / (0.221 * height^(1.983 - 0.103 * isPlantation))) # all data
  mutate(dbhWeight = treeCount / (0.632 * dbh^0.886), heightWeight = treeCount / (0.227 * height^(1.975 - 0.092 * isPlantation)), # mapped plots
         heightMax = 75, dbhMax = 250, heightDiameterRatioMax = 2 + 1000 * dbh^-0.52, heightDiameterRatioMin = 22)
abgr2022 = trees2022 %>% filter(species == "ABGR", is.na(height) == FALSE, isLiveUnbroken) %>%
  #mutate(dbhWeight = treeCount / (0.645 * dbh^0.930), heightWeight = treeCount / (1.913 * height^1.159)) # all data
  mutate(dbhWeight = treeCount / (0.590 * dbh^0.961), heightWeight = treeCount / (1.811 * height^1.158), # mapped plots
         heightMax = 75, dbhMax = 150, heightDiameterRatioMax = 5 + 1000 * dbh^-0.59, heightDiameterRatioMin = 27)
acma2022 = trees2022 %>% filter(species == "ACMA", is.na(dbh) == FALSE, is.na(height) == FALSE, isLiveUnbroken) %>%
  #mutate(dbhWeight = treeCount / ((1.783 - 0.201 * isPlantation) * dbh^0.754), heightWeight = treeCount / (4.291 * height^1.269)) # all data
  mutate(dbhWeight = treeCount / ((1.783 - 0.101 * isPlantation) * dbh^0.753), heightWeight = treeCount / (4.437 * height^1.267), # mapped plots
         heightMax = 52, dbhMax = 160, heightDiameterRatioMax = 5 + 1000 * dbh^-0.61, heightDiameterRatioMin = 11 * 0.01 * dbh)
other2022 = trees2022 %>% filter(speciesModel == "other", is.na(dbh) == FALSE, is.na(height) == FALSE, isLiveUnbroken) %>%
  # mutate(dbhWeight = treeCount / ((0.940 + 0.213 * isPlantation) * dbh^0.898), heightWeight = treeCount / (1.989 * height^1.552)) # all data
  mutate(dbhWeight = treeCount / ((0.915 + 0.229 * isPlantation) * dbh^0.906), heightWeight = treeCount / (2.047 * height^1.552),
         heightMax = 50, dbhMax = 125, heightDiameterRatioMax = 10 + 1000 * dbh^-0.66, heightDiameterRatioMin = 12)


## weighting for heteroskedasticity
if (htDiaOptions$includeSetup)
{
  # Douglas-fir weighting for heteroskedasticity
  # Update geom_smooth() calls below when changing smooths' k.
  psmeHeightFromDbh = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 8), data = psme2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8)
  psmeDbhFromHeight = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 16), data = psme2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8)
  #k.check(psmeHeightFromDbh)
  
  psmeHeightHeteroskedasticity = psme2022 %>% select(treeCount, dbh, isPlantation) %>% mutate(residual = psmeHeightFromDbh$residuals, squaredResidual = residual^2)
  psmeHeightVariance = gsl_nls(squaredResidual ~ a1 * dbh^b1, psmeHeightHeteroskedasticity, start = list(a1 = 0.5, b1 = 1), weight = treeCount / dbh^0.89)
  #psmeHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^b1, psmeHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1), weight = treeCount / dbh^0.88) # a1p not significant
  #psmeHeightVariance = gsl_nls(squaredResidual ~ a1 * dbh^(b1 + b1p * isPlantation), psmeHeightHeteroskedasticity, start = list(a1 = 0.5, b1 = 1, b1p = 0), weight = treeCount / dbh^0.90) # b1p not significant
  #psmeHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^(b1 + b1p * isPlantation), psmeHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / dbh^0.90) # b1p not significant
  #confint(psmeHeightVariance, level = 0.99)
  psmeDbhHeteroskedasticity = psme2022 %>% select(treeCount, height, isPlantation) %>% mutate(residual = psmeDbhFromHeight$residuals, squaredResidual = residual^2)
  psmeDbhVariance = gsl_nls(squaredResidual ~ a1 * height^(b1 + b1p * isPlantation), psmeDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1, b1p = 0), weight = treeCount / height^1.98)
  #psmeDbhVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * height^(b1 + b1p * isPlantation), psmeDbhHeteroskedasticity, start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / height^1.98) # a1p not signifcant
  #confint(psmeDbhVariance, level = 0.99)
  
  psmeHeightVarianceInterval = tibble(dbh = seq(0, 242)) %>% mutate(predictedHeight = as_tibble(predict(psmeHeightVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedHeight)
  psmeDbhVarianceInterval = crossing(isPlantation = c(TRUE, FALSE), height = seq(0, 73)) %>% mutate(predictedDbh = as_tibble(predict(psmeDbhVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedDbh)
  
  # grand fir weighting for heteroskedasticity
  abgrHeightFromDbh = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = abgr2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8)
  abgrDbhFromHeight = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 18), data = abgr2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8)
  #k.check(abgrDbhFromHeight)
  
  abgrHeightHeteroskedasticity = abgr2022 %>% select(treeCount, dbh, isPlantation) %>% mutate(residual = abgrHeightFromDbh$residuals, squaredResidual = residual^2)
  abgrHeightVariance = gsl_nls(squaredResidual ~ a1 * dbh^b1, abgrHeightHeteroskedasticity, start = list(a1 = 0.5, b1 = 1), weight = treeCount / dbh^0.96)
  #abgrHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^b1, abgrHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1), weight = treeCount / dbh^0.93) # a1p not significant
  #abgrHeightVariance = gsl_nls(squaredResidual ~ a1 * dbh^(b1 + b1p * isPlantation), abgrHeightHeteroskedasticity, start = list(a1 = 0.5, b1 = 1, b1p = 0), weight = treeCount / dbh^0.90) # b1p not significant
  #abgrHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^(b1 + b1p * isPlantation), abgrHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / dbh^0.90) # a1p, b1p not significant
  confint(abgrHeightVariance, level = 0.99)
  abgrDbhHeteroskedasticity = abgr2022 %>% select(treeCount, height, isPlantation) %>% mutate(residual = abgrDbhFromHeight$residuals, squaredResidual = residual^2)
  abgrDbhVariance = gsl_nls(squaredResidual ~ a1 * height^b1, abgrDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1), weight = treeCount / height^1.16)
  #abgrDbhVariance = gsl_nls(squaredResidual ~ a1 * height^(b1 + b1p * isPlantation), abgrDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1, b1p = 0), weight = treeCount / height^1.15)
  #confint(abgrDbhVariance, level = 0.99)
  
  abgrHeightVarianceInterval = tibble(dbh = seq(0, 132)) %>% mutate(predictedHeight = as_tibble(predict(abgrHeightVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedHeight)
  abgrDbhVarianceInterval = tibble(height = seq(0, 75)) %>% mutate(predictedDbh = as_tibble(predict(abgrDbhVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedDbh)
  
  # bigleaf maple weighting for heteroskedasticity
  acmaHeightFromDbh = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = acma2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8) # low k to avoid following limited measurements around 100 cm
  acmaDbhFromHeight = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 14), data = acma2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8)
  #k.check(acmaDbhFromHeight)
  
  acmaHeightHeteroskedasticity = acma2022 %>% select(treeCount, dbh, isPlantation) %>% mutate(residual = acmaHeightFromDbh$residuals, squaredResidual = residual^2)
  acmaHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^b1, acmaHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1), weight = treeCount / dbh^0.75)
  #acmaHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^(b1 + b1p * isPlantation), acmaHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / dbh^0.78) # b1p not significant
  #confint(acmaHeightVariance, level = 0.99)
  acmaDbhHeteroskedasticity = acma2022 %>% select(treeCount, height, isPlantation) %>% mutate(residual = acmaDbhFromHeight$residuals, squaredResidual = residual^2)
  acmaDbhVariance = gsl_nls(squaredResidual ~ a1 * height^b1, acmaDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1), weight = treeCount / height^1.27)
  #acmaDbhVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * height^b1, acmaDbhHeteroskedasticity, start = list(a1 = 1, a1p = 0, b1 = 1), weight = treeCount / height^1.01) # a1p not significant
  #acmaDbhVariance = gsl_nls(squaredResidual ~ a1 * height^(b1 + b1p * isPlantation), acmaDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1, b1p = 0), weight = treeCount / height^0.94) # b1p not significant
  #acmaDbhVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * height^(b1 + b1p * isPlantation), acmaDbhHeteroskedasticity, start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / height^0.94) # a1p, b1p not significant
  #confint(acmaDbhVariance, level = 0.99)
  
  acmaHeightVarianceInterval = crossing(isPlantation = c(TRUE, FALSE), dbh = seq(0, 147)) %>% mutate(predictedHeight = as_tibble(predict(acmaHeightVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedHeight)
  acmaDbhVarianceInterval = tibble(height = seq(0, 50)) %>% mutate(predictedDbh = as_tibble(predict(acmaDbhVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedDbh)
  
  # other species weighting for heteroskedasticity
  otherHeightFromDbh = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 10), data = other2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8) # low k to avoid following limited measurements around 100 cm
  otherDbhFromHeight = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 10), data = other2022, method = "REML", select = TRUE, weights = pmin(treeCount/dbh, treeCount), nthreads = 8)
  #k.check(otherDbhFromHeight)
  
  otherHeightHeteroskedasticity = other2022 %>% select(treeCount, dbh, isPlantation) %>% mutate(residual = otherHeightFromDbh$residuals, squaredResidual = residual^2)
  otherHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^b1, otherHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1), weight = treeCount / dbh^0.91) # a1p significant
  #otherHeightVariance = gsl_nls(squaredResidual ~ a1 * dbh^(b1 + b1p * isPlantation), otherHeightHeteroskedasticity, start = list(a1 = 0.5, b1 = 1, b1p = 0), weight = treeCount / dbh^(0.85 + 0.065 * isPlantation)) # b1p significant but higher error than retaining a1p
  #otherHeightVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * dbh^(b1 + b1p * isPlantation), otherHeightHeteroskedasticity, start = list(a1 = 0.5, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / dbh^0.78) # a1p, b1p not significant
  #confint(otherHeightVariance, level = 0.99)
  otherDbhHeteroskedasticity = other2022 %>% select(treeCount, height, isPlantation) %>% mutate(residual = otherDbhFromHeight$residuals, squaredResidual = residual^2)
  otherDbhVariance = gsl_nls(squaredResidual ~ a1 * height^b1, otherDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1), weight = treeCount / height^1.55)
  #otherDbhVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * height^b1, otherDbhHeteroskedasticity, start = list(a1 = 1, a1p = 0, b1 = 1), weight = treeCount / height^1.59) # a1p not significant
  #otherDbhVariance = gsl_nls(squaredResidual ~ a1 * height^(b1 + b1p * isPlantation), otherDbhHeteroskedasticity, start = list(a1 = 1, b1 = 1, b1p = 0), weight = treeCount / height^1.53) # b1p not significant
  #otherDbhVariance = gsl_nls(squaredResidual ~ (a1 + a1p * isPlantation) * height^(b1 + b1p * isPlantation), otherDbhHeteroskedasticity, start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0), weight = treeCount / height^1.49) # a1p, b1p not significant
  #confint(otherDbhVariance, level = 0.99)
  
  otherHeightVarianceInterval = crossing(isPlantation = c(TRUE, FALSE), dbh = seq(0, 117)) %>% mutate(predictedHeight = as_tibble(predict(otherHeightVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedHeight)
  otherDbhVarianceInterval = tibble(height = seq(0, 48)) %>% mutate(predictedDbh = as_tibble(predict(otherDbhVariance, ., interval = "confidence", level = 0.99))) %>% unnest_wider(predictedDbh)
  
  ggplot() +
    geom_point(aes(x = dbh, y = squaredResidual), psmeHeightHeteroskedasticity, alpha = 0.05, shape = 16, size = 0.5) +
    geom_smooth(aes(x = dbh, y = squaredResidual, color = "GAM"), psmeHeightHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 8), level = 0.99, alpha = 0.12, linewidth = 0.3) +
    geom_ribbon(aes(x = dbh, ymin = lwr, ymax = upr), psmeHeightVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = dbh, y = fit, color = "power"), psmeHeightVarianceInterval) +
    coord_transform(y = scales::pseudo_log_trans(), xlim = c(0, 245), ylim = c(0, 1500)) +
    labs(x = "DBH, cm", y = "squared residual, m²", color = NULL, title = paste(plotLetters[1], "Douglas-fir")) +
    scale_x_continuous(breaks = seq(0, 250, by = 100)) +
    scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 2000), minor_breaks = c(2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900)) +
  ggplot() +
    geom_point(aes(x = dbh, y = squaredResidual), abgrHeightHeteroskedasticity, alpha = 0.1, shape = 16, size = 0.5) +
    geom_smooth(aes(x = dbh, y = squaredResidual, color = "GAM"), abgrHeightHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 5), level = 0.99, alpha = 0.12, linewidth = 0.3) +
    geom_ribbon(aes(x = dbh, ymin = lwr, ymax = upr), abgrHeightVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = dbh, y = fit, color = "power"), abgrHeightVarianceInterval) +
    coord_transform(y = scales::pseudo_log_trans(), xlim = c(0, 245), ylim = c(0, 1500)) +
    labs(x = "DBH, cm", y = NULL, color = NULL, title = paste(plotLetters[2], "grand fir")) +
    scale_x_continuous(breaks = seq(0, 250, by = 100)) +
    scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 2000), minor_breaks = c(2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900)) +
  ggplot() +
    geom_point(aes(x = dbh, y = squaredResidual), acmaHeightHeteroskedasticity, alpha = 0.1, shape = 16, size = 0.5) +
    geom_smooth(aes(x = dbh, y = squaredResidual, color = "GAM"), acmaHeightHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 5), level = 0.99, alpha = 0.12, linewidth = 0.3) +
    geom_ribbon(aes(x = dbh, ymin = lwr, ymax = upr, group = isPlantation), acmaHeightVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = dbh, y = fit, color = "power", group = isPlantation, linetype = isPlantation), acmaHeightVarianceInterval) +
    coord_transform(y = scales::pseudo_log_trans(), xlim = c(0, 245), ylim = c(0, 1500)) +
    labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL, title = paste(plotLetters[3], "bigleaf maple")) +
    scale_x_continuous(breaks = seq(0, 250, by = 100)) +
    scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 2000), minor_breaks = c(2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900)) +
  ggplot() +
    geom_point(aes(x = dbh, y = squaredResidual), otherHeightHeteroskedasticity, alpha = 0.1, shape = 16, size = 0.5) +
    geom_smooth(aes(x = dbh, y = squaredResidual, color = "GAM"), otherHeightHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 10), level = 0.99, alpha = 0.12, linewidth = 0.3) +
    geom_ribbon(aes(x = dbh, ymin = lwr, ymax = upr, group = isPlantation), otherHeightVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = dbh, y = fit, color = "power", group = isPlantation, linetype = isPlantation), otherHeightVarianceInterval) +
    coord_transform(y = scales::pseudo_log_trans(), xlim = c(0, 245), ylim = c(0, 1500)) +
    labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL, title = paste(plotLetters[4], "other species height")) +
    scale_x_continuous(breaks = seq(0, 250, by = 100)) +
    scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 2000), minor_breaks = c(2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900)) +
  ggplot() +
    geom_point(aes(x = squaredResidual, y = height), psmeDbhHeteroskedasticity, alpha = 0.05, shape = 16, size = 0.5) +
    geom_smooth(aes(x = squaredResidual, y = height, color = "GAM"), psmeDbhHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 18), level = 0.99, orientation = "y", alpha = 0.2, linewidth = 0.3) + # by fails to find isPlantation, even when vector explictly provided
    geom_ribbon(aes(xmin = lwr, xmax = upr, y = height, group = isPlantation), psmeDbhVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = fit, y = height, color = "power", group = isPlantation, linetype = isPlantation), psmeDbhVarianceInterval) +
    coord_transform(x = scales::pseudo_log_trans(), xlim = c(0, 15000), ylim = c(0, 77)) +
    guides(color = "none", fill = "none") + # color and fill legends don't collect across height and DBH (but linetype does, oddly)
    labs(x = "squared residual, cm²", y = "height, m", color = NULL, linetype = NULL, title = paste(plotLetters[5], "Douglas-fir")) +
    scale_x_continuous(breaks = c(0, 100, 10000), minor_breaks = c(1, 10, 1000), labels = scales::comma) +
  ggplot() +
    geom_point(aes(x = squaredResidual, y = height), abgrDbhHeteroskedasticity, alpha = 0.1, shape = 16, size = 0.5) +
    geom_smooth(aes(x = squaredResidual, y = height, color = "GAM"), abgrDbhHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 16), level = 0.99, orientation = "y", alpha = 0.2, linewidth = 0.3) +
    geom_ribbon(aes(xmin = lwr, xmax = upr, y = height), abgrDbhVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = fit, y = height, color = "power"), abgrDbhVarianceInterval) +
    coord_transform(x = scales::pseudo_log_trans(), xlim = c(0, 15000), ylim = c(0, 77)) +
    guides(color = "none", fill = "none") +
    labs(x = "squared residual, cm²", y = NULL, color = NULL, title = paste(plotLetters[6], "grand fir")) +
    scale_x_continuous(breaks = c(0, 100, 10000), minor_breaks = c(1, 10, 1000), labels = scales::comma) +
  ggplot() +
    geom_point(aes(x = squaredResidual, y = height), acmaDbhHeteroskedasticity, alpha = 0.1, shape = 16, size = 0.5) +
    geom_smooth(aes(x = squaredResidual, y = height, color = "GAM"), acmaDbhHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 10), level = 0.99, orientation = "y", alpha = 0.2, linewidth = 0.3) +
    geom_ribbon(aes(xmin = lwr, xmax = upr, y = height), acmaDbhVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = fit, y = height, color = "power"), acmaDbhVarianceInterval) +
    coord_transform(x = scales::pseudo_log_trans(), xlim = c(0, 15000), ylim = c(0, 77)) +
    guides(color = "none", fill = "none") +
    labs(x = "squared residual, cm²", y = NULL, color = NULL, title = paste(plotLetters[7], "bigleaf maple")) +
    scale_x_continuous(breaks = c(0, 100, 10000), minor_breaks = c(1, 10, 1000), labels = scales::comma) +
  ggplot() +
    geom_point(aes(x = squaredResidual, y = height), otherDbhHeteroskedasticity, alpha = 0.1, shape = 16, size = 0.5) +
    geom_smooth(aes(x = squaredResidual, y = height, color = "GAM"), otherDbhHeteroskedasticity, method = "gam", formula = y ~ s(x, bs = "ts", k = 14), level = 0.99, orientation = "y", alpha = 0.2, linewidth = 0.3) +
    geom_ribbon(aes(xmin = lwr, xmax = upr, y = height), otherDbhVarianceInterval, alpha = 0.08) +
    geom_line(aes(x = fit, y = height, color = "power"), otherDbhVarianceInterval) +
    coord_transform(x = scales::pseudo_log_trans(), xlim = c(0, 15000), ylim = c(0, 77)) +
    guides(color = "none", fill = "none") +
    labs(x = "squared residual, cm²", y = NULL, color = NULL, title = paste(plotLetters[8], "other species DBH")) +
    scale_x_continuous(breaks = c(0, 100, 10000), minor_breaks = c(1, 10, 1000), labels = scales::comma) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, guides = "collect") &
    scale_color_manual(breaks = c("GAM", "power"), values = c("red", "cyan2")) &
    scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("unfactored\nor natural\nregeneration", "plantation"), values = c("solid", "longdash")) &
    theme(legend.justification = c(0.5, 0.9), legend.spacing.y = unit(0.2, "line"))
  #ggsave("trees/height-diameter/figures/Figure S variance models.png", height = 12, width = 16.5, units = "cm", dpi = 300)
}

## data checks and summaries
if (htDiaOptions$includeInvestigatory)
{
  table(trees2022$species, useNA = "ifany")
  table(trees2022$treeCount, useNA = "ifany")
  trees2022 %>% group_by(species) %>% summarize(records = n(), trees = sum(isLive * treeCount), snags = sum((isLive == FALSE) * treeCount), dbhLive = sum(isLive * (is.na(dbh) == FALSE) * treeCount), unbrokenHeightLive = sum(isLiveUnbroken * treeCount), brokenHeightLive = sum(isLive * (is.na(height) == FALSE) * (damage == "brokenTop") * treeCount), dbhSnag = sum((is.na(dbh) == FALSE) * (isLive == FALSE) * treeCount), heightSnag = sum((is.na(height) == FALSE) * (isLive == FALSE) * treeCount), mapped = sum(is.na(elevation) == FALSE)) %>% 
    mutate(pctTrees = 100 * trees / sum(trees)) %>% arrange(desc(pctTrees)) %>%
    relocate(species, records, trees, pctTrees, snags) %>%
    arrange(desc(pctTrees))
  trees2022 %>% filter(basalAreaLarger > standBasalAreaPerHectare) # should be empty, ok if BAL = stand BA
  trees2022 %>% filter(basalAreaTaller >= bootstrapStandBasalAreaPerHectare) # should be empty
  colSums(is.na(trees2022 %>% filter(((height < 1.37) & is.na(dbh)) != TRUE) %>% select(-height, -predictedHeight, -bootstrapDbh, -heightTopKill, -crownRatio, -originYear, -age, -BHAGE_2024, -harvestType, -harvestYear, -baf, -fixedRadiusExpansionFactor, -topHeightTph, -topHeightWeight, -topHeightMask, -inventoryAge, -INV_YR, -allocation, -severity))) # should be zero in all retained columns but 1 DBH and 1 relative diameter is expected
  trees2022 %>% filter(is.na(inventoryAge)) %>% group_by(stand) %>% summarize(records = n()) # two stands without age estimates
  
  # inventory coverage
  stands2022 %>% mutate(hasInventory = stand %in% unique(trees2022$stand)) %>% group_by(hasInventory) %>%
    summarize(stands = n(), area = sum(area), areaNet = sum(areaNet)) %>% relocate(stands) %>%
    bind_rows(., tibble(area = sum(.$area), areaNet = sum(.$areaNet)))
  
  # availability of physiographic variables
  physio2022 = trees2022 %>% filter(is.na(elevation), str_starts(stand, "11050") == FALSE) %>% 
    mutate(hasBreastHeight = height >= 1.37, isFittable = hasBreastHeight & (is.na(dbh) == FALSE)) %>% group_by(stand, hasBreastHeight) %>% 
    summarize(isFittable = sum(isFittable), naHeight = sum(is.na(height)), naElevation = sum(is.na(elevation)), naSlope = sum(is.na(slope)), naAspect = sum(is.na(aspect)), naTerrainRoughness = sum(is.na(terrainRoughness)), naWetness = sum(is.na(topographicWetnessFD8f)), .groups = "drop") %>%
    group_by(stand) %>%
    mutate(mixedHeightStand = (any(hasBreastHeight) == FALSE) & (any(hasBreastHeight) == TRUE), hasNAheights = sum(naHeight) > 0) %>%
    group_by(hasBreastHeight) %>%
    summarize(stands = length(unique(stand)), mixedHeightStands = sum(mixedHeightStand), fittableRecords = sum(isFittable), naHeightStands = sum(hasNAheights), naHeight = sum(naHeight), naElevation = sum(naElevation), naSlope = sum(naSlope), naAspect = sum(naAspect), naTerrainRoughness = sum(naTerrainRoughness), naWetness = sum(naWetness), standIDs = list(stand))
  print(trees2022 %>% filter(stand %in% physio2022$standIDs[2][[1]]) %>% group_by(stand) %>% # fails on converting inventoryYear to size 0 if no stands listed
          reframe(inventoryYear = inventoryYear[1], fittableRecords = sum((is.na(height) == FALSE) & (height >= 1.37) & (is.na(dbh) == FALSE)), quantiles = c(0, 0.8, 0.9, 1), height = quantile(height, probs = quantiles, na.rm = TRUE)) %>% 
          pivot_wider(id_cols = c("stand", "inventoryYear", "fittableRecords"), names_from = "quantiles", names_prefix = "heightQuantileQ", values_from = "height") %>% 
          arrange(desc(fittableRecords)), n = 50)

  standSummary2022 = trees2022 %>% group_by(stand) %>% summarize(standArea = standArea[1], standAreaNet = standAreaNet[1], isPlantation = isPlantation[1], standTreesPerHectare = standTreesPerHectare[1], topHeight = topHeight[1], standQmd = standQmd[1], standBasalAreaPerHectare = standBasalAreaPerHectare[1], baf = baf[1], plots = plotsInStand[1], bootstrapStandBasalAreaPerHectare = bootstrapStandBasalAreaPerHectare[1])
  colSums(is.na(standSummary2022)) # should all be zero except for young, non-BAF stands
  standSummary2022 %>% group_by(isPlantation) %>% summarize(stands = n(), area = sum(standArea), areaNet = sum(standAreaNet))

  # consistency between FVS and GIS plot counts
  plotDefinitions2022 = left_join(trees2022 %>% group_by(stand) %>% summarize(fvsPlots = length(unique(plot))),
                                  plots2022 %>% group_by(STAND) %>% summarize(gisPlots = length(unique(fvsPlotID))),
                                  by = join_by(stand == STAND))
  plotDefinitions2022 %>% filter(fvsPlots > gisPlots) # should be empty

  # inconsistencies between FVS and GIS plot numbering are most readily detected through missing physiographic variables
  unmappedPlots2022 = left_join(trees2022 %>% filter(is.na(elevation)) %>% group_by(stand, plot) %>% slice_head(n = 1) %>% select(stand, plot) %>% group_by(stand) %>% summarize(plots = list(plot)), 
                                plots2022 %>% select(STAND, PlotID, Year) %>% group_by(STAND) %>% summarize(year = Year[1], gisPlots = list(PlotID)), 
                                by = join_by(stand == STAND)) %>%
    relocate(stand, year, plots, gisPlots) %>%
    rowwise() %>%
    mutate(inventoryPlotCount = length(plots),
           gisPlotCount = length(gisPlots),
           unmatchedInventoryPlotIDs = length(setdiff(plots, gisPlots))) %>%
    ungroup() %>%
    mutate(firstInventoryPlotID = sapply(plots, function(x) { if (is.null(x)) { return(NA) } else { return(x[[1]]) } }),
           firstGisPlotID = sapply(gisPlots, function(x) { if (is.null(x)) { return(NA) } else { return(x[[1]]) } }))
  print(unmappedPlots2022, n = 160) # should be empty
  #writexl::write_xlsx(unmappedPlots2022, "trees/height-diameter/data/unmapped plots.xlsx")
  
  print(trees2022 %>% filter(stand %in% c("030206", "030306")) %>% select(stand, plot, tree, species, dbh, height, elevation), n = 50)
  plots2022 %>% filter(STAND == "030206")
  
  table(trees2022$plot)
  
  # nominal plot sizes
  # Fixed radius plots mostly 3.58 m with some 5.07 m.
  plotRadii2022 = trees2022 %>% filter(dbh > maxFixedDbh) %>% mutate(plotRadius = dbh * sqrt(0.25 / baf)) %>%
    group_by(stand, plot) %>%
    summarize(effectivePlotSize = sum(treeCount * 2/3 *plotRadius) / sum(treeCount), baf = baf[1], .groups = "drop")
  ggplot() +
    geom_histogram(aes(x = effectivePlotSize, fill = as.factor(round(baf, 2)), group = baf), plotRadii2022) +
    labs(x = "effective plot radius, m", y = "plots", fill = "BAF, m²/ha") +
  ggplot() +
    geom_segment(aes(x = 6.2, y = 0, xend = 6.2, yend = 1), color = "grey90") +
    geom_segment(aes(x = 8.1, y = 0, xend = 8.1, yend = 1), color = "grey90") +
    stat_ecdf(aes(x = effectivePlotSize), plotRadii2022) +
    labs(x = "effective plot radius, m", y = "cumulative distribution") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout() # looks better without collected guides

  # stand structure 
  # sdi = tph * (qmd/25)^1.605 -> tph = sdi / (qmd/25)^1.605, qmd = 25 (sdi/tph)^(1/1.605)
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 50, yend = 40), color = "grey90") +
    geom_point(aes(x = standAreaNet, y = plots, color = isPlantation), standSummary2022, alpha = 0.2, shape = 16) +
    labs(x = "net stand area, ha", y = "plots", color = NULL) +
    scale_y_continuous(expand = expansion(mult = 0.02, add = 0)) +
  ggplot() +
    geom_segment(aes(x = 108.067, y = 100, xend = 5000, yend = 9.1716), color = "grey90") + # SDI 1000
    geom_text(aes(x = 130, y = 95, label = "SDI = 1000"), color = "grey70", hjust = 0, size = 3) +
    geom_point(aes(x = standTreesPerHectare, y = standQmd, color = isPlantation, size = standAreaNet), standSummary2022 %>% filter(is.na(standQmd) == FALSE), alpha = 0.15, shape = 16) + # stands too young for trees to have reached DBH have undefined QMDs
    coord_cartesian(xlim = c(70, NA), ylim = c(0, 95)) +
    labs(x = bquote("trees, ha"^-1), y = "QMD, cm", color = NULL, size = "stand\narea,\nha") +
    scale_x_continuous(breaks = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000), minor_breaks = c(30, 40, 60, 70, 80, 90, 200, 300, 400, 600, 700, 800, 900, 3000, 4000, 6000, 7000, 8000, 9000), transform = scales::pseudo_log_trans(sigma = 5)) +
    scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100), minor_breaks = c(3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90), transform = scales::pseudo_log_trans(sigma = 5), expand = expansion(mult = 0.02, add = 0)) +
    guides(color = "none") + # collection doesn't work
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 75, yend = 50), color = "grey90") +
    geom_point(aes(x = standBasalAreaPerHectare, y = topHeight, color = isPlantation, size = standAreaNet), standSummary2022, alpha = 0.2, shape = 16) +
    labs(x = bquote("basal area, m"^2*"ha"^-1), y = "top height, m", color = NULL, size = "stand\narea,\nha") +
    guides(color = "none", size = "none") + # collection doesn't work
    scale_y_continuous(breaks = c(0, 1.37, seq(10, 50, by = 10)), minor_breaks = seq(5, 45, by = 10), expand = expansion(mult = 0.02, add = 0)) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    scale_color_manual(breaks = c(FALSE, TRUE), labels = c("natural regen", "plantation"), values = c("forestgreen", "blue"))
  
  # stand age
  ggplot() +
    geom_segment(aes(x = 1681, y = 2024 - 1681, xend = 2024, yend = 0), color = "grey90") +
    geom_point(aes(x = originYear, BHAGE_2024, color = isPlantation, size = area), stands2023gis, alpha = 0.2, shape = 16) +
    coord_cartesian(xlim = c(1681, NA)) +
    labs(x = "stand origin, CE", y = "stand age, years", color = NULL, size = "stand\narea,\nha") +
    scale_color_manual(breaks = c(FALSE, TRUE), labels = c("natural regen", "plantation"), values = c("forestgreen", "blue"))

  # data distribution and species limits on height-diameter ratios
  # Keep species limits in sync with species setup (psme2022, abgr2022, acma2022, other2022). See also Figure 1 in results.R.
  # With height in m and DBH in cm,
  #   h-d ratio = height/dbh = a0 + a1 * dbh^b1 => height = (a0 + a1 * dbh^b1) * 0.01 * dbh
  #                          = a0 => height = a0 * 0.01 * dbh
  liveUnbrokenSampling2022 = trees2022 %>% filter(isLiveUnbroken, is.na(dbh) == FALSE, is.na(height) == FALSE) %>%
    mutate(heightClass = round(height), dbhClass = 2.5 * round(dbh / 2.5)) %>%
    group_by(speciesModel, heightClass, dbhClass) %>%
    summarize(treesMeasured = sum(treeCount), .groups = "drop")
  
  ggplot() +
    geom_path(aes(x = seq(0, 250), y = 22 * 0.01 * seq(0, 250)), color = "grey80") + # height-diameter ratio lower bound
    geom_path(aes(x = seq(0.1, 250), y = (2 + 1000 * seq(0.1, 250)^-0.52) * 0.01 * seq(0.1, 250)), color = "grey80") + # height-diameter ratio upper bound
    geom_tile(aes(x = dbhClass, y = heightClass, fill = treesMeasured), liveUnbrokenSampling2022 %>% filter(speciesModel == "PSME")) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 75)) +
    labs(x = "DBH, cm", y = "height, m", fill = "trees\nmeasured", title = paste(plotLetters[1], "Douglas-fir")) +
  ggplot() +
    geom_path(aes(x = seq(0, 150), y = 27 * 0.01 * seq(0, 150)), color = "grey80") + # height-diameter ratio lower bound
    geom_path(aes(x = seq(0.1, 150), y = (5 + 1000 * seq(0.1, 150)^-0.59) * 0.01 * seq(0.1, 150)), color = "grey80") +
    geom_tile(aes(x = dbhClass, y = heightClass, fill = treesMeasured), liveUnbrokenSampling2022 %>% filter(speciesModel == "ABGR")) +
    coord_cartesian(xlim = c(0, 150), ylim = c(0, 75)) +
    labs(x = "DBH, cm", y = NULL, fill = "trees\nmeasured", title = paste(plotLetters[2], "grand fir")) +
  ggplot() +
    geom_path(aes(x = seq(0, 150), y = 11 * 0.01^2 * seq(0, 150)^2), color = "grey80") + # height-diameter ratio lower bound
    geom_path(aes(x = seq(0.1, 150), y = (5 + 1000 * seq(0.1, 150)^-0.61) * 0.01 * seq(0.1, 150)), color = "grey80") +
    geom_tile(aes(x = dbhClass, y = heightClass, fill = treesMeasured), liveUnbrokenSampling2022 %>% filter(speciesModel == "ACMA")) +
    coord_cartesian(xlim = c(0, 150), ylim = c(0, 75)) +
    labs(x = "DBH, cm", y = NULL, fill = "trees\nmeasured", title = paste(plotLetters[3], "bigleaf maple")) +
  ggplot() +
    geom_path(aes(x = seq(0, 150), y = 12 * 0.01 * seq(0, 150)), color = "grey80") + # height-diameter ratio lower bound
    geom_path(aes(x = seq(0.1, 150), y = (10 + 1000 * seq(0.1, 150)^-0.66) * 0.01 * seq(0.1, 150)), color = "grey80") +
    geom_tile(aes(x = dbhClass, y = heightClass, fill = treesMeasured), liveUnbrokenSampling2022 %>% filter(speciesModel == "other")) +
    coord_cartesian(xlim = c(0, 150), ylim = c(0, 75)) +
    labs(x = "DBH, cm", y = NULL, fill = "trees\nmeasured", title = paste(plotLetters[4], "other species")) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect", widths = c(250, 150, 150, 150)) &
    scale_fill_viridis_c(trans = "log10", breaks = c(1, 3, 10, 30, 100, 200), labels = c(1, 3, 10, 30, 100, "200+"), limits = c(1, 200), na.value = "yellow")
}

## variable importance
if (htDiaOptions$includeInvestigatory)
{
  library(ranger)
  
  heightCandidateVariables = c("height", "dbh", "standTreesPerHectare", "standBasalAreaPerHectare", "basalAreaLarger", "standQmd", "relativeDiameter", "topHeight", "isPlantation", "elevation", "slope", "aspect", "terrainRoughness", "topographicWetnessFD8f", "inventoryAge")
  dbhCandidateVariables  = c("dbh", "height", "topHeight", "relativeHeight", "predictedBasalAreaPerHectare", "basalAreaTaller", "isPlantation", "elevation", "slope", "aspect", "terrainRoughness", "topographicWetnessFD8f", "inventoryAge")
  variableImportanceTrees2022 = bind_rows(psme2022, abgr2022, acma2022, other2022)
  
  speciesModelWeights = tibble(psme = sum(psme2022$treeCount), abgr = sum(abgr2022$treeCount), acma = sum(acma2022$treeCount), other = sum(other2022$treeCount)) %>%
    mutate(total = sum(c_across(everything()))) %>%
    mutate(across(-total) / total)
  
  get_dbh_local_importance = function(data2022, mtry, minNodeSize, samplingFraction)
  {
    localImportance = ranger(dbh ~ ., data = data2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)),
                             importance = "permutation", local.importance = TRUE, # https://github.com/imbs-hl/ranger/issues/552
                             mtry = mtry, splitrule = 'variance', min.node.size = minNodeSize, sample.fraction = samplingFraction,
                             num.threads = htDiaOptions$rangerThreads)
    return(as_tibble(localImportance$variable.importance.local) %>%
             summarize(across(everything(), mean)) %>%
             rowwise() %>%
             mutate(across(where(is.numeric), ~100 * .x / max(across(where(is.numeric))))) %>%
             ungroup() %>%
             mutate(speciesModel = data2022$speciesModel[1]) %>%
             relocate(speciesModel))
  }
  
  get_height_local_importance = function(data2022, mtry, minNodeSize, samplingFraction)
  {
    localImportance = ranger(height ~ ., data = data2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)),
                             importance = "permutation", local.importance = TRUE, # https://github.com/imbs-hl/ranger/issues/552
                             mtry = mtry, splitrule = 'variance', min.node.size = minNodeSize, sample.fraction = samplingFraction,
                             num.threads = htDiaOptions$rangerThreads)
    return(as_tibble(localImportance$variable.importance.local) %>%
                       summarize(across(everything(), mean)) %>%
                       rowwise() %>%
                       mutate(across(where(is.numeric), ~100 * .x / max(across(where(is.numeric))))) %>%
                       ungroup() %>%
                       mutate(speciesModel = data2022$speciesModel[1]) %>%
                       relocate(speciesModel))
  }
  
  # ranger fits for variable importance (using tunings below)
  heightImportance = ranger(height ~ ., data = variableImportanceTrees2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)),
                            importance = "impurity_corrected", # https://github.com/imbs-hl/ranger/issues/664
                            mtry = 3, splitrule = 'variance', min.node.size = 2, sample.fraction = 0.753,
                            num.threads = htDiaOptions$rangerThreads)
  heightImportance = tibble(predictor = names(heightImportance$variable.importance), relativeImportance = 100 * heightImportance$variable.importance / max(heightImportance$variable.importance)) %>% 
                       arrange(desc(relativeImportance)) %>% 
                       mutate(predictor = factor(predictor, levels = predictor))
  
  psmeHeightLocalImportance = get_height_local_importance(psme2022, mtry = 3, minNodeSize = 2, samplingFraction = 0.761)
  abgrHeightLocalImportance = get_height_local_importance(abgr2022, mtry = 5, minNodeSize = 3, samplingFraction = 0.616)
  acmaHeightLocalImportance = get_height_local_importance(acma2022, mtry = 3, minNodeSize = 2, samplingFraction = 0.710)
  otherHeightLocalImportance = get_height_local_importance(other2022, mtry = 3, minNodeSize = 2, samplingFraction = 0.843)
  heightLocalImportance = bind_rows(psmeHeightLocalImportance, abgrHeightLocalImportance, acmaHeightLocalImportance, otherHeightLocalImportance) %>% 
    gather("predictor", "relativeImportance", -speciesModel) %>% spread("speciesModel", "relativeImportance") %>% # transpose
    mutate(predictor = factor(predictor, levels = levels(heightImportance$predictor))) %>%
    mutate(weightedMean = speciesModelWeights$psme * PSME + speciesModelWeights$abgr * ABGR + speciesModelWeights$acma * ACMA + speciesModelWeights$other * other) %>%
    arrange(desc(weightedMean))
  
  dbhImportance = ranger(dbh ~ ., data = variableImportanceTrees2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)),
                         importance = "impurity_corrected", # https://github.com/imbs-hl/ranger/issues/664
                         mtry = 6, splitrule = 'variance', min.node.size = 2, sample.fraction = 0.817,
                         num.threads = htDiaOptions$rangerThreads)
  dbhImportance = tibble(predictor = names(dbhImportance$variable.importance), relativeImportance = 100 * dbhImportance$variable.importance / max(dbhImportance$variable.importance)) %>% 
    arrange(desc(relativeImportance)) %>% 
    mutate(predictor = factor(predictor, levels = predictor))
  
  psmeDbhLocalImportance = get_dbh_local_importance(psme2022, mtry = 6, minNodeSize = 2, samplingFraction = 0.824)
  abgrDbhLocalImportance = get_dbh_local_importance(abgr2022, mtry = 11, minNodeSize = 2, samplingFraction = 0.532)
  acmaDbhLocalImportance = get_dbh_local_importance(acma2022, mtry = 10, minNodeSize = 3, samplingFraction = 0.698)
  otherDbhLocalImportance = get_dbh_local_importance(other2022, mtry = 8, minNodeSize = 3, samplingFraction = 0.743)
  dbhLocalImportance = bind_rows(psmeDbhLocalImportance, abgrDbhLocalImportance, acmaDbhLocalImportance, otherDbhLocalImportance) %>% 
    gather("predictor", "relativeImportance", -speciesModel) %>% spread("speciesModel", "relativeImportance") %>% # transpose
    mutate(predictor = factor(predictor, levels = levels(dbhImportance$predictor))) %>%
    mutate(weightedMean = speciesModelWeights$psme * PSME + speciesModelWeights$abgr * ABGR + speciesModelWeights$acma * ACMA + speciesModelWeights$other * other) %>%
    arrange(desc(weightedMean))
  
  ggplot() +
    geom_tile(aes(x = "all species", y = predictor, fill = relativeImportance), heightImportance) +
    labs(x = NULL, y = NULL, fill = "relative\nimportance", title = paste(plotLetters[1], "height prediction importance")) +
    scale_y_discrete(limits = rev) +
  ggplot() +
    geom_tile(aes(x = speciesModel, y = predictor, fill = relativeImportance), heightLocalImportance %>% select(-weightedMean) %>% pivot_longer(cols = -predictor, names_to = "speciesModel", values_to = "relativeImportance")) +
    labs(x = NULL, y = NULL, fill = "relative\nimportance") +
    scale_x_discrete(limits = c("PSME", "ACMA", "ABGR", "other"), labels = c("Douglas-fir", "bigleaf maple", "grand fir", "other species")) +
    scale_y_discrete(limits = rev, labels = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  ggplot() +
    geom_tile(aes(x = "all species", y = predictor, fill = relativeImportance), dbhImportance) +
    labs(x = NULL, y = NULL, fill = "relative\nimportance", title = paste(plotLetters[2], "DBH prediction importance")) +
    scale_y_discrete(limits = rev) +
  ggplot() +
    geom_tile(aes(x = speciesModel, y = predictor, fill = relativeImportance), dbhLocalImportance %>% select(-weightedMean) %>% pivot_longer(cols = -predictor, names_to = "speciesModel", values_to = "relativeImportance")) +
    labs(x = NULL, y = NULL, fill = "relative\nimportance") +
    scale_x_discrete(limits = c("PSME", "ACMA", "ABGR", "other"), labels = c("Douglas-fir", "bigleaf maple", "grand fir", "other species")) +
    scale_y_discrete(limits = rev, labels = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 1, widths = c(1, 4, 1, 4), guides = "collect") &
    scale_fill_viridis_c() &
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major = element_blank())
  #ggsave("trees/height-diameter/figures/Figure S random forest variable importance.png", height = 9, width = 16.5, units = "cm", dpi = 300)
  
  # ranger tuning for height prediction, 9950X timings
  #                 height                                        diameter
  # species         tune time   mtry  node size  sample fraction  tune time   mtry  node size  sample fraction  notes
  # all             1.00 min    3     2          0.753            1.28 min    6     2          0.817
  # Douglas-fir     42 s        3     2          0.761            49 s        6     2          0.824            36 rows (0.5%) excluded for lack of inventory age
  # grand fir       24 s        5     3          0.616            20 s        11    3          0.532            13 rows (0.6%) excluded for lack of inventory age
  # bigleaf maple   26 s        3     2          0.710            24 s        10    3          0.698            14 rows (0.6%) excluded for lack of inventory age
  # other           22 s        3     2          0.843            20 s        8     2          0.743            11 rows (0.8%) excluded for lack of inventory age
  # See https://mlr.mlr-org.com/articles/tutorial/measures.html#regression-1 for tuning targets. These include MAE and RMSE.
  library(mlr)
  library(tuneRanger)
  
  allHeightRangerTask = makeRegrTask(id = "all2022height", data = as.data.frame(variableImportanceTrees2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "height") # convert isPlanation to numeric as logicals aren't supported
  allHeightTuneStart = Sys.time()
  allHeightTuneResult = tuneRanger(allHeightRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70, # passing htDiaOptions$rangerThreads or a simple variable with the same value errors out
                                    build.final.model = FALSE) # tuneRanger()$model is a wrapped ranger fit which the tuneRanger be loaded as a package and a task be used, easier just to refit with ranger directly and save out that .Rds (tuneRanger()$model$learner.model is a ranger object which can be used directly, but is a probability tree rather than a classification tree)
  allHeightTuneTime = Sys.time() - allHeightTuneStart
  allHeightTuning = tibble(mtry = allHeightTuneResult$recommended.pars$mtry, minNodeSize = allHeightTuneResult$recommended.pars$min.node.size, samplingFraction = allHeightTuneResult$recommended.pars$sample.fraction)

  psmeHeightRangerTask = makeRegrTask(id = "psme2022height", data = as.data.frame(psme2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "height") # convert isPlanation to numeric as logicals aren't supported
  psmeHeightTuneStart = Sys.time()
  psmeHeightTuneResult = tuneRanger(psmeHeightRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70, # passing htDiaOptions$rangerThreads or a simple variable with the same value errors out
                                    build.final.model = FALSE) # tuneRanger()$model is a wrapped ranger fit which the tuneRanger be loaded as a package and a task be used, easier just to refit with ranger directly and save out that .Rds (tuneRanger()$model$learner.model is a ranger object which can be used directly, but is a probability tree rather than a classification tree)
  psmeHeightTuneTime = Sys.time() - psmeHeightTuneStart
  psmeHeightTuning = tibble(mtry = psmeHeightTuneResult$recommended.pars$mtry, minNodeSize = psmeHeightTuneResult$recommended.pars$min.node.size, samplingFraction = psmeHeightTuneResult$recommended.pars$sample.fraction)

  abgrHeightRangerTask = makeRegrTask(id = "abgr2022height", data = as.data.frame(abgr2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "height") # convert isPlanation to numeric as logicals aren't supported
  abgrHeightTuneStart = Sys.time()
  abgrHeightTuneResult = tuneRanger(abgrHeightRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                    build.final.model = FALSE)
  abgrHeightTuneTime = Sys.time() - abgrHeightTuneStart
  abgrHeightTuning = tibble(mtry = abgrHeightTuneResult$recommended.pars$mtry, minNodeSize = abgrHeightTuneResult$recommended.pars$min.node.size, samplingFraction = abgrHeightTuneResult$recommended.pars$sample.fraction)

  acmaHeightRangerTask = makeRegrTask(id = "acma2022height", data = as.data.frame(acma2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "height") # convert isPlanation to numeric as logicals aren't supported
  acmaHeightTuneStart = Sys.time()
  acmaHeightTuneResult = tuneRanger(acmaHeightRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                    build.final.model = FALSE)
  acmaHeightTuneTime = Sys.time() - acmaHeightTuneStart
  acmaHeightTuning = tibble(mtry = acmaHeightTuneResult$recommended.pars$mtry, minNodeSize = acmaHeightTuneResult$recommended.pars$min.node.size, samplingFraction = acmaHeightTuneResult$recommended.pars$sample.fraction)

  otherHeightRangerTask = makeRegrTask(id = "other2022height", data = as.data.frame(other2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(heightCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "height") # convert isPlanation to numeric as logicals aren't supported
  otherHeightTuneStart = Sys.time()
  otherHeightTuneResult = tuneRanger(otherHeightRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                     build.final.model = FALSE)
  otherHeightTuneTime = Sys.time() - otherHeightTuneStart
  otherHeightTuning = tibble(mtry = otherHeightTuneResult$recommended.pars$mtry, minNodeSize = otherHeightTuneResult$recommended.pars$min.node.size, samplingFraction = otherHeightTuneResult$recommended.pars$sample.fraction)

  # ranger tuning for diameter prediction
  allDbhRangerTask = makeRegrTask(id = "all2022dbh", data = as.data.frame(variableImportanceTrees2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "dbh") # convert isPlanation to numeric as logicals aren't supported
  allDbhTuneStart = Sys.time()
  allDbhTuneResult = tuneRanger(allDbhRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                 build.final.model = FALSE)
  allDbhTuneTime = Sys.time() - allDbhTuneStart
  allDbhTuning = tibble(mtry = allDbhTuneResult$recommended.pars$mtry, minNodeSize = allDbhTuneResult$recommended.pars$min.node.size, samplingFraction = allDbhTuneResult$recommended.pars$sample.fraction)
  
  psmeDbhRangerTask = makeRegrTask(id = "psme2022dbh", data = as.data.frame(psme2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "dbh") # convert isPlanation to numeric as logicals aren't supported
  psmeDbhTuneStart = Sys.time()
  psmeDbhTuneResult = tuneRanger(psmeDbhRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                 build.final.model = FALSE)
  psmeDbhTuneTime = Sys.time() - psmeDbhTuneStart
  psmeDbhTuning = tibble(mtry = psbhTuneResult$recommended.pars$mtry, minNodeSize = psmeDbhTuneResult$recommended.pars$min.node.size, samplingFraction = psmeDbhTuneResult$recommended.pars$sample.fraction)
  
  abgrDbhRangerTask = makeRegrTask(id = "abgr2022dbh", data = as.data.frame(abgr2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "dbh") # convert isPlanation to numeric as logicals aren't supported
  abgrDbhTuneStart = Sys.time()
  abgrDbhTuneResult = tuneRanger(abgrDbhRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                 build.final.model = FALSE)
  abgrDbhTuneTime = Sys.time() - abgrDbhTuneStart
  abgrDbhTuning = tibble(mtry = abgrDbhTuneResult$recommended.pars$mtry, minNodeSize = abgrDbhTuneResult$recommended.pars$min.node.size, samplingFraction = abgrDbhTuneResult$recommended.pars$sample.fraction)
  
  acmaDbhRangerTask = makeRegrTask(id = "acma2022dbh", data = as.data.frame(acma2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "dbh") # convert isPlanation to numeric as logicals aren't supported
  acmaDbhTuneStart = Sys.time()
  acmaDbhTuneResult = tuneRanger(acmaDbhRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                 build.final.model = FALSE)
  acmaDbhTuneTime = Sys.time() - acmaDbhTuneStart
  acmaDbhTuning = tibble(mtry = acmaDbhTuneResult$recommended.pars$mtry, minNodeSize = acmaDbhTuneResult$recommended.pars$min.node.size, samplingFraction = acmaDbhTuneResult$recommended.pars$sample.fraction)
  
  otherDbhRangerTask = makeRegrTask(id = "other2022dbh", data = as.data.frame(other2022 %>% filter(is.na(inventoryAge) == FALSE) %>% select(all_of(dbhCandidateVariables)) %>% mutate(isPlantation = as.numeric(isPlantation))), target = "dbh") # convert isPlanation to numeric as logicals aren't supported
  otherDbhTuneStart = Sys.time()
  otherDbhTuneResult = tuneRanger(otherDbhRangerTask, measure = list(mae), num.trees = 500, num.threads = 16, iters = 70,
                                  build.final.model = FALSE)
  otherDbhTuneTime = Sys.time() - otherDbhTuneStart
  otherDbhTuning = tibble(mtry = otherDbhTuneResult$recommended.pars$mtry, minNodeSize = otherDbhTuneResult$recommended.pars$min.node.size, samplingFraction = otherDbhTuneResult$recommended.pars$sample.fraction)
  
  detach("package:tuneRanger", unload = TRUE)
  detach("package:mlrMBO", unload = TRUE)
  detach("package:mlr", unload = TRUE)
  detach("package:smoof", unload = TRUE)
  detach("package:ParamHelpers", unload = TRUE)
  #detach("package:parallel", unload = TRUE)
  #detach("package:checkmate", unload = TRUE)
  detach("package:lhs", unload = TRUE)
}
