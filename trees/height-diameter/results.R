# setup for results analysis
# 1) get regressions and summaries from <species>.R by running <species> job.R for Douglas-fir, grand fir, bigleaf maple, and the other species group
# 2) load functions and forest inventory data by running the first 1000 lines or so of setup.R
# 3) load regression results and analyze them here
library(WeightedROC)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
future::plan(future::multisession, workers = 8)

figureDpi = 300
speciesGroupColors = c("forestgreen", "green3", "blue2", "grey65")

# rank model forms by estimated prediction ability (using AUC) for form selection
# ~24 s for PSME+ABGR+ACMA3+other initial fits, 9900X
get_model_aucs = function(heightDiameterFitStats)
{
  aucStart = Sys.time()
  progressr::with_progress({
    crossValidatedModelCount = heightDiameterFitStats %>% group_by(responseVariable, species, name, fitting) %>% 
      summarize(n = 1, .groups = "drop_last") %>%
      summarize(n = sum(n), .groups = "drop")
    progressBar = progressr::progressor(steps = sum(crossValidatedModelCount$n))
    
    heightDiameterModelAucs = heightDiameterFitStats %>%
      group_by(responseVariable, species, name, fitting) %>%
      group_split() %>%
      furrr::future_map_dfr(function(fitResults)
        #purrr::map_dfr(function(fitResults) # for debugging
      {
        #fitResults = heightDiameterFitStats %>% filter(responseVariable == "DBH", species == "other species", name == "REML GAM", fitting == "gam") # for debugging
        if ((nrow(fitResults) == 1) | all(is.na(fitResults$distinct)) | all(fitResults$distinct == FALSE))
        {
          # two cases
          # - no need to calculate AUCs because this model is not distinct
          # - no distribution to compare to since this model has only a fit failed result
          # Could calculate AUCs for indistinct models but this is not currently of interest and doing so roughly doubles runtime.
          progressBar(str_pad(paste(fitResults$responseVariable[1], fitResults$species[1], fitResults$name[1], fitResults$fitting[1]), 60, "right"))
          #print(paste("single self row or NA model efficiency:", fitResults$species[1], fitResults$responseVariable[1], fitResults$name[1], fitResults$fitting[1]))
          return(tibble(responseVariable = fitResults$responseVariable[1], species = fitResults$species[1], name = fitResults$name[1], fitting = fitResults$fitting[1], distinct = fitResults$distinct[1], isBaseForm = fitResults$isBaseForm[1],
                        otherModelName = NA_character_, otherModelFitting = NA_character_, otherModeldistinct = NA_integer_,
                        hasPhysio = fitResults$hasPhysio[1], hasStand = fitResults$hasStand[1], hasRelative = fitResults$hasRelative[1],
                        aucDeltaAicN = NA_real_, aucMae = NA_real_, aucNse = NA_real_, aucOutOfRange = NA_real_, aucRmse = NA_real_,
                        speciesFraction = fitResults$speciesFraction[1]))
        }
        
        # get all other cross-validations for distinct models for this response variable and species
        # Could include pairings between distinct and not distinct models but this is not currently of interest.
        matchingFitResults = heightDiameterFitStats %>% filter(distinct, responseVariable == fitResults$responseVariable[1], species == fitResults$species[1], ((name == fitResults$name[1]) & (fitting == fitResults$fitting[1])) == FALSE, isBaseForm == fitResults$isBaseForm[1])
        matchingModels = matchingFitResults %>% group_by(name, fitting) %>% summarize(name = name[1], fitting = fitting[1], .groups = "drop")
        #print(paste0(nrow(matchingModels), " AUC intercomparisons...")) # for debugging
        pairwiseAucs = bind_rows(purrr::map2(matchingModels$name, matchingModels$fitting, function(otherModelName, otherModelFitting)
        {
          otherFitResults = matchingFitResults %>% filter(name == otherModelName, fitting == otherModelFitting)
          #otherFitResults = matchingFitResults %>% filter(name == "REML GAM", fitting == "gamm") # for debugging
          
          if ((nrow(otherFitResults) == 1) | all(is.na(otherFitResults$nse)))
          {
            # also no distribution to compare to if the other model has only a no fit result or wasn't cross validated
            #print(paste("single other row or NA model efficiency:", fitResults$species[1], fitResults$responseVariable[1], fitResults$name[1], fitResults$fitting[1]))
            return(tibble(responseVariable = fitResults$responseVariable[1], species = fitResults$species[1], name = fitResults$name[1], fitting = fitResults$fitting[1], distinct = fitResults$distinct[1], isBaseForm = fitResults$isBaseForm[1],
                          otherModelName = otherModelName, otherModelFitting = otherModelFitting, otherModelDistinct = otherFitResults$distinct[1],
                          hasPhysio = fitResults$hasPhysio[1], hasStand = fitResults$hasStand[1], hasRelative = fitResults$hasRelative[1],
                          aucDeltaAicN = NA_real_, aucMae = NA_real_, aucNse = NA_real_, aucOutOfRange = NA_real_, aucRmse = NA_real_,
                          speciesFraction = fitResults$speciesFraction[1]))
          }
          if (nrow(fitResults) != nrow(otherFitResults))
          {
            stop(paste0(fitResults$species[1], " ", fitResults$responseVariable[1], " ", fitResults$name[1], " ", fitResults$fitting[1], " @", nrow(fitResults), " folds x ", otherModelName, " @ ", nrow(otherFitResults), " folds: number of cross validation folds mismatched between fits."), quote = FALSE)
          }
          
          # find AUCs: for species with few samples (e.g. THPL) mean absolute bias may be all NA, which case WeightedROC() errors
          # AUC is the probability a sample from one empirical distribution is greater than a sample from an another 
          # distribution (see AUC.R). If the labels are reversed, the probability changes to a sample from the first distribution
          # being less than a sample from the other distribution.
          # The labels here are thus set up to measure the probabilities goodness of fit metrics from fitResults are preferable
          # to those from otherFitResults.
          lowerIsBetter = factor(c(rep(0, nrow(fitResults)), rep(1, nrow(otherFitResults))), levels = c(0, 1))
          higherIsBetter = factor(c(rep(1, nrow(fitResults)), rep(0, nrow(otherFitResults))), levels = c(0, 1))
          #if ((length(lowerIsBetter) < 2) | (length(higherIsBetter) < 2) | (n_distinct(lowerIsBetter) < 2) | (n_distinct(higherIsBetter) < 2))
          #{
          #  stop(paste0("ROC label formation error with name = ", otherModelName, " for ", fitResults$species[1], " ", fitResults$responseVariable[1], ".  nrow(fitResults) = ", nrow(fitResults), ", nrow(otherFitResults) = ", nrow(otherFitResults), "."))
          #}
          
          # estimate AUC for bias based on what data is available, which is potentially none for species with few measurements
          if ((sum(is.na(fitResults$mae)) + sum(is.na(otherFitResults$mae)) + sum(is.na(fitResults$nse)) + sum(is.na(otherFitResults$nse)) + 
               sum(is.na(fitResults$rmse)) + sum(is.na(otherFitResults$rmse))) > 0)
          {
            # WeightedROC() errors on NAs in its arguments but can't do so informatively
            # Fail informatively instead so that investigation is possible.
            stop(paste0(fitResults$species[1], " ", fitResults$responseVariable[1], " ", fitResults$name[1], " @ ", nrow(fitResults), " folds x ", otherModelName, " @ ", nrow(otherFitResults), " folds:",
                        " range ", sum(is.na(fitResults$pctOutOfRange)), " ", sum(is.na(otherFitResults$pctOutOfRange)),
                        " mae ", sum(is.na(fitResults$mae)), " ", sum(is.na(otherFitResults$mae)),
                        " nse ", sum(is.na(fitResults$nse)), " ", sum(is.na(otherFitResults$nse)),
                        " rmse ", sum(is.na(fitResults$rmse)), " ", sum(is.na(otherFitResults$rmse))), quote = FALSE)
          }
          
          outOfRangeGuess = c(fitResults$pctOutOfRange, otherFitResults$pctOutOfRange)
          maeGuess = c(fitResults$mae, otherFitResults$mae)
          nseGuess = c(fitResults$nse, otherFitResults$nse)
          rmseGuess = c(fitResults$rmse, otherFitResults$rmse)
          expectedLength = nrow(fitResults) + nrow(otherFitResults)
          if ((length(maeGuess) != expectedLength) | (length(nseGuess) != expectedLength) | (length(rmseGuess) != expectedLength))
          {
            stop(paste0(fitResults$species[1], " ", fitResults$responseVariable[1], " ", fitResults$name[1], " @ ", nrow(fitResults), " folds x ", otherModelName, " @ ", nrow(otherFitResults), " folds: NULLs in data", 
                        " MAE ", length(maeGuess),
                        " NSE ", length(nseGuess),
                        " RMSE ", length(rmseGuess)))
          }
          
          aucAicN = NA_real_
          if ((fitResults$fitting[1] != "ranger") & (otherModelFitting != "ranger")) # AIC not defined for random forests, so no AUC for random forest pairings
          {
            deltaAicNguess = c(fitResults$deltaAicN, otherFitResults$deltaAicN)
            aucAicN = WeightedAUC(WeightedROC(guess = deltaAicNguess, label = lowerIsBetter))
          }
          
          #print(paste0("AUC: ", fitResults$species[1], " ", fitResults$responseVariable[1], " ", fitResults$name[1], " ", fitResults$fitting[1], " @ ", nrow(fitResults)," folds x ", otherModelName, " ", otherModelFitting, " @ ", nrow(otherFitResults), " folds."))
          return(tibble(responseVariable = fitResults$responseVariable[1], species = fitResults$species[1], name = fitResults$name[1], fitting = fitResults$fitting[1], distinct = fitResults$distinct[1], isBaseForm = fitResults$isBaseForm[1],
                        otherModelName = otherModelName, otherModelFitting = otherModelFitting, otherModelDistinct = otherFitResults$distinct[1],
                        hasPhysio = fitResults$hasPhysio[1], hasStand = fitResults$hasStand[1], hasRelative = fitResults$hasRelative[1],
                        aucDeltaAicN = aucAicN,
                        aucOutOfRange = WeightedAUC(WeightedROC(guess = outOfRangeGuess, label = lowerIsBetter)),
                        aucMae = WeightedAUC(WeightedROC(guess = maeGuess, label = lowerIsBetter)),
                        aucNse = WeightedAUC(WeightedROC(guess = nseGuess, label = higherIsBetter)),
                        aucRmse = WeightedAUC(WeightedROC(guess = rmseGuess, label = lowerIsBetter)),
                        speciesFraction = fitResults$speciesFraction[1]))
        }))
        
        progressBar(str_pad(paste(fitResults$responseVariable[1], fitResults$species[1], fitResults$name[1]), 60, "right"))
        return(pairwiseAucs)
      }) %>%
      mutate(folds = max(heightDiameterFitStats$fold, na.rm = TRUE),
             repetitions = max(heightDiameterFitStats$repetition, na.rm = TRUE))
    })
  Sys.time() - aucStart
  
  #table(heightDiameterModelAucs$distinct, useNA = "ifany")
  #table(heightDiameterModelAucs$otherModelDistinct, useNA = "ifany")
  #ggplot() +
  #  geom_histogram(aes(x = aucOutOfRange), heightDiameterModelAucs %>% filter(distinct, otherModelDistinct == 1), binwidth = 0.01) +
  #ggplot() +
  #  geom_histogram(aes(x = aucMae), heightDiameterModelAucs %>% filter(distinct, otherModelDistinct == 1), binwidth = 0.01) +
  #ggplot() +
  #  geom_histogram(aes(x = aucRmse), heightDiameterModelAucs %>% filter(distinct, otherModelDistinct == 1), binwidth = 0.01) +
  #ggplot() +
  #  geom_histogram(aes(x = aucDeltaAicN), heightDiameterModelAucs %>% filter(distinct, otherModelDistinct == 1), binwidth = 0.01) +
  #ggplot() +
  #  geom_histogram(aes(x = aucNse), heightDiameterModelAucs %>% filter(distinct, otherModelDistinct == 1), binwidth = 0.01) +
  #plot_annotation(theme = theme(plot.margin = margin())) +
  #plot_layout()
  
  return(heightDiameterModelAucs)
}

get_model_comparison = function(heightDiameterModelRanking, predictedTreeAttribute = "height", heightDiameterModelRankingForDisplaySort = heightDiameterModelRanking)
{
  heightDiameterModelDisplaySort = heightDiameterModelRankingForDisplaySort %>% 
    filter(responseVariable == predictedTreeAttribute, nameAndFit %in% heightDiameterModelRanking$nameAndFit) %>%
    group_by(responseVariable, name, fitting) %>%
    summarize(penalizedBlendedAuc = sum(speciesFraction * if_else(is.na(aucBlended), 0, aucBlended)),
              nameAndFit = nameAndFit[1],
              .groups = "drop") %>% # can also weight not distinct fits differently
    arrange(responseVariable, penalizedBlendedAuc)
  return(heightDiameterModelRanking %>% filter(responseVariable == predictedTreeAttribute) %>%
           mutate(isMixed = fitting %in% c("gamm", "nlme"),
                  nameAndFit = factor(paste0(name, if_else(isMixed, " mixed", "")), levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == predictedTreeAttribute))$nameAndFit)),
                  isMixedGeneralizedDisplay = isMixed & ((aucMaeRank <= 2) | (aucRmseRank <= 2) | (aucDeltaAicNRank <= 2) | (aucNseRank <= 2))) %>%
           group_by(responseVariable, name) %>%
           mutate(isMixedGeneralizedDisplay = any(isMixedGeneralizedDisplay)) %>% # if a mixed effects model is preferred for any species, show AUCs for all species
           ungroup())
}

# take median AUCs over otherModelName, excluding failed and nondistinct fits
get_model_ranking = function(heightDiameterModelAucs)
{
  # AUCs for pairs with not distinct members are (currently) mostly excluded above but filtering for significance is repeated here as some
  # rows for failed fits and not distinct models are present in the AUC tibble. Exclusion of unsuccessful fits is debatable. If a model could 
  # not be fit then its AUC could reasonably be taken to be zero rather than NA, implying all models which could be fit have an AUC of 1 in 
  # comparison.
  # Rankings are only meaningful between distinct model forms.
  #heightDiameterModelAucs %>% filter(responseVariable == "DBH", species == "other species", name == "REML GAM", fitting == "gam") # for debugging
  #heightDiameterModelAucs %>% filter(is.na(distinct)) # for debugging
  #heightDiameterModelRanking %>% filter(responseVariable == "height", name == "Chapman-Richards RelDbh") %>% group_by(species, name, fitting) %>% summarize(n = n())
  #ggplot() +
  #  geom_histogram(aes(x = aucOutOfRange, fill = fitting, group = fitting), heightDiameterModelAucs %>% filter(distinct, otherModelDistinct == 1, responseVariable == "DBH", species == "other species", name == "REML GAM"), binwidth = 0.01) +
  #  coord_cartesian(xlim = c(0, 1)) +
  #  labs(x = "AUC", y = "model fits", fill = NULL)
  heightDiameterModelRanking = heightDiameterModelAucs %>%
    group_by(responseVariable, species, name, fitting) %>%
    summarize(folds = folds[1], repetitions = repetitions[1],
              fitting = fitting[1], isBaseForm = isBaseForm[1], hasPhysio = hasPhysio[1], hasStand = hasStand[1], hasRelative = hasRelative[1],
              aucDeltaAicN = median(if_else((name != otherModelName) & distinct & (otherModelDistinct == 1), aucDeltaAicN, NA_real_), na.rm = TRUE),
              aucMae = median(if_else((name != otherModelName) & distinct & (otherModelDistinct == 1), aucMae, NA_real_), na.rm = TRUE),
              aucNse = median(if_else((name != otherModelName) & distinct & (otherModelDistinct == 1), aucNse, NA_real_), na.rm = TRUE),
              aucOutOfRange = median(if_else((name != otherModelName) & distinct & (otherModelDistinct == 1), aucOutOfRange, NA_real_), na.rm = TRUE),
              aucRmse = median(if_else((name != otherModelName) & distinct & (otherModelDistinct == 1), aucNse, NA_real_), na.rm = TRUE),
              distinct = distinct[1],
              speciesFraction = speciesFraction[1],
              .groups = "drop_last") %>%
    group_by(responseVariable, species, isBaseForm) %>%
    mutate(aucBlended = if_else(distinct, 0.1 * aucOutOfRange + 0.3 * if_else(responseVariable == "height", aucMae, aucRmse) + 0.1 * if_else(responseVariable == "DBH", aucMae, aucRmse) + 0.4 * if_else(fitting == "ranger", 0.5, aucDeltaAicN) + 0.1 * aucNse, 0.5),
           aucDeltaAicNRank = min_rank(desc(round(aucDeltaAicN * distinct, 3))), 
           aucOutOfRangeRank = min_rank(desc(round(aucOutOfRange * distinct, 3))), 
           aucMaeRank = min_rank(desc(round(aucMae * distinct, 3))), 
           aucNseRank = min_rank(desc(round(aucNse * distinct, 3))), 
           aucRmseRank = min_rank(desc(round(aucRmse * distinct, 3)))) %>%
    ungroup() %>%
    mutate(nameAndFit = paste0(name, if_else(fitting %in% c("gamm", "nlme"), " mixed", "")))
  #print(heightDiameterModelDisplaySort, n = 177)
  #heightDiameterModelRanking %>% group_by(responseVariable, species) %>%
  #  summarize(n = n(), outOfRangeN = sum(is.na(aucOutOfRange) == FALSE), maeN = sum(is.na(aucMae) == FALSE), rmseN = sum(is.na(aucRmse) == FALSE), nseN = sum(is.na(aucNse) == FALSE), .groups = "drop")
  return(heightDiameterModelRanking)
}

# summarize predictor variable selection
# TODO: report selection in GAMs and random forests
get_predictor_results = function(heightDiameterCoefficients)
{
  #tibble(parameter = names(heightDiameterCoefficients)) %>% filter(str_starts(parameter, "b1")) %>% arrange(parameter)
  #tibble(parameter = names(heightDiameterCoefficients)) %>% filter(str_ends(parameter, "tr")) %>% arrange(parameter)
  predictorVariableResults = heightDiameterCoefficients %>% 
    filter(fitSet == "primary", fitting %in% c("gsl_nls")) %>% # exclude GAMs and random forests for now, NLME fits have same fixed effects as GSL NLS fits
    mutate(hasPlantation = (is.na(a1p) == FALSE) | (is.na(a2p) == FALSE) | (is.na(a3p) == FALSE) | 
             (is.na(b1p) == FALSE) | (is.na(b2p) == FALSE) | (is.na(b2rdp) == FALSE) | (is.na(b4p) == FALSE) | (is.na(Hap) == FALSE) | (is.na(dp) == FALSE) |
             (is.na(b1abap) == FALSE) | (is.na(b1aatp) == FALSE),
           hasRelDbhOrRelHt = (is.na(a1rd) == FALSE) | (is.na(a1rh) == FALSE) |
             (is.na(a2rd) == FALSE) | 
             (is.na(b1rd) == FALSE) | (is.na(b1rh) == FALSE) | 
             (is.na(b2rd) == FALSE) | (is.na(b2rh) == FALSE) |
             (is.na(b3rd) == FALSE) | 
             (is.na(b4rd) == FALSE),
           hasDistinctRelative = hasRelDbhOrRelHt * distinct,
           nDistinctRelHtOrDbh = ((is.na(a1rd) == FALSE) + (is.na(a1rh) == FALSE) + 
                                    (is.na(a2rd) == FALSE) + 
                                    (is.na(b1rd) == FALSE) + (is.na(b1rh) == FALSE) + 
                                    (is.na(b2rd) == FALSE) + (is.na(b2rh) == FALSE) + 
                                    (is.na(b3rd) == FALSE) + 
                                    (is.na(b4rd) == FALSE)) * distinct,
           hasBasalArea = (is.na(b1ba) == FALSE) | (is.na(b2ba) == FALSE) | (is.na(b4ba) == FALSE) | 
             (is.na(a1bal) == FALSE) | (is.na(a2bal) == FALSE) | (is.na(b1bal) == FALSE) | (is.na(b2bal) == FALSE) | (is.na(b3bal) == FALSE) | (is.na(b4bal) == FALSE) |
             (is.na(a1aba) == FALSE) | (is.na(a1aat) == FALSE) | (is.na(b1aba) == FALSE) | (is.na(b1aat) == FALSE) | (is.na(b2aba) == FALSE) | (is.na(b2aat) == FALSE) |
             ((responseVariable == "height") & str_detect(name, "(Sharma-Parton|Sharma-Zhang)")),
           hasDistinctBasalArea = hasBasalArea * distinct,
           nDistinctBAorAba = ((is.na(b1ba) == FALSE) + (is.na(b2ba) == FALSE) + (is.na(b4ba) == FALSE) + 
                                 (is.na(a1aba) == FALSE) + (is.na(b1aba) == FALSE) + (is.na(b2aba) == FALSE)) * distinct,
           nDistinctBalOrAat = ((is.na(a1bal) == FALSE) + (is.na(a2bal) == FALSE) + (is.na(b1bal) == FALSE) + (is.na(b2bal) == FALSE) + (is.na(b3bal) == FALSE) + (is.na(b4bal) == FALSE) + # a1balr, b1balr, b2balr, ... need not counted as they're always copresent with a1bal, b1bal, ...
                                  (is.na(a1aat) == FALSE) + (is.na(b1aat) == FALSE) + (is.na(b2aat) == FALSE)) * distinct,
           nBasalAreaDistinct = ((is.na(b1ba) == FALSE) + (is.na(b2ba) == FALSE) + (is.na(b4ba) == FALSE) + # Sharma-Parton and Sharma-Zhang height forms are not counted here as there is no test for whether stand basal area is a distinct predictor (modified Sharma-Parton for DBH does not use basal area)
                                   (is.na(a1bal) == FALSE) + (is.na(a2bal) == FALSE) + (is.na(b1bal) == FALSE) + (is.na(b2bal) == FALSE) + (is.na(b3bal) == FALSE) + (is.na(b4bal) == FALSE) +
                                   (is.na(a1aba) == FALSE) + (is.na(a1aat) == FALSE) + (is.na(b1aba) == FALSE) + (is.na(b1aat) == FALSE) + (is.na(b2aba) == FALSE) + (is.na(b2aat) == FALSE)) * distinct,
           hasPhysio = (is.na(a1e) == FALSE) |  (is.na(b1e) == FALSE) |  (is.na(b2e) == FALSE) |                           (is.na(b4e) == FALSE) | 
             (is.na(a1s) == FALSE) |                           (is.na(b2s) == FALSE) |                           
             (is.na(a1as) == FALSE) | (is.na(b1as) == FALSE) | (is.na(b2as) == FALSE) |                          (is.na(b4as) == FALSE) | 
             (is.na(a1ac) == FALSE) | (is.na(b1ac) == FALSE) | (is.na(b2ac) == FALSE) |                          (is.na(b4ac) == FALSE) | 
             (is.na(a1tr) == FALSE) |                                                   (is.na(b3tr) == FALSE) | 
             (is.na(a1tw) == FALSE) | (is.na(b1tw) == FALSE) | (is.na(b2tw) == FALSE) |                          (is.na(b4tw) == FALSE),
           nPhysioDistinct = ((is.na(a1e) == FALSE) +  (is.na(b1e) == FALSE) +  (is.na(b2e) == FALSE) +                           (is.na(b4e) == FALSE) + 
                                (is.na(a1s) == FALSE) +                           (is.na(b2s) == FALSE) +                           
                                (is.na(a1as) == FALSE) + (is.na(b1as) == FALSE) + (is.na(b2as) == FALSE) +                          (is.na(b4as) == FALSE) + 
                                (is.na(a1ac) == FALSE) + (is.na(b1ac) == FALSE) + (is.na(b2ac) == FALSE) +                          (is.na(b4ac) == FALSE) + 
                                (is.na(a1tr) == FALSE) +                                                   (is.na(b3tr) == FALSE) + 
                                (is.na(a1tw) == FALSE) + (is.na(b1tw) == FALSE) + (is.na(b2tw) == FALSE) +                          (is.na(b4tw) == FALSE)) * distinct,
           hasDistinctPhysio = hasPhysio * distinct,
           hasDistinctSharmaB4d = is.na(b4d) * distinct,
           nDistinctElevation = ((is.na(a1e) == FALSE) + (is.na(b1e) == FALSE) + (is.na(b2e) == FALSE) + (is.na(b4e) == FALSE)) * distinct,
           nDistinctSlope = ((is.na(a1s) == FALSE) + (is.na(b2s) == FALSE)) * distinct,
           nDistinctSinAspect = ((is.na(a1as) == FALSE) + (is.na(b1as) == FALSE) + (is.na(b2as) == FALSE) + (is.na(b4as) == FALSE)) * distinct,
           nDistinctCosAspect = ((is.na(a1ac) == FALSE) + (is.na(b1ac) == FALSE) + (is.na(b2ac) == FALSE) + (is.na(b4ac) == FALSE)) * distinct,
           nDistinctTerrainRoughness = ((is.na(a1tr) == FALSE) + (is.na(b3tr) == FALSE)) * distinct,
           nDistinctTopographicWetness = ((is.na(a1tw) == FALSE) + (is.na(b1tw) == FALSE) + (is.na(b2tw) == FALSE) + (is.na(b4tw) == FALSE)) * distinct)
  # check for missing coefficients
  setdiff(names(heightDiameterCoefficients), c("responseVariable", "species", "fitSet", "name", "distinct", "isBaseForm", "hasRelative", "hasStand", "hasPhysio", "fitting", "repetition", "fold", "n", "fitTimeInS", "isConverged", "pctOutOfRange", "mae", "mape", "rmse", "rmspe", "aic", "deltaAicN", "nse", "meanAbsolutePlantationEffect", "meanAbsolutePercentPlantationEffect", "a0", "a1", "a2", "a3", "b1", "b2", "b3", "b4", "b4d", "Ha", "d", "kU", "randomEffect",
                                               "a1p", "a2p", "a3p", "b1p", "b1aatp", "b2p", "b2rdp", "b4p", "Hap", "dp",
                                               "a1rd", "a1rh", "a2rd", "b1rd", "b1rh", "b2rd", "b2rh", "b3rd", "b3rh", "b4rd", "b4rh", 
                                               "a1bal", "a1balr", "a2bal", "a2balr", "b1ba", "b1bal", "b1balr", "b2ba", "b2balr", "b2bal", "b3bal", "b3balr", "b4ba", "b4bal", "b4balr",
                                               "a1aba", "a1aat", "b1aba", "b1abap", "b1aat", "b2aba", "b2aat",
                                               "a1e", "b1e", "b2e", "b4e", "a1s", "b2s", "a1as", "b1as", "b2as", "b4as", "a1ac", "b1ac", "b2ac", "b4ac", "a1tr", "b3tr", "a1tw", "b1tw", "b2tw", "b4tw")) # should be empty
  return(predictorVariableResults)
}

# find preferred model forms
get_preferred_models = function(heightDiameterModelRanking, nPreferredModelForms = 4)
{
  preferredForms = full_join(full_join(heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucMae, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucMae)) %>% mutate(maeName = name, maeFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, maeName, maeFitting, aucMae),
                                       #exclude out of range due to many-way ties
                                       #full_join(heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucOutOfRange, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucOutOfRange)) %>% mutate(outOfRangeName = name, outOfRangeFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, outOfRangeName, outOfRangeFitting, aucOutOfRange),
                                       #          heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucMae, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucMae)) %>% mutate(maeName = name, maeFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, maeName, maeFitting, aucMae),
                                       #          by = join_by(responseVariable, species, isBaseForm, rank)),
                                       full_join(heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucRmse, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucRmse)) %>% mutate(rmseName = name, rmseFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, rmseName, rmseFitting, aucRmse),
                                                 heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucDeltaAicN, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucDeltaAicN)) %>% mutate(aicName = name, aicFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, aicName, aicFitting, aucDeltaAicN),
                                                 by = join_by(responseVariable, species, isBaseForm, rank)),
                                       by = join_by(responseVariable, species, isBaseForm, rank)),
                             full_join(heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucNse, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucNse)) %>% mutate(nseName = name, nseFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, nseName, nseFitting, aucNse),
                                       heightDiameterModelRanking %>% filter(distinct) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucBlended, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucBlended)) %>% mutate(blendedName = name, blendedFitting = fitting, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, blendedName, blendedFitting, aucBlended),
                                       by = join_by(responseVariable, species, isBaseForm, rank)),
                             by = join_by(responseVariable, species, isBaseForm, rank)) %>%
    arrange(desc(responseVariable), species, desc(isBaseForm), rank) %>%
    relocate(responseVariable, species, isBaseForm, rank) %>%
    ungroup()
  #preferredForms %>% group_by(responseVariable) %>% mutate(speciesPresent = n_distinct(species)) %>% group_by(responseVariable, species) %>% summarize(speciesPresent = speciesPresent[1], nPreferredForms = n(), .groups = "drop")
  return(preferredForms)
}
  
read_fit_stats = function(folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions)
{
  psmeFitStats = readRDS(paste0("trees/height-diameter/data/PSME fit stats ", folds, "x", repetitions, ".Rds")) 
  abgrFitStats = readRDS(paste0("trees/height-diameter/data/ABGR fit stats ", folds, "x", repetitions, ".Rds")) 
  acmaFitStats = readRDS(paste0("trees/height-diameter/data/ACMA3 fit stats ", folds, "x", repetitions, ".Rds")) 
  otherFitStats = readRDS(paste0("trees/height-diameter/data/other fit stats ", folds, "x", repetitions, ".Rds")) 

  heightDiameterFitStats = bind_rows(psmeFitStats, abgrFitStats, acmaFitStats, otherFitStats) %>%
    mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name)),
           species = factor(species, labels = c("Douglas-fir", "grand fir", "bigleaf maple", "other species"), levels = c("PSME", "ABGR", "ACMA3", "other")),
           speciesFraction = recode(species, "Douglas-fir" = 0.543, "bigleaf maple" = 0.183, "grand fir" = 0.144, "other species" = 0.130),
           isBaseForm = (str_detect(name, "Sharma-") == FALSE) & (str_detect(name, "ABA\\+T") == FALSE) & (str_detect(name, "BA\\+L") == FALSE) & (str_detect(name, "physio") == FALSE) & (str_detect(name, "RelDbh") == FALSE) & (str_detect(name, "RelHt") == FALSE),
           hasPhysio = str_detect(name, "physio"),
           hasStand = str_detect(name, "ABA\\+T") | str_detect(name, "BA\\+L"),
           hasRelative = str_detect(name, "RelDbh") | str_detect(name, "RelHt"),
           distinct = as.logical(distinct)) %>% # since R lacks NA_logical_ distinct can end up being either of type double (0/1/NA_real_) or logical (TRUE/FALSE), standardize back to logical (TRUE/FALSE/NA)
    group_by(fitSet, responseVariable, fitting, name, species, repetition, fold) %>% # calculate out of range fraction within each fold
    mutate(pctOutOfRange = 100 * nPredictionOutOfRange / n) %>%
    group_by(fitSet, responseVariable, species) %>% # ΔAICn within response variable and species, needed for AUCs and model selection
    mutate(deltaAicN = aic/n - min(aic/n, na.rm = TRUE)) %>%
    ungroup()
  return(heightDiameterFitStats)
}

read_model_coefficients = function(heightDiameterFitStats, folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions)
{
  psmeCoefficients = readRDS(paste0("trees/height-diameter/data/PSME coefficients ", folds, "x", repetitions, ".Rds")) 
  abgrCoefficients = readRDS(paste0("trees/height-diameter/data/ABGR coefficients ", folds, "x", repetitions, ".Rds")) 
  acmaCoefficients = readRDS(paste0("trees/height-diameter/data/ACMA3 coefficients ", folds, "x", repetitions, ".Rds")) 
  otherCoefficients = readRDS(paste0("trees/height-diameter/data/other coefficients ", folds, "x", repetitions, ".Rds")) 

  heightDiameterCoefficients = left_join(bind_rows(psmeCoefficients, abgrCoefficients, acmaCoefficients, otherCoefficients) %>%
                                           mutate(species = factor(species, labels = c("Douglas-fir", "grand fir", "bigleaf maple", "other species"), levels = c("PSME", "ABGR", "ACMA3", "other"))),
                                         heightDiameterFitStats %>% select(fitSet, responseVariable, species, name, fitting, repetition, fold, isBaseForm, hasRelative, hasStand, hasPhysio, n, fitTimeInS, pctOutOfRange, mae, mape, rmse, rmspe, deltaAicN, nse, meanAbsolutePlantationEffect,	meanAbsolutePercentPlantationEffect), # join in extracted model flags and supporting properties flowed through fit stats
                                         by = join_by(fitSet, responseVariable, species, name, fitting, repetition, fold)) %>%
    mutate(isConverged = as.logical(isConverged)) %>%
    select(-fixedEffects) %>%
    relocate(responseVariable, species, fitSet, name, distinct, isBaseForm, hasRelative, hasStand, hasPhysio, fitting, repetition, fold, n, fitTimeInS, isConverged, pctOutOfRange, mae, mape, rmse, rmspe, deltaAicN, nse, meanAbsolutePlantationEffect,	meanAbsolutePercentPlantationEffect, 
             a0, a1, a1p, a1rd, a1rh, a1bal, a1balr, a1aba, a1aat, a1e, a1s, a1as, a1ac, a1tr, a1tw, a2, a2rd, a2bal, a2balr, a2p, a3, a3p, b1, b1p, b1rd, b1ba, b1bal, b1balr, b1aba, b1abap, b1aat, b1aatp, b1rh, b1as, b1ac, b1e, b1tw, b2, b2p, b2rd, b2rdp, b2rh, b2ba, b2bal, b2balr, b2aba, b2aat, b2s, b2as, b2ac, b2e, b2tw, b3, b3rd, b3bal, b3balr, b3tr, b4, b4p, b4d, b4rd, b4ba, b4bal, b4balr, b4e, b4as, b4ac, b4tw, Ha, Hap, d, any_of(c("dp")), kU, randomEffect)
  return(heightDiameterCoefficients)
}

## assemble result tibbles from species results
heightDiameterFitStats2x50 = read_fit_stats(folds = 2, repetitions = 50)
heightDiameterCoefficients2x50 = read_model_coefficients(heightDiameterFitStats2x50, folds = 2, repetitions = 50)
heightDiameterFitStats10x10 = read_fit_stats(folds = 10, repetitions = 10)
heightDiameterCoefficients10x10 = read_model_coefficients(heightDiameterFitStats10x10, folds = 10, repetitions = 10)
heightDiameterFitStats10x25 = read_fit_stats(folds = 10, repetitions = 25)
heightDiameterCoefficients10x25 = read_model_coefficients(heightDiameterFitStats10x25, folds = 10, repetitions = 25)
heightDiameterFitStats10x50 = read_fit_stats(folds = 10, repetitions = 50)
heightDiameterCoefficients10x50 = read_model_coefficients(heightDiameterFitStats10x50, folds = 10, repetitions = 50)

# report duplicate naming and fit failures
heightDiameterFitStats10x10 %>% group_by(fitSet, fitting, responseVariable, species, name) %>% summarize(n = n(), distinct = any(distinct), .groups = "drop") %>% filter(n != htDiaOptions$folds * htDiaOptions$repetitions, ((n == 1) & (distinct == FALSE)) == FALSE)
# report unexpected NA AICs
#heightDiameterFitStats10x10 %>% filter(is.na(deltaAicN)) %>% group_by(fitSet, fitting, responseVariable, species, name) %>% summarize(n = n(), distinct = any(distinct), .groups = "drop") %>% filter(fitting != "ranger", ((fitting == "gam") & (n == 1) & (distinct == FALSE)) == FALSE)
#names(heightDiameterCoefficients10x10) # check for unordered coefficient columns at end of list
#writexl::write_xlsx(heightDiameterCoefficients10x10 %>% filter(fitSet == "primary") %>% select(-fitSet) %>% mutate(distict = as.logical(distinct)), "trees/height-diameter/data/height-diameter model coefficients 10x10.xlsx")

heightDiameterPredictors2x50 = get_predictor_results(heightDiameterCoefficients2x50)
heightDiameterModelAucs2x50 = get_model_aucs(heightDiameterFitStats2x50)
heightDiameterModelRanking2x50 = get_model_ranking(heightDiameterModelAucs2x50)
heightDiameterModelsPreferred2x50 = get_preferred_models(heightDiameterModelRanking2x50)

heightDiameterPredictors10x10 = get_predictor_results(heightDiameterCoefficients10x10)
heightDiameterModelAucs10x10 = get_model_aucs(heightDiameterFitStats10x10)
heightDiameterModelRanking10x10 = get_model_ranking(heightDiameterModelAucs10x10)
heightDiameterModelsPreferred10x10 = get_preferred_models(heightDiameterModelRanking10x10)

heightDiameterPredictors10x25 = get_predictor_results(heightDiameterCoefficients10x25)
heightDiameterModelAucs10x25 = get_model_aucs(heightDiameterFitStats10x25)
heightDiameterModelRanking10x25 = get_model_ranking(heightDiameterModelAucs10x25)
heightDiameterModelsPreferred10x25 = get_preferred_models(heightDiameterModelRanking10x25)

heightDiameterPredictors10x50 = get_predictor_results(heightDiameterCoefficients10x50)
heightDiameterModelAucs10x50 = get_model_aucs(heightDiameterFitStats10x50)
heightDiameterModelRanking10x50 = get_model_ranking(heightDiameterModelAucs10x50)
heightDiameterModelsPreferred10x50 = get_preferred_models(heightDiameterModelRanking10x50)

predictorVariableStats10x10 = heightDiameterPredictors10x10 %>% 
  group_by(responseVariable, species) %>%
  summarize(n = n(),
            hasBasalArea = sum(hasBasalArea),
            hasPhysio = sum(hasPhysio),
            hasPlantation = sum(hasPlantation),
            hasRelDbhOrRelHt = sum(hasRelDbhOrRelHt),
            hasDistinctBasalArea = sum(hasDistinctBasalArea),
            hasDistinctPhysio = sum(hasDistinctPhysio),
            nBasalAreaDistinct = sum(nBasalAreaDistinct),
            nPhysioDistinct = sum(nPhysioDistinct),
            nDistinctRelHtOrDbh = sum(nDistinctRelHtOrDbh),
            nDistinctBAorAba = sum(nDistinctBAorAba),
            nDistinctBalOrAat = sum(nDistinctBalOrAat),
            nDistinctElevation = sum(nDistinctElevation),
            nDistinctSlope = sum(nDistinctSlope),
            nDistinctSinAspect = sum(nDistinctSinAspect),
            nDistinctCosAspect = sum(nDistinctCosAspect),
            nDistinctTerrainRoughness = sum(nDistinctTerrainRoughness),
            nDistinctTopographicWetness = sum(nDistinctTopographicWetness),
            nDistinct = sum(distinct),
            plantationPct = 100 * hasPlantation / n, 
            nDistinctPct = 100 * nDistinct / n,
            nDistinctRelativePct = 100 * nDistinctRelHtOrDbh / hasRelDbhOrRelHt,
            .groups = "drop")


## summary for Abstract (plus AUC counting for Section 2.3)
# form counts, changes in model efficiency with inclusion of generalizing predictors
# 39 height forms + 38 DBH forms = 77 + 3 control forms => 539 * k * r fits
#
# number of AUCs per k x r cross validation combination (Section 2.3)
#   number of unique pairs in set = 1/2 * n * (n - 1) => 741 height pairs + 703 DBH pairs
#   number of AUCs = pairs * five goodness of fit stats * seven species = (TBD + TBD) * 5 * 7 = TBD max, TBD actual
#   but actual number is reduced
(modelCounts = heightDiameterFitStats10x10 %>% filter(fitSet == "primary") %>% 
  group_by(responseVariable, species) %>%
  summarize(nForms = length(unique(name)), nFits = n(), nAucsIfAllFitsSucceeded = 0.5 * nForms * (nForms - 1) * 5))
modelCounts %>% summarize(nAucsIfAllFitsSucceeded = sum(nAucsIfAllFitsSucceeded))

(actualAucs = heightDiameterModelAucs10x10 %>% filter(name != otherModelName) %>%
    group_by(responseVariable, species) %>% 
    mutate(otherBaseName = if_else(word(otherModelName) %in% c("REML", "modified", "unified"), paste(word(otherModelName, 1), word(otherModelName, 2)), word(otherModelName)),
           otherIsBaseForm = (str_detect(otherBaseName, "Sharma-") == FALSE) & (str_detect(otherBaseName, "ABA\\+T") == FALSE) & (str_detect(otherBaseName, "BA\\+L") == FALSE) & (str_detect(otherBaseName, "physio") == FALSE) & (str_detect(otherBaseName, "RelDbh") == FALSE) & (str_detect(otherBaseName, "RelHt") == FALSE),
           isBaseToBaseOrSharmaToSharma = (isBaseForm | (name %in% c("Sharma-Parton", "Sharma-Zhang", "modified Sharma-Parton"))) & (otherIsBaseForm | (name %in% c("Sharma-Parton", "Sharma-Zhang", "modified Sharma-Parton")))) %>%
    summarize(nPairs = n(), nBaseSharma = sum(isBaseToBaseOrSharmaToSharma, na.rm = TRUE)) %>% 
    mutate(nAucs = ceiling(0.5 * nPairs) * 5,
           nBaseSharmaAucs = ceiling(0.5 * nBaseSharma) * 5,
           nTotalAucs = sum(nAucs),
           nTotalBaseSharmaAucs = sum(nBaseSharmaAucs)))

# generalization effect sizes for stand and physiographic variables
heightDiameterFitStats10x10 %>% filter(distinct) %>% 
  mutate(deltaAicN = aic/n - min(aic/n, na.rm = TRUE)) %>% # change to unrestricted ΔAICn
  filter(baseName %in% c("Chapman-Richards", "Ruark", "Sharma-Parton", "Sharma-Zhang", "Sibbesen", "Weibull", "REML GAM")) %>% # remove base forms which weren't generalized
  group_by(responseVariable, species, baseName) %>%
  mutate(nseBase = median(if_else(isBaseForm | (baseName == "Sharma-Parton") | (baseName == "Sharma-Zhang"), nse, NA_real_), na.rm = TRUE),
         nseStandRelDelta = if_else((hasStand | hasRelative) & (hasPhysio == FALSE), nse, NA_real_) - nseBase,
         nsePhysioDelta = if_else(hasPhysio & (hasStand == FALSE) & (hasRelative == FALSE), nse, NA_real_) - nseBase,
         nseCombinedDelta = if_else(hasPhysio & (hasStand | hasRelative), nse, NA_real_) - nseBase) %>%
  #group_by(fitSet, responseVariable, species) %>%
  #slice_max(nse, n = 3 * 10 * 10) %>% # nmodels * nfolds * nrepetitions
  group_by(fitSet, responseVariable, baseName) %>%
  summarize(nForms = n_distinct(name),
            nFits = n(),
            mae = mean(mae),
            mape = mean(mape),
            rmse = mean(rmse),
            #deltaAicN = mean(deltaAicN),
            nseMean = mean(nse),
            nseBase = nseBase[1],
            nseStandRelDelta = median(nseStandRelDelta, na.rm = TRUE),
            nsePhysioDelta = median(nsePhysioDelta, na.rm = TRUE),
            nseCombinedDelta = median(nseCombinedDelta, na.rm = TRUE),
            .groups = "drop") %>%
  #filter(baseName %in% c("Chapman-Richards", "Sharma-Parton", "Ruark", "Sibbesen", "REML GAM"), (responseVariable != "height") | (baseName != "Sibbesen")) %>%
  arrange(desc(fitSet), desc(responseVariable))

# highest ranked base forms
# slice_max(auc) is likely to yield ties when a sufficiently large number of poor models are compared to good models that multiple good models have median AUCs of 
# 1. left_join() raises warnings about multiple rows matching in this case, which is obscure but indicates something needs to be done to disambiguate the tie and
# select a single most preferred model. Tie breaking is currently based on parsimony but other factors, such as the mean AUC, can be considered besides the median
# AUC.
get_top_base_form = function(preferredForms, statisticAucSuffix)
{
  if (statisticAucSuffix == "DeltaAicN")
  {
    aucPrefix = "aic"
  } else if (statisticAucSuffix == "OutOfRange") {
    aucPrefix = "outOfRange"
  } else {
    aucPrefix = str_to_lower(statisticAucSuffix)
  }
  return(preferredForms %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(.data[[str_c("auc", statisticAucSuffix)]], n = 1) %>%
           mutate(nTiedForms = n(), isMixed = maeFitting %in% c("nlme", "gamm")) %>% # for now, break ties by excluding more complex mixed models in favor of fixed effect models
           filter((nTiedForms == 1) | (isMixed == FALSE)) %>%
           select(responseVariable, species, all_of(str_c(aucPrefix, "Name")), all_of(str_c(aucPrefix, "Fitting"))))
}

#heightDiameterModelsPreferred10x10 %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(aucOutOfRange, n = 1) %>% select(responseVariable, species, outOfRangeName, outOfRangeFitting, aucOutOfRange)
#heightDiameterModelsPreferred10x10 %>% filter(responseVariable == "DBH", species == "other species", outOfRangeName == "REML GAM")
(topBaseForms = left_join(left_join(left_join(get_top_base_form(heightDiameterModelsPreferred10x10, "Mae"), get_top_base_form(heightDiameterModelsPreferred10x10, "Rmse"), by = join_by(responseVariable, species)), # left_join(get_top_base_form(heightDiameterModelsPreferred10x10, "OutOfRange"), 
                                    get_top_base_form(heightDiameterModelsPreferred10x10, "DeltaAicN"), by = join_by(responseVariable, species)),
                          get_top_base_form(heightDiameterModelsPreferred10x10, "Nse"), by = join_by(responseVariable, species)) %>%
    mutate(isGam = (maeName == "REML GAM") + (rmseName == "REML GAM") + (aicName == "REML GAM") + (nseName == "REML GAM"), # (outOfRangeName == "REML GAM") + 
           isRandomForest = (maeName == "random forest") + (rmseName == "random forest") + (aicName == "random forest") + (nseName == "random forest")) %>% # (outOfRangeName == "random forest") + 
    relocate(responseVariable, species, maeName, rmseName, aicName, nseName, isGam, isRandomForest) %>% # outOfRangeName, 
    arrange(desc(responseVariable)))
(pctGamPreference = 100 * sum(topBaseForms$isGam) / (5 * nrow(topBaseForms))) # maybe for abstract
(pctRandomForestPreference = 100 * sum(topBaseForms$isRandomForest) / (5 * nrow(topBaseForms))) # maybe for abstract
#print(topBaseForms, width = Inf)


## summaries for Results
# fraction of nonlinear base forms which are more accurate than a base GAM
improvementThresholdProbability = 0.5
heightDiameterModelAucs10x10 %>% filter(name != "REML GAM", otherModelName == "REML GAM", isBaseForm) %>%
  group_by(responseVariable, species) %>% # restriction to primary fit set and base forms is above
  summarize(deltaAic = 100 * sum(aucDeltaAicN > improvementThresholdProbability) / n(),
            mab = 100 * sum(aucOutOfRange > improvementThresholdProbability) / n(),
            mae = 100 * sum(aucMae > improvementThresholdProbability) / n(),
            nse = 100 * sum(aucNse > improvementThresholdProbability) / n(),
            rmse = 100 * sum(aucRmse > improvementThresholdProbability) / n(),
            aicMostPreferred = name[which.max(aucDeltaAicN)],
            .groups = "drop") %>%
  arrange(desc(responseVariable))

heightDiameterFitStats10x10 %>% group_by(responseVariable, species) %>% # plantation effect quantiles
  reframe(quantiles = c(0.05, 0.50, 0.95), plantationFx = quantile(meanAbsolutePlantationEffect, probs = quantiles, na.rm = TRUE), plantationFxPct = quantile(meanAbsolutePercentPlantationEffect, probs = quantiles, na.rm = TRUE)) %>%
  pivot_wider(id_cols = c("responseVariable", "species"), names_from = "quantiles", values_from = c("plantationFx", "plantationFxPct")) %>%
  arrange(desc(responseVariable))

# AIC selection of generalizing predictors based on ΔAICn AUC
improvementThresholdProbability = 0.5
heightDiameterModelAucs10x10 %>% 
  mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name))) %>%
  filter(isBaseForm == FALSE, fitting %in% c("gsl_nls", "nlrob"), distinct, baseName %in% c("Chapman-Richards", "Sharma-Parton", "Ruark", "Sibbesen"), (responseVariable == "height") | (baseName != "Sibbesen")) %>%
  group_by(responseVariable) %>%
  summarize(n = sum(is.na(aucDeltaAicN) == FALSE), 
            baseBetter = sum(aucDeltaAicN <= improvementThresholdProbability), # base form is ΔAICn preferable
            generalizedBetter = sum(aucDeltaAicN > improvementThresholdProbability), # a generalized form is preferable
            physioAicMedianAuc = median(if_else(hasPhysio & (hasStand == FALSE) & (hasRelative == FALSE), aucDeltaAicN, NA_real_), na.rm = TRUE), 
            standRelAicMedianAuc = median(if_else((hasPhysio == FALSE) & (hasStand | hasRelative), NA_real_, aucDeltaAicN), na.rm = TRUE),
            bothAicMedianAuc = median(if_else(hasPhysio & (hasStand | hasRelative), NA_real_, aucDeltaAicN), na.rm = TRUE)) %>%
  mutate(baseBetterPct = 100 * baseBetter / n, 
         genBetterPct = 100 * generalizedBetter / n) %>%
  arrange(desc(responseVariable))

# frequency of generalizing predictor selection
predictorVariableStats10x10 %>%
  mutate(balAat = nDistinctBalOrAat / hasDistinctBasalArea, baAba = nDistinctBAorAba / hasDistinctBasalArea,
         elev = nDistinctElevation / hasDistinctPhysio, slope = nDistinctSlope / hasDistinctPhysio, sinAspect = nDistinctSinAspect / hasDistinctPhysio, cosAspect = nDistinctCosAspect / hasDistinctPhysio, roughness = nDistinctTerrainRoughness / hasDistinctPhysio, wetness = nDistinctTopographicWetness / hasDistinctPhysio) %>%
  select(responseVariable, species, balAat, baAba, elev, slope, sinAspect, cosAspect, roughness, wetness) %>%
  arrange(desc(responseVariable)) # since relative DBH and height are singleton groups they are always selected with probability 1 if when distinct

# highest ranked model forms by goodness of fit metric
# Potential issue: while with_ties = FALSE avoids multijoin situations by forcing slice_max() to choose only one model the most parsimonious model is selected
# only by chance.
print(left_join(left_join(left_join(heightDiameterModelsPreferred10x10 %>% group_by(responseVariable, species) %>% slice_max(aucMae, n = 1, with_ties = FALSE) %>% select(responseVariable, species, maeName, maeFitting),
                                    #left_join(heightDiameterModelsPreferred10x10 %>% group_by(responseVariable, species) %>% slice_max(aucOutOfRange, n = 1, with_ties = FALSE) %>% select(responseVariable, species, outOfRangeName, outOfRangeFitting),
                                    #          heightDiameterModelsPreferred10x10 %>% group_by(responseVariable, species) %>% slice_max(aucMae, n = 1, with_ties = FALSE) %>% select(responseVariable, species, maeName, maeFitting),
                                    #          by = join_by(responseVariable, species)),
                                    heightDiameterModelsPreferred10x10 %>% group_by(responseVariable, species) %>% slice_max(aucRmse, n = 1, with_ties = FALSE) %>% select(responseVariable, species, rmseName, rmseFitting),
                                    by = join_by(responseVariable, species)),
                          heightDiameterModelsPreferred10x10 %>% group_by(responseVariable, species) %>% slice_max(aucDeltaAicN, n = 1, with_ties = FALSE) %>% select(responseVariable, species, aicName, aicFitting),
                          by = join_by(responseVariable, species)),
                heightDiameterModelsPreferred10x10 %>% group_by(responseVariable, species) %>% slice_max(aucNse, n = 1, with_ties = FALSE) %>% select(responseVariable, species, nseName, nseFitting),
                by = join_by(responseVariable, species)) %>%
        pivot_longer(cols = ends_with("Name"), names_to = "statistic", names_pattern = "(.*)Name", values_to = "name") %>%
        mutate(fitting = case_match(statistic, "mae" ~ maeFitting, "rmse" ~ rmseFitting, "aic" ~ aicFitting, "nse" ~ nseFitting)) %>% # "outOfRange" ~ outOfRangeFitting, 
        select(-maeFitting, -rmseFitting, -aicFitting, -nseFitting) %>% # -outOfRangeFitting, 
        relocate(responseVariable, species, statistic, name) %>%
        arrange(desc(responseVariable)), # %>%
      #group_by(responseVariable, species) %>%
      #summarize(modelsSelected = length(unique(name)))
      n = 70)

# summarize fitting success
# Show results for all fit sets here, whether or not primary.
heightDiameterFitStats10x10 %>% filter(is.na(distinct) | (distinct == TRUE)) %>%
  group_by(fitSet, responseVariable, species, name, fitting) %>%
  # reduce cross validated fits to a single record
  summarize(gslNls = fitting[1] == "gsl_nls", 
            lm = fitting[1] == "lm", 
            gam = fitting[1] == "gam",
            nlme = fitting[1] == "nlme",
            fail = is.na(distinct[1]),
            .groups = "drop") %>%
  # summarize across fitting attempts at equal weight, regardless of cross validation settings
  group_by(fitSet) %>%
  summarize(n = n(),
            gslNls = sum(gslNls, na.rm = TRUE), 
            lm = sum(lm, na.rm = TRUE), 
            gam = sum(gam, na.rm = TRUE),
            nlme = sum(nlme, na.rm = TRUE),
            fail = sum(fail),
            uniqueNonlinear = n() - lm - gam) %>%
  mutate(gslNlsPct = 100 * gslNls / uniqueNonlinear, 
         nlmePct = 100 * nlme / uniqueNonlinear,
         totalFailPct = 100 * fail / uniqueNonlinear) %>%
  arrange(desc(fitSet))
# fittings which failed
heightDiameterFitStats10x10 %>% filter(is.na(distinct)) %>% select(fitSet, responseVariable, species, name, fitting) %>% arrange(desc(fitSet))

# summaries for supplemental material
standVariableSelection10x10 = heightDiameterPredictors10x10 %>% filter(distinct == 1, isBaseForm == FALSE) %>% 
  group_by(responseVariable, species, name, fitting) %>%
  summarize(hasBasalArea = hasBasalArea[1], hasRelative = hasRelative[1],
            .groups = "drop")
standVariableSelection10x10 %>% group_by(responseVariable, species) %>%
  summarize(n = n(), pBasalArea = sum(hasBasalArea) / n(), pRelative = sum(hasRelative) / n())


## dataset summary figure
# ggplot() +
#   stat_summary_2d(aes(x = dbhClass, y = heightClass, z = treeCount), psme2022 %>% mutate(dbhClass = 2.5 * round(dbh/2.5, 0), heightClass = round(height)), binwidth = c(2.5, 1), fun = sum) +
#   coord_cartesian(xlim = c(0, 250), ylim = c(0, 75)) +
#   labs(x = "DBH, cm", y = "height, m", fill = "trees\nmeasured", title = paste(plotLetters[1], "Douglas-fir")) +
# ggplot() +
#   stat_summary_2d(aes(x = dbhClass, y = heightClass, z = treeCount), abgr2022 %>% mutate(dbhClass = 2.5 * round(dbh/2.5, 0), heightClass = round(height)), binwidth = c(2.5, 1), fun = sum) +
#   coord_cartesian(xlim = c(0, 150), ylim = c(0, 75)) +
#   labs(x = "DBH, cm", y = NULL, fill = "trees\nmeasured", title = paste(plotLetters[2], "grand fir")) +
# ggplot() +
#   stat_summary_2d(aes(x = dbhClass, y = heightClass, z = treeCount), acma2022 %>% mutate(dbhClass = 2.5 * round(dbh/2.5, 0), heightClass = round(height)), binwidth = c(2.5, 1), fun = sum) +
#   coord_cartesian(xlim = c(0, 150), ylim = c(0, 75)) +
#   labs(x = "DBH, cm", y = NULL, fill = "trees\nmeasured", title = paste(plotLetters[3], "bigleaf maple")) +
# ggplot() +
#   stat_summary_2d(aes(x = dbhClass, y = heightClass, z = treeCount), other2022 %>% mutate(dbhClass = 2.5 * round(dbh/2.5, 0), heightClass = round(height)), binwidth = c(2.5, 1), fun = sum) +
#   coord_cartesian(xlim = c(0, 150), ylim = c(0, 75)) +
#   labs(x = "DBH, cm", y = NULL, fill = "trees\nmeasured", title = paste(plotLetters[4], "other species")) +
# plot_annotation(theme = theme(plot.margin = margin())) +
# plot_layout(guides = "collect", widths = c(250, 150, 150, 150)) &
#   scale_fill_viridis_c(trans = "log10", breaks = c(1, 3, 10, 30, 100, 300, 700), labels = c(1, 3, 10, 30, 100, 300, "700+"), limits = c(1, 700), na.value = "yellow")
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.png", height = 13, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.tif", height = 13, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.pdf", height = 13, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figures 1 and 2: height-diameter AUCs
heightFromDiameterModelComparison2x50 = get_model_comparison(heightDiameterModelRanking2x50, heightDiameterModelRankingForDisplaySort = heightDiameterModelRanking10x10)
heightFromDiameterModelComparison10x10 = get_model_comparison(heightDiameterModelRanking10x10)
heightFromDiameterModelComparison10x25 = get_model_comparison(heightDiameterModelRanking10x25, heightDiameterModelRankingForDisplaySort = heightDiameterModelRanking10x10)
heightFromDiameterModelComparison10x50 = get_model_comparison(heightDiameterModelRanking10x50, heightDiameterModelRankingForDisplaySort = heightDiameterModelRanking10x10)

plot_auc_bank(heightFromDiameterModelComparison2x50, heightFromDiameterModelComparison10x10)
plot_auc_bank(heightFromDiameterModelComparison10x25, heightFromDiameterModelComparison10x10)
plot_auc_bank(heightFromDiameterModelComparison10x50, heightFromDiameterModelComparison10x25)
#ggsave("trees/height-diameter/figures/Figure 01 height AUCs base.png", height = 20, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 01 height AUCs base.tif", height = 20, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 01 height AUCs base.pdf", height = 20, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)
plot_auc_bank(heightFromDiameterModelComparison10x50, heightFromDiameterModelComparison10x25, generalized = TRUE, legendHjustification = 1.05, plotRightMargin = 10)
#ggsave("trees/height-diameter/figures/Figure 02 height AUCs generalized.png", height = 17, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 02 height AUCs generalized.tif", height = 20, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 02 height AUCs generalized.pdf", height = 20, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)

## Figure 3: diameter AUCs
diameterFromHeightModelComparison2x50 = get_model_comparison(heightDiameterModelRanking2x50, "DBH", heightDiameterModelRanking10x10)
diameterFromHeightModelComparison10x10 = get_model_comparison(heightDiameterModelRanking10x10, "DBH")
diameterFromHeightModelComparison10x25 = get_model_comparison(heightDiameterModelRanking10x25, "DBH", heightDiameterModelRanking10x10)
diameterFromHeightModelComparison10x50 = get_model_comparison(heightDiameterModelRanking10x50, "DBH", heightDiameterModelRanking10x10)

plot_auc_bank(diameterFromHeightModelComparison10x50 %>% filter(isBaseForm), diameterFromHeightModelComparison10x25 %>% filter(isBaseForm), legendHjustification = 1.02)
#ggsave("trees/height-diameter/figures/Figure 03 diameter AUCs base.png", height = 20.5, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 03 diameter AUCs base.tif", height = 20.5, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 03 diameter AUCs base", height = 20.5, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)
plot_auc_bank(diameterFromHeightModelComparison10x50, diameterFromHeightModelComparison10x25, generalized = TRUE, legendHjustification = 1.02)
#ggsave("trees/height-diameter/figures/Figure 04 diameter AUCs generalized.png", height = 17.5, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 04 diameter AUCs generalized.tif", height = 20.5, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 04 diameter AUCs generalized", height = 20.5, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 5: model efficiency
# Long tail of negative model efficiencies is suppressed as histogram plotting becomes very slow even if bins
# are not within the display window set by coord_cartesian().
ggplot(heightDiameterFitStats10x10 %>% filter(responseVariable == "height", isBaseForm)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = NULL, y = "fraction of fits, %", fill = NULL, title = bquote(.(plotLetters[1])~"base form height prediction")) +
  #labs(x = NULL, y = "fraction of fits, %", fill = NULL, title = bquote(bold(.(plotLetters[1]))~"base form height prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(heightDiameterFitStats10x10 %>% filter(responseVariable == "DBH", isBaseForm)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[2])~"base form DBH prediction")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[2]))~"base form DBH prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(heightDiameterFitStats10x10 %>% filter(responseVariable == "height", isBaseForm == FALSE)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = "model efficiency", y = "fraction of fits, %", fill = NULL, title = bquote(.(plotLetters[3])~"generalized height prediction")) +
  #labs(x = "model efficiency", y = "fraction of fits, %", fill = NULL, title = bquote(bold(.(plotLetters[3]))~"generalized height prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(heightDiameterFitStats10x10 %>% filter(responseVariable == "DBH", isBaseForm == FALSE)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = "model efficiency", y = NULL, fill = NULL, title = bquote(.(plotLetters[4])~"generalized DBH prediction")) +
  #labs(x = "model efficiency", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[4]))~"generalized DBH prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 2, ncol = 2, guides = "collect") &
  guides(color = "none", fill = guide_legend(ncol = 7)) &
  scale_color_manual(breaks = levels(heightDiameterFitStats10x10$species), limits = levels(heightDiameterFitStats10x10$species), values = speciesGroupColors) &
  scale_fill_manual(breaks = levels(heightDiameterFitStats10x10$species), limits = levels(heightDiameterFitStats10x10$species), values = speciesGroupColors) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not distinct"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not distinct"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not distinct"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.6, "line"), legend.position = "bottom", legend.text = element_text(margin = margin(l = 1, r = 6.5)))
#ggsave("trees/height-diameter/figures/Figure 05 model efficiency.png", height = 10, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 05 model efficiency.tif", height = 10, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 05 model efficiency.pdf", height = 10, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 6: fitting success, model significance, and predictor variable use
modelFitStats = heightDiameterFitStats10x10 %>% group_by(responseVariable, species, name, fitting) %>%
  summarize(modelFit = is.na(distinct[1]) == FALSE, .groups = "drop_last") %>% 
  summarize(fitPct = 100 * mean(modelFit), .groups = "drop") %>%
  mutate(species = factor(species, levels = levels(predictorVariableResults$species)))
  
ggplot(modelFitStats %>% filter(responseVariable == "height")) +
  geom_col(aes(x = fitPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[1])~"height fits")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[1]))~"height fits")) +
  scale_y_discrete(limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableStats10x10 %>% filter(responseVariable == "height")) +
  geom_col(aes(x = distinctPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[2])~"height forms")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[2]))~"height forms")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(heightDiameterFitStats10x10 %>% filter(responseVariable == "height", distinct)) +
  geom_violin(aes(x = meanAbsolutePercentPlantationEffect, y = species, color = species, group = species), quantiles = c(0.25, 0.5, 0.75), quantile.linetype = "solid", na.rm = TRUE, width = 1.2) +
  coord_cartesian(xlim = c(0, 125)) +
  labs(x = NULL, y = NULL, color = NULL, title = bquote(.(plotLetters[3])~"height")) +
  #labs(x = NULL, y = NULL, color = NULL, title = bquote(bold(.(plotLetters[3]))~"height")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "height", hasBasalArea)) +
  geom_count(aes(x = nBasalAreaDistinct, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.2, 2.2)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[4])~"BA+L")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[4]))~"BA+L")) +
  scale_x_continuous(breaks = seq(0, 2), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "height", hasPhysio)) +
  geom_count(aes(x = nPhysioDistinct, y = species, color = species), shape = 18) + 
  coord_cartesian(xlim = c(-0.3, 5.3)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[5])~"height physiographic")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[5]))~"height physiographic")) +
  scale_x_continuous(breaks = seq(0, 5), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "height", hasRelDbhOrRelHt)) +
  geom_count(aes(x = distinctRelHtOrDbh, y = species, color = species)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[6])~"RelDbh")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[6]))~"RelDbh")) +
  coord_cartesian(xlim = c(-0.3, 1.3)) +
  scale_x_continuous(breaks = seq(0, 1), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(modelFitStats %>% filter(responseVariable == "DBH")) +
  geom_col(aes(x = fitPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "model forms\nfit, %", y = NULL, fill = NULL, title = bquote(.(plotLetters[7])~"DBH fits")) +
  #labs(x = "model forms\nfit, %", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[7]))~"DBH fits")) +
  scale_y_discrete(limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableStats10x10 %>% filter(responseVariable == "DBH")) +
  geom_col(aes(x = distinctPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "distinct\nforms, %", y = NULL, fill = NULL, title = bquote(.(plotLetters[8])~"DBH forms")) +
  #labs(x = "distinct\nforms, %", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[8]))~"DBH forms")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(heightDiameterFitStats10x10 %>% filter(responseVariable == "DBH", distinct)) +
  geom_violin(aes(x = meanAbsolutePercentPlantationEffect, y = species, color = species, group = species), quantiles = c(0.25, 0.5, 0.75), quantile.linetype = "solid", na.rm = TRUE, width = 1.2) +
  coord_cartesian(xlim = c(0, 125)) +
  labs(x = "mean absolute\nplantation effect, %", y = NULL, color = NULL, title = bquote(.(plotLetters[9])~"DBH")) +
  #labs(x = "mean absolute\nplantation effect, %", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[9]))~"DBH")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "DBH", hasBasalArea)) +
  geom_count(aes(x = nBasalAreaDistinct, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.2, 2.2)) +
  labs(x = "distinct\nvariables", y = NULL, fill = NULL, title = bquote(.(plotLetters[10])~"ABA+T")) +
  #labs(x = "distinct\nvariables", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[10]))~"ABA+T")) +
  scale_x_continuous(breaks = seq(0, 2), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "DBH", hasPhysio)) +
  geom_count(aes(x = nPhysioDistinct, y = species, color = species), shape = 18) + 
  coord_cartesian(xlim = c(-0.3, 5.3)) +
  labs(x = "distinct\nvariables", y = NULL, fill = NULL, title = bquote(.(plotLetters[11])~"DBH physiographic")) +
  #labs(x = "distinct\nvariables", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[11]))~"DBH physiographic")) +
  scale_x_continuous(breaks = seq(0, 5), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "DBH", hasRelDbhOrRelHt)) +
  geom_count(aes(x = distinctRelHtOrDbh, y = species, color = species)) + 
  labs(x = "distinct\nvariables", y = NULL, fill = NULL, title = bquote(.(plotLetters[12])~"RelHt")) +
  #labs(x = "distinct\nvariables", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[12]))~"RelHt")) +
  coord_cartesian(xlim = c(-0.3, 1.3)) +
  scale_x_continuous(breaks = c(0, 1), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 6, widths = c(1.7, 1.7, 2.2, 1.7, 2.7, 1.2), guides = "collect") &
  scale_color_manual(breaks = levels(predictorVariableStats10x10$species), limits = levels(predictorVariableResults$species), values = speciesGroupColors) &
  scale_fill_manual(breaks = levels(predictorVariableStats10x10$species), limits = levels(predictorVariableResults), values = speciesGroupColors) &
  scale_size_area(max_size = 4.5) & # 4.5 for Cairo .pdf, 5.0 for .tif
  theme(legend.position = "none")
#ggsave("trees/height-diameter/figures/Figure 06 predictor variables.png", height = 10, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 06 predictor variables.tif", height = 10, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 06 predictor variables.pdf", height = 10, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 7: Douglas-fir and grand fir preferred models
# preferred forms from 10x10 cross validation
#print(heightDiameterModelsPreferred10x10 %>% filter(species %in% c("Douglas-fir", "red alder")) %>% select(-outOfRangeName, -aucOutOfRange, -nseName, -aucNse) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "ACMA3", "ABGR", "other"), levels = c("Douglas-fir", "bigleaf maple", "grand fir", "other species")), maeName = str_trunc(maeName, 28, ellipsis = ""), rmseName = str_trunc(rmseName, 28, ellipsis = ""), aicName = str_trunc(aicName, 28, ellipsis = "")), n = 33)
print(heightDiameterModelRanking10x10 %>% filter(distinct, isBaseForm) %>% select(-aucBlended, -starts_with("has"), -speciesFraction) %>%
        pivot_longer(starts_with("auc"), names_to = "statistic", values_to = "auc") %>% mutate(statistic = factor(statistic, levels = c("aucOutOfRange", "aucMae", "aucRmse", "aucDeltaAicN", "aucNse"))) %>%
        group_by(species, responseVariable, statistic) %>%
        slice_max(auc, n = 1) %>% arrange(species, desc(responseVariable), statistic), n = 70)

if (exists("psmeHeightFromDiameterPreferred") == FALSE) 
{ 
  psmeHeightFromDiameterPreferred = readRDS("trees/height-diameter/data/PSME preferred height models.Rds")
}
if (exists("psmeDiameterFromHeightPreferred") == FALSE) 
{ 
  psmeDiameterFromHeightPreferred = readRDS("trees/height-diameter/data/PSME preferred diameter models.Rds")
}
if (exists("abgrHeightFromDiameterPreferred") == FALSE)
{ 
  abgrHeightFromDiameterPreferred = readRDS("trees/height-diameter/data/ABGR preferred height models.Rds")
}
if (exists("abgrDiameterFromHeightPreferred") == FALSE)
{ 
  abgrDiameterFromHeightPreferred = readRDS("trees/height-diameter/data/ABGR preferred diameter models.Rds")
}

# Temesgen et al. 2007 height = 1.3 + exp(b1 - b2 * DBH^b3) => b1 - b2 * DBH^b3 = ln(height - 1.3) => DBH^b3 = 1/b2 * (b1 - ln(height - 1.3))
#                      DBH = (1/b2 * (b1 - ln(height - 1.3)))^(1/b3)
psmeReference = bind_rows(bind_rows(psmeHeightFromDiameterPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    psmeHeightFromDiameterPreferred$weibull$stats %>% mutate(model = "base form 2"),
                                    psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Temesgen et al. 2007", "height", psme2016, 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), 4, psme2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(psmeDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    psmeDiameterFromHeightPreferred$sibbesenReplace$stats %>% mutate(model = "base form 2"),
                                    psmeDiameterFromHeightPreferred$gamAbatPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Temesgen et al. 2007", "DBH", psme2016, (1/6.7792 * (5.7567 - log(psme2016$TotalHt - 1.3)))^(1/-0.2795), 4, psme2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "PSME", adaptiveWeightFraction = 0)
psmeReferenceDbh = seq(0, ceiling(max(psme2016$DBH))) # Temesgen et al. 2007 dataset limit of 190 cm, extended of necessity

# Larsen et al. 1987 https://ir.library.oregonstate.edu/concern/technical_reports/h702q764 (but see also Tinkham et al. 2016 http://dx.doi.org/10.1080/07038992.2016.1232587)
#   H = 1.37 + 0.3048 * exp(6.74974 - 5.49823 * (dbh/2.54)^-0.327093) (eq. 3)
#     = 1.37 + 0.3048 * exp(6.08798 - 5.02407 * (dbh/2.54)^-0.375209 + 0.00100968 * 3.28084^2/2.47105 * standBasalArea) (eq 4, m²/ha = 3.28084^2/2.47 ft²/ac)
#   H = a0 + 0.3048 * exp(b0 + b1 dbh^b2) => exp(b0 + b1 dbh^b2) = 3.28084 * (H - a0) => dbh^b2 = 1/b1 (log(3.28084 * (H - a0)) - b0)
# dbh = (1/b1 (log(3.28084 * (H - a0)) - b0))^(1/b2)
abgrReference = bind_rows(bind_rows(abgrHeightFromDiameterPreferred$hossfeld$stats %>% mutate(model = "base form 1"),
                                    abgrHeightFromDiameterPreferred$chapmanRichards$stats %>% mutate(model = "base form 2"),
                                    abgrHeightFromDiameterPreferred$sharmaPartonPhysio$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Larsen et al. 1999", "height", abgr2016, 1.37 + 0.3048 * exp(6.74974 - 5.49823 * (abgr2016$dbh / 2.54)^-0.327093), 4, abgr2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(abgrDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    abgrDiameterFromHeightPreferred$parabolic$stats %>% mutate(model = "base form 2"),
                                    abgrDiameterFromHeightPreferred$gamPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Larsen et al. 1999", "DBH", abgr2016, (-1/5.49823 * (log(3.28084 * (abgr2016$height - 1.37)) + 5.49823))^(1/-0.327093), 4, abgr2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "THPL", adaptiveWeightFraction = 0)
abgrReferenceDbh = seq(0, 213) # Larsen et al. 1987


ggplot() +
  geom_point(aes(x = abgr2016$DBH, y = abgr2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = abgr2016$DBH, y = predict(abgrHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = abgr2016$isPlantation, linetype = abgr2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = abgr2016$DBH, y = predict(abgrHeightFromDiameterPreferred$hossfeld), color = "Hossfeld IV (MAE, RMSE, ME)", group = abgr2016$isPlantation, linetype = abgr2016$isPlantation), alpha = 0.8) + # Hossfeld and Michaelis-Menten fits are visually indistinguishable
  geom_line(aes(x = abgr2016$DBH, y = predict(abgrHeightFromDiameterPreferred$chapmanRichards), color = "Chapman-Richards (ΔAICn)", group = abgr2016$isPlantation, linetype = abgr2016$isPlantation), alpha = 0.8) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = abgrReferenceDbh, y = 1.37 + 0.3048 * exp(7.232880669 - 5.746899904 * (0.393701 * abgrReferenceDbh)^-0.271564741), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(.(plotLetters[1])~"grand fir height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[1]))~"grand fir height")) +
  scale_color_manual(breaks = c("Hossfeld IV (MAE, RMSE, ME)", "Chapman-Richards (ΔAICn)", "Sharma-Parton physio", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = abgr2016$DBH, y = abgr2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(abgrDiameterFromHeightPreferred$gamPhysio), y = abgr2016$TotalHt, color = "REML GAM physio", group = abgr2016$isPlantation, linetype = abgr2016$isPlantation)) +
  geom_line(aes(x = predict(abgrDiameterFromHeightPreferred$gam), y = abgr2016$TotalHt, color = "REML GAM (all but MAB)", group = abgr2016$isPlantation, linetype = abgr2016$isPlantation)) +
  geom_line(aes(x = predict(abgrDiameterFromHeightPreferred$parabolic), y = abgr2016$TotalHt, color = "parabolic (MAB)", group = abgr2016$isPlantation, linetype = abgr2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = abgrReferenceDbh, y = 1.37 + 0.3048 * exp(7.232880669 - 5.746899904 * (0.393701 * abgrReferenceDbh)^-0.271564741), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(.(plotLetters[2])~"grand fir DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[2]))~"grand fir DBH")) +
  scale_color_manual(breaks = c("REML GAM (all but MAB)", "parabolic (MAB)", "REML GAM physio", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh), color = "Sharma-Parton BA+L RelDbh physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$gam), color = "REML GAM (MAB, MAE, RMSE, ME)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$ratkowsky), color = "Ratkowsky (ΔAICn)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) + # or Prodan
  geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL, title = bquote(.(plotLetters[1])~"Douglas-fir height")) +
  #labs(x = NULL, y = "height, m", color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[1]))~"Douglas-fir height")) +
  scale_color_manual(breaks = c("REML GAM (MAB, MAE, RMSE, ME)", "Ratkowsky (ΔAICn)", "Sharma-Parton BA+L RelDbh physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey60")) +
  theme(legend.key = element_rect(fill = alpha("white", 0.5)), legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gamAbatPhysio), y = psme2016physio$TotalHt, color = "REML GAM ABA+T physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gam), y = psme2016$TotalHt, color = "REML GAM (MAB, MAE, RMSE, ME)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$sibbesenReplace), y = psme2016$TotalHt, color = "Sibbesen replace (ΔAICn)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = bquote(.(plotLetters[2])~"Douglas-fir DBH")) +
  #labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[2]))~"Douglas-fir DBH")) +
  scale_color_manual(breaks = c("REML GAM (MAB, MAE, RMSE, ME)", "Sibbesen replace (ΔAICn)", "REML GAM ABA+T physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "previous model"), labels = c("natural regeneration", "plantation", "previous model"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 07 PSME-ABGR fits.png", height = 12, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/height-diameter/figures/Figure 07 PSME-ABGR fits.tif", height = 12, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 07 PSME-ABGR fits 1000.pdf", height = 12, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 8: bigleaf maple and other species preferred models
#print(heightDiameterModelsPreferred10x10 %>% filter(speciesModel %in% c("grand fir", "other")) %>% select(-outOfRangeName, -aucOutOfRange, -nseName, -aucNse) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "ACMA3", "ABGR", "other"), levels = c("Douglas-fir", "bigleaf maple", "grand fir", "other species")), maeName = str_trunc(maeName, 28, ellipsis = ""), rmseName = str_trunc(rmseName, 28, ellipsis = ""), aicName = str_trunc(aicName, 28, ellipsis = "")), n = 32)
if (exists("acmaHeightFromDiameterPreferred") == FALSE)
{ 
  acmaHeightFromDiameterPreferred = readRDS("trees/height-diameter/data/ACMA3 preferred height models.Rds")
}
if (exists("acmaDiameterFromHeightPreferred") == FALSE)
{ 
  acmaDiameterFromHeightPreferred = readRDS("trees/height-diameter/data/ACMA3 preferred diameter models.Rds")
}
if (exists("otherHeightFromDiameterPreferred") == FALSE) 
{ 
  otherHeightFromDiameterPreferred = readRDS("trees/height-diameter/data/other preferred height models.Rds") 
}
if (exists("otherDiameterFromHeightPreferred") == FALSE) 
{ 
  otherDiameterFromHeightPreferred = readRDS("trees/height-diameter/data/other preferred diameter models.Rds")
}

acmaReference = bind_rows(bind_rows(acmaHeightFromDiameterPreferred$chapmanRichards$stats %>% mutate(model = "base form 1"),
                                    acmaHeightFromDiameterPreferred$michaelisMenten$stats %>% mutate(model = "base form 2"),
                                    acmaHeightFromDiameterPreferred$weibullBal$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Wang and Hann 1988", "height", acma2016, 1.37 + 0.3048 * exp(5.21462 - 2.70252 * (0.393701 * acma2016$DBH)^-0.354756), 4, acma2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(acmaDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    acmaDiameterFromHeightPreferred$ruark$stats %>% mutate(model = "base form 2"),
                                    acmaDiameterFromHeightPreferred$ruarkAbatPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Wang and Hann 1988", "DBH", acma2016, 1/0.393701 * (1/2.70252 * (5.21462 - log(1/1.37 * (acma2016$TotalHt - 1.37))))^(1/-0.354756), 4, acma2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "ACMA3", adaptiveWeightFraction = 0)
acmaReferenceDbh = seq(0, 182) # Wang and Hann 1988 dataset limit of 130 cm, extended of necessity plus for plot visibility

# no prior height model for other species, so no reference or reference DBH
otherReference = bind_rows(bind_rows(otherHeightFromDiameterPreferred$gam$stats %>% mutate(model = "base form 1"),
                                     otherHeightFromDiameterPreferred$linear$stats %>% mutate(model = "base form 2"),
                                     otherHeightFromDiameterPreferred$gamRelDbhPhysio$stats %>% mutate(model = "generalized height")) %>% 
                             mutate(responseVariable = "height"),
                           bind_rows(otherDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                     otherDiameterFromHeightPreferred$sibbesenReplace$stats %>% mutate(model = "base form 2"),
                                     otherDiameterFromHeightPreferred$gamRelHt$stats %>% mutate(model = "generalized DBH")) %>% 
                             mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "other", adaptiveWeightFraction = 0)

ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$weibullBal), color = "Weibull BA+L", group = acma2016$isPlantation, linetype = acma2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$michaelisMenten), color = "Michaelis-Menten (ΔAICn)", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$ratkowsky), color = "Ratkowsky (RMSE, ME)", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = acmaReferenceDbh, y = 1.37 + 0.3048 * exp(5.21462 - 2.70252 * (0.393701 * acmaReferenceDbh)^-0.354756), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(.(plotLetters[3])~"bigleaf maple height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[3]))~"bigleaf maple height")) +
  scale_color_manual(breaks = c("Ratkowsky (RMSE, ME)", "Michaelis-Menten (ΔAICn)", "Weibull BA+L", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  #geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$gamAbatPhysioRelHt), y = acma2016$TotalHt, color = "REML GAM ABA+T RelHt physio", group = acma2016$isPlantation, linetype = acma2016$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$ruarkAbatPhysio), y = acma2016$TotalHt, color = "Ruark ABA+T physio", group = acma2016$isPlantation, linetype = acma2016$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$gam), y = acma2016$TotalHt, color = "REML GAM (all five statistics)", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$ruark), y = acma2016$TotalHt, color = "Ruark", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = acmaReferenceDbh, y = 1.37 + 0.3048 * exp(5.21462 - 2.70252 * (0.393701 * acmaReferenceDbh)^-0.354756), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(.(plotLetters[4])~"bigleaf maple DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[4]))~"bigleaf maple DBH")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "Ruark", "Ruark ABA+T physio", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = other2016physio$DBH, y = predict(otherHeightFromDiameterPreferred$gamRelDbhPhysio), color = "REML GAM RelDbh physio", group = other2016physio$isPlantation, linetype = other2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$gam), color = "REML GAM (MAE, RMSE, ME)", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$linear), color = "linear (ΔAICn)", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(.(plotLetters[3])~"other species height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[3]))~"other species height")) +
  scale_color_manual(breaks = c("REML GAM (MAE, RMSE, ME)", "linear (ΔAICn)", "REML GAM RelDbh physio", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$gamRelHt), y = other2016$TotalHt, color = "REML GAM RelHt", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$gam), y = other2016$TotalHt, color = "REML GAM (all five statistics)", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$sibbesenReplace), y = other2016$TotalHt, color = "Sibbesen replace", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  #geom_line(aes(x = predict(otherDiameterFromHeightPreferred$parabolic), y = other2016$TotalHt, color = "parabolic", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(.(plotLetters[4])~"other species DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[4]))~"other species DBH")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "Sibbesen replace", "REML GAM RelHt", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "previous model"), labels = c("natural regeneration", "plantation", "previous model"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 08 ACMA3-other fits.png", height = 12, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 08 ACMA3-other fits.tif", height = 12, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 08 ACMA3-other fits.pdf", height = 12, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 9: comparison of goodness of fit statistics between reference and revised models
selectedModels = bind_rows(psmeReference, abgrReference, acmaReference, otherReference) %>% 
  select(-coefficients) %>%
  mutate(deltaAicN = aic/n - min(aic/n, na.rm = TRUE),
         model = as.factor(model),
         species = factor(species, labels = c("Douglas-fir", "bigleaf maple", "grand fir", "other species"), levels = c("PSME", "ACMA3", "ABGR", "other")))
#writexl::write_xlsx(list(preferredModelStats = selectedModels %>% select(species, responseVariable, name, pctOutOfRange, mae, mape, rmse, rmpse, deltaAicN, nse), "trees/height-diameter/figures/selected model stats.xlsx")
# selectedModels %>% filter(species == "bigleaf maple") %>% select(responseVariable, name, pctOutOfRange, mae, rmse, deltaAicN, nse)

ggplot() +
  geom_point(aes(x = pctOutOfRange, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "out of range, %", y = NULL, color = NULL, shape = NULL, size = NULL, title = paste0(plotLetters[1], " height prediction")) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(limits = rev) +
ggplot() +
  geom_point(aes(x = mae, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "MAE, m", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = rmse, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "RMSE, m", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = deltaAicN, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "ΔAICn", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 6)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = nse, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0.1, 1)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "model efficiency", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1.0), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = pctOutOfRange, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "out of range, %", y = NULL, color = NULL, shape = NULL, size = NULL, title = paste0(plotLetters[2], " DBH prediction")) +
  scale_y_discrete(limits = rev) +
ggplot() +
  geom_point(aes(x = mae, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "MAE, cm", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = rmse, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "RMSE, cm", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = deltaAicN, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "ΔAICn", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 6)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = nse, y = species, color = model, shape = model, size = model), selectedModels %>% filter((responseVariable == "DBH") | ((species == "bigleaf maple") & (model == "generalized height"))) %>% mutate(nse = if_else(model != "generalized height", nse, -1)), alpha = 0.6) + # flow fake height data point off plot to get ggplot2 3.5.1 to draw the generalized height legend key, for some reason inserting a fake row with add_row() causes ggplot to draw additional red alder points
  annotate("text", x = 0.06, y = "bigleaf maple", label = "🡄 −1.47", alpha = 0.8, color = "grey60", hjust = 0, size = 2.0, vjust = 0.4) +
  coord_cartesian(xlim = c(0.1, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 0.8))) +
  labs(x = "model efficiency", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1.0), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL, limits = rev) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, guides = "collect") &
  scale_color_manual(breaks = c("base form 1", "base form 2", "generalized height", "generalized DBH", "previous model"), values = c("dodgerblue2", "red2", "green2", "cyan", "grey60"), drop = FALSE) &
  scale_shape_manual(breaks = c("base form 1", "base form 2", "generalized height", "generalized DBH", "previous model"), values = c(15, 18, 16, 16, 17), drop = FALSE) &
  scale_size_manual(breaks = c("base form 1", "base form 2", "generalized height", "generalized DBH", "previous model"), values = c(1.9, 2.5, 2.2, 2.2, 1.7), drop = FALSE) &
  theme(plot.title.position = "plot")
#ggsave("trees/height-diameter/figures/Figure 09 model accuracy.png", height = 9, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 09 model accuracy.tif", height = 9, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure 09 model accuracy.pdf", height = 9, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure S1: comparison of accuracy metrics
accuracyCorrelation = bind_rows(as.data.frame(cor(heightDiameterFitStats10x10 %>% filter(responseVariable == "height") %>% 
                                                    select(pctOutOfRange, mape, rmse, deltaAicN, nse), use = "pairwise.complete.obs")) %>%
                                  tibble::rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "height"),
                                as.data.frame(cor(heightDiameterFitStats10x10 %>% filter(responseVariable == "DBH") %>% 
                                                    select(pctOutOfRange, mape, rmse, deltaAicN, nse), use = "pairwise.complete.obs")) %>%
                                  tibble::rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "DBH")) %>%
  mutate(metricX = factor(metricX, levels = c("pctOutOfRange", "mape", "rmse", "deltaAicN", "nse"), labels = c("MAB", "MAE", "RMSE", "ΔAICn", "model\nefficiency")),
         metricY = factor(metricY, levels = c("pctOutOfRange", "mape", "rmse", "deltaAicN", "nse"), labels = c("MAB", "MAE", "RMSE", "ΔAICn", "model efficiency")))

ggplot(accuracyCorrelation %>% filter(responseVariable == "height")) + 
  coord_equal() +
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  labs(x = NULL, y = NULL, fill = "correlation", title = bquote(.(plotLetters[1])~"height prediction")) +
  #labs(x = NULL, y = NULL, fill = "correlation", title = bquote(bold(.(plotLetters[1]))~"height prediction")) +
  scale_y_discrete(limits = rev) +
ggplot(accuracyCorrelation %>% filter(responseVariable == "DBH")) + 
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "correlation", title = bquote(.(plotLetters[2])~"DBH prediction")) +
  #labs(x = NULL, y = NULL, fill = "correlation", title = bquote(bold(.(plotLetters[2]))~"DBH prediction")) +
  scale_y_discrete(labels = NULL, limits = rev) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  scico::scale_fill_scico(palette = "vik", limits = c(-1, 1)) &
  theme(legend.spacing.y = unit(0.5, "line"))
#ggsave("trees/height-diameter/figures/Figure S01 goodness of fit correlation.png", height = 8.5, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/height-diameter/figures/Figure S01 goodness of fit correlation.tif", height = 8.5, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure S01 goodness of fit correlation.pdf", height = 8.5, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


# Figure S2: height prediction accuracy
# If axis limits, vertical dodges, or aesthetics are changed Figure S6 should be updated.
heightFromDiameterAccuracyLevels = heightDiameterFitStats10x10 %>% 
  filter(responseVariable == "height") %>%
  group_by(name, species) %>%
  summarize(speciesFraction = speciesFraction[1],
            meanPenalizedMae = mean(if_else(is.na(mae), 100, mae)),
            meanPenalizedMape = mean(if_else(is.na(mape), 100, mape)),
            meanPenalizedRmspe = mean(if_else(is.na(rmspe), 100, rmse)),
            .groups = "drop_last") %>% # change to grouping only by name
  summarize(weightedMae = sum(speciesFraction * meanPenalizedMae),
            .groups = "drop") %>%
  arrange(weightedMae)
heightFromDiameterResults = heightDiameterFitStats10x10 %>% 
  filter(responseVariable == "height", str_detect(name, "GNLS") == FALSE, str_detect(name, "RelHt") == FALSE) %>%
  mutate(name = factor(name, levels = heightFromDiameterAccuracyLevels$name)) %>%
  group_by(species, name) %>%
  reframe(quantiles = c(0.025, 0.5, 0.975),
          pctOutOfRange = quantile(pctOutOfRange, probs = quantiles, na.rm = TRUE), 
          mape = quantile(mape, probs = quantiles, na.rm = TRUE), 
          rmspe = quantile(rmspe, probs = quantiles, na.rm = TRUE), 
          deltaAicN = quantile(deltaAicN, probs = quantiles, na.rm = TRUE), 
          nse = quantile(nse, probs = quantiles, na.rm = TRUE),
          distinct = distinct[1]) %>%
  pivot_wider(names_from = quantiles, names_sep = "_q", values_from = c("pctOutOfRange", "mape", "rmspe", "deltaAicN", "nse")) %>%
  mutate(verticalDodge = recode(species, "Douglas-fir" = 0.15, "bigleaf maple" = 0.05, "grand fir" = -0.05, "other species" = -0.15))

ggplot(heightFromDiameterResults) +
  geom_errorbar(aes(xmin = pctOutOfRange_q0.025, xmax = pctOutOfRange_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = pctOutOfRange_q0.5, y = name, color = species, alpha = distinct, size = distinct, shape = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "MAB, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbar(aes(xmin = mape_q0.025, xmax = mape_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = mape_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbar(aes(xmin = rmspe_q0.025, xmax = rmspe_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = rmspe_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "RMSE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  #scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  #geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM", yend = "REML GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM BA+L", yend = "REML GAM BA+L", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM physio", yend = "REML GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(xmin = deltaAicN_q0.025, xmax = deltaAicN_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = deltaAicN_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  labs(x = "ΔAICn", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  coord_cartesian(xlim = c(0, 10)) + # exclude high AIC of linear, parabolic, and Douglas-fir power+Curtis fits to avoid squashing of nonlinear differences
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbar(aes(xmin = nse_q0.025, xmax = nse_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = nse_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0.1, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 0, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 1), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_alpha_manual(breaks = c("fixed weights", "not distinct"), labels = c("form distinct", "not distinct"), values = c(0.75, 0.3)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = speciesGroupColors) &
  scale_shape_manual(breaks = c("fixed weights", "not distinct"), labels = c("form distinct", "not distinct"), values = c(16, 3)) &
  scale_size_manual(breaks = c("fixed weights", "not distinct"), labels = c("form distinct", "not distinct"), values = c(1.9, 1.4)) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure S02 height accuracy.png", height = 16, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/height-diameter/figures/Figure S02 height accuracy.tif", height = 22, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure S02 height accuracy.pdf", height = 22, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


# Figure S3: DBH prediction accuracy
# If axis limits, vertical dodges, or aesthetics are changed Figure S4 should be updated.
diameterFromHeightAccuracyLevels = heightDiameterFitStats10x10 %>% 
  filter(responseVariable == "DBH") %>%
  group_by(name, species) %>%
  summarize(speciesFraction = speciesFraction[1],
            meanPenalizedMae = mean(if_else(is.na(mae), 100, mae)),
            meanPenalizedMape = mean(if_else(is.na(mape), 100, mape)),
            meanPenalizedRmspe = mean(if_else(is.na(rmspe), 100, rmspe)),
            .groups = "drop_last") %>% # change to grouping only by name
  summarize(weightedRmse = sum(speciesFraction * meanPenalizedRmspe), 
            .groups = "drop") %>%
  arrange(weightedRmse)
diameterFromHeightResults = heightDiameterFitStats10x10 %>% 
  filter(responseVariable == "DBH", str_detect(name, "BA\\+L") == FALSE, str_detect(name, "GNLS") == FALSE, name != "REML GAM ABA+T physio") %>%
  mutate(name = factor(name, levels = diameterFromHeightAccuracyLevels$name)) %>%
  group_by(species, name) %>%
  reframe(quantiles = c(0.025, 0.5, 0.975),
          pctOutOfRange = quantile(pctOutOfRange, probs = quantiles, na.rm = TRUE), 
          mape = quantile(mape, probs = quantiles, na.rm = TRUE), 
          rmspe = quantile(rmspe, probs = quantiles, na.rm = TRUE), 
          deltaAicN = quantile(deltaAicN, probs = quantiles, na.rm = TRUE), 
          nse = quantile(nse, probs = quantiles, na.rm = TRUE), 
          distinct = distinct[1]) %>%
  pivot_wider(names_from = quantiles, names_sep = "_q", values_from = c("pctOutOfRange", "mape", "rmspe", "deltaAicN", "nse")) %>%
  mutate(verticalDodge = recode(species, "Douglas-fir" = 0.15, "bigleaf maple" = 0.05, "grand fir" = -0.05, "other species" = -0.15))

ggplot(diameterFromHeightResults) +
  geom_errorbar(aes(xmin = pctOutOfRange_q0.025, xmax = pctOutOfRange_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = pctOutOfRange_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "MAB, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbar(aes(xmin = mape_q0.025, xmax = mape_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = mape_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbar(aes(xmin = rmspe_q0.025, xmax = rmspe_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = rmspe_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "RMSE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  #geom_segment(x = 0.255, xend = 0.275, y = "Schnute inverse", yend = "Schnute inverse", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM", yend = "REML GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM ABA+T", yend = "REML GAM ABA+T", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM physio", yend = "REML GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(xmin = deltaAicN_q0.025, xmax = deltaAicN_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = deltaAicN_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "ΔAICn", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), trans = scales::pseudo_log_trans()) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbar(aes(xmin = nse_q0.025, xmax = nse_q0.975, y = name, color = species, alpha = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), linewidth = 0.5, orientation = "y", width = 0.1) +
  geom_point(aes(x = nse_q0.5, y = name, color = species, alpha = distinct, shape = distinct, size = distinct), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0.1, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 0, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 1), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not distinct"), labels = c("NLME/nlrob", "form distinct", "not distinct"), values = c(0.75, 0.75, 0.3)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = speciesGroupColors) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not distinct"), labels = c("NLME/nlrob", "form distinct", "not distinct"), values = c(16, 18, 3)) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not distinct"), labels = c("NLME/nlrob", "form distinct", "not distinct"), values = c(1.5, 1.9, 1.4)) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure S03 DBH accuracy.png", height = 16, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/height-diameter/figures/Figure S03 DBH accuracy.tif", height = 22, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
#ggsave("trees/height-diameter/figures/Figure S03 DBH accuracy.pdf", height = 22, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure S4 in setup.R


## supplementary spreadsheet: model forms by species and use 
preferredModelForms = heightDiameterModelRanking10x10 %>% filter(distinct) %>% 
  group_by(responseVariable, species) %>% 
  arrange(desc(aucBlended)) %>% 
  mutate(rank = row_number()) %>% 
  arrange(desc(responseVariable), species, isBaseForm, rank) %>%
  select(-nameAndFit, -starts_with("has"), -speciesFraction) %>%
  relocate(responseVariable, species, rank, name, fitting, isBaseForm, distinct, aucOutOfRange, aucMae, aucRmse, aucDeltaAicN, aucNse)
preferredModelFormsHeightNonPhysio = heightDiameterModelRanking10x10 %>% filter(distinct, hasPhysio == FALSE) %>% 
  group_by(responseVariable, species) %>% 
  arrange(desc(aucBlended)) %>% 
  mutate(rank = row_number()) %>% 
  arrange(desc(responseVariable), species, isBaseForm, rank) %>%
  select(-nameAndFit, -starts_with("has"), -speciesFraction) %>%
  relocate(responseVariable, species, rank, name, fitting, isBaseForm, distinct, aucOutOfRange, aucMae, aucRmse, aucDeltaAicN, aucNse)
preferredModelFormsInitialDbh = heightDiameterModelRanking10x10 %>% filter(responseVariable == "DBH", distinct, hasStand == FALSE) %>% 
  group_by(responseVariable, species) %>% 
  arrange(desc(aucBlended)) %>% 
  mutate(rank = row_number()) %>% 
  arrange(desc(responseVariable), species, isBaseForm, rank) %>%
  select(-nameAndFit, -starts_with("has"), -speciesFraction) %>%
  relocate(responseVariable, species, rank, name, fitting, isBaseForm, distinct, aucOutOfRange, aucMae, aucRmse, aucDeltaAicN, aucNse)
#writexl::write_xlsx(list(preferred = preferredModelForms, heightNonPhysio = preferredModelFormsHeightNonPhysio, initialDbh = preferredModelFormsInitialDbh), "trees/height-diameter/figures/height-diameter model AUCs.xlsx")


## additional info
if (htDiaOptions$includeInvestigatory)
{
  # validation fold size consistency
  bind_rows(heightDiameterCoefficients2x50 %>% filter(distinct == 1, is.na(n) == FALSE) %>% group_by(species) %>%
              summarize(crossValidation = "2x50", minN = min(n), meanN = mean(n), maxN = max(n), sdN = sd(n)),
            heightDiameterCoefficients10x10 %>% filter(distinct == 1, is.na(n) == FALSE) %>% group_by(species) %>%
              summarize(crossValidation = "10x10", minN = min(n), meanN = mean(n), maxN = max(n), sdN = sd(n)),
            heightDiameterCoefficients10x25 %>% filter(distinct == 1, is.na(n) == FALSE) %>% group_by(species) %>%
              summarize(crossValidation = "10x25", minN = min(n), meanN = mean(n), maxN = max(n), sdN = sd(n))) %>%
    arrange(crossValidation, species)

  # nonphysical predictions
  outOfRange = heightDiameterCoefficients10x10 %>% filter(fitSet == "primary", distinct) %>% 
    group_by(responseVariable, species, name) %>%
    summarize(baseName = if_else(word(name[1]) %in% c("REML", "modified", "unified"), paste(word(name[1], 1), word(name[1], 2)), word(name[1])),
              isBaseForm = isBaseForm[1], pctOutOfRange = mean(pctOutOfRange), .groups = "drop")
  
  outOfRange %>% group_by(responseVariable, species) %>% 
    summarize(pctBaseNonPhysical = sum(isBaseForm * pctOutOfRange > 0) / sum(isBaseForm),
              pctGeneralizedNonPhysical = sum((isBaseForm == FALSE) * pctOutOfRange > 0) / sum(isBaseForm == FALSE))
  
  print(outOfRange %>% group_by(responseVariable) %>%
    slice_max(pctOutOfRange, n = 20) %>%
    arrange(desc(responseVariable), desc(pctOutOfRange)), n = 40)
  print(outOfRange %>% filter(pctOutOfRange > 0) %>% group_by(responseVariable) %>%
          slice_min(pctOutOfRange, n = 20) %>%
          arrange(desc(responseVariable), desc(pctOutOfRange)), n = 40)
  print(outOfRange %>% filter(isBaseForm, pctOutOfRange == 0), n = 250)
  print(outOfRange %>% filter(baseName == "Ruark"), n = 30)
  
  outOfRange$pctOutOfRange
  ggplot() +
    geom_histogram(aes(x = pctOutOfRange, fill = responseVariable, group = responseVariable), outOfRange, binwidth = 1) +
    coord_trans(y = scales::pseudo_log_trans()) +
    labs(x = "fraction of fits with non-physical trees, %", y = "model fits", fill = NULL) +
    scale_y_continuous(breaks = c(0, 1, 10, 100, 1000)) +
  ggplot() +
    stat_ecdf(aes(x = pctOutOfRange, color = responseVariable, group = responseVariable), outOfRange) +
    labs(x = "non-physical fraction of trees, %", y = "fraction of fits", color = NULL) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    guides(fill = "none")
}