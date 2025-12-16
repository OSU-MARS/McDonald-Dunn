# load libraries, functions, and other2022 from setup.R

otherOptions = tibble(fitHeight = TRUE, 
                      fitHeightMixed = TRUE,
                      fitDbh = TRUE,
                      fitDbhMixed = TRUE,
                      recalcPreferredModels = FALSE)

## other species height regressions
if (otherOptions$fitHeight)
{
  otherHeightFromDiameter = list(linear = fit_lm("linear", height ~ dbh, other2022)) # a1p not significant in 50% of fits
  otherHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ dbh + I(dbh^2), other2022) # a1p, a2p not significant in 50% of fits

  otherHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, b1 = -0.02, b2 = 0.84)) # a1p, b1p, b2p (25%) not significant
  otherHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a1bal = 0, b1 = -0.02, b2 = 0.84), significant = FALSE) # a1p, a1ba, a2p, a1bal, a1balp, b1p, b2p, exponent BAL not significant
  otherHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1tw * topographicWetnessFD8f) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a1bal = 0, a1tw = 0, b1 = -0.02, b2 = 0.8), significant = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  otherHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 80, a1bal = 0.3, a1rd = 0, a1e = 0, b1 = -0.01, b2 = 0.8), significant = FALSE) # by propagation
  otherHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a1bal = 0, a1rd = 0, b1 = -0.02, b2 = 0.84), significant = FALSE) # singular gradient with a1bal, a1rd, a1rdp not significant, exponents NaN-Inf or singular gradient
  otherHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp(b1 * dbh))^(b2 + b2ac * cos(pi/180 * aspect)), other2022, start = list(a1 = 30, b1 = -0.02, b2 = 0.81, b2ac = -0.08), significant = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant, though b2as + b2ac intermittently test significant
  otherHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a1rd = 0, b1 = -0.02, b2 = 0.8), significant = FALSE) # a1rd, a1rdp, b1rd, b2rd not significant, signs of a1-b1 evaporation
  otherHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a1rd = 0, a1e = 0, b1 = -0.02, b2 = 0.81), significant = FALSE) # by propagation
  otherHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1*dbh / (1 + dbh)^b1, other2022, start = list(a1 = 2, b1 = 0.4))
  otherHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), other2022, start = list(a1 = 44, a2 = 33, b1 = -0.9)) # prone to singular gradient
  otherHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1 * dbh^b2), other2022, start = list(a1 = 300, b1 = -5, b2 = -0.1)) # a1p, b1p, b2p not significant, form not well posed (should probably be a1 * (exp() - 1))
  otherHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), other2022, start = list(a1 = 50, a2 = 40, b1 = 0.9)) # a1p, a2p, b1p not significant
  otherHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3), other2022, start = list(a1 = -0.01, a2 = 2.5, a3 = -5)) # a1p, a2p, a3p not significant
  otherHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, other2022, start = list(a1 = 1.5, b1 = 0.7)) # a1p, b1p not significant
  otherHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2)), other2022, start = list(a1 = 30, b1 = -11, b2 = 3)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$richardsW = tryCatch({ fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), other2022, start = list(Ha = 28, d = 0.2, kU = 0.025)) }, # always NaN-inf @ 2x50
                                               error = function(e) { return(create_fit_statistics("unified Richards", fittingMethod = "gsl_nls")) })
  otherHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^b4, other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant
  otherHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant
  otherHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0.1)) # { a1, b1, b2, b3, b4 } x { e, s, ac, tr, tw }, { a1, b1, b2, b3 } x as not significant
  otherHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, start = list(a1 = 9, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.06, b4 = 1.0, b4as = -0.1)) # a1rd, b1rd, b4rd not significant, b2rd often significant
  otherHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, start = list(a1 = 10, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.07, b4 = 1.0)) # a1rd, b1rd, b2rd, b4rd not significant
  otherHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a1tw * topographicWetnessFD8f)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare))^(b3)*dbh))^(b4), other2022, start = list(a1 = 5, a1tw = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0), significant = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw } not significant
  otherHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + (a1)*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, start = list(a1 = 8, b1 = 0.6, b2 = -0.01, b3 = 0.4, b3rd = -0.05, b4 = 1.0), significant = FALSE) # a1rd not significant
  otherHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1e * elevation)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, start = list(a1 = 10, a1rd = 0, a1e = 0, b1 = 0.5, b2 = -0.007, b3 = 0.4, b4 = 1.1), significant = FALSE) # by propagation
  otherHeightFromDiameter$sharmaZhang = tryCatch({ fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, other2022, start = list(a1 = 25, b1 = 0.05, b2 = -0.01, b3 = 0.2, b4 = 0.9)) }, # a1p, b1, b1p, b2p, b3, b3p, b4p not significant, NaN-inf @ 2x50
                                                 error = function(e) { return(create_fit_statistics("Sharma-Zhang", fittingMethod = "gsl_nls")) })
  otherHeightFromDiameter$sharmaZhangBal = tryCatch({ fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2bal * basalAreaLarger)*standTreesPerHectare^b3*dbh))^b4, other2022, start = list(a1 = 25, b1 = 0.05, b2 = -0.01, b2bal = 0, b3 = 0.2, b4 = 0.9), significant = FALSE) }, # a1ba, a2p, a1bal, b1bal, b2bal, b3bal, b4bal not significant, NaN-inf @ 2x50
                                                    error = function(e) { return(create_fit_statistics("Sharma-Zhang BA+L", fittingMethod = "gsl_nls")) })
  otherHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1 * dbh^b2), other2022, start = list(a1 = 1.2, b1 = 1.0, b2 = -0.08)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1 * (1 - exp(b1 * dbh^b2)), other2022, start = list(a1 = 30, b1 = -0.05, b2 = 0.9)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + a1 * (1 - exp(b1 * dbh^(b2 + b2bal * basalAreaLarger))), other2022, start = list(a1 = 40, b1 = -0.04, b2 = 0.9, b2bal = 0), significant = FALSE) # a1ba, a1bal, b1ba, b1bal, b2ba, b2bal not significant
  #print(to_parameter_confidence_intervals(otherHeightFromDiameter$weibullBal), n = 16)
  #to_fixed_coeffficients(otherHeightFromDiameter$sharmaPartonRelDbhPhysio)
  #bind_rows(lapply(otherHeightFromDiameter$weibullBal$fit, function(fit) { return(fit$stats$validation) }))
  #ggplot() +
  #  geom_point(aes(x = dbh, y = height), other2022, alpha = 0.1, shape = 16) +
  #  geom_line(aes(x = dbh, y = 1.37 + 35*exp(-6*dbh^-0.6)), other2022) # Korf
  #  geom_line(aes(x = dbh, y = 1.37 + 70 / (1 + 60 * dbh^-0.9)), other2022) # Hossfeld IV
  #  geom_line(aes(x = dbh, y = 1.37 + 75*dbh^0.75 / (75 + dbh^0.75)), other2022) # Michaelis-Menten
  #  geom_line(aes(x = dbh, y = 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d))), other2022) # Richards
  #  geom_line(aes(x = dbh, y = 1.37 + 150*topHeight^-0.15*(1 - exp(-0.01*(standTreesPerHectare/standBasalAreaPerHectare)^-0.3*dbh))^0.77), other2022) # Sharma-Parton
  #ggplot() +
  #  geom_point(aes(x = dbh, y = height), other2022, alpha = 0.1, shape = 16) +
  #  geom_line(aes(x = dbh, y = predict(otherHeightFromDiameter$hossfeld$fit[[1]], other2022)), other2022)
  
  # all significant GAM forms NaN with Γ link
  otherHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I(1.523 * dbh^0.697) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022) # plantation significant, 10x10 MAE 2.25, AIC 649
  #otherHeightFromDiameter$gamKorf = fit_gam("REML GAM", height ~ I(332.5834 * exp(-5.5934 * dbh^-0.1847)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022) # 10x10 MAE 2.23, AIC 766, marginally avoids collapse to linear
  #otherHeightFromDiameter$gamSibbesen = fit_gam("REML GAM", height ~ I(1.24408 * dbh^(1.01570 * dbh^-0.08286)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022, significant = FALSE) # 10x10 MAE 2.25, AIC 673 but collapses to linear
  #otherHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, bs = "ts", k = 4), data = other2022, significant = FALSE) # k = 4 minimum collapses to linear, plantation not significant
  #otherHeightFromDiameter$gamBal = fit_gam("REML GAM BAL", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, bs = "ts", k = 4), data = other2022, significant = FALSE) # k = 4 minimum collapses to linear, plantation not significant
  otherHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", k = 11), data = other2022, significant = FALSE) # k = 11 minimum > ~0, plantation not significant, 10x0 AIC 
  otherHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, bs = "ts", by = isPlantation, k = 11), data = other2022, significant = FALSE) # by propagation
  otherHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), data = other2022, significant = FALSE) # by propagation
  otherHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 11), data = other2022) # k = 16 minimum possible -> drop basal area -> k = 11 minimum, significant but not by much, plantation not reliably significant, 10x10 AIC 790
  # greatest accuracy increase from slope and roughness, all but wetness AIC disadvantageous
  #otherHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I(1.523 * dbh^0.697) + s(dbh, elevation, bs = "ts", by = as.factor(isPlantation), k = 6), data = other2022) # 10x10 AIC 755
  #otherHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I(1.523 * dbh^0.697) + s(dbh, slope, bs = "ts", by = as.factor(isPlantation), k = 9), data = other2022) # 10x10 AIC 744
  #otherHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I(1.523 * dbh^0.697) + s(dbh, sin(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 4), data = other2022) # k = 4 minimum, 10x10 AIC 822
  #otherHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I(1.523 * dbh^0.697) + s(dbh, cos(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022) # 10x10 AIC 761
  #otherHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 8), data = other2022) # 10x10 AIC 759
  #otherHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I(1.523 * dbh^0.697) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022) # 10x10 AIC 720
  otherHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, bs = "ts", k = 4), data = other2022, threads = 4) # k = 16 minimum with slope, roughness, wetness -> drop slope k = 11 minimum -> drop wetness -> k = 4 minimum, plantation not significant, 10x10 AIC 817
  otherHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(1.523 * dbh^0.697) + s(dbh, relativeDiameter, bs = "ts", k = 6), data = other2022) # plantation not reliably significant, 10x10 AIC 812
  otherHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, relativeDiameter, bs = "ts", k = 11, by = as.factor(isPlantation)), data = other2022, threads = 4) # k = 11 minimum > ~7, debatable if should be marked significant
  #otherHeightFromDiameter$gamBalRelDbh$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #lapply(otherHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(otherHeightFromDiameter$gamRelDbhPhysio$fit, summary)

  otherHeightFromDiameter$randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"),
                                                    other2022, mtry = 2, minNodeSize = 76, sampleFraction = 0.869) # from setup.R tuneRanger @ 500 trees
  otherHeightFromDiameter$randomForestBalPhysioRelDbh = fit_ranger("random forest BA+L RelDbh physio", c("height", "dbh", "relativeDiameter", "standQmd", "topHeight", "standBasalAreaPerHectare", "basalAreaLarger", "standTreesPerHectare", "slope", "elevation", "terrainRoughness", "topographicWetnessFD8f", "aspect", "isPlantation"), # from setup.R VSURF
                                                                   other2022, mtry = 3, minNodeSize = 2, sampleFraction = 0.839) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(otherHeightFromDiameter, paste0("trees/height-diameter/data/other height ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}

if (otherOptions$fitHeightMixed)
{
  #otherHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1 * dbh))^(b2 + b2r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve
  #                                                               start = list(fixed = c(a1 = 32, b1 = -0.5, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 1E-4)))
  otherHeightFromDiameterMixed = list(chapmanRichards = create_fit_statistics(name = "Chapman-Richards", fitting = "nlme"))
  #otherHeightFromDiameterMixed = list(chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp((b1) * dbh))^(b2 + b2r), other2022, # a1r, b1r, b24 singularity in backsolve
  #                                                                  fixedFormula = a1 + a1bal + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot, # |stand singularity in backsolve
  #                                                                  start = list(fixed = c(a1 = 30, a1bal = 0, b1 = -0.02, b2 = 0.84))))
  otherHeightFromDiameterMixed$chapmanRichardsBal = create_fit_statistics(name = "Chapman-Richards BA+L", fitting = "nlme", significant = FALSE)
  #otherHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + a1bal * basalAreaLarger + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                                 fixedFormula = a1 + a1bal + a1e + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                                 start = list(fixed = c(a1 = 80, a1bal = 0.3, a1e = 0.005, b1 = -0.005, b2 = 0.8)))
  otherHeightFromDiameterMixed$chapmanRichardsBalPhysio = create_fit_statistics(name = "Chapman-Richards BA+L physio", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = create_fit_statistics("Chapman-Richards BA+L RelDbh physio", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsBalRelDbh = create_fit_statistics("Chapman-Richards BA+L RelDbh", fitting = "nlme", significant = FALSE)
  #otherHeightFromDiameterMixed = list(chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                                     fixedFormula = a1 + a1e + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                                     start = list(fixed = c(a1 = 70, a1e = -0.01, b1 = -0.006, b2 = 0.81))))
  otherHeightFromDiameterMixed$chapmanRichardsPhysio = create_fit_statistics(name = "Chapman-Richards physio", fitting = "nlme")
  otherHeightFromDiameterMixed$chapmanRichardsRelDbh = create_fit_statistics("Chapman-Richards RelDbh", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsRelDbhPhysio = create_fit_statistics("Chapman-Richards RelDbh physio", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1)*dbh / (1 + dbh)^(b1 + b1r), other2022, # a1r less accurate
                                                 fixedFormula = a1 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 1.0, b1 = 0.18)))
  #otherHeightFromDiameterMixed$curtisA1 = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1r)*dbh / (1 + dbh)^(b1), other2022,
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.0, b1 = 0.18)))
  otherHeightFromDiameterMixed$hossfeld = tryCatch({ fit_nlme("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^(b1 + b1r)), other2022, # a1r, a2r max iterations, b1r 2x50 job computationally singular
                                                              fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 34, b1 = 25, b2 = -0.9))) },
                                                     error = function(e) { return(create_fit_statistics("Hossfeld IV", fittingMethod = "nlme")) })
  #otherHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1)*exp((b1 + b1r) * dbh^(b2)), other2022, # |stand/plot a1r job singularity in backsolve, b1r singularity in backsolve, b2r step halving, |plot b2r step halving, a1r also tractable
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot, # |stand, |uniquePlotID step halving
  #                                             start = list(fixed = c(a1 = 100, b1 = -6, b2 = -0.2)))
  otherHeightFromDiameterMixed$korf = create_fit_statistics(name = "Korf", fitting = "nlme")
  #otherHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + a1*dbh^(b1 + b1r) / (a2 + dbh^(b1 + b1r)), other2022, # a1r max iterations, a2r singularity in backsolve, potential condition number and job step halving with b1r|stand/plot -> reduced to stand -> job max iterations
  #                                                        fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1r ~ 1|stand, # |uniquePlotID step halving
  #                                                        start = list(fixed = c(a1 = 100, a2 = 100, b1 = 0.8)))
  otherHeightFromDiameterMixed$michaelisMenten = create_fit_statistics(name = "Michaelis-Menten", fitting = "nlme")
  otherHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + (a2 + a2r)*dbh + a3), other2022, # a1r typically less accurate than a2r, a3r max iterations
                                                 fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 0.025, a2 = 1.08, a3 = -0.01)))
  otherHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1)*dbh^(b1 + b1r), other2022, # a1r less accurate
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.0, b1 = 0.85)))
  #otherHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1 + a1r)*dbh^(b1), other2022,
  #                                              fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 1.0, b1 = 0.85)))
  otherHeightFromDiameterMixed$ratkowsky = tryCatch({ fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1 + b1r)/(dbh + b2)), other2022, # a1r less accurate, b1r 2x50 job step halving, b2r step halving
                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 20, b1 = -10, b2 = 2))) },
                                                    error = function(e) { return(create_fit_statistics("Ratkowsky", fittingMethod = "nlme")) })
  #otherHeightFromDiameterMixed$ratkowskyA1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1r)*exp((b1)/(dbh + b2)), other2022,
  #                                                    fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 20, b1 = -10, b2 = 2)))
  #otherHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU + kUr) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), other2022, # Har, kUr, dr singularity in backsolve
  #                                                  fixedFormula = Ha + d + kU ~ 1, randomFormula = kUr ~ 1|stand/plot,
  #                                                  start = list(fixed = c(Ha = 50, d = 0.2, kU = 0.01)))
  otherHeightFromDiameterMixed$richardsW = create_fit_statistics(name = "unified Richards", fitting = "nlme")
  otherHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1r)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, other2022, # a1r, b2r, b4r singularity in backsolve, b3r somewhat less accurate
                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  #otherHeightFromDiameterMixed$sharmaPartonB3 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3r)*dbh))^b4, other2022,
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  otherHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, other2022, # a1r, b2r, b4r singularity in backsolve, b3r job max iterations
                                                          fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 12, b1 = 0.3, b2 = -0.02, b3 = 0.3, b4 = 0.9)))
  otherHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, # a1r, b2r, b4r singularity in backsolve, b3r less accurate
                                                                fixedFormula = a1 + b1 + b2 + b3 + b4 + b4as ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0)))
  #otherHeightFromDiameterMixed$sharmaPartonBalPhysioB3 = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3r)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022,
  #                                                                fixedFormula = a1 + b1 + b2 + b3 + b4 + b4as ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0)))
  otherHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, # a1r, b4r singularity in backsolve, b2r, b3r step halving
                                                                      fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 + b4as ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                                      start = list(fixed = c(a1 = 9, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.06, b4 = 1.0, b4as = -0.1)))
  otherHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton BA+L RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, # a1r, b2r, b4r singularity in backsolve, b3r less accurate
                                                                fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 30, b1 = 0.1, b2 = -0.02, b3 = 0.2, b3rd = -0.04, b4 = 1.0)))
  #otherHeightFromDiameterMixed$sharmaPartonBalRelDbhB3 = fit_nlme("Sharma-Parton BA+L RelDbh", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter + b3r)*dbh))^b4, other2022,
  #                                                                fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 30, b1 = 0.1, b2 = -0.02, b3 = 0.2, b3rd = -0.04, b4 = 1.0)))
  #otherHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1r + a1tw * topographicWetnessFD8f)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, 
  #                                                           fixedFormula = a1 + a1tw + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 5, a1tw = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  otherHeightFromDiameterMixed$sharmaPartonPhysio = create_fit_statistics("Sharma-Parton physio", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$sharmaPartonRelDbh = create_fit_statistics("Sharma-Parton RelDbh", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$sharmaPartonRelDbhPhysio = create_fit_statistics("Sharma-Parton RelDbh physio", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1) * (1 - exp(b2*standTreesPerHectare^(b3 + b3r)*dbh))^b4, other2022, # a1r step halving, b2r, b4r singularity in backsolve, b1r of similar accuracy with notably conflicting signals between MAE and AIC
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 20, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  #otherHeightFromDiameterMixed$sharmaZhangB1 = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1r) * (1 - exp(b2*standTreesPerHectare^(b3)*dbh))^b4, other2022,
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 20, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  #otherHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^(b3 + a1r)*dbh))^b4, other2022, # a1r, b1r, b2r, b3r, b4r singularity in backsolve, presumably due to a1bal lacking significance
  #                                                       fixedFormula = a1 + a1bal + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 20, a1bal = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  otherHeightFromDiameterMixed$sharmaZhangBal = create_fit_statistics("Sharma-Zhang BA+L", fitting = "nlme", significant = FALSE)
  otherHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1r) * dbh^b2), other2022, # a1r, b2r max iterations
                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 1.2, b1 = 0.8, b2 = -0.05)))
  otherHeightFromDiameterMixed$weibull = tryCatch({ fit_nlme("Weibull", height ~ 1.37 + a1*(1 - exp(b1 * dbh^(b2 + b2r))), other2022, # a1r max iterations, b1r step halving, b2r 2x50 job step halving
                                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 90, b1 = -0.015, b2 = 0.8))) },
                                                  error = function(e) { return(create_fit_statistics("Weibull", fittingMethod = "nlme")) })
  otherHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp(b1 * dbh^(b2 + b2r))), other2022, # a1r max iterations, b1r singular, a1bal not significant
                                                     fixedFormula = a1 + a1ba + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 90, a1ba = 0, b1 = -0.015, b2 = 0.8)), significant = FALSE)
  #to_fixed_coeffficients(otherHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(otherHeightFromDiameterMixed$weibullBal), n = 16)
  #otherHeightFromDiameterMixed$sharmaZhang$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  otherHeightFromDiameterMixed$gam = fit_gam("REML GAM", height ~ I(1.523 * dbh^0.697) + s(dbh, bs = "ts", k = 5), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE) # plantation not significant, collapses to linear
  otherHeightFromDiameterMixed$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", k = 11), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE)
  otherHeightFromDiameterMixed$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, bs = "ts", by = isPlantation, k = 11), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE)
  otherHeightFromDiameterMixed$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE)
  otherHeightFromDiameterMixed$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 11), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE) # k = 11 minimum> ~6
  otherHeightFromDiameterMixed$gamPhysio = fit_gam("REML GAM physio", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, bs = "ts", k = 4), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE) # k = 4 minimum, collapses to linear
  otherHeightFromDiameterMixed$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(1.523 * dbh^0.697) + s(dbh, relativeDiameter, bs = "ts", k = 4), random = list(uniquePlotID = ~1), data = other2022) # k = 4 minimum, only marginally significant
  otherHeightFromDiameterMixed$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, relativeDiameter, bs = "ts", k = 11, by = as.factor(isPlantation)), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE) # k = 11 minimum > ~5
  #lapply(otherHeightFromDiameterMixed$gamRelDbhPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(otherHeightFromDiameterMixed$gamRelDbhPhysio$fit, function(fit) { return(summary(fit$gam)) })
  
  saveRDS(otherHeightFromDiameterMixed, paste0("trees/height-diameter/data/other height mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## other species diameter regressions
if (otherOptions$fitDbh)
{
  otherDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37), other2022))
  otherDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I((height - 1.37)^2), other2022)

  otherDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^b2, other2022, start = list(a1 = 4, b1 = 0.3, b2 = 1.0), control = gsl_nls_control(scale = "levenberg")) # a1-b1 evaporation, a1p, b1p, b2p not significant
  otherDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(exp((b1)*(height - 1.37)) - 1)^(b2), other2022, start = list(a1 = 4, a1aba = 0.0, b1 = 0.3, b2 = 1.0), significant = FALSE) # a1-b1 evaporation, { a1, b1, b2 } x { aba, bat } not significant
  otherDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^(b2 + b2rh * relativeHeight), other2022, start = list(a1 = 4, b1 = 0.3, b2 = 1.0, b2rh = 0), significant = FALSE) # a1rh, b1rh, b2rh not significant
  otherDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, start = list(a1 = -6.5, a1p = 2, b1 = 0.3, b2 = 0.4, b2p = 0.08)) # b1p not significant
  otherDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*log(1 - pmin((b1)*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, start = list(a1 = -7, a1aba = 0.0, b1 = 0.35, b2 = 0.37, b2p = -0.025), significant = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  otherDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2as * sin(pi/180 * aspect) + b2e * elevation), 0.9999)), other2022, start = list(a1 = -5, b1 = 0.3, b2 = 0.4, b2as = 0.0, b2e = 0.0001)) # { a1, b1, b2 } x { s, ac, tr, tw }, a1as, b1as not significant, b1e often significant
  otherDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1p * isPlantation + a1rh * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, start = list(a1 = -6.5, a1p = 2, a1rh = 0, b1 = 0.3, b2 = 0.4, b2p = 0.08), significant = FALSE) # a1rh, b1rh, b2rh not significant
  otherDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, start = list(a1 = -100, a2 = -200, b1 = 1.5)) # a1p, a2p, b1p not significant, form not well posed
  otherDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + a2*sqrt(height - 1.37)), other2022, start = list(a1 = 2.3, a2 = -0.13)) # a1p, a2p not significant
  otherDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, other2022, start = list(a1 = 0.8, b1 = 1.3)) # a1p, b1p not significant
  otherDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1), other2022, start = list(a1 = 0.8, a1aba = 0, b1 = 1.3), significant = FALSE) # { a1, b1 } x { aba, bat } not significant
  otherDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1s * slope)*(height - 1.37)^b1, other2022, start = list(a1 = 2.7, a1s = 0, b1 = 0.8), significant= FALSE) # { a1, b1 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^b1, other2022, start = list(a1 = 0.8, a1rh = 0, b1 = 1.3), significant = FALSE) # a1rh, b1rh not significant
  otherDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, b1 = 1.5, b2 = 0), significant = FALSE) # a1p, b1p, b2, b2p not significant -> collapses to power
  otherDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), other2022, start = list(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0), significant = FALSE) # { a1, b1, b2 } x { aba, bat }, b2 not significant
  otherDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, start = list(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0, b2s = 0), significant = FALSE) # by propagation
  otherDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, start = list(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0, b2s = 0), significant = FALSE) # by propagation
  otherDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0), significant = FALSE) # by propagation
  otherDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, start = list(a1 = 1, b2s = 0, b1 = 1.5, b2 = 0), significant = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, a1rh = 0, b1 = 1.5, b2 = 0), significant = FALSE) # a1rh, a1rhp, b1rh, b1rhp, b2rh, b2rhp not significant
  otherDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a1rh * relativeHeight + a1e * elevation)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, a1rh = 0, a1e = 0, b1 = 1.5, b2 = 0), significant = FALSE)
  otherDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), other2022, start = list(a1 = 0.00003, a2 = 0.01, b1 = 1.4, Ha = 70)) # a1p, a2p, b1p not significant
  otherDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1), other2022, start = list(a1 = 10, b1 = 0.2, b2 = 0.1, b3 = -0.1), control = gsl_nls_control(scale = "levenberg")) # a1p, b1p, b2p, b3p not significant, a1-b2 evaporation, nan-inf with b4
  otherDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, b1 = 1.5, b2 = 0)) # a1p, b1p, b2p not significant
  otherDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), other2022, start = list(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0), significant = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  otherDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, a1aba = 0, a1tr = 0, b1 = 1.5, b2 = 0), significant = FALSE) # by propagation
  otherDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, a1aba = 0, a1tr = 0, a1rh = 0, b1 = 1.5, b2 = 0), significant = FALSE) # by propagation
  otherDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0), significant = FALSE)
  otherDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2tw * topographicWetnessFD8f)), other2022, start = list(a1 = 1, b1 = 1.5, b2 = 0, b2tw = 0), significant = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, a1rh = 0, b1 = 1.5, b2 = 0), significant = FALSE) # a1rh, b1rh, b2rh not significant
  otherDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a1rh * relativeHeight + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 3.5, a1rh = 0, a1tr = 0, b1 = 0.3, b2 = 0.33), significant = FALSE) # by propagation
  otherDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, other2022, start = list(a1 = -10000, b1 = 0.0001, b2 = 1.3))
  #to_fixed_coeffficients(otherDiameterFromHeight$weibull)
  #print(to_parameter_confidence_intervals(otherDiameterFromHeight$weibull), n = 16)
  #ggplot() +
  #  geom_point(aes(x = dbh, y = height), other2022, alpha = 0.1, shape = 16) +
  #  geom_line(aes(x = 6*(exp(0.3*(height - 1.37)) - 1)^0.3, y = height), other2022) + # Chapman-Richards replace
  #  geom_line(aes(x = -160 * (height - 1.37)^1.2 / (-250 - (height - 1.37)^1.5), y = height), other2022) # Michaelis-Menten replace
  #  geom_line(aes(x = (-10000*log(1 - pmin(0.0001*(height - 1.37), 0.9999)))^1.2, y = height), other2022) # Weibull

  # all significant GAM forms NaN with Γ link
  otherDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(0.834 * dbh^1.280) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022) # 10x10 RMSE 1.65, AIC 486
  #otherDiameterFromHeight$gamMichaelis = fit_gam("REML GAM", dbh ~ I(-202.573 * (height - 1.37)^1.479/(-339.456 - (height - 1.37)^1.479)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022) # 10x10 RMSE 9.27, AIC 963
  #otherDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.834 * dbh^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 15), data = other2022)
  #otherDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(0.834 * dbh^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, slope, cos(pi/180*aspect), sin(pi/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 21), data = other2022)
  #otherDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.834 * dbh^1.280) + s(height, basalAreaTaller, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 19), data = other2022)
  # all physiographic predictors increase AIC with only slope providing a slight accuracy increase -> physio only marginally significant
  #otherDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I(0.834 * dbh^1.280) + s(height, elevation, bs = "ts", by = as.factor(isPlantation), k = 6), data = other2022) # 10x10 AIC 586
  #otherDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I(0.834 * dbh^1.280) + s(height, slope, bs = "ts", by = as.factor(isPlantation), k = 11), data = other2022) # 10x10 AIC 562
  #otherDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I(0.834 * dbh^1.280) + s(height, sin(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022) # 10x10 AIC 531
  #otherDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I(0.834 * dbh^1.280) + s(height, cos(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 9), data = other2022) # 10x10 AIC 550
  #otherDiameterFromHeight$gamRoughness = fit_gam("REML GAM roughness", dbh ~ I(0.834 * dbh^1.280) + s(height, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 8), data = other2022) # 10x10 AIC 593
  #otherDiameterFromHeight$gamWetness = fit_gam("REML GAM westness", dbh ~ I(0.834 * dbh^1.280) + s(height, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 9), data = other2022) # 10x10 AIC 529
  otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * dbh^1.280) + s(height, slope, bs = "ts", by = as.factor(isPlantation), k = 11), data = other2022)
  otherDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.834 * dbh^1.280) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022) # 10x10 AIC 562, only marginally significant
  otherDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 11), data = other2022, significant = FALSE) # 10x10 AIC 996, 5.9x increase in median RMSE
  #otherDiameterFromHeight$gamRelHtPhysio$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #lapply(otherDiameterFromHeight$gam$fit, k.check)
  #lapply(otherDiameterFromHeight$gam$fit, summary)

  otherDiameterFromHeight$randomForest = fit_ranger("random forest", c("dbh", "height"),
                                                    other2022, mtry = 1, minNodeSize = 75, sampleFraction = 0.688) # from setup.R tuneRanger @ 500 trees
  #otherDiameterFromHeight$randomForestAbat = fit_ranger("random forest ABA+T", c("dbh", "height", "bootstrapStandBasalAreaPerHectare"), # from setup.R VSURF
  #                                                      other2022, mtry = 2, minNodeSize = 2, sampleFraction = 0.592) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(otherDiameterFromHeight, paste0("trees/height-diameter/data/other DBH ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


if (otherOptions$fitDbhMixed)
{
  #otherDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1)*(exp((b1)*(height - 1.37)) - 1)^(b2 + b2r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve
  #                                                              start = list(fixed = c(a1 = 4, b1 = 0.3, b2 = 1.0))))
  otherDiameterFromHeightMixed = list(chapmanReplace = create_fit_statistics("Chapman-Richards replace", fitting = "nlme", significant = FALSE))
  #otherDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(exp((b1)*(height - 1.37)) - 1)^(b2), other2022,
  #                                                           fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 4, a1aba = 0.0, b1 = 0.3, b2 = 1.0)), significant = FALSE)
  otherDiameterFromHeightMixed$chapmanReplaceAbat = create_fit_statistics("Chapman-Richards replace ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(exp((b1)*(height - 1.37)) - 1)^(b2 + b2r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                            fixedFormula = a1 + a1rh + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 1.0, a1rh = 0, b1 = 1.3, b2 = 0.35)), significant = FALSE)
  otherDiameterFromHeightMixed$chapmanReplaceRelHt = create_fit_statistics("Chapman-Richards replace RelHt", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2r), 0.9999)), other2022, # a1r max iterations, b1r, b2r step halving
  #                                                        fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = -6.5, a1p = 2, b1 = 0.3, b2 = 0.4, b2p = 0.08))))
  otherDiameterFromHeightMixed$chapmanRichards = create_fit_statistics("Chapman-Richards inverse", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*log(1 - pmin((b1)*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, 
  #                                                            fixedFormula = a1 + a1aba + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = -7, a1aba = 0.0, b1 = 0.35, b2 = 0.37, b2p = -0.025)), significant = FALSE)
  otherDiameterFromHeightMixed$chapmanRichardsAbat = create_fit_statistics("Chapman-Richards inverse ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2as * sin(pi/180 * aspect) + b2e * elevation + b2r), 0.9999)), other2022, # a1r, b1r, b2r step halving
  #                                                                     fixedFormula = a1 + b1 + b2 + b2as + b2e ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                                     start = list(fixed = c(a1 = -5, b1 = 0.3, b2 = 0.4, b2as = 0.0, b2e = 0.0001))))
  otherDiameterFromHeightMixed$chapmanRichardsPhysio = create_fit_statistics("Chapman-Richards inverse physio", fitting = "nlme")
  #otherDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2rh * relativeHeight), 0.9999)), other2022, # a1r job step halving, b1r, b2r step halving
  #                                                                    fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                    start = list(fixed = c(a1 = -2.5, b1 = 0.5, b2 = 0.4, b2rh = 2))))
  otherDiameterFromHeightMixed$chapmanRichardsRelHt = create_fit_statistics("Chapman-Richards inverse RelHt", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation + a1r) * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, # a1r, a2r, b1r step halving
  #                                                               fixedFormula = a1 + a1p + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot, # |stand max iterations, |uniquePlotID singularity in backsolve
  #                                                               start = list(fixed = c(a1 = -100, a1p = 0, a2 = -200, b1 = 1.5)))
  otherDiameterFromHeightMixed$michaelisMentenReplace = create_fit_statistics("Michaelis-Menten replace", fitting = "nlme")
  otherDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + (a2 + a2r)*sqrt(height - 1.37)), other2022, # a1r max iterations
                                                  fixedFormula = a1 + a2 ~ 1, randomFormula = a2r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.8, a2 = -0.2)))
  otherDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1)*(height - 1.37)^(b1 + a1r), other2022, # a1r step halving
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot, 
                                                start = list(fixed = c(a1 = 1.3, b1 = 1.1)))
  #otherDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^b1, other2022, 
  #                                                  fixedFormula = a1 + a1aba + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 0.8, a1aba = 0, b1 = 1.3)), significant = FALSE)
  otherDiameterFromHeightMixed$powerAbat = create_fit_statistics("power ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1r + a1e * elevation)*(height - 1.37)^b1, other2022, # fixed physiographic effects not significant
  #                                                    fixedFormula = a1 + a1e + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2.7, a1e = -0.001, b1 = 0.8)))
  otherDiameterFromHeightMixed$powerPhysio = create_fit_statistics("power physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight)*(height - 1.37)^b1, other2022, # relative height not significant
  #                                                   fixedFormula = a1 + a1rh + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 3.0, a1rh = 2.3, b1 = 0.4)/0))
  otherDiameterFromHeightMixed$powerRelHt = create_fit_statistics("power RelHt", fitting = "nlme", significant = FALSE)
  otherDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ a1*(height - 1.37)^((b1 + b1r) * exp((b2) * (height - 1.37))), other2022, # a1r step halving, b2r job max iterations
                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.3, b1 = 1.1, b2 = 0.002)))
  #otherDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), other2022, 
  #                                                  fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$ruarkAbat = create_fit_statistics("Ruark ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
  #                                                        fixedFormula = a1 + a1aba + b1 + b2 + b2s ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0, b2s = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$ruarkAbatPhysio = create_fit_statistics("Ruark ABA+T physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
  #                                                             fixedFormula = a1 + a1aba + a1e + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0, b2s = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$ruarkAbatPhysioRelHt = create_fit_statistics("Ruark RelHt ABA+T physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
  #                                                       fixedFormula = a1 + a1aba + a1rh + b1 + b2 + b2s ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$ruarkAbatRelHt = create_fit_statistics("Ruark RelHt ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1e * elevation + a1r)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed physiographic effects not significant
  #                                                    fixedFormula = a1 + a1e + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2.67, a1e = 0, b1 = 0.813, b2 = 0.0067)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5))
  otherDiameterFromHeightMixed$ruarkPhysio = create_fit_statistics("Ruark physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight + a1r)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed effects not significant
  #                                                   fixedFormula = a1 + a1rh + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1, a1rh = 0, b1 = 1.5, b2 = 0)))
  otherDiameterFromHeightMixed$ruarkRelHt = create_fit_statistics("Ruark RelHt", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a1e * elevation + a1rh * relativeHeight + a1r)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed effects not significant 
  #                                                         fixedFormula = a1 + a1e + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 3.4, a1e = -0.002, a1rh = 1.2, b1 = 0.5, b2 = 0.0035)))
  otherDiameterFromHeightMixed$ruarkRelHtPhysio = create_fit_statistics("Ruark RelHt physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/(a1 + a1r) * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), other2022, # a1r, a2r, Har singularity in backsolve
  #                                                fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.00003, a1bat = 0.01, b1 = 1.4, Ha = 70)))
  otherDiameterFromHeightMixed$schnute = create_fit_statistics("Schnute inverse", fitting = "nlme")
  otherDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + a1r)*(height - 1.37)) - 1), other2022, # a1r, b1r, b2r step halving
                                                       fixedFormula = a1 + b1 + b2 + b3 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 35, b1 = -0.2, b2 = 0.1, b3 = -0.1)))
  otherDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ a1*(height - 1.37)^((b1 + b1r)*(height - 1.37)^(b2)), other2022, # a1r step halving, b2r singular precision, a1-b1 evaporation risk
                                                          fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 1, b1 = 1.5, b2 = 0)))
  #otherDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), other2022, 
  #                                                            fixedFormula = a1 + a1bat + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0)))
  otherDiameterFromHeightMixed$sibbesenReplaceAbat = create_fit_statistics("Sibbesen replace ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022,
  #                                                                  fixedFormula = a1 + a1aba + a1tr + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 1, a1aba = 0, a1tr = 0, b1 = 1.5, b2 = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = create_fit_statistics("Sibbesen replace ABA+T physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, #
  #                                                                       fixedFormula = a1 + a1aba + a1tr + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                       start = list(fixed = c(a1 = 1, a1aba = 0, a1tr = 0, a1rh = 0, b1 = 1.5, b2 = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = create_fit_statistics("Sibbesen replace RelHt ABA+T physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
  #                                                                 fixedFormula = a1 + a1bat + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0)), significant = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = create_fit_statistics("Sibbesen replace RelHt ABA+T", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed physiographic effects not significant 
  #                                                              fixedFormula = a1 + a1tr + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = 3.5, a1tr = -0.01, b1 = 0.3, b2 = 0.33)))
  otherDiameterFromHeightMixed$sibbesenReplacePhysio = create_fit_statistics("Sibbesen replace physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed relative height not significant 
  #                                                             fixedFormula = a1 + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 3.0, a1rh = -0.7, b1 = 0.3, b2 = 0.348)))
  otherDiameterFromHeightMixed$sibbesenReplaceRelHt = create_fit_statistics("Sibbesen replace RelHt", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a1tr * terrainRoughness + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed effects not significant 
  #                                                                   fixedFormula = a1 + a1tr + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                   start = list(fixed = c(a1 = 3.5, a1tr = -0.01, a1rh = 0, b1 = 0.3, b2 = 0.33)))
  otherDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = create_fit_statistics("Sibbesen replace RelHt physio", fitting = "nlme", significant = FALSE)
  #otherDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2 + b2r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = -10000, b1 = 0.0001, b2 = 1.3)))
  otherDiameterFromHeightMixed$weibull = create_fit_statistics("Weibull inverse", fitting = "nlme")
  #to_fixed_coeffficients(otherDiameterFromHeightMixed$sibbesenReplace)

  otherDiameterFromHeightMixed$gam = fit_gam("REML GAM", dbh ~ I(0.834 * dbh^1.280) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), random = list(uniquePlotID = ~1), data = other2022)
  #otherDiameterFromHeightMixed$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.834 * dbh^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 15), random = list(uniquePlotID = ~1), data = other2022)
  #otherDiameterFromHeightMixed$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(0.834 * dbh^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, slope, cos(pi/180*aspect), sin(pi/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 21), random = list(uniquePlotID = ~1), data = other2022)
  #otherDiameterFromHeightMixed$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.834 * dbh^1.280) + s(height, basalAreaTaller, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 19), random = list(uniquePlotID = ~1), data = other2022)
  otherDiameterFromHeightMixed$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * dbh^1.280) + s(height, slope, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(uniquePlotID = ~1), data = other2022)
  otherDiameterFromHeightMixed$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.834 * dbh^1.280) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 8), random = list(uniquePlotID = ~1), data = other2022)
  otherDiameterFromHeightMixed$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(uniquePlotID = ~1), data = other2022, significant = FALSE)
  #lapply(otherDiameterFromHeightMixed$gamRelHt$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(otherDiameterFromHeightMixed$gamRelHt$fit, function(fit) { return(summary(fit$gam)) })
  
  saveRDS(otherDiameterFromHeightMixed, paste0("trees/height-diameter/data/other DBH mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## collect model parameters
if (otherOptions$fitHeight & otherOptions$fitHeightMixed & otherOptions$fitDbh & otherOptions$fitDbhMixed)
{
  if (exists("otherHeightFromDiameter") == FALSE)
  { 
    otherHeightFromDiameter = readRDS(paste0("trees/height-diameter/data/other height ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  if (exists("otherHeightFromDiameterMixed") == FALSE)
  { 
    otherHeightFromDiameterMixed = readRDS(paste0("trees/height-diameter/data/other height mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  }
  if (exists("otherDiameterFromHeight") == FALSE)
  { 
    otherDiameterFromHeight = readRDS(paste0("trees/height-diameter/data/other DBH ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  if (exists("otherDiameterFromHeightMixed") == FALSE)
  { 
    otherDiameterFromHeightMixed = readRDS(paste0("trees/height-diameter/data/other DBH mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  
  otherFixedEffects = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, unnest_fixed_effects)),
                                          bind_rows(lapply(otherHeightFromDiameterMixed, unnest_fixed_effects))) %>%
                                  mutate(responseVariable = "height"),
                                bind_rows(bind_rows(lapply(otherDiameterFromHeight, unnest_fixed_effects)),
                                          bind_rows(lapply(otherDiameterFromHeightMixed, unnest_fixed_effects))) %>%
                                  mutate(responseVariable = "DBH")) %>%
    mutate(species = "other")
  otherFitStats = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, unnest_validation_statistics)),
                                      bind_rows(lapply(otherHeightFromDiameterMixed, unnest_validation_statistics))) %>%
                              mutate(responseVariable = "height"),
                            bind_rows(bind_rows(lapply(otherDiameterFromHeight, unnest_validation_statistics)),
                                      bind_rows(lapply(otherDiameterFromHeightMixed, unnest_validation_statistics))) %>%
                              mutate(responseVariable = "DBH")) %>%
    mutate(species = "other")
  
  check_plot_fit_stats(otherFitStats)
  saveRDS(otherFixedEffects, paste0("trees/height-diameter/data/other coefficients ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  saveRDS(otherFitStats, paste0("trees/height-diameter/data/other fit stats ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## preferred forms identified (results.R, Figure 9)
if (otherOptions$recalcPreferredModels)
{
  otherHeightFromDiameterPreferred = list(randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"), other2022, mtry = 2, minNodeSize = 76, sampleFraction = 0.869, folds = 1, repetitions = 1))
  otherHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0.1), folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1 * dbh^b2), other2022, start = list(a1 = 1.2, b1 = 1.0, b2 = -0.08), folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$sibbesenMixed = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1r) * dbh^b2), other2022,
                                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 1.2, b1 = 0.8, b2 = -0.05)), folds = 1, repetitions = 1)
  saveRDS(otherHeightFromDiameterPreferred, "trees/height-diameter/data/other preferred height models.Rds")

  otherDiameterFromHeightPreferred = list(gam = fit_gam("REML GAM", dbh ~ I(0.834 * dbh^1.280) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022, folds = 1, repetitions = 1))
  saveRDS(otherDiameterFromHeightPreferred, "trees/height-diameter/data/other preferred diameter models.Rds")
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  # height
  #detach(package:gam, unload = TRUE)
  #library(mgcv)
  otherHeightGam = gam(height ~ I(1.523 * dbh^0.697) + 
                                s(dbh, bs = "ts", k = 3) + # plantation not significant
                                s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 7) + # near linear, arguably no effect < relativeDiameter = 3
                                s(standBasalAreaPerHectare, bs = "ts", k = 3) + # plantation not significant, linear
                                s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 5) + # hyperbolic-like; linear + additional effect 0-30 m²/ha
                                #s(elevation, bs = "ts", k = 6) + # not significant
                                #s(slope, bs = "ts", k = 6) + # not significant
                                #s(aspect, bs = "ts", k = 6) + # not significant
                                #s(terrainRoughness, bs = "ts", k = 6) + # not significant
                                s(topographicWetnessFD8f, bs = "ts", k = 4), # hyperbolic
                        data = other2022, method = "REML", select = TRUE, weights = dbhWeight)
  k.check(otherHeightGam)
  summary(otherHeightGam)
  plot_gam_smooths(otherHeightGam)

  otherDbhGam = gam(dbh ~ I(0.834 * dbh^1.280) + 
                          s(height, bs = "ts", by = as.factor(isPlantation), k = 8) +
                          s(relativeHeight, bs = "ts", k = 3) + # plantation not significant, linear
                          #s(bootstrapStandBasalAreaPerHectare, bs = "ts", k = 3) + # not significant
                          #s(basalAreaTaller, bs = "ts", k = 3) + # not significant
                          s(elevation, bs = "ts", k = 3) + # curved
                          s(slope, bs = "ts", k = 5) + # curved
                          s(aspect, bs = "ts", k = 3) + # curved
                          #s(terrainRoughness, bs = "ts", k = 3) + # not significant
                          s(topographicWetnessFD8f, bs = "ts", k = 6), # squiggly, near no effect if collapsed to linear
                    data = other2022, folds = 1, repetitions = 1, method = "REML", select = TRUE, , weights = heightWeight)
  k.check(otherDbhGam)
  summary(otherDbhGam)
  plot_gam_smooths(otherDbhGam)

  otherHeightGam = gam(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11), data = other2022, family = Gamma, method = "REML", select = TRUE, weights = dbhWeight) # quasi(inverse, mu^2) (not as stable as gamma, major outlier)
  otherHeightGammRidge = gam(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11) + s(uniquePlotID, bs = "re"), data = other2022, control = gam.control(nthreads = htDiaOptions$rangerThreads), method = "REML", select = TRUE, weights = dbhWeight)
  sqrt(mean((other2022$height - otherHeightGammRidge$fitted.values)^2))
  otherHeightGammRidgePopulationPrediction = predict(otherHeightGammRidge, other2022, exclude = "s(uniquePlotID)", newdata.guaranteed = TRUE) # drop random effect to make population predictions
  #otherHeightGamm = gamm(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11), data = other2022, family = Gamma, method = "REML", random = list(uniquePlotID = ~1), weights = dbhWeight)
  #sqrt(mean((other2022$height - otherHeightGamm$gam$fitted.values)^2))
  #otherStandIndices = other2022 %>% select(stand) %>% mutate(stand = as.factor(stand), rowNumber = row_number()) %>% group_by(stand) %>% nest() %>% tibble::deframe()
  #otherHeightGammMarkov = gam(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11) + s(uniquePlotID, bs = "mrf", xt = list(nb = otherStandIndices)), data = other2022, control = gam.control(nthreads = htDiaOptions$rangerThreads), method = "REML", select = TRUE, weights = dbhWeight)
  otherHeightGamDbh0.6 = gam(height ~ I(1.523 * dbh^0.697) + s(dbh, bs = "ts", k = 8), data = other2022, family = Gamma, method = "REML", select = TRUE, weights = dbhWeight) # unweighted Gaussian RMSE 3.774 @ k = 8, gamma 3.880 @ k = 5
  k.check(otherHeightGamDbh0.6)
  sqrt(mean((other2022$height - otherHeightGamDbh0.6$fitted.values)^2))
  otherHeightGamDbh1.0 = gam(height ~ s(dbh, bs = "ts", k = 6), data = other2022, method = "REML", select = TRUE, weights = dbhWeight) # unweighted RMSE 3.928
  #otherHeightGam = gam(height ~ s(dbh, by = standBasalAreaPerHectare, bs = "ts", k = 7), data = other2022, method = "REML", select = TRUE, weights = dbhWeight)
  #otherHeightGam = gam(height ~ s(dbh, bs = "ts", k = 10) + s(standBasalAreaPerHectare, bs = "ts", k = 4), data = other2022, method = "REML", select = TRUE, weights = dbhWeight)
  #otherHeightGam = gam(height ~ s(dbh, bs = "ts", k = 10) + s(standBasalAreaPerHectare, bs = "ts", k = 4) + ti(dbh, standBasalAreaPerHectare, bs = "ts", k = 2), data = other2022, method = "REML", select = TRUE, weights = dbhWeight)
  
  ggplot() +
    geom_point(aes(x = dbh, y = height), other2022, alpha = 0.2, shape = 16) +
    geom_line(aes(x = dbh, y = otherHeightGammRidge$fitted.values, color = "gam(dbh^0.70, BA, BAL + ridge)|plantation", linetype = isPlantation, group = isPlantation), other2022, alpha = 0.5) +
    geom_line(aes(x = dbh, y = otherHeightGammRidgePopulationPrediction, color = "gam(dbh^0.70, BA, BAL - ridge)|plantation", linetype = isPlantation, group = isPlantation), other2022, alpha = 0.5) +
    geom_line(aes(x = dbh, y = otherHeightGam$fitted.values, color = "gam(dbh^0.70, BA, BAL)|plantation", linetype = isPlantation, group = isPlantation), other2022, alpha = 0.5) +
    #geom_line(aes(x = dbh, y = otherHeightGamDbh0.6$fitted.values, color = "gam(dbh^0.70)"), other2022) +
    labs(x = "DBH, cm", y = "height, m", color = NULL, linetype = NULL) +
    scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural regen", "plantation"), values = c("solid", "longdash"))
  
  # mgcv versus gam
  #detach(package:mgcv, unload = TRUE)
  #library(gam)
  #otherHeightGam = gam(height ~ I(dbh^0.6) + s(dbh, df = 4), data = other2022, weights = dbhWeight) # unweighted RMSE 3.900 with linear DBH base, 3.790 with dbh^0.6
  #sqrt(mean(otherHeightGam$residuals^2))
  #summary(otherHeightGam)
}