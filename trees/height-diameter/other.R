# load libraries, functions, and other2022 from setup.R

otherOptions = tibble(fitHeight = FALSE, 
                      fitHeightMixed = TRUE,
                      fitDbh = FALSE,
                      fitDbhMixed = FALSE)

## other species height regressions
if (otherOptions$fitHeight)
{
  otherHeightFromDiameter = list(linear = fit_lm("linear", height ~ 0 + dbh, other2022)) # a1p not significant in 50% of fits
  otherHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ 0 + dbh + I(dbh^2), other2022) # a1p, a2p not significant in 50% of fits

  otherHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, b1 = -0.02, b2 = 0.84)) # a1p, b1p, b2p (25%) not significant
  otherHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a2 = 0, b1 = -0.02, b2 = 0.84), significant = FALSE) # a1p, a2, a2p, a3, a3p, b1p, b2p, exponent BAL not significant
  otherHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a5 * topographicWetnessFD8f) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a2 = 0, a5 = 0, b1 = -0.02, b2 = 0.8), signficant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant
  otherHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * relativeDiameter + a5 * elevation) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 80, a2 = 0.3, a4 = 0, a5 = 0, b1 = -0.01, b2 = 0.8), significant = FALSE) # by propagation
  otherHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * relativeDiameter) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a2 = 0, a4 = 0, b1 = -0.02, b2 = 0.84), significant = FALSE) # singular gradient with a3, a4, a4p not significant, exponents NaN-Inf or singular gradient
  otherHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp(b1 * dbh))^(b2 + b2a * cos(pi/180 * aspect)), other2022, start = list(a1 = 30, b1 = -0.02, b2 = 0.81, b2a = -0.08), signficant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant, though b2as + b2ac intermittently test significant
  otherHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a4 * relativeDiameter)*(1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a4 = 0, b1 = -0.02, b2 = 0.8), significant = FALSE) # a4, a4p, b1rd, b2rd not significant, signs of a1-b1 evaporation
  otherHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a4 * relativeDiameter + a5 * elevation) * (1 - exp(b1 * dbh))^b2, other2022, start = list(a1 = 30, a4 = 0, a5 = 0, b1 = -0.02, b2 = 0.81), significant = FALSE) # by propagation
  otherHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1*dbh / (1 + dbh)^b1, other2022, start = list(a1 = 2, b1 = 0.4))
  otherHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), other2022, start = list(a1 = 44, a2 = 33, b1 = -0.9)) # prone to singular gradient
  otherHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1 * dbh^b2), other2022, start = list(a1 = 300, b1 = -5, b2 = -0.1)) # a1p, b1p, b2p not significant, form not well posed (should probably be a1 * (exp() - 1))
  otherHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), other2022, start = list(a1 = 50, a2 = 40, b1 = 0.9)) # a1p, a2p, b1p not significant
  otherHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3), other2022, start = list(a1 = -0.01, a2 = 2.5, a3 = -5)) # a1p, a2p, a3p not significant
  otherHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, other2022, start = list(a1 = 1.5, b1 = 0.7)) # a1p, b1p not significant
  otherHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2)), other2022, start = list(a1 = 30, b1 = -11, b2 = 3)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), other2022, start = list(Ha = 28, d = 0.2, kU = 0.025))
  otherHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^b4, other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant
  otherHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant
  otherHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3sa * sin(pi/180 * aspect))*dbh))^(b4), other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b3sa = 0.04, b4 = 1.0)) # { a5, b1, b2, b3, b4 } x { e, s, ac, tr, tw } not significant, b4sa also significant, b1sa, b2sa often also signficant
  otherHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b4rd * relativeDiameter + b3sa * sin(pi/180 * aspect))*dbh))^b4, other2022, start = list(a1 = 7, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3sa = 0.04, b4 = 1.0, b4rd = -0.06)) # a4, b1rd, b4rd not significant, b2rd often significant
  otherHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, start = list(a1 = 10, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.05, b4 = 1.0)) # a4, b1rd, b2rd, b4rd not significant
  otherHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a5 * topographicWetnessFD8f)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare))^(b3)*dbh))^(b4), other2022, start = list(a1 = 5, a5 = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0), significant = FALSE) # { a5, b1, b2, b3, b4 } x { e, s, a, tr, tw } not significant
  otherHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + (a1)*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, start = list(a1 = 8, b1 = 0.6, b2 = -0.01, b3 = 0.4, b3rd = -0.05, b4 = 1.0), significant = FALSE) # a4 not significant
  otherHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a4 * relativeDiameter + a5 * elevation)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, start = list(a1 = 140, a4 = -0.01, a4 = -10, b1 = -0.33, b2 = -0.007, b3 = -0.17, b4 = 0.81), signficant = FALSE) # by propagation
  otherHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, other2022, start = list(a1 = 30, b1 = 0, b2 = -0.015, b3 = 0.1, b4 = 0.9)) # a1p, b1p, b2p, b3p, b4p not significant, signs of a1-b2 evaporation
  otherHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, other2022, start = list(a1 = 30, a2 = 0, b1 = -0.3, b2 = -0.01, b3 = 0.2, b4 = 0.9), significant = FALSE) # a2, a2p not significant
  otherHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1 * dbh^b2), other2022, start = list(a1 = 1.2, b1 = 1.0, b2 = -0.08)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1 * (1 - exp(b1 * dbh^b2)), other2022, start = list(a1 = 30, b1 = -0.05, b2 = 0.9)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1 * dbh^b2)), other2022, start = list(a1 = 40, a2 = 0, b1 = -0.04, b2 = 0.9), significant = FALSE) # a2, a3, b1ba, b1bal, b2ba not significant
  #print(to_parameter_confidence_intervals(otherHeightFromDiameter$weibullBal), n = 16)
  #to_fixed_coeffficients(otherHeightFromDiameter$chapmanRichards)
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
  
  otherHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I(dbh^0.697) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022) # plantation significant, 2x25 AIC 4092 (3137-5402)
  #otherHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, bs = "ts", k = 4), data = other2022, significant = FALSE) # k = 4 minimum, plantation not significant, 2x25 AIC 4344 (3343-6743)
  #otherHeightFromDiameter$gamBal = fit_gam("REML GAM BAL", height ~ I(dbh^0.697) + s(dbh, basalAreaLarger, bs = "ts", k = 4), data = other2022, significant = FALSE) # k = 4 minimum, plantation not significant, 2x25 AIC 4216 (2913-5535)
  otherHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", k = 11), data = other2022, significant = FALSE) # k = 11 minimum possible with BA and BAL, plantation not significant, 2x25 AIC 5874 (2892-12662)
  otherHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(pi/180*aspect), cos(pi/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = isPlantation, k = 6), data = other2022, significant = FALSE) # insufficient covariants for all predictors
  otherHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, terrainRoughness, topographicWetnessFD8f, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 57), data = other2022)
  otherHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", k = 16), data = other2022) # k = 16 minimum possible, plantation not reliably significant, AIC 5076 (2813-10597)
  # physiographic predictors are all AIC disadvantageous
  #otherHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I(dbh^0.697) + s(dbh, elevation, bs = "ts", by = as.factor(isPlantation), k = 6), data = other2022) # 2x25 AIC 4471 (2603-9034)
  #otherHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I(dbh^0.697) + s(dbh, slope, bs = "ts", by = as.factor(isPlantation), k = 9), data = other2022) # 2x25 AIC 4167 (2602-5730)
  #otherHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I(dbh^0.697) + s(dbh, sin(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 4), data = other2022) # k = 4 minimum, 2x25 AIC 4420 (3182-6706)
  #otherHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I(dbh^0.697) + s(dbh, cos(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022) # 2x25 AIC 4141 (2870-5837)
  #otherHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I(dbh^0.697) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 8), data = other2022) # 2x25 AIC 4162 (2555-5174)
  #otherHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I(dbh^0.697) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022) # 2x25 AIC 4151 (2657-5622)
  #otherHeightFromDiameter$gamWetness$validation %>% summarize(minAic = min(aic), meanAic = mean(aic), maxAic = max(aic))
  otherHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(dbh^0.697) + s(dbh, elevation, slope, sin(3.141593/180*aspect), cos(3.141593/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", k = 331), data = other2022, significant = FALSE, threads = 4) # k = 331 minimum with all predictors (23-72s/fit @ 9950X 4 threads), plantation not significant
  otherHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(dbh^0.697) + s(dbh, relativeDiameter, bs = "ts", k = 6), data = other2022) # plantation not reliably significant, 2x25 AIC 4054 (2367-5819)
  otherHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(dbh^0.697) + s(dbh, elevation, slope, sin(3.141593/180*aspect), cos(3.141593/180*aspect), terrainRoughness, topographicWetnessFD8f, relativeDiameter, bs = "ts", k = 496, by = as.factor(isPlantation)), data = other2022, significant = FALSE, threads = 4) # by propagation, k = 496 minimum, AIC ~4700
  #otherHeightFromDiameter$gam$validation %>% summarize(maeMin = min(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMean = mean(aic), aicMax = max(aic))
  #bind_rows(lapply(otherHeightFromDiameter$gamRelDbhPhysio$fit, function(fit) { return(fit$stats$validation) })) %>% summarize(aic = mean(aic))
  #lapply(otherHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(otherHeightFromDiameter$gamRelDbhPhysio$fit, summary)
  
  saveRDS(otherHeightFromDiameter, "trees/height-diameter/data/other height.Rds")
}

if (otherOptions$fitHeightMixed)
{
  #otherHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1 * dbh))^(b2 + b2r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = 32, b1 = -0.5, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 1E-4)))
  #otherHeightFromDiameterMixed = list(chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1r + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1 * dbh))^b2, other2022, # singularity in backsolve
  #                                                                  fixedFormula = a1 + a2 + a3 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 70, a2 = 0, a3 = 0, b1 = -0.006, b2 = 0.84))))
  #otherHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + a2 * basalAreaLarger + a4 * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                                 fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                                 start = list(fixed = c(a1 = 80, a2 = 0.3, a4 = 0.005, b1 = -0.005, b2 = 0.8)))
  #otherHeightFromDiameterMixed = list(chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r + a4 * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                                     fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                                     start = list(fixed = c(a1 = 70, a4 = -0.01, b1 = -0.006, b2 = 0.81))))
  otherHeightFromDiameterMixed = list(curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1r)*dbh / (1 + dbh)^b1, other2022,
                                                        fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 1.0, b1 = 0.18))))
  otherHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^(b1 + b1r)), other2022,
                                                   fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 34, b1 = 25, b2 = -0.9)))
  #otherHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp(b1 * dbh^b2), other2022, # a1r job singularity in backsolve, b1r singularity in backsolve, b2r step halving 
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 100, b1 = -6, b2 = -0.2)))
  otherHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + a1*dbh^(b1 + b1r) / (a2 + dbh^(b1 + b1r)), other2022, # max iterations with a1r, potential condition number with b1r
                                                          fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 100, a2 = 100, b1 = 0.8)))
  otherHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + (a2 + a2r)*dbh + a3), other2022,  # a1r maybe less accurate than a2r, a3r max iterations
                                                 fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 0.025, a2 = 1.08, a3 = -0.01)))
  #bind_rows(lapply(otherHeightFromDiameterMixed$prodan$fit, function(fit) { return(fit$stats$validation) }))
  otherHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1r)*dbh^b1, other2022, # b1r also tractable
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.0, b1 = 0.85)))
  otherHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1r)*exp(b1/(dbh + b2)), other2022, 
                                                    fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 20, b1 = -10, b2 = 2)))
  #otherHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU + kUr) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), other2022, # Har, kUr, dr singularity in backsolve
  #                                                  fixedFormula = Ha + d + kU ~ 1, randomFormula = kUr ~ 1|stand/plot,
  #                                                  start = list(fixed = c(Ha = 50, d = 0.2, kU = 0.01)))
  otherHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3r)*dbh))^b4, other2022, # a1r, b2r, b4r singularity in backsolve, b1r also tractable
                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  otherHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3r)*dbh))^b4, other2022, 
                                                          fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  #otherHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1r + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, other2022,
  #                                                              fixedFormula = a1 + a4 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = 200, a4 = 0, b1 = -0.3, b2 = -0.006, b3 = -0.2, b4 = 0.8)))
  #otherHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1r + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, 
  #                                                           fixedFormula = a1 + a4 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 140, a4 = 0, b1 = -0.14, b2 = -0.01, b3 = -0.3, b4 = 0.78)))
  otherHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^(b3 + b3r)*dbh))^b4, other2022, # a1r step halving, b2r, b4r singularity in backsolve, b1r also tractable
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 20, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  #otherHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^(b3 + a1r)*dbh))^b4, other2022, # a1r, b1r, b2r, b3r, b4r singularity in backsolve, presumably due to a2 lacking significance
  #                                                       fixedFormula = a1 + a2 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 20, a2 = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  otherHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + a1r) * dbh^b2), other2022, # a1r, b2r max iterations
                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 1.2, b1 = 0.8, b2 = -0.05)))
  otherHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + a1*(1 - exp(b1 * dbh^(b2 + b2r))), other2022, # a1r max iterations, b1r step halving
                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 90, b1 = -0.015, b2 = 0.8)))
  otherHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a3 * standBasalAreaPerHectare) * (1 - exp(b1 * dbh^(b2 + b2r))), other2022, # a2 step halving, a3 not signficant
                                                     fixedFormula = a1 + a3 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 90, a3 = 0, b1 = -0.015, b2 = 0.8)), significant = FALSE)
  #to_fixed_coeffficients(otherHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(otherHeightFromDiameterMixed$weibullBal), n = 16)

  saveRDS(otherHeightFromDiameterMixed, "trees/height-diameter/data/other height mixed.Rds")
}


## other species diameter regressions
if (otherOptions$fitDbh)
{
  otherDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37), other2022))
  otherDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I((height - 1.37)^2), other2022)

  otherDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^b2, other2022, start = list(a1 = 4, b1 = 0.3, b2 = 1.0), control = gsl_nls_control(scale = "levenberg")) # a1-b1 evaporation, a1p, b1p, b2p not significant
  #otherDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(exp(b1*(height - 1.37)) - 1)^b2, other2022, start = list(a1 = 10, a2 = 0.05, b1 = 0.10, b2 = 0.5))
  #otherDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a2 * basalAreaTaller) * (exp(b1*(height - 1.37)^b2) - 1), other2022, start = list(a1 = 1.0, a2 = 0, b1 = 0.3, b2 = 1.0), significant = FALSE)
  #otherDiameterFromHeight$chapmanReplaceAbatRelHt = fit_gsl_nls("Chapman-Richards replace ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight) * (exp(b1*(height - 1.37)^b2) - 1), other2022, start = list(a1 = 1.0, a2 = 0, a4 = 0, b1 = 0.3, b2 = 1.0), significant = FALSE)
  otherDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ a1*(exp(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)) - 1), other2022, start = list(a1 = 1.0, b1 = 1.3, b2 = 0.35, b2rh = 0), control = gsl_nls_control(scale = "levenberg"), signficant = FALSE) # a4, b1rh, b2rh not significant
  otherDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, start = list(a1 = -6.5, a1p = 2, b1 = 0.3, b2 = 0.4, b2p = 0.08)) # b1p not significant
  #otherDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, start = list(a1 = -7, a2 = -0.03, b1 = 0.35, b2 = 0.37, b2p = -0.025))
  otherDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2a * sin(pi/180 * aspect) + b2e * elevation), 0.9999)), other2022, start = list(a1 = -5, b1 = 0.3, b2 = 0.4, b2a = 0.0, b2e = 0.0001)) # { a5, b1, b2 } x { s, ac, tr, tw }, a5as, b1as not significant, b1e often significant
  print(to_parameter_confidence_intervals(otherDiameterFromHeight$chapmanRichardsPhysio), n = 16)
  to_fixed_coeffficients(otherDiameterFromHeight$chapmanRichardsPhysio)
  
  otherDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1p * isPlantation + a4 * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, start = list(a1 = -6.5, a1p = 2, a4 = 0, b1 = 0.3, b2 = 0.4, b2p = 0.08), significant = FALSE) # a4, b1rh, b2rh not significant
  otherDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, start = list(a1 = -100, a1p = 0, a2 = -200, b1 = 1.5)) # a1p, a2p, b1p not significant, form not well posed
  otherDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + a2*sqrt(height - 1.37)), other2022, start = list(a1 = 2.3, a2 = -0.13)) # a1p, a2p not significant
  otherDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, other2022, start = list(a1 = 0.8, b1 = 1.3)) # a1p, b1p not signficant
  #otherDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^b1, other2022, start = list(a1 = 2.8, a2 = 0.02, b1 = 0.4))
  otherDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a5 * slope)*(height - 1.37)^b1, other2022, start = list(a1 = 2.7, a5 = 0, b1 = 0.8), significant= FALSE) # { a5, b1 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a4 * relativeHeight)*(height - 1.37)^b1, other2022, start = list(a1 = 0.8, a4 = 0, b1 = 1.3), significant = FALSE) # a4, b1rh not significant
  otherDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, b1 = 1.5, b2 = 0)) # a1p, b1p, b2p not significant
  #otherDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 2.4, a2 = 0.03, b1 = 0.4, b2 = 0.05))
  #otherDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a2 * basalAreaTaller + a4 * elevation)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 3.0, a2 = 0.02, a4 = -0.002, b1 = 0.5, b2 = 0.04))
  #otherDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a4 * elevation + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 3.0, a3 = 0.01, a4 = -0.002, a4 = 2, b1 = 0.42, b2 = 0.0033))
  #otherDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 2.4, a2 = 0.04, a4 = 3, b1 = 0.4, b2 = 0.03))
  otherDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, start = list(a1 = 1, b2s = 0, b1 = 1.5, b2 = 0), significant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, a4 = 0, b1 = 1.5, b2 = 0), significant = FALSE) # a4, a4p, b1rh, b1rhp, b2rh, b2rhp not significant
  otherDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a4 * relativeHeight + a5 * elevation)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, start = list(a1 = 1, a4 = 0, a5 = 0, b1 = 1.5, b2 = 0), significant = FALSE)
  otherDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2022, start = list(a1 = 0.00003, a2 = 0.01, b1 = 1.4, Ha = 70)) # a1p, a2p, b1p not significant
  otherDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1), other2022, start = list(a1 = 10, b1 = 0.2, b2 = 0.1, b3 = -0.1), control = gsl_nls_control(scale = "levenberg")) # a1p, b1p, b2p, b3p not significant, a1-b2 evaporation
  otherDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, b1 = 1.5, b2 = 0)) # a1p, b1p, b2p not significant
  #otherDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 2.7, a2 = 0.01, b1 = 0.39, b2 = 0.24))
  #otherDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a8 * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 3.1, a3 = 0.002, a8 = 0, b1 = 0.4, b2 = 0.25))
  #otherDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a8 * terrainRoughness + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 3.4, a3 = -0.005, a8 = 0, a4 = -1, b1 = 0.4, b2 = 0.29))
  #otherDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 2.9, a2 = 0, a4 = 0, b1 = 0.4, b2 = 0.25))
  otherDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2tw * topographicWetnessFD8f)), other2022, start = list(a1 = 1, b1 = 1.5, b2 = 0, b2tw = 0), significant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 1, a4 = 0, b1 = 1.5, b2 = 0), significant = FALSE) # a4, b1rh, b2rh not significant
  otherDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a4 * relativeHeight + a8 * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, start = list(a1 = 3.5, a4 = 0, a8 = 0, b1 = 0.3, b2 = 0.33), significant = FALSE) # by propagation
  otherDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, other2022, start = list(a1 = -10000, b1 = 0.0001, b2 = 1.3))
  #to_fixed_coeffficients(otherDiameterFromHeight$weibull)
  #print(to_parameter_confidence_intervals(otherDiameterFromHeight$weibull), n = 16)
  #ggplot() +
  #  geom_point(aes(x = dbh, y = height), other2022, alpha = 0.1, shape = 16) +
  #  geom_line(aes(x = 6*(exp(0.3*(height - 1.37)) - 1)^0.3, y = height), other2022) + # Chapman-Richards replace
  #  geom_line(aes(x = -160 * (height - 1.37)^1.2 / (-250 - (height - 1.37)^1.5), y = height), other2022) # Michaelis-Menten replace
  #  geom_line(aes(x = (-10000*log(1 - pmin(0.0001*(height - 1.37), 0.9999)))^1.2, y = height), other2022) # Weibull

  # individual term selection: height by = isPlantation, ABA, slope, elevation, sin(aspect), tsi retained by AIC but not significant (p >= 0.09)
  otherDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(dbh^1.280) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022)
  #otherDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(dbh^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 15), data = other2022)
  #otherDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(dbh^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, slope, cos(pi/180*aspect), sin(pi/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 21), data = other2022)
  #otherDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ I(dbh^1.280) + s(height, basalAreaTaller, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 19), data = other2022)
  #otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(dbh^1.280) + s(height, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), data = other2022)
  otherDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(dbh^1.280) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022)
  #otherDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 15), data = other2022)
  #lapply(otherDiameterFromHeight$gam$fit, k.check)
  #lapply(otherDiameterFromHeight$gam$fit, summary)
  
  saveRDS(otherDiameterFromHeight, "trees/height-diameter/data/other DBH.Rdata")
}


if (otherOptions$fitDbhMixed)
{
  #otherDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1)*(exp(b1*(height - 1.37)) - 1)^(b2 + a1r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = 4, b1 = 0.3, b2 = 1.0))))
  #otherDiameterFromHeightMixed = list(chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(exp(b1*(height - 1.37)) - 1)^b2, other2022,
  #                                                                  fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 10, a2 = 0.05, b1 = 0.10, b2 = 0.5))))
  #otherDiameterFromHeightMixed = list(chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller) * (exp(b1*(height - 1.37)^b2) - 1), other2022,
  #                                                                  fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 1.0, a2 = 0, b1 = 0.3, b2 = 1.0))))
  #otherDiameterFromHeightMixed$chapmanReplaceAbatRelHt = fit_nlme("Chapman-Richards replace ABA+T RelHt", dbh ~ (a1 + a1r + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(height - 1.37)^b2) - 1), other2022,
  #                                                                 fixedFormula = a1 + a3 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 0.9, a3 = 0.003, a4 = 0, b1 = 1.3, b2 = 0.3)))
  #otherDiameterFromHeightMixed = list(chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a4 * relativeHeight)*(exp(b1*(height - 1.37)^(b2 + b2r)) - 1), other2022, # a1r singularity in backsolve, b1r, b2r step halving
  #                                                                   fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                                   start = list(fixed = c(a1 = 1.0, a4 = 0, b1 = 1.3, b2 = 0.35))))
  #otherDiameterFromHeightMixed = list(chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2r), 0.9999)), other2022, # a1r max iterations, b1r, b2r step halving
  #                                                               fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = -6.5, a1p = 2, b1 = 0.3, b2 = 0.4, b2p = 0.08))))
  #otherDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, 
  #                                                            fixedFormula = a1 + a2 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = -7, a2 = -0.03, b1 = 0.35, b2 = 0.37, b2p = -0.025)))
  #otherDiameterFromHeightMixed = list(chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2a * sin(pi/180 * aspect) + b2e * elevation + b2r), 0.9999)), other2022, # a1r, b1r, b2r step halving
  #                                                                     fixedFormula = a1 + b1 + b2 + b2a + b2e ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                                     start = list(fixed = c(a1 = -5, b1 = 0.3, b2 = 0.4, b2a = 0.0, b2e = 0.0001))))
  otherDiameterFromHeightMixed = list(chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2rh * relativeHeight), 0.9999)), other2022, # b1r, b2r step halving
                                                                      fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                      start = list(fixed = c(a1 = -2.5, b1 = 0.5, b2 = 0.4, b2rh = 2))))
  #otherDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation + a1r) * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, # a1r, a2r, b1r step halving
  #                                                               fixedFormula = a1 + a1p + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = -100, a1p = 0, a2 = -200, b1 = 1.5)))
  otherDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + (a2 + a2r)*sqrt(height - 1.37)), other2022, # a1r max iterations
                                                  fixedFormula = a1 + a2 ~ 1, randomFormula = a2r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.8, a2 = -0.2)))
  otherDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1)*(height - 1.37)^(b1 + a1r), other2022, # a1r step halving
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot, 
                                                start = list(fixed = c(a1 = 1.3, b1 = 1.1)))
  #otherDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^b1, other2022, 
  #                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 2.8, a2 = 0.02, b1 = 0.4)))
  #otherDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1r + a4 * elevation)*(height - 1.37)^b1, other2022, # fixed physiographic effects not significant
  #                                                    fixedFormula = a1 + a4 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2.7, a4 = -0.001, b1 = 0.8)))
  #otherDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^b1, other2022, # relative height not significant
  #                                                   fixedFormula = a1 + a4 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 3.0, a4 = 2.3, b1 = 0.4)/0))
  otherDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ a1*(height - 1.37)^((b1 + b1r) * exp(b2 * (height - 1.37))), other2022, # a1r step halving, b2r also tractable
                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.3, b1 = 1.1, b2 = 0.002)))
  otherDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
                                                    fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2.6, a2 = 0.03, b1 = 0.5, b2 = 0.04)))
  #otherDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * elevation)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
  #                                                        fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 3.0, a2 = 0.02, a4 = -0.002, b1 = 0.5, b2 = 0.04)))
  #otherDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark ABA+T RelHt physio", dbh ~ (a1 + a1r + a3 * bootstrapStandBasalAreaPerHectare + a4 * elevation + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
  #                                                             fixedFormula = a1 + a3 + a4 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 3.0, a3 = 0.01, a4 = -0.002, a4 = 2, b1 = 0.42, b2 = 0.0033)))
  #otherDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
  #                                                       fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 2.4, a2 = 0.04, a4 = 3, b1 = 0.4, b2 = 0.03)))
  #otherDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1r + a4 * elevation)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed physiographic effects not significant
  #                                                    fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2.67, a4 = 0, b1 = 0.813, b2 = 0.0067)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5))
  #otherDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed effects not significant
  #                                                   fixedFormula = a1 + a4 + a4p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1, a4 = 0, b1 = 1.5, b2 = 0)))
  #otherDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a4 * elevation + a4 * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed effects not significant 
  #                                                         fixedFormula = a1 + a4 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 3.4, a4 = -0.002, a4 = 1.2, b1 = 0.5, b2 = 0.0035)))
  #otherDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/(a1 + a1r) * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2022, # a1r, a2r, Har singularity in backsolve
  #                                                fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.00003, a2 = 0.01, b1 = 1.4, Ha = 70)))
  otherDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + a1r)*(height - 1.37)) - 1), other2022, # a1r, b1r, b2r step halving
                                                       fixedFormula = a1 + b1 + b2 + b3 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 35, b1 = -0.2, b2 = 0.1, b3 = -0.1)))
  otherDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ a1*(height - 1.37)^((b1 + b1r)*(height - 1.37)^(b2)), other2022, # a1r step halving, b2r singular precision, a1-b1 evaporation risk
                                                          fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 1, b1 = 1.5, b2 = 0)))
  #otherDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
  #                                                            fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 2.7, a2 = 0.01, b1 = 0.39, b2 = 0.24)))
  #otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r + a3 * bootstrapStandBasalAreaPerHectare + a8 * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022,
  #                                                                  fixedFormula = a1 + a3 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 3.1, a3 = 0.002, a8 = 0, b1 = 0.4, b2 = 0.25)))
  #otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1r + a3 * bootstrapStandBasalAreaPerHectare + a8 * terrainRoughness + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, #
  #                                                                       fixedFormula = a1 + a3 + a8 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                       start = list(fixed = c(a1 = 3.4, a3 = -0.005, a8 = 0, a4 = -1, b1 = 0.4, b2 = 0.29)))
  #otherDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
  #                                                                 fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 2.9, a2 = 0, a4 = 0, b1 = 0.4, b2 = 0.25)))
  #otherDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a8 * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed physiographic effects not significant 
  #                                                              fixedFormula = a1 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = 3.5, a8 = -0.01, b1 = 0.3, b2 = 0.33)))
  #otherDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed relative height not significant 
  #                                                             fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 3.0, a4 = -0.7, b1 = 0.3, b2 = 0.348)))
  #otherDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a8 * terrainRoughness + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed effects not significant 
  #                                                                   fixedFormula = a1 + a8 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                   start = list(fixed = c(a1 = 3.5, a8 = -0.01, a4 = 0, b1 = 0.3, b2 = 0.33)))
  #otherDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2 + b2r), other2022, # a1r, b1r, b2r singularity in backsolve
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = -10000, b1 = 0.0001, b2 = 1.3)))
  #to_fixed_coeffficients(otherDiameterFromHeightMixed$sibbesenReplace)

  saveRDS(otherDiameterFromHeightMixed, "trees/height-diameter/data/other DBH mixed.Rdata")
}


## collect model parameters
if (otherOptions$fitHeight & otherOptions$fitHeightMixed & otherOptions$fitDbh & otherOptions$fitDbhMixed)
{
  if (exists("otherHeightFromDiameter") == FALSE) { readRDS("trees/height-diameter/data/other height.Rds") }
  if (exists("otherHeightFromDiameterMixed") == FALSE) { readRDS("trees/height-diameter/data/other height mixed.Rds") }
  if (exists("otherDiameterFromHeight") == FALSE) { readRDS("trees/height-diameter/data/other DBH.Rds") }
  if (exists("otherDiameterFromHeightMixed") == FALSE) { readRDS("trees/height-diameter/data/other DBH mixed.Rds") }
  
  otherCoefficients = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, get_list_coefficients)),
                                          bind_rows(lapply(otherHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed"))) %>%
                                  mutate(responseVariable = "height"),
                                bind_rows(bind_rows(lapply(otherDiameterFromHeight, get_list_coefficients)),
                                          bind_rows(lapply(otherDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed"))) %>%
                                  mutate(responseVariable = "DBH")) %>%
    mutate(species = "other")
  otherResults = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, get_list_stats)),
                                     bind_rows(lapply(otherHeightFromDiameterMixed, get_list_stats, fitSet = "mixed"))) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(bind_rows(lapply(otherDiameterFromHeight, get_list_stats)),
                                     bind_rows(lapply(otherDiameterFromHeightMixed, get_list_stats, fitSet = "mixed"))) %>%
                             mutate(responseVariable = "DBH")) %>%
    mutate(species = "other")
  
  check_plot_results(otherResults)
  saveRDS(otherCoefficients, "trees/height-diameter/data/other coefficients.Rdata")
  saveRDS(otherResults, "trees/height-diameter/data/other results.Rdata")
}

## preferred forms identified (results.R, Figure 9)
if (otherOptions$fitHeight & otherOptions$fitDbh)
{
  otherHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022, folds = 1, repetitions = 1))
  #otherHeightFromDiameterPreferred$gamBal = fit_gam("REML GAM BA+L", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15), data = other2022, folds = 1, repetitions = 1)
  #otherHeightFromDiameterPreferred$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 23), data = other2022, folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ s(dbh, elevation, terrainRoughness, relativeDiameter, bs = "ts", k = 22, by = as.factor(isPlantation)), data = other2022, folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$linear = fit_lm("linear", height ~ 0 + dbh + I(isPlantation*dbh), other2022, folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$parabolic = fit_lm("parabolic", height ~ 0 + dbh + I(dbh^2) + I(isPlantation*dbh) + I(isPlantation*dbh^2), other2022, folds = 1, repetitions = 1)
  
  otherDiameterFromHeightPreferred = list(gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 9), data = other2022, folds = 1, repetitions = 1))
  otherDiameterFromHeightPreferred$gamPhysio = fit_gam("REML GAM physio", dbh ~ s(height, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), data = other2022, folds = 1, repetitions = 1)
  otherDiameterFromHeightPreferred$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9), data = other2022, folds = 1, repetitions = 1)
  otherDiameterFromHeightPreferred$naslund = fit_gsl_nls("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + a2*sqrt(height - 1.37)), other2022, start = list(a1 = 2.3, a2 = -0.13), folds = 1, repetitions = 1)
  otherDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)) + I((height - 1.37)^2), other2022, folds = 1, repetitions = 1)
  otherDiameterFromHeightPreferred$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1p * isPlantation)*(height - 1.37)^(b2 + b2p * isPlantation)), other2022, start = list(a1 = 0.67, a1p = 1.84, b1 = 1.70, b1p = -1.26, b2 = -0.063, b2p = 0.30), folds = 1, repetitions = 1)
  otherDiameterFromHeightPreferred$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation + a8 * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), other2022, start = list(a1 = 3.47, a1p = -0.43, a8 = -0.001, b1 = 0.34, b2 = 0.28, b2p = 0.019), folds = 1, repetitions = 1)
  
  saveRDS(otherHeightFromDiameterPreferred, "trees/height-diameter/data/other preferred height models.Rdata")
  saveRDS(otherDiameterFromHeightPreferred, "trees/height-diameter/data/other preferred diameter models.Rdata")
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  # height
  #detach(package:gam, unload = TRUE)
  #library(mgcv)
  otherHeightGam = gam(height ~ I(dbh^0.697) + s(dbh, bs = "ts", k = 8) + 
                                s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 6) +
                                s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 8) + # not significant (p = 0.04)
                                s(elevation, bs = "ts", k = 6) + 
                                s(slope, bs = "ts", k = 6) +
                                #s(aspect, bs = "ts", k = 6) + # significant but collapses to a linear effect of a few centimeters
                                s(terrainRoughness, bs = "ts", k = 6) +
                                s(relativeDiameter, bs = "ts", k = 6),
                        data = other2022, method = "REML", select = TRUE, weights = dbhWeight)
  k.check(otherHeightGam)
  summary(otherHeightGam)
  par(mfrow = c(2, 5), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(otherHeightGam, all.terms = TRUE)
  
  otherHeightGam = gam(height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11), data = other2022, family = Gamma, method = "REML", select = TRUE, weights = dbhWeight) # unweighted Gamma , quasi(inverse, mu^2) RMSE  (not as stable as gamma, major outlier), Gaussian RMSE 
  otherHeightGammRidge = gam(height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11) + s(uniquePlotID, bs = "re"), data = other2022, control = gam.control(nthreads = htDiaOptions$rangerThreads), method = "REML", select = TRUE, weights = dbhWeight) # unweighted gamma RMSE  ( if linear in DBH) but predicts only submeter population heights, Gaussian RMSE  but usable for prediction without the random smooths
  sqrt(mean((other2022$height - otherHeightGammRidge$fitted.values)^2))
  otherHeightGammRidgePopulationPrediction = predict(otherHeightGammRidge, other2022, exclude = "s(uniquePlotID)", newdata.guaranteed = TRUE) # drop random effect to make population predictions
  #otherHeightGamm = gamm(height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11), data = other2022, family = Gamma, method = "REML", random = list(uniquePlotID = ~1), weights = dbhWeight) # unweighted Gamma RMSE 
  #sqrt(mean((other2022$height - otherHeightGamm$gam$fitted.values)^2))
  #otherStandIndices = other2022 %>% select(stand) %>% mutate(stand = as.factor(stand), rowNumber = row_number()) %>% group_by(stand) %>% nest() %>% tibble::deframe()
  #otherHeightGammMarkov = gam(height ~ I(dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11) + s(uniquePlotID, bs = "mrf", xt = list(nb = otherStandIndices)), data = other2022, control = gam.control(nthreads = htDiaOptions$rangerThreads), method = "REML", select = TRUE, weights = dbhWeight)
  otherHeightGamDbh0.6 = gam(height ~ I(dbh^0.697) + s(dbh, bs = "ts", k = 8), data = other2022, family = Gamma, method = "REML", select = TRUE, weights = dbhWeight) # unweighted Gaussian RMSE 3.774 @ k = 8, gamma 3.880 @ k = 5
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

  # diameter
  otherDbhGam = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 8) +
                          #s(bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 3) + # not significant
                          s(basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 4),
                          #s(elevation, bs = "ts", k = 3) + # not significant
                          #s(slope, bs = "ts", k = 4), # not significant
                          #s(aspect, bs = "ts", k = 3) + # not significant
                          #s(terrainRoughness, bs = "ts", k = 3), # not significant
                          #s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 3), # not significant
                    data = other2022, folds = 1, repetitions = 1, method = "REML", select = TRUE, , weights = heightWeight)
  k.check(otherDbhGam)
  summary(otherDbhGam)
  par(mfrow = c(1, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(otherDbhGam, scale = 0)

  # mgcv versus gam
  #detach(package:mgcv, unload = TRUE)
  #library(gam)
  #otherHeightGam = gam(height ~ I(dbh^0.6) + s(dbh, df = 4), data = other2022, weights = dbhWeight) # unweighted RMSE 3.900 with linear DBH base, 3.790 with dbh^0.6
  #sqrt(mean(otherHeightGam$residuals^2))
  #summary(otherHeightGam)
}


## random forest regression
if (htDiaOptions$includeInvestigatory)
{
  otherHeightForest = train(height ~ dbh + standBasalAreaPerHectare + basalAreaLarger + elevation + slope + aspect + terrainRoughness + relativeDiameter, data = other2022, method = "ranger", trControl = repeatedCrossValidation, 
                            importance = "impurity_corrected",
                            tuneGrid = expand.grid(mtry = c(6, 8),
                                                   splitrule = "variance",
                                                   min.node.size = c(1, 2)))
  otherHeightForest
  varImp(otherHeightForest)
  
  otherDbhForest = train(dbh ~ height + bootstrapStandBasalAreaPerHectare + basalAreaTaller + elevation + slope + aspect + terrainRoughness + relativeHeight, data = other2022, method = "ranger", trControl = repeatedCrossValidation, 
                         importance = "impurity_corrected",
                         tuneGrid = expand.grid(mtry = c(7, 8),
                                                splitrule = "variance",
                                                min.node.size = c(1, 2)))
  otherDbhForest
  varImp(otherDbhForest)
}
