# load libraries, functions, and acma2022 from setup.R

acmaOptions = tibble(fitHeight = TRUE, 
                     fitHeightMixed = FALSE,
                     fitDbh = FALSE,
                     fitDbhMixed = FALSE)

## bigleaf maple height regressions
if (acmaOptions$fitHeight)
{
  acmaHeightFromDiameter = list(linear = fit_lm("linear", height ~ 0 + dbh + I(isPlantation*dbh), acma2022))
  acmaHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ 0 + dbh + I(dbh^2) + I(isPlantation*dbh) + I(isPlantation*dbh^2), acma2022)
  
  acmaHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, b1 = -0.06, b2 = 1.1)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a2 * basalAreaTaller) * (1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, a2 = 0, b1 = -0.06, b2 = 1.1), significant = FALSE) # a2, a3, b1ba, b1bal, b2ba, b2bal not significant, a2, b1bal significant
  acmaHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a2 * basalAreaTaller + a5 * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 27, a2 = 0, a5 = 0, b1 = -0.06, b2 = 1.1, b2p = -0.2), significant = FALSE) # by propagation
  acmaHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1 + a2 * basalAreaTaller + a4 * relativeDiameter + a5 * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 30, a2 = 0, a4 = 0, a5 = 0, b1 = -0.06, b2 = 1.1, b2p = -0.2), significant = FALSE) # by propagation
  acmaHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a2 * basalAreaTaller + a4 * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 27, a2 = 0, a4 = 0, b1 = -0.027, b2 = 1.08, b2p = -0.2), significant = FALSE) # a2, a4 not significant
  acmaHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp((b1)*dbh))^(b2 + b2p * isPlantation + a5 * topographicWetnessFD8f), acma2022, start = list(a1 = 27, a5 = 0, b1 = -0.06, b2 = 1.1, b2p = -0.2), significant = FALSE) # { a5, b1, b2 } x { s, e, a, tr, tw } not significant
  acmaHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a4 * relativeDiameter)*(1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, a4 = 0, b1 = -0.06, b2 = 1.1), significant = FALSE) # a4, b2p, b2rd not significant, b1rd nan-inf
  acmaHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a4 * relativeDiameter + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, a4 = 0, a5 = 0, b1 = -0.06, b2 = 1.1), significant = FALSE) # by propagation
  acmaHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1 * dbh / (1 + dbh)^b1, acma2022, start = list(a1 = 2.4, b1 = 0.4)) # a1p, b1p not significant
  acmaHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), acma2022, start = list(a1 = 32, a2 = 33, b1 = -1.2)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1*dbh^b2), acma2022, start = list(a1 = 70, b1 = -5.0, b2 = -0.35)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 33, a2 = 30, b1 = 1.1)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1 * dbh^2 + (a2 + a2p * isPlantation)*dbh + a3), acma2022, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19))
  acmaHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, acma2022, start = list(a1 = 2.0, b1 = 0.66)) # a1p, b1p not significant, all data power = 0.666
  acmaHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2)), acma2022, start = list(a1 = 30, b1 = -13, b2 = 3.5)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), acma2022, start = list(Ha = 25, d = 0.7, kU = 0.035)) # Hap, kUp, dp not significant
  acmaHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^b4, acma2022, start = list(a1 = 12.4, b1 = 0.2, b2 = -0.05, b3 = 0.1, b4 = 1.1)) # a1p, b1p, b2p, b3p, b4p not significant
  acmaHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaTaller))^b3*dbh))^b4, acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b4 = 1.1))
  acmaHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaTaller))^b3*dbh))^(b4 + b4tw * topographicWetnessFD8f), acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b4 = 1.1, b4tw = 0), significant = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw } not significant
  acmaHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaTaller))^(b3 + b3rd * relativeDiameter + b3tr * terrainRoughness)*dbh))^b4, acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b3rd = -0.04, b3tr = 0, b4 = 1.1), significant = FALSE) # by propagation, b3rd loses significance
  acmaHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaTaller))^(b3 + b3rd * relativeDiameter)*dbh))^b4, acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b3rd = -0.04, b4 = 1.1)) # a4, b4rd not significant, b2rd nan-inf
  acmaHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^(b4 + b4s * slope), acma2022, start = list(a1 = 12.4, b1 = 0.2, b2 = -0.05, b3 = 0.1, b4 = 1.1, b4s = 0), significant = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw } not significant
  acmaHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^b4, acma2022, start = list(a1 = 17, b1 = 0.14, b2 = -0.04, b3 = 0.12, b3rd = -0.04, b4 = 1.2)) # a4, b1rd, b4rd not significant, b2rd nan-inf
  acmaHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^(b3 + b3rd * relativeDiameter + b3as * sin(3.14159/180 * aspect))*dbh))^b4, acma2022, start = list(a1 = 17, b1 = 0.14, b2 = -0.04, b3 = 0.12, b3as = 0, b3rd = -0.04, b4 = 1.2), significant = FALSE) # by propagation
  acmaHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, acma2022, start = list(a1 = 17, b1 = 0.1, b2 = -0.05, b3 = -0.1, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant
  acmaHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^b3*dbh))^b4, acma2022, start = list(a1 = 10, b1 = 0.3, b2 = -0.06, b2ba = 0.0005, b3 = 0.04, b4 = 1.1)) # a2, b1bal, b2bal, b3ba, b3bal, b4ba, b4bal not significant, { a3, b1ba, b2ba } all significant
  acmaHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1*dbh^b2), acma2022, start = list(a1 = 0.7, b1 = 1.5, b2 = -0.16)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + (a1)*(1 - exp(b1*dbh^b2)), acma2022, start = list(a1 = 27, b1 = -0.04, b2 = 1.0)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a2 * basalAreaTaller) * (1 - exp(b1*dbh^b2)), acma2022, start = list(a1 = 27, a2 = 0, b1 = -0.04, b2 = 1.0), significant = FALSE) # a2, a3, b2ba, b2bal, b3ba, b3bal not significant but a2 significant (p < 0.01) in ~80% of fits
  #print(to_parameter_confidence_intervals(acmaHeightFromDiameter$weibullBal), n = 12)
  #to_fixed_coeffficients(acmaHeightFromDiameter$weibullBal)
  
  acmaHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I(dbh^0.666) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022) # 10x10 MAE 2.51, AIC 1259 (810-2019), gamma NaN
  #acmaHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(dbh^0.666) + s(dbh, standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022, significant = FALSE) # 10x10 AIC 1267 (746-1957)
  #acmaHeightFromDiameter$gamBal = fit_gam("REML GAM BAL", height ~ I(dbh^0.666) + s(dbh, basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 10), data = acma2022) # 10x10 AIC 1268 (682-2066), gamma NaN
  acmaHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(dbh^0.666) + s(dbh, standBasalAreaPerHectare, basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 11), data = acma2022, significant = FALSE) # k = 11 minimum, 10x10 AIC 1300 (822-1993)
  acmaHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(dbh^0.666) + s(dbh, basalAreaTaller, slope, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022, significant = FALSE) # by propagation
  acmaHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I(dbh^0.666) + s(dbh, basalAreaTaller, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  acmaHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I(dbh^0.666) + s(dbh, standBasalAreaPerHectare, basalAreaTaller, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), data = acma2022, significant = FALSE) # k = 16 is minimum, 10x10 AIC 1315 (798-2016)
  # slight MAE reductions with slope, sin(aspect), roughness, and FD8f, AIC neutral -> drop sin(aspect) on highest AIC + k = 57 > ~30 -> drop wetness on AIC + k = 16 > 10
  #acmaHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I(dbh^0.666) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022) # 10x10 AIC 1242 (781-1912)
  #acmaHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I(dbh^0.666) + s(dbh, slope, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1268 (814-1938)
  #acmaHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I(dbh^0.666) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 10), data = acma2022) # 10x10 AIC 1263 (806-1880)
  #acmaHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I(dbh^0.666) + s(dbh, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1273 (753-2007)
  #acmaHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I(dbh^0.666) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1264 (742-1984)
  #acmaHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I(dbh^0.666) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1264 (799-1985)
  acmaHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(dbh^0.666) + s(dbh, slope, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 11), data = acma2022) # k = 11 minimum
  acmaHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(dbh^0.666) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022) # 10x10 AIC 1258 (785-2087)
  acmaHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(dbh^0.666) + s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", k = 16, by = as.factor(isPlantation)), data = acma2022, significant = FALSE) # k = 16 minimum, 10x10 AIC 1485 (1028-2767)
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, summary)
  #acmaHeightFromDiameter$gamRelDbhPhysio$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(acmaHeightFromDiameter, "trees/height-diameter/data/ACMA3 height.Rds")
}

if (acmaOptions$fitHeightMixed)
{
  acmaHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1r)*(1 - exp(b1*dbh))^b2, acma2022, # b1r, b2r singularity in backsolve
                                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 24, b1 = -0.06, b2 = 1.1))))
  acmaHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1r + a2 * basalAreaTaller) * (1 - exp(b1*dbh))^b2, acma2022, # a2 not significant
                                                            fixedFormula = a1 + a2+ b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 24, a2 = 0, b1 = -0.06, b2 = 1.1)), significant = FALSE)
  #acmaHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + a2 * basalAreaTaller + a5 * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022,
  #                                                                fixedFormula = a1 + a2 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 30, a2 = 0.06, a5 = -0.1, b1 = -0.03, b2 = 1.1, b2p = -0.19)))
  #acmaHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, 
  #                                                             fixedFormula = a1 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 30, a5 = -7, b1 = -0.034, b2 = 1.1, b2p = -0.28)))
  acmaHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1) * dbh / (1 + dbh)^(b1 + b1r), acma2022, # a1r and b1r both viable, b1r maybe more tractable
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.5, b1 = 0.25)))
  acmaHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1r) / (1 + a2 * dbh^b1), acma2022, # a2r, b1r max iterations
                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 27, b1 = 30, b2 = -1.3)))
  #acmaHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp(b1*dbh^b2), acma2022, # a1r max iterations, b1r, b2r step halving
  #                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                            start = list(fixed = c(a1 = 70, b1 = -5.0, b2 = -0.35)), control = nlmeControl(maxIter = 500))
  acmaHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1r)*dbh^b1 / (a2 + dbh^b1), acma2022, # a1r, a2r NaN, b1r max iterations
                                                         fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 33, a2 = 30,  b1 = 1.3)))
  acmaHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1r) * dbh^2 + a2*dbh + a3), acma2022, # a2r max iterations, a3r also viable (a1r + a3r singular)
                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.024, a2 = 1.27, a3 = -0.19)))
  acmaHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1r)*dbh^b1, acma2022, # b1r also viable (a1r + b1r nlminb fails to converge)
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 1.3, b1 = 0.8)))
  acmaHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2 + b2r)), acma2022, # a1r, b1r, b2r all tractable
                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 28, b1 = -12, b2 = 2.8)))
  acmaHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Har) * (1 + ((1.37/(Ha + Har))^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), acma2022, # kUr singularity in backsolve, dr step halving
                                                   fixedFormula = Ha + d + kU ~ 1, randomFormula = Har ~ 1|stand/plot,
                                                   start = list(fixed = c(Ha = 24, d = 0.7, kU = 0.04)))
  acmaHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1r)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b3r max iterations
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 12, b1 = 0.25, b2 = -0.04, b3 = 0.1, b4 = 1.1)))
  acmaHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaTaller))^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b3r max iterations
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 11, b1 = 0.25, b2 = -0.04, b3 = 0.1, b4 = 1.0)))
  #acmaHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1r + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaTaller))^b3*dbh))^b4, acma2022, 
  #                                                             fixedFormula = a1 + a6 + b1 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 12, a6 = 0.4, b1 = 0.2, b1p = 0.04, b2 = -0.037, b3 = -0.046, b4 = 1.06)))
  #acmaHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1r + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, acma2022, 
  #                                                          fixedFormula = a1 + a6 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 11, a6 = 0.3, b1 = 0.2, b1p = 0.06, b2 = -0.03, b3 = 0, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1r)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b3r max iterations
                                                     fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 17, b1 = 0.1, b2 = -0.04, b3 = 0.6, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1r) * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b3r tractable but iterations may be an issue
                                                        fixedFormula = a1 + b1 + b2 + b2ba + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 10, b1 = 0.3, b2 = -0.05, b2ba = 0.0005, b3 = 0.04, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^(b1*dbh^(b2 + b2r)), acma2022, # a1r singularity in backsolve, b1r also tractable
                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.5, b1 =2.0, b2 = -0.18)))
  acmaHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1r)*(1 - exp(b1*dbh^b2)), acma2022, # b1r step halving, b2r max iterations
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 24, b1 = -0.05, b2 = 1.1)))
  acmaHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1r + a2 * basalAreaTaller) * (1 - exp(b1*dbh^b2)), acma2022, # a2 not significant
                                                    fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 24, a2 = 0, b1 = -0.05, b2 = 1.1)), significant = FALSE)
  #to_fixed_coeffficients(acmaHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(acmaHeightFromDiameterMixed$weibullBal), n = 16)
  
  saveRDS(acmaHeightFromDiameterMixed, "trees/height-diameter/data/ACMA3 height mixed.Rds")
}


## bigleaf maple diameter regressions
if (acmaOptions$fitDbh)
{
  acmaDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37)), acma2022) # plantation not significant
  acmaDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I((height - 1.37)^2), acma2022) # plantation, plantation^2 not significant

  acmaDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, start = list(a1 = 40, b1 = 0.02, b2 = 1.1)) # a1p, b1p, b2p not significant, a1-b1 evaporation
  #acmaDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, start = list(a1 = 50, a2 = 15, b1 = 0.001, b2 = 1.0))
  #acmaDiameterFromHeight$chapmanReplaceAbatRelHt = fit_gsl_nls("Chapman-Richards replace ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(height - 1.37)^b2) - 1), acma2022, start = list(a1 = 50, a2 = -15, a3 = 4, a4 = -200, b1 = 0.001, b2 = 1.00), significant = FALSE) # by propagation
  acmaDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ a1*(exp((b1 + b1rh * relativeHeight)*(height - 1.37)^b2) - 1), acma2022, start = list(a1 = 50, b1 = 0.02, b1p = 0, b2 = 1.1), significant = FALSE) # a4, b1rh, b2rh not significant
  acmaDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), acma2022, start = list(a1 = 40, b1 = -0.02, b2 = 1.3)) # a1p, b1p, b2p not significant, a1-b1 evaporation
  #acmaDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 40, a2 = 0, b1 = -0.021, b2 = 1.6))
  acmaDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), acma2022, start = list(a1 = 40, b1 = -0.023, b2 = 1.6, b2tw = 0), significant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant
  acmaDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2rh * relativeHeight), 0.9999)), acma2022, start = list(a1 = 40, b1 = -0.02, b2 = 1.3, b2rh = 0), signficant = FALSE) # a4, b1rh, b2rh not significant
  acmaDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), acma2022, start = list(a1 = -116, a2 = -128, b1 = 1.30)) # a1p, a2p, b1p not significant
  acmaDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), acma2022, start = list(a1 = 3, a1p = -1, a2 = -0.12, a2p = -0.02))
  acmaDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, acma2022, start = list(a1 = 0.8, b1 = 1.2)) # a1p, b1p not significant
  #acmaDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, start = list(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272))
  acmaDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ a1*(height - 1.37)^(b1 + b1tw * topographicWetnessFD8f), acma2022, start = list(a1 = 0.8, b1 = 1.1, b1tw = 0.02)) # { a5, b1 } x { e, s, a, tr }, a5tw not significant
  acmaDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ a1*(height - 1.37)^(b1 + b1rh * relativeHeight), acma2022, start = list(a1 = 0.8, b1 = 1.0, b1rh = 0), significant = FALSE) # a4, b1rh not significant
  acmaDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0), significant = FALSE) # a1p, b1p, b2, b2p not significant -> collapses to power form
  #acmaDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.3, a2 = -0.012, b1 = 1.45, b1p = -0.1, b2 = -0.03))
  #acmaDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.5, a3 = -0.007, a7 = -0.05, b1 = 1.5, b1p = -0.11, b2 = -0.03))
  #acmaDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a7 * sin(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.4, a3 = -0.007, a7 = -0.06, a4 = -0.2, b1 = 1.45, b1p = -0.1, b2 = 0.028))
  #acmaDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.3, a2 = -0.01, a4 = 0, b1 = 1.52, b1p = -0.09, b2 = -0.033))
  acmaDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2e * elevation) * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0, b2e = 0), significant = FALSE) # { a5, b1 } x { e, s, a, tr, tw }, b2 not significant
  acmaDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0, b2rh = 0), significant = FALSE) # a4, b1rh, b2rh not significant
  acmaDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2e * elevation + b2rh * relativeHeight) * (height - 1.37)), acma2022, start = list(a1 = 0.9, a4 = 0, b1 = 1.3, b2 = 0, b2e = 0, b2rh = 0), significant = FALSE) # by propagation
  acmaDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), acma2022, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.3, Ha = 161)) # a1p, b1p, Hap not significant, Ha-a1 evaporation increases with Levenberg
  acmaDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^b4, acma2022, start = list(a1 = 0.23, b1 = 1.3, b2 = 0.00001, b3 = 0.7, b4 = -0.1)) # a1p, b1p, b2p, b3p, b4p all NaN-inf, a1-b2 evaporation
  acmaDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b2 = 0), significant = FALSE) # a1p, b1p, b2, b2p not significant
  #acmaDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.5, a2 = 0, b1 = 2.2, b2 = -0.13))
  #acmaDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a2 * basalAreaTaller + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022, start = list(a1 = 0.4, a2 = -0.002, a5 = -0.15, b1 = 3.2, b2 = -0.18, b2p = -0.02))
  #acmaDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a5 * sin(3.14159/180 * slope) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.5, a3 = 0, a5 = -0.2, a4 = -0.1, b1 = 2.6, b2 = -0.14))
  #acmaDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.8, a2 = -0.003, a4 = -0.4, b1 = 2.1, b2 = -0.10))
  acmaDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ a1*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b1tw = 0.02, b2 = 0), significant = FALSE) # b2 not significant
  acmaDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), acma2022, start = list(a1 = 1.3, b1 = 1, b2 = 0, b2rh = 0), significant = FALSE) # a4, b1rh, b2rh not significant
  acmaDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ a1*(height - 1.37)^((b1 + b1rh * relativeHeight + b1tw * topographicWetnessFD8f)*(height - 1.37)^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b1rh = 0, b1tw = 0.02, b2 = 0), significant = FALSE) # b1rh not significant, also by propagation
  acmaDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, acma2022, start = list(a1 = -30, b1 = 0.07, b2 = 0.63)) # a1p, b1p, b2p not significant
  #print(to_parameter_confidence_intervals(acmaDiameterFromHeight$weibull), n = 12)
  #to_fixed_coeffficients(acmaDiameterFromHeight$weibull)
  
  # individual term selection: height + AAT + RelHt by = isPlantation, slope + elevation + cos(aspect) retained by AIC but not significant (p > 0.06)
  acmaDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(height^1.263) + s(height, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022) # 10x10 median AIC 1665, RMSE 10.8
  #acmaDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(height^1.263) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 13), data = acma2022)
  #acmaDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(height^1.263) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, slope, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  #acmaDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ I(height^1.263) + s(height, bootstrapStandBasalAreaPerHectare, basalAreaTaller, slope, terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  # modest RMSE reductions with sin(aspect), cos(aspect), roughness, wetness -> drop elevation, slope -> drop cos(aspect) on k = 56 > ~23
  #acmaDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I(height^1.263) + s(height, elevation, bs = "ts", k = 4), data = acma2022, significant = FALSE) # typically collapses to linear, plantation not significant
  #acmaDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I(height^1.263) + s(height, slope, bs = "ts", k = 8), data = acma2022, significant = FALSE) # plantation not significant, 10x10 AIC 1694
  #acmaDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I(height^1.263) + s(height, sin(3.14159/180 * aspect), bs = "ts", k = 6), data = acma2022) # plantation not significant, 10x10 AIC 1631
  #acmaDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I(height^1.263) + s(height, cos(3.14159/180 * aspect), bs = "ts", k = 7), data = acma2022) # plantation not significant, 10x10 AIC 1694
  #acmaDiameterFromHeight$gamRoughness = fit_gam("REML GAM physio", dbh ~ I(height^1.263) + s(height, terrainRoughness, bs = "ts", k = 7), data = acma2022) # plantation not significant, 10x10 AIC 1626
  #acmaDiameterFromHeight$gamWetness = fit_gam("REML GAM physio", dbh ~ I(height^1.263) + s(height, topographicWetnessFD8f, bs = "ts", k = 8), data = acma2022) # plantation not significant, 10x10 AIC 1676
  acmaDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(height^1.263) + s(height, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 16), data = acma2022, threads = 8) # k = 16 minimum
  acmaDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(height^1.263) + s(height, relativeHeight, bs = "ts", k = 4), data = acma2022) # plantation not significant
  acmaDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(height^1.263) + s(height, relativeHeight, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022, significant = FALSE) # k = 57 minimum > ~22
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, summary)

  saveRDS(acmaDiameterFromHeight, "trees/height-diameter/data/ACMA3 DBH.Rds")
}

if (acmaOptions$fitDbhMixed)
{
  #acmaDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ a1 * (exp(b1*(height - 1.37)) - 1)^(b2 + b2r), acma2022, # a1r, b1r, b2r step halving
  #                                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 40, b1 = 0.001, b2 = 1.1))))
  #acmaDiameterFromHeightMixed = list(chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, 
  #                                                                 fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 50, a2 = 15, b1 = 0.001, b2 = 1.0))))
  #acmaDiameterFromHeightMixed = list(chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller) * (exp(b1*(height - 1.37)^b2) - 1), acma2022,
  #                                                                fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 50, a2 = -15, b1 = 0.001, b2 = 1.0))))
  #acmaDiameterFromHeightMixed = list(chapmanReplaceAbatRelHt = fit_nlme("Chapman-Richards replace ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(height - 1.37)^b2) - 1), acma2022, 
  #                                                                     fixedFormula = a1 + a2 + a3 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                     start = list(fixed = c(a1 = 50, a2 = -15, a3 = 4, a4 = -200, b1 = 0.001, b2 = 1.00))))
  #acmaDiameterFromHeightMixed = list(chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(exp(b1*(height - 1.37)^b2) - 1), acma2022, # relative height not significant
  #                                                                  fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 50, a4 = -15, b1 = 0.1, b2 = 1.0))))
  #acmaDiameterFromHeightMixed = list(chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ a1 * log(1 - pmin(b1*(height - 1.37)^(b2 + b1r), 0.9999)), acma2022, # a1r, b1r, b2r step halving
  #                                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = 40, b1 = -0.02, b2 = 1.3))))
  #acmaDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, 
  #                                                           fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 40, a2 = 0, b1 = -0.021, b2 = 1.6)))
  #acmaDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1r + a8 * terrainRoughness)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, # physio not significant
  #                                                             fixedFormula = a1 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 40, a8 = 0, b1 = -0.023, b2 = 1.6)))
  #acmaDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, # relative height not significant
  #                                                            fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 80, a4 = -30, b1 = -0.012, b2 = 1.6)))
  #acmaDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1r) / (a2 - (height - 1.37)^(b1 + b1r)), acma2022, # a1r, b1r, b2r step halving
  #                                                              fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = -116, a2 = -128, b1 = 1.3)))
  acmaDiameterFromHeightMixed = list(naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1p * isPlantation + a1r) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), acma2022, # a2r max iterations
                                                        fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 3, a1p = -1, a2 = -0.12, a2p = -0.02))))  acmaDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, acma2022, start = list(a1 = 0.8, b1 = 1.2)) # a1p, b1p not significant
  acmaDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ a1*(height - 1.37)^(b1 + b1r), acma2022, # a1r max iterations
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 0.8, b1 = 1.2)))
  #acmaDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, 
  #                                                 fixedFormula = a1 + a1p + a2 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272)))
  acmaDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ a1*(height - 1.37)^(b1 + b1r + b1tw * topographicWetnessFD8f), acma2022, # a1r viable but noticeably slower to converge than b1r
                                                     fixedFormula = a1 + b1 + b1tw ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 0.8, b1 = 1.1, b1tw = 0.02)))
  #acmaDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^b1, acma2022, # relative height not significant
  #                                                  fixedFormula = a1 + a4 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 2.5, a4 = -0.6, b1 = 1.0)))
  #acmaDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1 + a1r) * (height - 1.37)^b1 * exp(b2 * (height - 1.37)), acma2022, # b2 not significant -> collapses to power
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 1.3, b1 = 1.45, b1p = -0.3, b2 = -0.033, b2p = 0.03)))
  #acmaDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
  #                                                 fixedFormula = a1 + a2 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 1.3, a2 = -0.012, b1 = 1.45, b1p = -0.1, b2 = -0.03)))
  #acmaDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1r + a3 * bootstrapStandBasalAreaPerHectare + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
  #                                                       fixedFormula = a1 + a3 + a7 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1.5, a3 = -0.007, a7 = -0.05, b1 = 1.5, b1p = -0.11, b2 = -0.03)))
  #acmaDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark ABA+T RelHt physio", dbh ~ (a1 + a1r + a3 * bootstrapStandBasalAreaPerHectare + a7 * sin(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022,
  #                                                            fixedFormula = a1 + a3 + a7 + a4 + b1 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 1.4, a3 = -0.007, a7 = -0.06, a4 = -0.2, b1 = 1.45, b1p = -0.1, b2 = 0.028)))
  #acmaDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
  #                                                      fixedFormula = a1 + a2 + a4 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 1.3, a2 = -0.01, a4 = 0, b1 = 1.52, b1p = -0.09, b2 = -0.033)))
  #acmaDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1r + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, # fixed effects not significant
  #                                                   fixedFormula = a1 + a7 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1.2, a7 = -0.055, b1 = 1.45, b1p = -0.064, b2 = 0.03)))
  #acmaDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, # fixed effects not significant
  #                                                  fixedFormula = a1 + a4 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 1.4, a4 = 0, b1 = 1.45, b1p = -0.3, b2 = -0.03)))
  #acmaDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a7 * sin(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, # fixed effects not significant
  #                                                        fixedFormula = a1 + a7 + a4 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 1.2, a7 = -0.06, a4 = 0, b1 = 1.45, b1p = -0.064, b2 = 0.03)))
  #acmaDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/(a1 + a1r) * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), acma2022, # a1r, b1r, Har step halving
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.000003, a2 = 0.002, b1 = 1.3, Ha = 161)))
  #acmaDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1 * (height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^(b4 + b4r), acma2022, # a1r, b1r, b2r, b3r, b4r singularity in backsolve
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b4r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 0.23, b1 = 1.14, b2 = 0.00001, b3 = 0.7, b4 = -0.1)))
  #acmaDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1r) * (height - 1.37)^(b1*(height - 1.37)^b2), acma2022, # b2 not significant -> collapses to power
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1.3, b1 = 1, b2 = 0)))
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
  #                                                           fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 0.5, a2 = 0, b1 = 2.2, b2 = -0.13)))
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022,
  #                                                                 fixedFormula = a1 + a2 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 0.4, a2 = -0.002, a5 = -0.15, b1 = 3.2, b2 = -0.18, b2p = -0.02)))
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1r + a3 * bootstrapStandBasalAreaPerHectare + a5 * sin(3.14159/180 * slope) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
  #                                                                      fixedFormula = a1 + a3 + a5 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                      start = list(fixed = c(a1 = 0.5, a3 = 0, a5 = -0.2, a4 = -0.1, b1 = 2.6, b2 = -0.14)))
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
  #                                                                fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 0.8, a2 = -0.003, a4 = -0.4, b1 = 2.1, b2 = -0.10)))
  #acmaDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022, # fixed effects not significant
  #                                                             fixedFormula = a1 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 0.5, a5 = -0.2, b1 = 2.6, b2 = -0.18, b2p = 0.02)))
  #acmaDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, # fixed effects not significant
  #                                                            fixedFormula = a1 + a4 +b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 0.65, a4 = -0.25, b1 = 2.0, b2 = -0.1)))
  #acmaDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a5 * sin(3.14159/180 * slope) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, # fixed effects not significant
  #                                                                  fixedFormula = a1 + a5 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 0.5, a5 = -0.2, a4 = -0.1, b1 = 2.6, b2 = -0.14)))
  acmaDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1 * log(1 - pmin((b1 + b1r)*(height - 1.37), 0.9999)))^b2, acma2022, # a1r singularity in backsolve, b2r step halving
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = -30, b1 = 0.07, b2 = 0.63)))

  saveRDS(acmaDiameterFromHeightMixed, "trees/height-diameter/data/ACMA3 DBH mixed.Rds")
}


## collect model parameters
if (acmaOptions$fitHeight & acmaOptions$fitHeightMixed & acmaOptions$fitDbh & acmaOptions$fitDbhMixed)
{
  if (exists("acmaHeightFromDiameter") == FALSE) { readRDS("trees/height-diameter/data/ACMA3 height.Rds") }
  if (exists("acmaHeightFromDiameterMixed") == FALSE) { readRDS("trees/height-diameter/data/ACMA3 height mixed.Rds") }
  if (exists("acmaDiameterFromHeight") == FALSE) { readRDS("trees/height-diameter/data/ACMA3 dbh.Rds") }
  if (exists("acmaDiameterFromHeightMixed") == FALSE) { readRDS("trees/height-diameter/data/ACMA3 dbh mixed.Rds") }
  
  acmaCoefficients = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, get_list_coefficients)),
                                         bind_rows(lapply(acmaHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed"))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(acmaDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(acmaDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed"))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "ACMA3")
  acmaResults = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, get_list_stats)),
                                    bind_rows(lapply(acmaHeightFromDiameterMixed, get_list_stats, fitSet = "mixed"))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(acmaDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Chapman-Richards replace ABA+T", fitSet = "primary", fittingMethod = "gsl_nls"),
                                    bind_rows(lapply(acmaDiameterFromHeightMixed, get_list_stats, fitSet = "mixed"))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "ACMA3")
  
  check_plot_results(acmaResults)
  saveRDS(acmaCoefficients, "trees/height-diameter/data/ACMA3 coefficients.Rds")
  saveRDS(acmaResults, "trees/height-diameter/data/ACMA3 results.Rds")
}

## preferred forms identified (results.R, Figure 7)
if (acmaOptions$fitHeight & acmaOptions$fitDbh)
{
  acmaHeightFromDiameterPreferred = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1*dbh))^(b2 + b2p*isPlantation), acma2022, start = list(a1 = 27, b1 = -0.03, b2 = 1.1, b2p = -0.2), folds = 1, repetitions = 1))
  #acmaHeightFromDiameterPreferred$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + (a1 + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect)) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 32.5, a5 = -8.37, a7 = 1.557, b1 = -0.031, b2 = 1.01, b2p = -0.1-4), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), acma2022, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986), folds = 1, repetitions = 1)
  #acmaHeightFromDiameterPreferred$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1 * dbh^2 + (a2 + a2p * isPlantation)*dbh + a3), acma2022, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1/(dbh + b2 + b2p * isPlantation)), acma2022, start = list(a1 = 31.8, a1p = 3.00, b1 = -25.7, b2 = 6.70, b2p = 0.62), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), acma2022, start = list(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/bootstrapStandBasalAreaPerHectare)^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 19, b1 = 0.1, b2 = -0.03, b3 = 0.06, b4 = 1.1, b4p = -0.25), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*dbh^b2)), acma2022, start = list(a1 = 26, a2 = 0.1, b1 = -0.027, b1p = -0.01, b2 = 0.99), folds = 1, repetitions = 1)
  #AIC(acmaHeightFromDiameterPreferred$chapmanRichards, acmaHeightFromDiameterPreferred$michaelisMenten, acmaHeightFromDiameterPreferred$prodan, acmaHeightFromDiameterPreferred$richardsW)
  
  acmaDiameterFromHeightPreferred = list(chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a8 * terrainRoughness)*log(1 - pmin((b1 + b1p * isPlantation)*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 20, a8 = 0.2, b1 = -0.03, b1p = 0.01, b2 = 1.8), folds = 1, repetitions = 1))
  acmaDiameterFromHeightPreferred$gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022, folds = 1, repetitions = 1)
  #acmaDiameterFromHeightPreferred$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ s(height, bootstrapStandBasalAreaPerHectare, basalAreaTaller, slope, terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022, folds = 1, repetitions = 1)
  #acmaDiameterFromHeightPreferred$power = fit_gsl_nls("power", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282), folds = 1, repetitions = 1)
  #acmaDiameterFromHeightPreferred$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*(height - 1.37)^b1, acma2022, start = list(a1 = 3.7, a1p = -0.9, a5 = -1.2, a8 = 0.02, b1 = 0.9), folds = 1, repetitions = 1)
  acmaDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), acma2022, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), folds = 1, repetitions = 1)
  acmaDiameterFromHeightPreferred$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a3 * bootstrapStandBasalAreaPerHectare + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.5, a3 = -0.007, a7 = -0.05, b1 = 1.5, b1p = -0.11, b2 = -0.03), folds = 1, repetitions = 1)
  
  saveRDS(acmaHeightFromDiameterPreferred, "trees/height-diameter/data/ACMA3 preferred height models.Rds")
  saveRDS(acmaDiameterFromHeightPreferred, "trees/height-diameter/data/ACMA3 preferred diameter models.Rds")
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  acmaHeightGam = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 9) + 
                               s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 5) + 
                               #s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 3) + # not significant
                               s(elevation, bs = "ts", k = 5) + 
                               s(slope, bs = "ts", k = 4) +
                               #s(aspect, bs = "ts", k = 3) + # not significant
                               #s(terrainRoughness, bs = "ts", k = 4) + # not significant
                               s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 4), 
                      data = acma2016physio, method = "REML", select = TRUE, weights = dbhWeight)
  k.check(acmaHeightGam)
  summary(acmaHeightGam)
  par(mfrow = c(2, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(acmaHeightGam, scale = 0)
  
  acmaDbhGam = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 9) +
                         s(bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 4),
                         #s(basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 3) + # not significant
                         #s(elevation, bs = "ts", k = 3) + # not significant
                         #s(slope, bs = "ts", k = 4), # not significant
                         #s(aspect, bs = "ts", k = 3), # not significant
                         #s(terrainRoughness, bs = "ts", k = 3), # not significant
                         #s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 3), # not significant
                    data = acma2016physio, method = "REML", select = TRUE, weights = heightWeight)
  k.check(acmaDbhGam)
  summary(acmaDbhGam)
  par(mfrow = c(1, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(acmaDbhGam, scale = 0)
}


## random forest regression
if (htDiaOptions$includeInvestigatory)
{
  acmaHeightForest = train(height ~ dbh + standBasalAreaPerHectare + basalAreaLarger + elevation + slope + aspect + terrainRoughness + relativeDiameter, data = acma2022, method = "ranger", trControl = repeatedCrossValidation, 
                           importance = "impurity_corrected",
                           tuneGrid = expand.grid(mtry = c(4, 5, 6),
                                                  splitrule = "variance",
                                                  min.node.size = c(1, 2, 3)))
  acmaHeightForest
  varImp(acmaHeightForest)
  
  acmaDbhForest = train(dbh ~ height + bootstrapStandBasalAreaPerHectare + basalAreaTaller + elevation + slope + aspect + terrainRoughness + relativeHeight, data = acma2022, method = "ranger", trControl = repeatedCrossValidation, 
                        importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = c(6, 7, 8),
                                               splitrule = "variance",
                                               min.node.size = c(1, 2, 3)))
  acmaDbhForest
  varImp(acmaDbhForest)
}
