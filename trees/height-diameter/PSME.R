# load libraries, functions, and psme2022 from setup.R

psmeOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightMixed = FALSE,
                     fitDbhPrimary = FALSE,
                     fitDbhMixed = FALSE)

## Douglas-fir height regressions
if (psmeOptions$fitHeightPrimary)
{
  # linear regressions
  psmeHeightFromDiameter = list(linear = fit_lm("linear", height ~ 0 + dbh + I(isPlantation*dbh), psme2022))
  psmeHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ 0 + dbh + I(dbh^2) + I(isPlantation*dbh^2), psme2022) # plantation, plantation^2 not mutually significant
  
  # nonlinear regressions
  psmeHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1 * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2p = 0.1)) # a1p not significant
  psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + b2p * isPlantation + (b2bal + b2balp * isPlantation) * basalAreaLarger), psme2022, start = list(a1 = 60, b1 = -0.02, b2 = 1.3, b2p = 0.1, b2bal = -0.01, b2balp = -0.005)) # a3, a3p, b1ba, b1bap, b2ba, b2bap not significant, a2, a2p, b1bal also significant
  psmeHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + a1 * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^(b2 + b2p * isPlantation + (b2bal + b2balp * isPlantation) * basalAreaLarger), psme2022, start = list(a1 = 63, b1 = -0.02, b1rd = 0.0005, b2 = 1.2, b2p = 0.06, b2bal = -0.008, b2balp = -0.004)) # a4p, b1rdp, b2rd not significant
  psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + a1 * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation + (b2bal + b2balp * isPlantation) * basalAreaLarger), psme2022, start = list(a1 = 60, b1 = -0.01, b1ac = -0.001, b2 = 1.2, b2p = 0.1, b2bal = -0.008, b2balp = -0.004))
  psmeHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + a1 * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2bal * basalAreaLarger + b2rd * relativeDiameter), psme2022, start = list(a1 = 62, b1 = -0.018, b1ac = -0.001, b2 = 1.3, b2bal = -0.01, b2rd = 0.06)) # b2p, b2balp, b2rdp not significant
  psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 60, b1 = -0.02, b1ac = -0.001, b2 = 1.1, b2p = 0.1)) # { a5, b1, b2 } x { e, s, tr, tw }, { a5, b2 } x a not significant, b1as not reliably significant
  psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psme2022, start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2rd = 0.1, b2rdp = 0.05)) # a4p, b2p not significant, b1rd NaN-inf, a4 also significant
  psmeHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psme2022, start = list(a1 = 58, b1 = -0.02, b1ac = -0.001, b2 = 1.1, b2rd = 0.1, b2rdp = 0.06))
  psmeHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 3, a1p = -1.4, b1 = 0.4, b1p = -0.16))
  psmeHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 73, a2 = 150, b1 = -1.3, b1p = 0.03)) # a1p not significant, a2p-b1p not mutually significant
  psmeHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 200, b1 = -6.8, b1p = -1.2, b2 = -0.3, b2p = -0.04)) # a1p-b1p not mutually significant
  psmeHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psme2022, start = list(a1 = 75, a2 = 140, a2p = 20, b1 = 1.2)) # a1p, b1p not significant
  psmeHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3 + a3p* isPlantation), psme2022, start = list(a1 = 0.012, a2 = 0.9, a3 = 1.4, a3p = 3.5)) # a1p, a2p not significant
  psmeHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation), psme2022, start = list(a1 = 2.5, a1p = -1, b1 = 0.65, b1p = 0.14))
  psmeHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(dbh + b2)), psme2022, start = list(a1 = 73, b1 = -46, b1p = -3, b2 = 10)) # a1p, b2p not significant
  psmeHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-kU * dbh)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2022, start = list(Ha = 55, Hap = -3, d = 0.3, dp = 0.2, kU = 0.013)) # kUp not significant
  psmeHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psme2022, start = list(a1 = 13, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4 = 1.0)) # a1p, b1p, b2p, b3, b3p, b4p not significant
  psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^b4, psme2022, start = list(a1 = 18, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4 = 1.0)) # a1p, b1p, b2p, b3, b3p, b4p not significant
  psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = -0.02, b4 = 1.0, b4e = 0.0001, b4as = -0.035, b4ac = -0.04)) # { a5, b1 } x { e, s, tr, tw }, a5ac, b1as, { b2, b3, b4 } x { s, tr, tw }, b3 not significant, { b2e, b3e}-b4e, { a4a, b1a, b2a, b3a}-b4a not mutually significant -> what is most accurate a placement?
  psmeHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a4 * relativeDiameter)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 18, a4 = 0.4, b1 = 0.3, b2 = -0.02, b3 = -0.05, b4 = 1.0, b4e = 0.0001, b4as = -0.04, b4ac = -0.04))
  psmeHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1 + a4 * relativeDiameter)*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^b4, psme2022, start = list(a1 = 19, a4 = 0.4, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0)) # a4p, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rdp not significant, b1rd, b4rd not reliably significant
  psmeHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)) # { a5, b1, b2, b3, b4 } x { e, s, tr, tw }, { a1, b1, b2, b3 } x as, b3 not significant, { a1ac, b1ac, b2ac, b3ac }-b4ac not mutually significant -> what is most accurate ac placement?
  psmeHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psme2022, start = list(a1 = 13, b1 = 0.3, b1rd = 0, b2 = -0.02, b3 = 0.02, b4 = 1.0), significant = FALSE) # a4, a4p, b1rd, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rd, b4rdp not significant
  psmeHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1rd * relativeDiameter) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 18, b1 = 0.3, b1rd = 0.004, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04), significant = FALSE) # b1rd, b3 not significant
  psmeHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, psme2022, start = list(a1 = 45, a1p = -12, b1 = 0.08, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)) # b2p, b3p, b4p not significant
  psmeHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp((b2 + (b2ba + b2bap * isPlantation) * standBasalAreaPerHectare)*standTreesPerHectare^(b3) * dbh))^(b4 + (b4bal + b4balp * isPlantation) * basalAreaLarger), psme2022, start = list(a1 = 40, a1p = -10, b1 = 0.1, b1p = 0.1, b2 = -0.03, b2ba = 0.0001, b2bap = 0.0001, b3 = -0.05, b4 = 1.2, b4bal = -0.005, b4balp = -0.004)) # a2p, a3p, { b1, b3 } x { bal, balp }, b4bap, b4balp not significant, b2bal, b2balp NaN-inf, b3 not reliably significant
  psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 0.2, b1 = 2.5, b1p = -0.2, b2 = -0.14, b2p = 0.02)) # a1p not significant
  psmeHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psme2022, start = list(a1 = 56, b1 = -0.01, b1p = 0.001, b2 = 1.1)) # a1p, b2p not significant
  psmeHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + a1 * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger))), psme2022, start = list(a1 = 60, b1 = -0.01, b2 = 1.0, b2bal = 0.003)) # a3, a3p, b1ba, b1bap, b1bal, b2ba, b2balp, b2bal, b2balp not significant
  #print(to_parameter_confidence_intervals(psmeHeightFromDiameter$weibullBal), n = 40)
  #to_fixed_coeffficients(psmeHeightFromDiameter$weibullBal)
  
  # GAMs
  psmeHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = psme2022, threads = 4) # 10x10 AIC 3881
  #psmeHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, threads = 4) # 10x10 AIC 3874
  psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 16), data = psme2022, threads = 4) # 10x10 AIC 3839
  #psmeHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 19), data = psme2022, threads = 4) # 10x10 AIC 3854
  psmeHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, significant = FALSE, threads = 4) # 10x10 AIC 3781
  psmeHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 4) # 10x10 AIC 3898
  if (psmeOptions$fitPhysioGams)
  {
    psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, threads = 4)
    psmeHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), terrainRoughness, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, threads = 4)
    # modest reductions in MAE from slope, aspect, roughness, and wetness, including elevation is detrimental
    #psmeHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, elevation, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, threads = 8, significant = FALSE) # 10x10 AIC 4125, support at k ~= 33 but clearly overfit, k = 10 minimum for curvature
    #psmeHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, threads = 8) # 10x10 AIC 3856
    #psmeHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 8) # 10x10 AIC 3937
    #psmeHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 8) # 10x10 AIC3948
    #psmeHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 8) # 10x10 AIC 3940
    #psmeHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 8) # 10x10 AIC 3941
    #psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 85), data = psme2022, threads = 4) # 10x10 AIC 4220 + higher MAE than slope or wetness alone @ k ~= 61 -> arguably overfit, fails on dplyr errors with threads = 8, warnings at threads = 4
    #psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 4) # 10x10 AIC 3945 + higher MAE than slope or wetness alone
    psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 8) # 10x10 AIC 3941
    psmeHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, topographicWetnessFD8f, bs = "ts", k = 15, by = as.factor(isPlantation)), data = psme2022, threads = 4) # 10x10 AIC 3925
    #lapply(psmeHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
    #lapply(psmeHeightFromDiameter$gamRelDbhPhysio$fit, summary)
    #psmeHeightFromDiameter$gamRelDbhPhysio$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
    
    psmeHeightFromDiameterGamBalPhysio = psmeHeightFromDiameter$gamBalPhysio
    psmeHeightFromDiameterGamBalPhysioRelDbh = psmeHeightFromDiameter$gamBalPhysioRelDbh
    psmeHeightFromDiameterGamPhysio = psmeHeightFromDiameter$gamPhysio
    psmeHeightFromDiameterGamRelDbhPhysio = psmeHeightFromDiameter$gamRelDbhPhysio
    saveRDS(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamBalPhysioRelDbh, psmeHeightFromDiameterGamPhysio, psmeHeightFromDiameterGamRelDbhPhysio, "trees/height-diameter/data/PSME height primary GAMs.Rds")
    rm(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamBalPhysioRelDbh, psmeHeightFromDiameterGamPhysio, psmeHeightFromDiameterGamRelDbhPhysio)
  } else {
    readRDS("trees/height-diameter/data/PSME height primary GAMs.Rds")
    psmeHeightFromDiameter$gamBalPhysio = psmeHeightFromDiameterGamBalPhysio
    psmeHeightFromDiameter$gamBalPhysioRelDbh = psmeHeightFromDiameter$gamBalPhysioRelDbh
    psmeHeightFromDiameter$gamPhysio = psmeHeightFromDiameterGamPhysio
    psmeHeightFromDiameter$gamRelDbhPhysio = psmeHeightFromDiameterGamRelDbhPhysio
    rm(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamPhysio)
  }
  saveRDS(psmeHeightFromDiameter, "trees/height-diameter/data/PSME height.Rds")
}


## Douglas-fir height mixed regressions
if (psmeOptions$fitHeightMixed)
{
  psmeHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1r) *(1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, # b1r, b2r singularity in backsolve
                                                                fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 45, b1 = -0.02, b2 = 1.2, b2p = -0.1))))
  psmeHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1r) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation + b2bal * basalAreaLarger), psme2022, # b1r, b2r singularity in backsolve, b2balp not significant
                                                            fixedFormula = a1 + b1 + b2 + b2p + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 45, b1 = -0.02, b2 = 1.3, b2p = -0.1, b2bal = -0.008)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation + b2bal * basalAreaLarger), psme2022, # b1r, b2r singularity in backsolve, b2balp not significant
                                                                  fixedFormula = a1 + b1 + b1ac + b2 + b2p + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 45, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2p = -0.1, b2bal = -0.006)))
  psmeHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), psme2022, 
                                                               fixedFormula = a1 + b1 + b1ac + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 45, b1 = -0.02, b1ac = -0.002, b2 = 1.2, b2p = 0.1)))
  psmeHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation + b1r), psme2022, # a1r also viable but slower to converge
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 3, a1p = -0.8, b1 = 0.4, b1p = -0.1)))
  psmeHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1r) / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psme2022, # a2r, b1r also tractable
                                                  fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 50, a2 = 90, b1 = -1.3, b1p = 0.01)))
  psmeHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, # b1r max iterations, b2r step halving
                                              fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                              start = list(fixed = c(a1 = 90, b1 = -8, b1p = 2, b2 = -0.5, b2p = 0.1)))
  psmeHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + a1r) / (a2 + a2p * isPlantation + dbh^(b1 + a1r)), psme2022, # a1r, a2r also tractable
                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 50, a2 = 100, a2p = -20, b1 = 1.3)))
  psmeHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1r)*dbh^2 + (a2)*dbh + a3), psme2022,  # a2r also tractable, a3r max iterations, a3p not significant
                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.018, a2 = 0.7, a3 = 2.7)))
  psmeHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation + b1r), psme2022, # b1r also tractable
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2.2, a1p = -0.7, b1 = 0.67, b1p = 0.08)))
  psmeHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1r)*exp((b1 + b1p * isPlantation)/(dbh + b2)), psme2022, # b1r, b2r also tractable
                                                   fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 53, b1 = -29, b1p = 3, b2 = 5)))
  psmeHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), psme2022, # kUr max iterations, dr step halving, dp not significant
                                                   fixedFormula = Ha + Hap + d + kU ~ 1, randomFormula = Har ~ 1|stand/plot,
                                                   start = list(fixed = c(Ha = 47, Hap = -5, d = 0.4, kU = 0.018)))
  psmeHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1r)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, psme2022, # a1r, b2r, b4r singularity in backsolve, b3r also tractable, b3 not significant
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^b4, psme2022, # a1r, b2r, b4r singularity in backsolve, b3r also tractable
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 32, b1 = 0.1, b2 = -0.02, b3 = -0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, # a1r, b24, b4r singularity in backsolve, b3r also tractable
                                                               fixedFormula = a1 + b1 + b2 + b3 + b4 + b4e + b4as + b4ac ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 30, b1 = 0.1, b2 = -0.02, b3 = -0.1, b4 = 1.0, b4e = 0.0002, b4as = -0.035, b4ac = -0.04)))
  psmeHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, # a1r, b4r singularity in backsolve, b2r step halving, b3r also tractable, b3 not significant
                                                            fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac + b4as ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  psmeHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1r)*(1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4), psme2022, # a1r, b4r singularity in backsolve, b2r step halving, b3r also tractable
                                                     fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 35, a1p = -10, b1 = 0.1, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation) * standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1r) * (1 - exp(b2 * standTreesPerHectare^(b3)*dbh))^(b4 + (b4bal + b4balp * isPlantation) * basalAreaLarger), psme2022, # a1r, b4r singularity in backsolve, b2r step halving, b2ba, b2bap not significant
                                                        fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 + b4bal + b4balp ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 40, a1p = -10, b1 = 0.1, b1p = 0.1, b2 = -0.03, b3 = -0.05, b4 = 1.2, b4bal = -0.005, b4balp = -0.004)))
  psmeHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + (a1)*dbh^((b1 + b1p * isPlantation + b1r)*dbh^(b2 + b2p * isPlantation)), psme2022, # a1r, b2r also tractable
                                                  fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.2, b1 = 2.3, b1p = -0.2, b2 = -0.14, b2p = 0.02)))
  psmeHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1r)*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psme2022, # b1r, b2r max iterations
                                                 fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 45, b1 = -0.01, b1p = 0.001, b2 = 1.1)))
  psmeHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger + b2r))), psme2022, # b1r step halving, a1r also tractable
                                                    fixedFormula = a1 + b1 + b2 + b2bal ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 50, b1 = -0.01, b2 = 1.0, b2bal = 0.003)))
  #to_fixed_coeffficients(psmeHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(psmeHeightFromDiameterMixed$sharmaZhang), n = 40)

  saveRDS(psmeHeightFromDiameterMixed, "trees/height-diameter/data/PSME height mixed.Rds")
}


## Douglas-fir diameter regressions
if (psmeOptions$fitDbhPrimary)
{
  # linear regressions
  psmeDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)), psme2022))
  psmeDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I((height - 1.37)^2) + I(isPlantation*(height - 1.37)^2), psme2022) # plantation not significant
  
  # nonlinear regressions
  psmeDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation)*(exp(b1*(height - 1.37)) - 1)^b2, psme2022, start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8)) # b1p, b2p not significant
  #psmeDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 75, a1p = -35, a2 = 0.5, a3 = -0.07, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07))
  psmeDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1)*(exp((b1)*(height - 1.37)^(b2 + b2rh * relativeHeight)) - 1), psme2022, start = list(a1 = 45, b1 = 0.04, b2 = 0.8, b2rh = 0.02)) # a1p-b1rh, b2rh not mutually significant
  psmeDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -100, b1 = 0.012, b2 = 0.9)) # a1p, b1p, b2p not significant
  #psmeDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -175, a1p = 100, a2 = 1.2, a3 = 0.14, b1 = 0.007, b1p = 0.0057, b2 = 0.79))
  psmeDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), psme2022, start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8, b2tw = 0), significant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant
  psmeDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), psme2022, start = list(a1 = -90, b1 = 0.012, b2 = 0.9, b2rh = -0.07)) # a4, b1rh also significant
  psmeDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 130, a2 = 70, b1 = 0.9, b1p = -0.008)) # a1p, a2p not significant
  psmeDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psme2022, start = list(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006))
  psmeDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2))
  #psmeDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053))
  psmeDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1e * elevation), psme2022, start = list(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2, a6 = 0), significant = FALSE) # { a5, b1 } x { e, s, a, tr, tw } not significant
  psmeDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight), psme2022, start = list(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)) # a1p-a4 not mutually significant
  psmeDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2p = -0.001)) # a1p not significant, b1p-b2p not mutually significant
  #psmeDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013))
  #psmeDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012))
  #psmeDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.02, a4 = 0.2, b1 = 0.91, b1p = -0.16, b2 = 0, b2p = 0.009))
  #psmeDiameterFromHeight$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark ABA+T RelHt physio", dbh ~ (a1 + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect) + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.4, a2 = -0.013, a3 = -0.0037, a6 = -0.05, a4 = 0.1, b1 = 0.94, b1p = -0.16, b2 = 0.001, b2p = 0.01))
  psmeDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1 + a6 * slope) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2, a6 = 0, b1 = 0.8, b2 = 0.01, b2p = -0.001)) # { a6, b1 } x { e, s, tr, tw } not significant
  print(to_parameter_confidence_intervals(psmeDiameterFromHeight$ruarkPhysio), n = 40)
  to_fixed_coeffficients(psmeDiameterFromHeight$ruarkPhysio)
  psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.2, a4 = 0.37, b1 = 0.85, b1p = -0.11, b2 = 0.003, b2p = 0.007))
  psmeDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a6 * cos(3.14159/180 * aspect) + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.2, a6 = -0.06, a4 = 0.35, b1 = 0.86, b1p = -0.11, b2 = 0.0033, b2p = 0.007))
  #psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1), 0.9999)), psme2022, start = list(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90))
  psmeDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + b3p * isPlantation)*(height - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13))
  psmeDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017))
  #psmeDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, b1 = 0.61, b2 = 0.071, b2p = 0.021))
  #psmeDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077))
  #psmeDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.01, a3 = -0.006, a4 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021))
  #psmeDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.1, a1p = -0.5, a2 = -0.01, a3 = -0.005, a6 = -0.06, a4 = 0.1, b1 = 0.65, b2 = 0.07))
  psmeDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.0, a1p = -0.3, a6 = -0.06, b1 = 0.60, b2 = 0.09))
  psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1p * isPlantation + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.2, a1p = -0.65, a4 = 0.43, b1 = 0.62, b2 = 0.07, b2p = 0.04))
  psmeDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.0, a1p = -0.5, a6 = -0.05, a4 = 0.4, b1 = 0.62, b2 = 0.07))
  psmeDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(height - 1.37), 0.9999)))^b2, psme2022, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81))
  # GAMs
  # individual term selection: height + ABA + AAT + RelHt by = isPlantation + slope + elevation + sin(aspect) + TSI + RelHt
  #   primary effects NSE 0.861
  psmeDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, threads = 2)
  #psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 28), data = psme2022, threads = 2)
  psmeDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 2)
  #to_fixed_coeffficients(psmeDiameterFromHeight$weibullBal)
  #print(to_parameter_confidence_intervals(psmeDiameterFromHeight$sharmaZhang), n = 40)

  if (psmeOptions$fitPhysioGams)
  {
    #psmeDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, threads = 2)
    psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), data = psme2022, threads = 2)
    psmeDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, threads = 2)
    
    #psmeDiameterFromHeightGamAbatPhysio = psmeDiameterFromHeight$gamAbatPhysio
    psmeDiameterFromHeightGamPhysio = psmeDiameterFromHeight$gamPhysio
    psmeDiameterFromHeightGamRelHtPhysio = psmeDiameterFromHeight$gamRelHtPhysio
    saveRDS(psmeDiameterFromHeightGamAbatPhysio, psmeDiameterFromHeightGamPhysio, psmeDiameterFromHeightGamRelHtPhysio, "trees/height-diameter/data/PSME dbh primary GAMs.Rds")
    rm(psmeDiameterFromHeightGamAbatPhysio, psmeDiameterFromHeightGamPhysio, psmeDiameterFromHeightGamRelHtPhysio)
  } else {
    readRDS("trees/height-diameter/data/PSME dbh primary GAMs.Rds")
    psmeDiameterFromHeight$gamAbatPhysio = psmeDiameterFromHeightGamAbatPhysio
    psmeDiameterFromHeight$gamPhysio = psmeDiameterFromHeightGamPhysio
    psmeDiameterFromHeight$gamRelHtPhysio = psmeDiameterFromHeightGamRelHtPhysio
    rm(psmeDiameterFromHeightGamAbatPhysio, psmeDiameterFromHeightGamPhysio, psmeDiameterFromHeightGamRelHtPhysio)
  }
  
  if (psmeOptions$fitAbatRelHtPhysioGam)
  {
    psmeDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ I(height^) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496), data = psme2022, threads = 2)
    psmeDiameterFromHeightGamAbatPhysioRelHt = psmeDiameterFromHeight$gamAbatPhysioRelHt
    saveRDS(psmeDiameterFromHeightGamAbatPhysioRelHt, "trees/height-diameter/data/PSME DBH primary GAM ABA+T RelHt physio.Rds")
    rm(psmeDiameterFromHeightGamAbatPhysioRelHt)
  } else {
    readRDS("trees/height-diameter/data/PSME DBH primary GAM ABA+T RelHt physio.Rds")
    psmeDiameterFromHeight$gamAbatPhysioRelHt = psmeDiameterFromHeightGamAbatPhysioRelHt
    rm(psmeDiameterFromHeightGamAbatPhysioRelHt)
  }
  
  saveRDS(psmeDiameterFromHeight, "trees/height-diameter/data/PSME DBH primary.Rds")
}

if (psmeOptions$fitDbhMixed)
{
  psmeDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1 + a1r) * (exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, 
                                                               fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 40, b1 = 0.03, b1p = -0.005, b2 = 0.61, b2p = 0.15))))
  psmeDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, 
                                                            fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 75, a1p = -35, a2 = 0.5, a3 = -0.07, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07)))
  psmeDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a4 * relativeHeight)*(exp(b1*(height - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2022, 
                                                             fixedFormula = a1 + a1p + a4 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 30, a1p = -10, a4 = 5, b1 = 0.12, b2 = 0.60, b2p = 0.035)))
  psmeDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1r + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, 
                                                         fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = -123, a1p = 53.1, b1 = 0.0085, b1p = 0.0041, b2 = 0.77)))
  psmeDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, 
                                                             fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = -175, a1p = 100, a2 = 1.2, a3 = 0.14, b1 = 0.007, b1p = 0.0057, b2 = 0.79)))
  psmeDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1r + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022,
                                                               fixedFormula = a1 + a1p + a5 + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = -12, a1p = -3.9, a5 = -2.2, a8 = 0.04, b1 = 0.020, b1p = 0.0054, b2 = 0.45)))
  psmeDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2022, 
                                                              fixedFormula = a1 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = -9.0, a4 = -12.6, b1 = 0.013, b1p = 0.009, b2 = 0.004, b2p = 0.31)))
                                                               control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.001, msTol = 1E-5)
  psmeDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1r + a1p * isPlantation) * (height - 1.37)^b1 / (a2 + a2p * isPlantation - (height - 1.37)^b1), psme2022, 
                                                                fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 150, a1p = -77, a2 = 50, a2p = -15, b1 = 0.72)))
  psmeDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1r + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psme2022, 
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018)))
  psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1r + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, 
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108)))
  psmeDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, 
                                                   fixedFormula = a1 + a1p + a2 + a2p + a3 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053)))
  psmeDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1r + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation), psme2022,
                                                     fixedFormula = a1 + a1p + a4 + a5 + a6 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102)))
  psmeDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation), psme2022, 
                                                    fixedFormula = a1 + a4 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 1.95, a4 = 0.361, b1 = 0.943, b1p = -0.068)))
  psmeDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1 + a1r + a1r) * (height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                               fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096)))
  psmeDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                   fixedFormula = a1 + a2 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 2.5, a2 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013)))
  psmeDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                         fixedFormula = a1 + a2 + a3 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012)))
  psmeDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                        fixedFormula = a1 + a2 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 2.5, a2 = -0.027, a4 = 0.3, b1 = 0.91, b1p = -0.2, b2 = 0, b2p = 0.013)))
  psmeDiameterFromHeightMixed$ruarkAbatRelHtPhysio = fit_nlme("Ruark ABA+T RelHt physio", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect) + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                              fixedFormula = a1 + a2 + a3 + a6 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2.4, a2 = -0.017, a3 = -0.004, a6 = -0.05, a4 = 0.3, b1 = 0.92, b1p = -0.19, b2 = 0.001, b2p = 0.012)))
  psmeDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1r + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                     fixedFormula = a1 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 2.5, a6 = -0.05, b1 = 0.84, b1p = -0.11, b2 = 0.005, b2p = 0.008)))
  psmeDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                    fixedFormula = a1 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2.4, a4 = 0.45, b1 = 0.83, b1p = -0.13, b2 = 0.004, b2p = 0.008)))
  psmeDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a6 * cos(3.14159/180 * aspect) + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                          fixedFormula = a1 + a6 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 2.4, a6 = -0.05, a4 = 0.4, b1 = 0.83, b1p = -0.13, b2 = 0.0042, b2p = 0.008)))
  psmeDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1), 0.9999)), psme2022, 
                                                 fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90)))
  psmeDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ (a1 + a1r + a1r) * (height - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + b3p * isPlantation)*(height - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2022, 
                                                      fixedFormula = a1 + b1 + b1p + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13)))
  psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1r + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                             fixedFormula = a1 + a1p + a2 + a3 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, b1 = 0.61, b2 = 0.071, b2p = 0.021)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                                   fixedFormula = a1 + a1p + a2 + a3 + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                   start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                                  fixedFormula = a1 + a1p + a2 + a3 + a4 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, a4 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysio = fit_nlme("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                                        fixedFormula = a1 + a1p + a2 + a3 + a6 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = 0, a3 = 0, a6 = -0.06, a4 = 0.5, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                               fixedFormula = a1 + a1p + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 3.0, a1p = -0.3, a6 = -0.06, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                              fixedFormula = a1 + a1p + a4 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 3.90, a1p = -0.95, a4 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                                    fixedFormula = a1 + a1p + a6 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 3.2, a1p = -0.6, a6 = -0.06, a4 = 0.5, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ ((a1 + a1r + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(height - 1.37), 0.9999)))^b2, psme2022, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81)))
  
  saveRDS(psmeDiameterFromHeightMixed, "trees/height-diameter/data/PSME DBH mixed.Rds")
}

if (psmeOptions$fitHeightPrimary & psmeOptions$fitHeightMixed & psmeOptions$fitHeightNlrobAndFixedWeight & 
    psmeOptions$fitDbhPrimary & psmeOptions$fitDbhMixed & psmeOptions$fitDbhNlrobAndFixedWeight)
{
  # assemble parameter and results tibbles using incremental load and reduce to fit in memory
  # height primary + nlrob + gsl_nls + dbh primary + nlrob = 120 GB DDR at 10x10 cross validation with models dropped
  # Therefore, batch formation height and diameter results.
  if (exists("psmeHeightFromDiameter") == FALSE) { readRDS("trees/height-diameter/data/PSME height primary.Rds") }
  if (exists("psmeHeightFromDiameterMixed") == FALSE) { readRDS("trees/height-diameter/data/PSME height mixed.Rds") }
  if (exists("psmeDiameterFromHeight") == FALSE) { readRDS("trees/height-diameter/data/PSME dbh primary.Rds") }
  if (exists("psmeDiameterFromHeightMixed") == FALSE) { readRDS("trees/height-diameter/data/PSME dbh mixed.Rds") }
  
  psmeCoefficients = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_list_coefficients)),
                                         bind_rows(lapply(psmeHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed"))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(psmeDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed"))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")
  psmeResults = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_list_stats)),
                                    bind_rows(lapply(psmeHeightFromDiameterMixed, get_list_stats, fitSet = "mixed"))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitting = "gsl_nls", fitSet = "primary"),
                                    bind_rows(lapply(psmeDiameterFromHeightMixed, get_list_stats, fitSet = "mixed"))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")
  
  check_plot_results(psmeResults)
  saveRDS("trees/height-diameter/data/PSME results.Rds", psmeCoefficients, psmeResults)
} else if (psmeOptions$fitHeightPrimary & psmeOptions$fitDbhPrimary)
{
  if (exists("psmeHeightFromDiameter") == FALSE) { readRDS("trees/height-diameter/data/PSME height primary.Rds") }
  if (exists("psmeDiameterFromHeight") == FALSE) { readRDS("trees/height-diameter/data/PSME DBH primary.Rds") }
  
  psmeCoefficients = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_list_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_list_coefficients))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")
  psmeResults = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_list_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitting = "gsl_nls", fitSet = "primary")) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")
  
  check_plot_results(psmeResults)
  saveRDS("trees/height-diameter/data/PSME results.Rds", psmeCoefficients, psmeResults)
}


## preferred forms identified (results.R, Figure 8 and trees.R)
if (psmeOptions$fitHeightPrimary & psmeOptions$fitHeightNlrobAndFixedWeight & psmeOptions$fitDbhPrimary & psmeOptions$fitDbhNlrobAndFixedWeight)
{
  psmeHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, folds = 1, repetitions = 1))
  psmeHeightFromDiameterPreferred$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + (a2 + a2p * isPlantation)*dbh + a3 + a3p* isPlantation), psme2022, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(dbh + b2 + b2p * isPlantation)), psme2022, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness + (a4 + a4p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 22, a1p = -7.6, a4 = -0.0020, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.06, a4 = -0.35, a4p = 0.79, b1 = 0.28, b2 = -0.021, b2p = -0.026, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.40), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation + (a4 + a4p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 4.2, a1p = 14.5, a4 = 0.16, a4p = 0.27, b1 = 0.63, b1p = -0.43, b2 = -0.025, b3 = -0.09, b4 = 1.73, b4p = -0.66), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050), folds = 1, repetitions = 1)
  #AIC(psmeHeightFromDiameterPreferred$prodan, psmeHeightFromDiameterPreferred$sibbesen)
  
  psmeDiameterFromHeightPreferred = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), folds = 1, repetitions = 1))
  psmeDiameterFromHeightPreferred$gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, threads = 2, folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, folds = 1, repetitions = 1, threads = 4)
  psmeDiameterFromHeightPreferred$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496), data = psme2022, folds = 1, repetitions = 1, threads = 4)
  psmeDiameterFromHeightPreferred$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017), folds = 1, repetitions = 1)
  
  saveRDS(psmeHeightFromDiameterPreferred, "trees/height-diameter/data/PSME preferred height models.Rds")
  saveRDS(psmeDiameterFromHeightPreferred, "trees/height-diameter/data/PSME preferred diameter models.Rds")
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  psmeHeightGam = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 11) + 
                               s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 8) + 
                               s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 8) + 
                               s(elevation, bs = "ts", k = 6) + 
                               s(slope, bs = "ts", k = 6) + 
                               s(aspect, bs = "ts", k = 8) + 
                               s(terrainRoughness, bs = "ts", k = 6) + 
                               s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 6), 
                      data = psme2022, method = "REML", select = TRUE, weights = dbhWeight)
  k.check(psmeHeightGam)
  summary(psmeHeightGam) # all signficant
  par(mfrow = c(3, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(psmeHeightGam, scale = 0)
  
  psmeDbhGam = gam(dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 10) +
                         s(bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 8) +
                         s(basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 15) +
                         s(elevation, bs = "ts", k = 9) +
                         s(slope, bs = "ts", k = 11) +
                         s(aspect, bs = "ts", k = 25) + # high complexity but only a few mm range, likely effectively noise
                         s(terrainRoughness, bs = "ts", k = 15) +
                         s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 14), 
                   data = psme2022, method = "REML", select = TRUE, weights = heightWeight)
  k.check(psmeDbhGam)
  summary(psmeDbhGam) # all significant
  plot.gam(psmeDbhGam, scale = 0) #, ylim = c(-15, 15))
}

## random forest regression
if (htDiaOptions$includeInvestigatory)
{
  startTime = Sys.time()
  psmeHeightForest = train(height ~ dbh + standBasalAreaPerHectare + basalAreaLarger + elevation + slope + aspect + terrainRoughness + relativeDiameter, data = psme2016physio, method = "ranger", trControl = repeatedCrossValidation, 
                           importance = "impurity_corrected",
                           tuneGrid = expand.grid(mtry = c(4, 5, 6),
                                                  splitrule = "variance",
                                                  min.node.size = c(2, 3, 4)))
  Sys.time() - startTime
  psmeHeightForest
  varImp(psmeHeightForest)
  
  startTime = Sys.time()
  psmeDbhForest = train(dbh ~ height + bootstrapStandBasalAreaPerHectare + basalAreaTaller + elevation + slope + aspect + terrainRoughness + relativeHeight, data = psme2016physio, method = "ranger", trControl = repeatedCrossValidation, 
                        importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = c(5, 6, 7),
                                               splitrule = "variance",
                                               min.node.size = c(4, 5, 6)))
  Sys.time() - startTime
  psmeDbhForest
  varImp(psmeDbhForest)
}
