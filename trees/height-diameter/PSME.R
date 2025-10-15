# load libraries, functions, and psme2022 from setup.R

psmeOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightMixed = FALSE,
                     fitDbhPrimary = FALSE,
                     fitDbhGslNlsAndGams = TRUE,
                     fitDbhMixed = FALSE)

## Douglas-fir height regressions
if (psmeOptions$fitHeightPrimary)
{
  # linear regressions
  psmeHeightFromDiameter = list(linear = fit_lm("linear", height ~ 0 + dbh + I(isPlantation*dbh), psme2022))
  psmeHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ 0 + dbh + I(dbh^2) + I(isPlantation*dbh) + I(isPlantation*dbh^2), psme2022)
  # nonlinear regressions
  psmeHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))
  psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 65, a1p = -9, a2 = 0.07, a2p = 0.6, a3 = 0.07, a3p = -0.05, b1 = -0.016, b2 = 1.25, b2p = -0.07))
  psmeHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 64, a1p = -23, a2 = 0, a2p = 0.48, a3 = 0.07, a3p = 0.09, a9 = -0.77, a9p = 2.1, b1 = -0.020, b2 = 1.44, b2p = -0.28)) # a2 not significant
  psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness) * (1 - exp((b1 + b1p * isPlantation)*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 63, a2 = 0.03, a2p = 0.67, a3 = 0.084, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, b1 = -0.021, b1p = 0.009, b2 = 1.5, b2p = -0.4))
  psmeHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness + (a9 + a9p * isPlantation) * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 61, a2 = 0, a2p = 0.72, a3 = 0.11, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, a9 = -0.4, a9p = 1.2, b1 = -0.023, b1p = 0.010, b2 = 1.54, b2p = -0.49))
  psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * terrainRoughness) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31))
  psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 64, a1p = -22, a2 = 0, a2p = 0.5, a3 = -0.07, a3p = 0.08, a9 = -0.77, a9p = 2.1, b1 = -0.020, b2 = 1.43, b2p = -0.28))
  psmeHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * terrainRoughness + (a9 + a9p * isPlantation) * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 72, a1p = -18, a4 = -0.0045, a5 = -8.6, a6 = 0.7, a7 = 0.8, a8 = 0.23, a9 = -1.2, a9p = 1.6, b1 = -0.022, b2 = 1.5, b2p = -0.36))
  psmeHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156))
  psmeHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (a2 + a2p * isPlantation) * dbh^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 75.4, a1p = -11.4, a2 = 462, a2p = -322, b1 = -1.54, b1p = 0.28))
  psmeHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = gsl_nls_control())
  psmeHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1 + a1p*isPlantation)*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psme2022, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30))
  psmeHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + (a2 + a2p * isPlantation)*dbh + a3 + a3p* isPlantation), psme2022, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6))
  psmeHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation), psme2022, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14))
  psmeHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(dbh + b2 + b2p * isPlantation)), psme2022, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52))
  psmeHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * dbh)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2022, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126))
  psmeHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3p * isPlantation) * dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 15, b1 = 0.35, b1p = -0.047, b2 = -0.023, b2p = -0.009, b3 = 0.002, b3p = -0.09, b4 = 1.52, b4p = -0.42))
  psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 6, a1p = 13, b1 = 0.6, b1p = -0.33, b2 = -0.025, b3 = 0.0, b4 = 1.53, b4p = -0.15))
  psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 20, a1p = -2.8, a4 = -0.0016, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.07, b1 = 0.30, b2 = -0.035, b2p = -0.007, b3 = -0.003, b3p = -0.07, b4 = 1.57, b4p = -0.51))
  psmeHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness + (a9 + a9p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 22, a1p = -7.6, a4 = -0.0020, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.06, a9 = -0.35, a9p = 0.79, b1 = 0.28, b2 = -0.021, b2p = -0.026, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.40))
  psmeHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 4.2, a1p = 14.5, a9 = 0.16, a9p = 0.27, b1 = 0.63, b1p = -0.43, b2 = -0.025, b3 = -0.09, b4 = 1.73, b4p = -0.66))
  psmeHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 20, a4 = -0.0008, a5 = -0.04, a6 = 0.15, a7 = 0.23, a8 = 0.095, b1 = 0.30, b2 = -0.025, b3 = -0.012, b3p = -0.12, b4 = 1.6, b4p = -0.65))
  psmeHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + (a1 + (a9 + a9p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3p * isPlantation) * dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 13.3, a9 = -0.2, a9p = 0.73, b1 = 0.40, b1p = -0.12, b2 = -0.020, b2p = -0.0233, b3 = 0.032, b3p = -0.20, b4 = 1.5, b4p = -0.39))
  psmeHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness + (a9 + a9p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 15.5, a4 = -0.0016, a5 = -0.028, a7 = 0.17, a8 = 0.06, a9 = 0.48, a9p = -0.3, b1 = 0.33, b2 = -0.027, b3 = -0.043, b3p = -0.074, b4 = 1.52, b4p = -0.55))
  psmeHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*standTreesPerHectare^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 54, a1p = -33, b1 = 0.05, b1p = 0.2, b2 = -0.03, b2p = -0.05, b3 = -0.04, b3p = -0.16, b4 = 1.56, b4p = -0.48))
  psmeHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare)*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 60, a2 = 0.026, a2p = 0.65, a3 = 0.09, a3p = -0.07, b1 = 0.05, b2 = -0.017, b3 = 0.02, b3p = -0.05, b4 = 1.35, b4p = -0.22))
  psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050))
  psmeHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation))), psme2022, start = list(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1))
  psmeHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation))), psme2022, start = list(a1 = 60, a2 = 0, a2p = 0.68, a3 = 0.066, b1 = -0.0053, b1p = -0.0042, b2 = 1.3, b2p = -0.20))
  
  # GAMs
  psmeHeightFromDiameter$gam = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, nthreads = 2)
  psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 26), data = psme2022, nthreads = 2)
  psmeHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 37), data = psme2022)
  psmeHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 40), data = psme2022)
  if (psmeOptions$fitPhysioGams)
  {
    psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, nthreads = 2)
    psmeHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ s(dbh, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), terrainRoughness, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022)
    psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), data = psme2022, nthreads = 2)
    psmeHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeDiameter, bs = "ts", k = 331, by = as.factor(isPlantation)), data = psme2022)
    
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
  psmeHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation + a1r) *(1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, 
                                                                fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))))
  psmeHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, 
                                                            fixedFormula = a1 + a1p + a2 + a2p + a3 + a3p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 65, a1p = -9, a2 = 0.07, a2p = 0.6, a3 = -0.07, a3p = -0.05, b1 = -0.016, b2 = 1.25, b2p = -0.07)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness) * (1 - exp((b1 + b1p * isPlantation)*dbh))^(b2 + b2p * isPlantation), psme2022, 
                                                                  fixedFormula = a1 + a2 + a2p + a3 + a4 + a5 + a6 + a7 + a8 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 63, a2 = 0.03, a2p = 0.67, a3 = 0.084, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, b1 = -0.021, b1p = 0.009, b2 = 1.5, b2p = -0.4)), control = nlmeControl())
  psmeHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * terrainRoughness) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, 
                                                               fixedFormula = a1 + a1p + a4 + a5 + a6 + a7 + a8 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31)))
  psmeHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation + a1r) * dbh / (1 + dbh)^(b1 + b1p*isPlantation), psme2022, 
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156)))
  psmeHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1p * isPlantation + a1r) / (1 + (a2 + a2p * isPlantation) * dbh^(b1 + b1p * isPlantation)), psme2022, 
                                                  fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 75.4, a1p = -11.4, a2 = 462, a2p = -322, b2 = -1.54, b2p = 0.28)))
  psmeHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp(b1*dbh^(b2 + b2p * isPlantation)), psme2022, 
                                              fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                              start = list(fixed = c(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084)), control = nlmeControl())
  psmeHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1p*isPlantation + a1r)*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psme2022, 
                                                         fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30)))
  psmeHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + (a2 + a2p * isPlantation)*dbh + a3 + a3p* isPlantation + a3r), psme2022, 
                                                fixedFormula = a1 + a2 + a2p + a3 + a3p ~ 1, randomFormula = a3r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6)))
  psmeHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*dbh^(b1 + b1p * isPlantation), psme2022, 
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14)))
  psmeHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*exp((b1 + b1p * isPlantation)/(dbh + b2 + b2p * isPlantation)), psme2022, 
                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52)))
  psmeHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * dbh)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2022, 
                                                   fixedFormula = Ha + Hap + d + dp + kU + kUp ~ 1, randomFormula = Har ~ 1|stand/plot,
                                                   start = list(fixed = c(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1 + a1r)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, 
                                                      fixedFormula = a1 + b1 + b1p + b2 + b2p + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4p * isPlantation), psme2022, 
                                                         fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 6, a1p = 13, b1 = 0.6, b1p = -0.33, b2 = -0.025, b3 = 0.0, b4 = 1.53, b4p = -0.15)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, 
                                                               fixedFormula = a1 + a1p + a4 + a5 + a6 + a7 + a8 + b1 + b2 + b2p + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 20, a1p = -2.8, a4 = -0.0016, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.07, b1 = 0.30, b2 = -0.035, b2p = -0.007, b3 = -0.003, b3p = -0.07, b4 = 1.57, b4p = -0.51)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1r + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, 
                                                            fixedFormula = a1 + a4 + a5 + a6 + a7 + a8 + b1 + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 20, a4 = -0.0008, a5 = -0.04, a6 = 0.15, a7 = 0.23, a8 = 0.095, b1 = 0.30, b2 = -0.025, b3 = -0.012, b3p = -0.12, b4 = 1.6, b4p = -0.65)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*standTreesPerHectare^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, 
                                                     fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 54, a1p = -33, b1 = 0.05, b1p = 0.2, b2 = -0.03, b2p = -0.05, b3 = -0.04, b3p = -0.16, b4 = 1.56, b4p = -0.48)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare)*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, 
                                                        fixedFormula = a1 + a2 + a2p + a3 + a3p + b1 + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 60, a2 = 0.026, a2p = 0.65, a3 = 0.09, a3p = -0.07, b1 = 0.05, b2 = -0.017, b3 = 0.02, b3p = -0.05, b4 = 1.35, b4p = -0.22)), control = nlmeControl())
  psmeHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, 
                                                  fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050)), control = nlmeControl())
  psmeHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*(1 - exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation))), psme2022, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1)))
  psmeHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation))), psme2022, 
                                                    fixedFormula = a1 + a2+ a2p + a3 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 60, a2 = 0, a2p = 0.68, a3 = 0.066, b1 = -0.0053, b1p = -0.0042, b2 = 1.3, b2p = -0.20)), control = nlmeControl())
  
  psmeHeightFromDiameterMixed$gamm = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 15) + s(StandID, bs = "re"), data = psme2022, mixed = TRUE)
  psmeHeightFromDiameterMixed$gammBal = fit_gam("REML GAM BA+L", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 26) + s(StandID, bs = "re"), data = psme2022, mixed = TRUE)
  
  saveRDS(psmeHeightFromDiameterMixed, "trees/height-diameter/data/PSME height mixed.Rds")
}


## Douglas-fir diameter regressions
if (psmeOptions$fitDbhPrimary)
{
  # linear regressions
  psmeDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)), psme2022))
  psmeDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)) + I((height - 1.37)^2) + I(isPlantation*(height - 1.37)^2), psme2022)
  
  # nonlinear regressions
  if (psmeOptions$fitDbhGslNlsAndGams)
  {
    psmeDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 40, b1 = 0.03, b1p = -0.005, b2 = 0.61, b2p = 0.15))
    psmeDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 75, a1p = -35, a2 = 0.5, a3 = -0.07, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07))
    psmeDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(height - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2022, start = list(a1 = 40, a1p = -13, a9 = 5, b1 = 0.7, b2 = 0.7, b2p = 0.03))
    psmeDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -123, a1p = 53.1, b1 = 0.0085, b1p = 0.0041, b2 = 0.77))
    psmeDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -175, a1p = 100, a2 = 1.2, a3 = 0.14, b1 = 0.007, b1p = 0.0057, b2 = 0.79))
    psmeDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -12, a1p = -3.9, a5 = -2.2, a8 = 0.04, b1 = 0.020, b1p = 0.0054, b2 = 0.45))
    psmeDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2022, start = list(a1 = -18, a9 = 0, b1 = 0.016, b1p = 0.0084, b2 = 0.08, b2p = 0.4))
    psmeDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation) * (height - 1.37)^b1 / (a2 + a2p * isPlantation - (height - 1.37)^b1), psme2022, start = list(a1 = 150, a1p = -77, a2 = 50, a2p = -15, b1 = 0.72))
    psmeDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psme2022, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018))
    psmeDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108))
    psmeDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053))
    psmeDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation), psme2022, start = list(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102))
    psmeDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation), psme2022, start = list(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068))
    psmeDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096))
    psmeDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013))
    psmeDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012))
    psmeDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", dbh ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.02, a9 = 0.2, b1 = 0.91, b1p = -0.16, b2 = 0, b2p = 0.009))
    psmeDiameterFromHeight$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark ABA+T RelHt physio", dbh ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.4, a2 = -0.013, a3 = -0.0037, a6 = -0.05, a9 = 0.1, b1 = 0.94, b1p = -0.16, b2 = 0.001, b2p = 0.01))
    psmeDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ (a1 + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a6 = -0.05, b1 = 0.84, b1p = -0.11, b2 = 0.005, b2p = 0.008))
    psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a9*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.2, a9 = 0.37, b1 = 0.85, b1p = -0.11, b2 = 0.003, b2p = 0.007))
    psmeDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.2, a6 = -0.06, a9 = 0.35, b1 = 0.86, b1p = -0.11, b2 = 0.0033, b2p = 0.007))
    #psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1), 0.9999)), psme2022, start = list(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90))
    psmeDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + b3p * isPlantation)*(height - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13), control = nlmeControl())
    psmeDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017))
    psmeDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, b1 = 0.61, b2 = 0.071, b2p = 0.021))
    psmeDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077))
    psmeDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.01, a3 = -0.006, a9 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021))
    psmeDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.1, a1p = -0.5, a2 = -0.01, a3 = -0.005, a6 = -0.06, a9 = 0.1, b1 = 0.65, b2 = 0.07))
    psmeDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.0, a1p = -0.3, a6 = -0.06, b1 = 0.60, b2 = 0.09))
    psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.2, a1p = -0.65, a9 = 0.43, b1 = 0.62, b2 = 0.07, b2p = 0.04))
    psmeDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.0, a1p = -0.5, a6 = -0.05, a9 = 0.4, b1 = 0.62, b2 = 0.07))
    psmeDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(height - 1.37), 0.9999)))^b2, psme2022, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81))
    # GAMs
    # individual term selection: height + ABA + AAT + RelHt by = isPlantation + slope + elevation + sin(aspect) + TSI + RelHt
    #   primary effects NSE 0.861
    psmeDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, nthreads = 2)
    psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28), data = psme2022, nthreads = 2)
    psmeDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, nthreads = 2)
    #lapply(psmeDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio$fit, confint2, level = 0.99)
    #lapply(psmeDiameterFromHeight$sibbesenReplaceAbat$fit, get_model_coefficients)
  }
  
  if (psmeOptions$fitPhysioGams)
  {
    psmeDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, nthreads = 2)
    psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ s(height, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), data = psme2022, nthreads = 2)
    psmeDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, nthreads = 2)
    
    psmeDiameterFromHeightGamAbatPhysio = psmeDiameterFromHeight$gamAbatPhysio
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
    psmeDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496), data = psme2022, nthreads = 2)
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
                                                               start = list(fixed = c(a1 = 40, b1 = 0.03, b1p = -0.005, b2 = 0.61, b2p = 0.15)), control = nlmeControl()))
  psmeDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, 
                                                            fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 75, a1p = -35, a2 = 0.5, a3 = -0.07, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07)), control = nlmeControl())
  psmeDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(height - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2022, 
                                                             fixedFormula = a1 + a1p + a9 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 30, a1p = -10, a9 = 5, b1 = 0.12, b2 = 0.60, b2p = 0.035)), control = nlmeControl())
  psmeDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1r + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, 
                                                         fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = -123, a1p = 53.1, b1 = 0.0085, b1p = 0.0041, b2 = 0.77)), control = nlmeControl())
  psmeDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, 
                                                             fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = -175, a1p = 100, a2 = 1.2, a3 = 0.14, b1 = 0.007, b1p = 0.0057, b2 = 0.79)), control = nlmeControl())
  psmeDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1r + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022,
                                                               fixedFormula = a1 + a1p + a5 + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = -12, a1p = -3.9, a5 = -2.2, a8 = 0.04, b1 = 0.020, b1p = 0.0054, b2 = 0.45)), control = nlmeControl())
  psmeDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2022, 
                                                              fixedFormula = a1 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = -9.0, a9 = -12.6, b1 = 0.013, b1p = 0.009, b2 = 0.004, b2p = 0.31)), control = nlmeControl())
                                                               control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.001, msTol = 1E-5)
  psmeDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1r + a1p * isPlantation) * (height - 1.37)^b1 / (a2 + a2p * isPlantation - (height - 1.37)^b1), psme2022, 
                                                                fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 150, a1p = -77, a2 = 50, a2p = -15, b1 = 0.72)), control = nlmeControl())
  psmeDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1r + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psme2022, 
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018)), control = nlmeControl())
  psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1r + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, 
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108)))
  psmeDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, 
                                                   fixedFormula = a1 + a1p + a2 + a2p + a3 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053)))
  psmeDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1r + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation), psme2022,
                                                     fixedFormula = a1 + a1p + a4 + a5 + a6 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102)))
  psmeDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation), psme2022, 
                                                    fixedFormula = a1 + a9 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068)))
  psmeDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1 + a1r + a1r) * (height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                               fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096)))
  psmeDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                   fixedFormula = a1 + a2 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 2.5, a2 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013)))
  psmeDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                         fixedFormula = a1 + a2 + a3 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012)))
  psmeDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                        fixedFormula = a1 + a2 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 2.5, a2 = -0.027, a9 = 0.3, b1 = 0.91, b1p = -0.2, b2 = 0, b2p = 0.013)), control = nlmeControl())
  psmeDiameterFromHeightMixed$ruarkAbatRelHtPhysio = fit_nlme("Ruark ABA+T RelHt physio", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                              fixedFormula = a1 + a2 + a3 + a6 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2.4, a2 = -0.017, a3 = -0.004, a6 = -0.05, a9 = 0.3, b1 = 0.92, b1p = -0.19, b2 = 0.001, b2p = 0.012)))
  psmeDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1r + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                     fixedFormula = a1 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 2.5, a6 = -0.05, b1 = 0.84, b1p = -0.11, b2 = 0.005, b2p = 0.008)))
  psmeDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a9*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                    fixedFormula = a1 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2.4, a9 = 0.45, b1 = 0.83, b1p = -0.13, b2 = 0.004, b2p = 0.008)))
  psmeDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
                                                          fixedFormula = a1 + a6 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 2.4, a6 = -0.05, a9 = 0.4, b1 = 0.83, b1p = -0.13, b2 = 0.0042, b2p = 0.008)))
  psmeDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1), 0.9999)), psme2022, 
                                                 fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90)), control = nlmeControl())
  psmeDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ (a1 + a1r + a1r) * (height - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + b3p * isPlantation)*(height - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2022, 
                                                      fixedFormula = a1 + b1 + b1p + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13)), control = nlmeControl())
  psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1r + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                             fixedFormula = a1 + a1p + a2 + a3 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, b1 = 0.61, b2 = 0.071, b2p = 0.021)), control = nlmeControl())
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                                   fixedFormula = a1 + a1p + a2 + a3 + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                   start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                                  fixedFormula = a1 + a1p + a2 + a3 + a9 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, a9 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysio = fit_nlme("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                                        fixedFormula = a1 + a1p + a2 + a3 + a6 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = 0, a3 = 0, a6 = -0.06, a9 = 0.5, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                               fixedFormula = a1 + a1p + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 3.0, a1p = -0.3, a6 = -0.06, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
                                                              fixedFormula = a1 + a1p + a9 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 3.90, a1p = -0.95, a9 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016)), control = nlmeControl())
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
                                                                    fixedFormula = a1 + a1p + a6 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 3.2, a1p = -0.6, a6 = -0.06, a9 = 0.5, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ ((a1 + a1r + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(height - 1.37), 0.9999)))^b2, psme2022, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81)))
  
  psmeDiameterFromHeightMixed = list(gamm = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 10) + s(StandID, bs = "re"), data = psme2022, mixed = TRUE))
  psmeDiameterFromHeightMixed$gammAbat = fit_gam("REML GAM ABA+T", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28) + s(StandID, bs = "re"), data = psme2022, mixed = TRUE)
  psmeDiameterFromHeightMixed$gammRelHt = fit_gam("REML GAM RelHt", dbh ~ s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 13) + s(StandID, bs = "re"), data = psme2022, mixed = TRUE)
  
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
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * terrainRoughness + (a9 + a9p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 22, a1p = -7.6, a4 = -0.0020, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.06, a9 = -0.35, a9p = 0.79, b1 = 0.28, b2 = -0.021, b2p = -0.026, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.40), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^(b4 + b4p * isPlantation), psme2022, start = list(a1 = 4.2, a1p = 14.5, a9 = 0.16, a9p = 0.27, b1 = 0.63, b1p = -0.43, b2 = -0.025, b3 = -0.09, b4 = 1.73, b4p = -0.66), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050), folds = 1, repetitions = 1)
  #AIC(psmeHeightFromDiameterPreferred$prodan, psmeHeightFromDiameterPreferred$sibbesen)
  
  psmeDiameterFromHeightPreferred = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), folds = 1, repetitions = 1))
  psmeDiameterFromHeightPreferred$gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, nthreads = 2, folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, folds = 1, repetitions = 1, nthreads = 4)
  psmeDiameterFromHeightPreferred$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496), data = psme2022, folds = 1, repetitions = 1, nthreads = 4)
  psmeDiameterFromHeightPreferred$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012), folds = 1, repetitions = 1)
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
                         s(standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 8) +
                         s(tallerApproxBasalArea, bs = "ts", by = as.factor(isPlantation), k = 15) +
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
  psmeDbhForest = train(dbh ~ height + standBasalAreaApprox + tallerApproxBasalArea + elevation + slope + aspect + terrainRoughness + relativeHeight, data = psme2016physio, method = "ranger", trControl = repeatedCrossValidation, 
                        importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = c(5, 6, 7),
                                               splitrule = "variance",
                                               min.node.size = c(4, 5, 6)))
  Sys.time() - startTime
  psmeDbhForest
  varImp(psmeDbhForest)
}
