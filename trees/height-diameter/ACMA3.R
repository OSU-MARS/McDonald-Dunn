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
  
  acmaHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1*dbh))^(b2 + b2p*isPlantation), acma2022, start = list(a1 = 27, b1 = -0.03, b2 = 1.1, b2p = -0.2))
  acmaHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 26, a2 = 0.05, b1 = -0.025, b2 = 1.08, b2p = -0.2))
  acmaHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a5 * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 30, a2 = 0.06, a5 = -0.1, b1 = -0.03, b2 = 1.1, b2p = -0.19))
  acmaHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a5 * slope + a9 * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 30, a2 = 0.06, a5 = -0.1, a9 = -0.2, b1 = -0.03, b2 = 1.1, b2p = -0.21))
  acmaHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a9 * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 27, a2 = 0.05, a9 = 0, b1 = -0.027, b2 = 1.08, b2p = -0.2))
  acmaHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", height ~ 1.37 + (a1 + a2 * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 1.0, a2 = 0.08, a9 = 53, a9p = -16, b1 = -0.04, b2 = 0.50))
  acmaHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + (a1 + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 30, a5 = -7, b1 = -0.034, b2 = 1.1, b2p = -0.28))
  acmaHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a9 * relativeDiameter)*(1 - exp(b1*dbh))^(b2 + b2p*isPlantation), acma2022, start = list(a1 = 28, a9 = -0.5, b1 = -0.03, b2 = 1.1, b2p = -0.2))
  acmaHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 31, a5 = -6.7, a9 = -0.6, b1 = -0.032, b2 = 1.1, b2p = -0.25))
  acmaHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + (a1 + a1p * isPlantation) * dbh / (1 + dbh)^b1, acma2022, start = list(a1 = 1.24, a1p = 0.23, b1 = 0.28))
  acmaHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1 * dbh^b2), acma2022, start = list(a1 = 40.1, a1p = 6.30, b1 = 44.9, b2 = -0.97))
  acmaHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1*dbh^(b2 + b2p * isPlantation)), acma2022, start = list(a1 = 70, b1 = -5.0, b2 = -0.35, b2p = -0.02))
  acmaHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), acma2022, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986))
  acmaHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1 * dbh^2 + (a2 + a2p * isPlantation)*dbh + a3), acma2022, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19))
  acmaHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^b1, acma2022, start = list(a1 = 1.11, a1p = 0.21, b1 = 0.75))
  acmaHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1/(dbh + b2 + b2p * isPlantation)), acma2022, start = list(a1 = 31.8, a1p = 3.00, b1 = -25.7, b2 = 6.70, b2p = 0.62))
  acmaHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), acma2022, start = list(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064))
  acmaHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 19, b1 = 0.1, b2 = -0.03, b3 = 0.06, b4 = 1.1, b4p = -0.25))
  acmaHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 15, b1 = 0.14, b2 = -0.03, b3 = -0.06, b4 = 1.1, b4p = -0.24))
  acmaHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, acma2022, start = list(a1 = 12, a6 = 0.4, b1 = 0.2, b1p = 0.04, b2 = -0.037, b3 = -0.046, b4 = 1.06), control = gsl_nls_control())
  acmaHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a6 * sin(3.14159/180 * aspect) + a9 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, acma2022, start = list(a1 = 13, a6 = 0, a9 = -0.5, b1 = 0.2, b1p = 0.05, b2 = -0.024, b3 = 0.06, b4 = 1.0))
  acmaHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1 + a9 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 16, a9 = -0.5, b1 = 0.14, b2 = -0.026, b3 = 0.06, b4 = 1.12, b4p = -0.24))
  acmaHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*dbh))^b4, acma2022, start = list(a1 = 11, a6 = 0.3, b1 = 0.2, b1p = 0.06, b2 = -0.03, b3 = 0, b4 = 1.0))
  acmaHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + (a1 + a9 * relativeDiameter)*topHeight^b1*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 17, a9 = -0.6, b1 = 0.14, b2 = -0.03, b3 = 0.08, b4 = 1.1, b4p = -0.23))
  acmaHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a6 * sin(3.14159/180 * aspect) + a9 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*dbh))^b4, acma2022, start = list(a1 = 11, a6 = 0.3, a9 = -0.6, b1 = 0.23, b1p = 0.05, b2 = -0.02, b3 = 0.10, b4 = 0.96))
  acmaHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*dbh))^b4, acma2022, start = list(a1 = 17, a1p = 2.0, b1 = 0.1, b2 = -0.05, b3 = -0.1, b4 = 1.0))
  acmaHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 22, a2 = 0.05, b1 = 0.05, b2 = -0.02, b3 = 0.04, b4 = 1.09, b4p = -0.2))
  acmaHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1*dbh^b2), acma2022, start = list(a1 = 0.752, a1p = 0.120, b1 = 1.180, b2 = -0.087))
  acmaHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*dbh^b2)), acma2022, start = list(a1 = 27, a1p = 3, b1 = -0.033, b2 = 0.94))
  acmaHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*dbh^b2)), acma2022, start = list(a1 = 26, a2 = 0.1, b1 = -0.026, b1p = -0.008, b2 = 0.99))
  acmaHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", height ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * pmin(relativeHeight, 1.5)) * (1 - exp(b1*dbh^b2)), acma2022, start = list(a1 = 10, a2 = -2, a4 = -10, b1 = 0.3, b2 = 0.43), control = gsl_nls_control())
  
  acmaHeightFromDiameter$gam = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022)
  acmaHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 14), data = acma2022)
  acmaHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ s(dbh, basalAreaLarger, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022)
  acmaHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ s(dbh, basalAreaLarger, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  acmaHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 25), data = acma2022)
  acmaHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022)
  acmaHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 15), data = acma2022)
  acmaHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", k = 85, by = as.factor(isPlantation)), data = acma2022)
  
  save(file = "trees/height-diameter/data/ACMA3 height.Rdata", acmaHeightFromDiameter)
}

if (acmaOptions$fitHeightMixed)
{
  acmaHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1r)*(1 - exp(b1*dbh))^(b2 + b2p*isPlantation), acma2022,
                                                                fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                start = list(fixed = c(a1 = 27, b1 = -0.03, b2 = 1.1, b2p = -0.2))))
  acmaHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1r + a2 * basalAreaLarger) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, 
                                                            fixedFormula = a1 + a2+ b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 26, a2 = 0.05, b1 = -0.025, b2 = 1.08, b2p = -0.2)), control = nlmeControl())
  acmaHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + a2 * basalAreaLarger + a5 * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022,
                                                                  fixedFormula = a1 + a2 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 30, a2 = 0.06, a5 = -0.1, b1 = -0.03, b2 = 1.1, b2p = -0.19)))
  acmaHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, 
                                                               fixedFormula = a1 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 30, a5 = -7, b1 = -0.034, b2 = 1.1, b2p = -0.28)))
  acmaHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1p * isPlantation + a1r) * dbh / (1 + dbh)^b1, acma2022, 
                                                fixedFormula = a1 + a1p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                start = list(fixed = c(a1 = 1.24, a1p = 0.23, b1 = 0.28)), control = nlmeControl())
  acmaHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1p * isPlantation + a1r) / (1 + b1 * dbh^b2), acma2022, 
                                                  fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                  start = list(fixed = c(a1 = 40.1, a1p = 6.30, b1 = 44.9, b2 = -0.97)), control = nlmeControl())
  acmaHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp(b1*dbh^(b2 + b2p * isPlantation)), acma2022, 
                                              fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                              start = list(fixed = c(a1 = 70, b1 = -5.0, b2 = -0.35, b2p = -0.02)), control = nlmeControl())
  acmaHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1r)*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), acma2022, 
                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986)), control = nlmeControl())
  acmaHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / (a1 * dbh^2 + (a2 + a2p * isPlantation)*dbh + a3 + a3r), acma2022, 
                                                fixedFormula = a1 + a2 + a2p + a3 ~ 1, randomFormula = a3r ~ 1,
                                                start = list(fixed = c(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19)), control = nlmeControl())
  acmaHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*dbh^b1, acma2022, 
                                               fixedFormula = a1 + a1p + b1 ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 1.11, a1p = 0.21, b1 = 0.75)), control = nlmeControl())
  acmaHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*exp(b1/(dbh + b2 + b2p * isPlantation)), acma2022, 
                                                   fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 31.8, a1p = 3.00, b1 = -25.7, b2 = 6.70, b2p = 0.62)), control = nlmeControl())
  acmaHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Har) * (1 + ((1.37/(Ha + Har))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), acma2022, 
                                                   fixedFormula = Ha + d + kU + kUp ~ 1, randomFormula = Har ~ 1,
                                                   start = list(fixed = c(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064)))
  acmaHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1 + a1r)*topHeight^b1*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4p * isPlantation), acma2022, 
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
                                                      start = list(fixed = c(a1 = 19, b1 = 0.1, b2 = -0.03, b3 = 0.06, b4 = 1.1, b4p = -0.25)), control = nlmeControl())
  acmaHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1 + a1r)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4p * isPlantation), acma2022,
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 15, b1 = 0.14, b2 = -0.03, b3 = -0.06, b4 = 1.1, b4p = -0.24)), control = nlmeControl())
  acmaHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1r + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, acma2022, 
                                                               fixedFormula = a1 + a6 + b1 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 12, a6 = 0.4, b1 = 0.2, b1p = 0.04, b2 = -0.037, b3 = -0.046, b4 = 1.06)), control = nlmeControl())
  acmaHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1r + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*dbh))^b4, acma2022, 
                                                            fixedFormula = a1 + a6 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 11, a6 = 0.3, b1 = 0.2, b1p = 0.06, b2 = -0.03, b3 = 0, b4 = 1.0)), control = nlmeControl())
  acmaHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*dbh))^b4, acma2022,
                                                     fixedFormula = a1 + a1p + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 17, a1p = 2.0, b1 = 0.1, b2 = -0.05, b3 = -0.1, b4 = 1.0)), control = nlmeControl())
  acmaHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1r + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*dbh))^(b4 + b4p * isPlantation), acma2022, 
                                                        fixedFormula = a1 + a2 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 22, a2 = 0.05, b1 = 0.05, b2 = -0.02, b3 = 0.04, b4 = 1.09, b4p = -0.2)), control = nlmeControl())
  acmaHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*dbh^(b1*dbh^b2), acma2022, 
                                                  fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                  start = list(fixed = c(a1 = 0.752, a1p = 0.120, b1 = 1.180, b2 = -0.087)), control = nlmeControl())
  acmaHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*(1 - exp(b1*dbh^b2)), acma2022, 
                                                 fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 27, a1p = 3, b1 = -0.033, b2 = 0.94)))
  acmaHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1r + a2 * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*dbh^b2)), acma2022, 
                                                    fixedFormula = a1 + a2 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 26, a2 = 0.1, b1 = -0.026, b1p = -0.008, b2 = 0.99)), control = nlmeControl())
  
  acmaHeightFromDiameterMixed$gamm = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7) + s(StandID, bs = "re"), data = acma2022, mixed = TRUE)
  acmaHeightFromDiameterMixed$gammBal = fit_gam("REML GAM BA+L", height ~ s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 14) + s(StandID, bs = "re"), data = acma2022, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/ACMA3 height mixed.Rdata", acmaHeightFromDiameterMixed)
}


## bigleaf maple diameter-height regressions
if (acmaOptions$fitDbh)
{
  acmaDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)), acma2022))
  acmaDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)) + I(isPlantation*(height - 1.37)^2), acma2022)
  
  acmaDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, start = list(a1 = 50, b1 = 0.001, b2 = 1.0), control = gsl_nls_control())
  acmaDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, start = list(a1 = 50, a2 = 15, b1 = 0.001, b2 = 1.0), control = gsl_nls_control())
  acmaDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", dbh ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(height - 1.37)^b2) - 1), acma2022, start = list(a1 = 50, a2 = -15, b1 = 0.001, b2 = 1.0), control = gsl_nls_control())
  acmaDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", dbh ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(height - 1.37)^b2) - 1), acma2022, start = list(a1 = 50, a2 = -15, a3 = 4, a9 = -200, b1 = 0.001, b2 = 1.00), control = gsl_nls_control())
  acmaDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1 + a9 * relativeHeight)*(exp(b1*(height - 1.37)^b2) - 1), acma2022, start = list(a1 = 50, a9 = -15, b1 = 0.1, b2 = 1.0), control = gsl_nls_control())
  acmaDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 39, b1 = -0.022, b2 = 1.6))
  acmaDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 40, a2 = 0, b1 = -0.021, b2 = 1.6), control = gsl_nls_control())
  acmaDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a8 * terrainRoughness)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 40, a8 = 0, b1 = -0.023, b2 = 1.6))
  acmaDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 80, a9 = -30, b1 = -0.012, b2 = 1.6))
  acmaDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), acma2022, start = list(a1 = -116, a2 = -128, b1 = 1.50))
  acmaDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), acma2022, start = list(a1 = 5.1, a1p = -2.0, a2 = -0.12, a2p = -0.023))
  acmaDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282))
  acmaDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, start = list(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272))
  acmaDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*(height - 1.37)^b1, acma2022, start = list(a1 = 3.7, a1p = -0.9, a5 = -1.2, a8 = 0.02, b1 = 0.9))
  acmaDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a9 * relativeHeight)*(height - 1.37)^b1, acma2022, start = list(a1 = 2.5, a9 = -0.6, b1 = 1.0))
  acmaDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), acma2022, start = list(a1 = 1.3, b1 = 1.45, b1p = -0.3, b2 = -0.033, b2p = 0.03))
  acmaDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.3, a2 = -0.012, b1 = 1.45, b1p = -0.1, b2 = -0.03))
  acmaDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a3 * standBasalAreaApprox + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.5, a3 = -0.007, a7 = -0.05, b1 = 1.5, b1p = -0.11, b2 = -0.03))
  acmaDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", dbh ~ (a1 + a3 * standBasalAreaApprox + a7 * sin(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.4, a3 = -0.007, a7 = -0.06, a9 = -0.2, b1 = 1.45, b1p = -0.1, b2 = 0.028))
  acmaDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", dbh ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.3, a2 = -0.01, a9 = 0, b1 = 1.52, b1p = -0.09, b2 = -0.033))
  acmaDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ (a1 + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.2, a7 = -0.055, b1 = 1.45, b1p = -0.064, b2 = 0.03))
  acmaDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.4, a9 = 0, b1 = 1.45, b1p = -0.3, b2 = -0.03))
  acmaDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a7 * sin(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.2, a7 = -0.06, a9 = 0, b1 = 1.45, b1p = -0.064, b2 = 0.03))
  acmaDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), acma2022, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 161))
  acmaDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(height - 1.37)) - 1)^b4, acma2022, start = list(a1 = 0.29, b1 = 1.14, b2 = 0.00001, b3 = 0.8, b4 = -0.2), control = gsl_nls_control())
  acmaDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 6, b1 = 2, b2 = -0.13))
  acmaDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.5, a2 = 0, b1 = 2.2, b2 = -0.13))
  acmaDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a2 * tallerApproxBasalArea + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022, start = list(a1 = 0.4, a2 = -0.002, a5 = -0.15, b1 = 3.2, b2 = -0.18, b2p = -0.02))
  acmaDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.5, a3 = 0, a5 = -0.2, a9 = -0.1, b1 = 2.6, b2 = -0.14))
  acmaDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.8, a2 = -0.003, a9 = -0.4, b1 = 2.1, b2 = -0.10))
  acmaDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022, start = list(a1 = 0.5, a5 = -0.2, b1 = 2.6, b2 = -0.18, b2p = 0.02))
  acmaDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.65, a9 = -0.25, b1 = 2.0, b2 = -0.1))
  acmaDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 0.5, a5 = -0.2, a9 = -0.1, b1 = 2.6, b2 = -0.14))
  acmaDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, acma2022, start = list(a1 = -40, b1 = 0.1, b2 = 0.65), control = gsl_nls_control())
  #lapply(acmaDiameterFromHeight$sibbesenReplaceAbatRelHt$fit, confint2, level = 0.99)
  #lapply(acmaDiameterFromHeight$sibbesenReplaceAbat$fit, get_model_coefficients)
  
  # individual term selection: height + AAT + RelHt by = isPlantation, slope + elevation + cos(aspect) retained by AIC but not significant (p > 0.06)
  acmaDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022)
  acmaDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 13), data = acma2022)
  acmaDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, slope, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  acmaDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ s(height, standBasalAreaApprox, tallerApproxBasalArea, slope, terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  acmaDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ s(height, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022)
  acmaDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022)
  acmaDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, elevation, slope, sin(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022)
  
  save(file = "trees/height-diameter/data/ACMA3 dbh.Rdata", acmaDiameterFromHeight, acmaDiameterFromHeightNlrob, acmaDiameterFromHeightGslNlsDefault)
}

if (acmaOptions$fitDbhMixed)
{
  acmaDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1 + a1r) * (exp(b1*(height - 1.37)) - 1)^b2, acma2022,
                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 50, b1 = 0.001, b2 = 1.0)), control = nlmeControl()))
  acmaDiameterFromHeightMixed = list(chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, 
                                                                   fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                   start = list(fixed = c(a1 = 50, a2 = 15, b1 = 0.001, b2 = 1.0)), control = nlmeControl()))
  acmaDiameterFromHeightMixed = list(chapmanReplaceBal = fit_nlme("Chapman-Richards replace BA+L", dbh ~ (a1 + a1r + a2 * basalAreaLarger) * (exp(b1*(height - 1.37)^b2) - 1), acma2022,
                                                                  fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 50, a2 = -15, b1 = 0.001, b2 = 1.0)), control = nlmeControl()))
  acmaDiameterFromHeightMixed = list(chapmanReplaceBalRelHt = fit_nlme("Chapman-Richards replace BA+L RelHt", dbh ~ (a1 + a1r + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(height - 1.37)^b2) - 1), acma2022, 
                                                                       fixedFormula = a1 + a2 + a3 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                       start = list(fixed = c(a1 = 50, a2 = -15, a3 = 4, a9 = -200, b1 = 0.001, b2 = 1.00)), control = nlmeControl()))
  acmaDiameterFromHeightMixed = list(chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*(exp(b1*(height - 1.37)^b2) - 1), acma2022, 
                                                                    fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                    start = list(fixed = c(a1 = 50, a9 = -15, b1 = 0.1, b2 = 1.0)), control = nlmeControl()))
  acmaDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1r) * log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 39, b1 = -0.022, b2 = 1.6)), control = nlmeControl())
  acmaDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                             start = list(fixed = c(a1 = 40, a2 = 0, b1 = -0.021, b2 = 1.6)), control = nlmeControl())
  acmaDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1r + a8 * terrainRoughness)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022,
                                                               fixedFormula = a1 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 40, a8 = 0, b1 = -0.023, b2 = 1.6)), control = nlmeControl())
  acmaDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, 
                                                              fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 80, a9 = -30, b1 = -0.012, b2 = 1.6)), control = nlmeControl())
  acmaDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1r) * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), acma2022, 
                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1,
                                                                start = list(fixed = c(a1 = -116, a2 = -128, b1 = 1.50)), control = nlmeControl())
  acmaDiameterFromHeightMixed = list(naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1r + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), acma2022, 
                                                        fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 5.1, a1p = -2.0, a2 = -0.12, a2p = -0.023))))
  acmaDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1r + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022,
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 3.2, a1p = -1.7, b1 = 0.94, b1p = 0.25)), control = nlmeControl())
  acmaDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, 
                                                   fixedFormula = a1 + a1p + a2 + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272)))
  acmaDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1r + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*(height - 1.37)^b1, acma2022,
                                                     fixedFormula = a1 + a1p + a5 + a8 + b1 ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 3.7, a1p = -0.9, a5 = -1.2, a8 = 0.02, b1 = 0.9)))
  acmaDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*(height - 1.37)^b1, acma2022, 
                                                    fixedFormula = a1 + a9 + b1 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 2.5, a9 = -0.6, b1 = 1.0)))
  acmaDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1 + a1r) * (height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), acma2022, 
                                               fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 1.3, b1 = 1.45, b1p = -0.3, b2 = -0.033, b2p = 0.03)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
                                                   fixedFormula = a1 + a2 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 1.3, a2 = -0.012, b1 = 1.45, b1p = -0.1, b2 = -0.03)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1r + a3 * standBasalAreaApprox + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
                                                         fixedFormula = a1 + a3 + a7 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 1.5, a3 = -0.007, a7 = -0.05, b1 = 1.5, b1p = -0.11, b2 = -0.03)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark ABA+T RelHt physio", dbh ~ (a1 + a1r + a3 * standBasalAreaApprox + a7 * sin(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022,
                                                              fixedFormula = a1 + a3 + a7 + a9 + b1 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 1.4, a3 = -0.007, a7 = -0.06, a9 = -0.2, b1 = 1.45, b1p = -0.1, b2 = 0.028)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
                                                        fixedFormula = a1 + a2 + a9 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 1.3, a2 = -0.01, a9 = 0, b1 = 1.52, b1p = -0.09, b2 = -0.033)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1r + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022,
                                                     fixedFormula = a1 + a7 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 1.2, a7 = -0.055, b1 = 1.45, b1p = -0.064, b2 = 0.03)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, 
                                                    fixedFormula = a1 + a9 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 1.4, a9 = 0, b1 = 1.45, b1p = -0.3, b2 = -0.03)), control = nlmeControl())
  acmaDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a7 * sin(3.14159/180 * aspect) + a9 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022,
                                                          fixedFormula = a1 + a7 + a9 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                          start = list(fixed = c(a1 = 1.2, a7 = -0.06, a9 = 0, b1 = 1.45, b1p = -0.064, b2 = 0.03)), control = nlmeControl())
  acmaDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1)), acma2022, 
                                                 fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1,
                                                 start = list(fixed = c(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 161)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ (a1 + a1r) * (height - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(height - 1.37)) - 1)^b4, acma2022, 
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                      start = list(fixed = c(a1 = 0.29, b1 = 1.14, b2 = 0.00001, b3 = 0.8, b4 = -0.2)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1r) * (height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 6, b1 = 2, b2 = -0.13)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                             start = list(fixed = c(a1 = 0.5, a2 = 0, b1 = 2.2, b2 = -0.13)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022,
                                                                   fixedFormula = a1 + a2 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                   start = list(fixed = c(a1 = 0.4, a2 = -0.002, a5 = -0.15, b1 = 3.2, b2 = -0.18, b2p = -0.02)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1r + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
                                                                        fixedFormula = a1 + a3 + a5 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                        start = list(fixed = c(a1 = 0.5, a3 = 0, a5 = -0.2, a9 = -0.1, b1 = 2.6, b2 = -0.14)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1r + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
                                                                  fixedFormula = a1 + a2 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 0.8, a2 = -0.003, a9 = -0.4, b1 = 2.1, b2 = -0.10)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a5 * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022,
                                                               fixedFormula = a1 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 0.5, a5 = -0.2, b1 = 2.6, b2 = -0.18, b2p = 0.02)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, 
                                                              fixedFormula = a1 + a9 +b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 0.65, a9 = -0.25, b1 = 2.0, b2 = -0.1)), control = nlmeControl())
  acmaDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022,
                                                                    fixedFormula = a1 + a5 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                    start = list(fixed = c(a1 = 0.5, a5 = -0.2, a9 = -0.1, b1 = 2.6, b2 = -0.14)), control = nlmeControl())
  acmaDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ ((a1 + a1r) * log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, acma2022, 
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = -40, b1 = 0.1, b2 = 0.65)), control = nlmeControl())
  
  acmaDiameterFromHeightMixed$gamm = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 8) + s(StandID, bs = "re"), data = acma2022, mixed = TRUE)
  acmaDiameterFromHeightMixed$gammAbat = fit_gam("REML GAM ABA+T", dbh ~ s(height, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 13) + s(StandID, bs = "re"), data = acma2022, mixed = TRUE)
  acmaDiameterFromHeightMixed$gammRelHt = fit_gam("REML GAM RelHt", dbh ~ s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 7) + s(StandID, bs = "re"), data = acma2022, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/ACMA3 dbh mixed.Rdata", acmaDiameterFromHeightMixed)
}


## collect model parameters
if (acmaOptions$fitHeight & acmaOptions$fitHeightMixed & acmaOptions$fitDbh & acmaOptions$fitDbhMixed)
{
  if (exists("acmaHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/ACMA3 height.Rdata") }
  if (exists("acmaHeightFromDiameterMixed") == FALSE) { load("trees/height-diameter/data/ACMA3 height mixed.Rdata") }
  if (exists("acmaDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/ACMA3 dbh.Rdata") }
  if (exists("acmaDiameterFromHeightMixed") == FALSE) { load("trees/height-diameter/data/ACMA3 dbh mixed.Rdata") }
  
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
  save(file = "trees/height-diameter/data/ACMA3 results.Rdata", acmaCoefficients, acmaResults)
} else if (acmaOptions$fitHeight & acmaOptions$fitHeightMixed & acmaOptions$fitDbh & acmaOptions$fitDbhMixed)
{
  if (exists("acmaHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/ACMA3 height.Rdata") }
  if (exists("acmaDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/ACMA3 dbh.Rdata") }
  
  acmaCoefficients = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, get_list_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(acmaDiameterFromHeight, get_list_coefficients))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "ACMA3")
  acmaResults = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, get_list_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(acmaDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Chapman-Richards replace ABA+T", fitSet = "primary", fittingMethod = "gsl_nls")) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "ACMA3")
  
  check_plot_results(acmaResults)
  save(file = "trees/height-diameter/data/ACMA3 results.Rdata", acmaCoefficients, acmaResults)
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
  acmaHeightFromDiameterPreferred$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4p * isPlantation), acma2022, start = list(a1 = 19, b1 = 0.1, b2 = -0.03, b3 = 0.06, b4 = 1.1, b4p = -0.25), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*dbh^b2)), acma2022, start = list(a1 = 26, a2 = 0.1, b1 = -0.027, b1p = -0.01, b2 = 0.99), folds = 1, repetitions = 1)
  #AIC(acmaHeightFromDiameterPreferred$chapmanRichards, acmaHeightFromDiameterPreferred$michaelisMenten, acmaHeightFromDiameterPreferred$prodan, acmaHeightFromDiameterPreferred$richardsW)
  
  acmaDiameterFromHeightPreferred = list(chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a8 * terrainRoughness)*log(1 - pmin((b1 + b1p * isPlantation)*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 20, a8 = 0.2, b1 = -0.03, b1p = 0.01, b2 = 1.8), folds = 1, repetitions = 1))
  acmaDiameterFromHeightPreferred$gam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022, folds = 1, repetitions = 1)
  #acmaDiameterFromHeightPreferred$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ s(height, standBasalAreaApprox, tallerApproxBasalArea, slope, terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022, folds = 1, repetitions = 1)
  #acmaDiameterFromHeightPreferred$power = fit_gsl_nls("power", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation), acma2022, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282), folds = 1, repetitions = 1)
  #acmaDiameterFromHeightPreferred$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * terrainRoughness)*(height - 1.37)^b1, acma2022, start = list(a1 = 3.7, a1p = -0.9, a5 = -1.2, a8 = 0.02, b1 = 0.9), folds = 1, repetitions = 1)
  acmaDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), acma2022, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), folds = 1, repetitions = 1)
  acmaDiameterFromHeightPreferred$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a3 * standBasalAreaApprox + a7 * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 1.5, a3 = -0.007, a7 = -0.05, b1 = 1.5, b1p = -0.11, b2 = -0.03), folds = 1, repetitions = 1)
  
  save(file = "trees/height-diameter/data/ACMA3 preferred models.Rdata", acmaHeightFromDiameterPreferred, acmaDiameterFromHeightPreferred)
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  acmaHeightGam = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 9) + 
                            s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 5) + 
                            #s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 3) + # not significant
                            s(elevation, bs = "ts", k = 5) + 
                            s(slope, bs = "ts", k = 4) +
                            #s(aspect, bs = "ts", k = 3) + # not significant
                            #s(terrainRoughness, bs = "ts", k = 4) + # not significant
                            s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 4), 
                          data = acma2016physio, constraint = acma2016gamConstraint, folds = 1, repetitions = 1)
  k.check(acmaHeightGam)
  summary(acmaHeightGam)
  par(mfrow = c(2, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(acmaHeightGam, scale = 0)
  
  acmaDbhGam = fit_gam("REML GAM", dbh ~ s(height, bs = "ts", by = as.factor(isPlantation), k = 9) +
                         s(standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 4),
                       #s(tallerApproxBasalArea, bs = "ts", by = as.factor(isPlantation), k = 3) + # not significant
                       #s(elevation, bs = "ts", k = 3) + # not significant
                       #s(slope, bs = "ts", k = 4), # not significant
                       #s(aspect, bs = "ts", k = 3), # not significant
                       #s(terrainRoughness, bs = "ts", k = 3), # not significant
                       #s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 3), # not significant
                       data = acma2016physio, constraint = acma2016gamConstraint, folds = 1, repetitions = 1)
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
  
  acmaDbhForest = train(dbh ~ height + standBasalAreaApprox + tallerApproxBasalArea + elevation + slope + aspect + terrainRoughness + relativeHeight, data = acma2022, method = "ranger", trControl = repeatedCrossValidation, 
                        importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = c(6, 7, 8),
                                               splitrule = "variance",
                                               min.node.size = c(1, 2, 3)))
  acmaDbhForest
  varImp(acmaDbhForest)
}
