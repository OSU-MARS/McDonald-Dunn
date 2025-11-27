# load libraries, functions, and psme2022 from setup.R

psmeOptions = tibble(fitHeightPrimary = FALSE,
                     fitHeightMixed = TRUE,
                     fitDbhPrimary = TRUE,
                     fitDbhMixed = TRUE,
                     fitPhysioGams = TRUE)

## Douglas-fir height regressions
if (psmeOptions$fitHeightPrimary)
{
  # linear regressions
  psmeHeightFromDiameter = list(linear = fit_lm("linear", height ~ 0 + dbh + I(isPlantation*dbh), psme2022))
  psmeHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ 0 + dbh + I(dbh^2) + I(isPlantation*dbh^2), psme2022) # plantation, plantation^2 not mutually significant
  
  # nonlinear regressions
  psmeHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1 * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2p = 0.1)) # a1p not significant
  psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a2 * standBasalAreaPerHectare) * (1 - exp((b1)*dbh))^(b2 + b2bal * basalAreaLarger), psme2022, start = list(a1 = 57, a2 = 0.1, b1 = -0.02, b2 = 1.3, b2bal = -0.005)) # a2p, a3p, b1ba, b1bap, b1balp, b2p, b2ba, b2bap, b2balp not significant, a3, b1bal also significant
  #psmeHeightFromDiameter$chapmanRichardsBalA3 = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a2 * standBasalAreaPerHectare + a3 * basalAreaLarger) * (1 - exp((b1)*dbh))^(b2), psme2022, start = list(a1 = 57, a2 = 0.1, a3 = 0, a3p = 0, b1 = -0.02, b2 = 1.3)) # job singular gradient
  #psmeHeightFromDiameter$chapmanRichardsBalB1 = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a2 * standBasalAreaPerHectare) * (1 - exp((b1 + b1bal * basalAreaLarger)*dbh))^(b2), psme2022, start = list(a1 = 57, a2 = 0.1, b1bal = 0, b1balp = 0, b1 = -0.02, b2 = 1.3)) # job singular gradient
  psmeHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a4 * relativeDiameter) * (1 - exp((b1)*dbh))^(b2 + b2bal * basalAreaLarger), psme2022, start = list(a1 = 70, a4 = -1.7, b1 = -0.015, b2 = 1.3, b2bal = -0.005)) # a2, a4p, b1rdp, b2rdp not significant, a4, b1rd also significant
  #psmeHeightFromDiameter$chapmanRichardsBalRelDbhB2 = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a2 * standBasalAreaPerHectare) * (1 - exp((b1)*dbh))^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger), psme2022, start = list(a1 = 57, a2 = 0.05, b1 = -0.02, b2 = 1.25, b2rd = 0.1, b2bal = -0.005))
  #psmeHeightFromDiameter$chapmanRichardsBalRelDbhB1 = fit_gsl_nls("Chapman-Richards BA+L RelDbh", height ~ 1.37 + a1 * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^(b2 + b2bal * basalAreaLarger), psme2022, start = list(a1 = 65, b1 = -0.02, b1rd = 0.001, b2 = 1.3, b2bal = -0.005))
  psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a2 * standBasalAreaPerHectare) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation + b2bal * basalAreaLarger), psme2022, start = list(a1 = 55, a2 = 0.1, b1 = -0.01, b1ac = -0.001, b2 = 1.3, b2p = 0.1, b2bal = -0.008))
  psmeHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger), psme2022, start = list(a1 = 61, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2rd = 0.1, b2bal = -0.005)) # a1, b1ba, b2ba not significant
  psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 60, b1 = -0.02, b1ac = -0.001, b2 = 1.1, b2p = 0.1)) # { a5, b1, b2 } x { e, s, tr, tw }, { a5, b2 } x a not significant, b1as not reliably significant
  psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psme2022, start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2rd = 0.1, b2rdp = 0.05)) # a4p, b2p not significant, b1rd NaN-inf, a4 also significant
  #psmeHeightFromDiameter$chapmanRichardsRelDbhA4 = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a4 * relativeDiameter) * (1 - exp((b1)*dbh))^(b2), psme2022, start = list(a1 = 72, a4 = -2.5, b1 = -0.02, b2 = 1.1))
  psmeHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psme2022, start = list(a1 = 58, b1 = -0.02, b1ac = -0.001, b2 = 1.1, b2rd = 0.1, b2rdp = 0.06))
  psmeHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 3, a1p = -1.4, b1 = 0.4, b1p = -0.16))
  psmeHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 73, a2 = 150, b1 = -1.3, b1p = 0.03)) # a1p not significant, a2p-b1p not mutually significant
  psmeHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 200, b1 = -6.8, b1p = -1.2, b2 = -0.3, b2p = -0.04)) # a1p-b1p not mutually significant
  psmeHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psme2022, start = list(a1 = 75, a2 = 140, a2p = 20, b1 = 1.2)) # a1p, b1p not significant
  psmeHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3 + a3p* isPlantation), psme2022, start = list(a1 = 0.012, a2 = 0.9, a3 = 1.4, a3p = 3.5)) # a1p, a2p not significant
  psmeHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation), psme2022, start = list(a1 = 2.5, a1p = -1, b1 = 0.65, b1p = 0.14))
  psmeHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(dbh + b2)), psme2022, start = list(a1 = 73, b1 = -46, b1p = -3, b2 = 10)) # a1p, b2p not significant
  psmeHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-kU * dbh)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2022, start = list(Ha = 55, Hap = -3, d = 0.3, dp = 0.2, kU = 0.013)) # kUp not significant
  psmeHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psme2022, start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0)) # a1p, b1p, b2p, b3, b3p, b4p not significant
  psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^b4, psme2022, start = list(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = -0.02, b4 = 1.0)) # a1p, b1p, b2p, b3, b3p, b4p not significant
  # A5, B1, B2, B3, b4 forms all quite close with preference to A5 over B1 only on AIC parsimony
  psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a5ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation), psme2022, start = list(a1 = 25, a5ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0, b4e = 0.0002))  # { a5, b1 } x { e, s, tr, tw }, a5 x { as, ac }, b1as, { b2, b3, b4 } x { s, tr, tw }, b3 not significant, { b2e, b3e }-b4e, { a4a, b1a, b2a, b3a }-b4a not mutually significant
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB1 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation), psme2022, start = list(a1 = 25, b1 = 0.2, b1as = 0.004, b1ac = 0.01, b2 = -0.02, b3 = -0.04, b4 = 1.0, b4e = 0.0002))
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB2 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2 + b2e * elevation + b2as * sin(3.14159/180 * aspect) + b2ac * cos(3.14159/180 * aspect))*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4), psme2022, start = list(a1 = 25, b1 = 0.3, b2 = -0.02, b2e = 0.000001, b2as = -0.001, b2ac = -0.001, b3 = -0.05, b4 = 1.0))
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB3 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3e * elevation + b3as * sin(3.14159/180 * aspect) + b3ac * cos(3.14159/180 * aspect))*dbh))^(b4), psme2022, start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3e = -0.0001, b3as = 0.01, b3ac = 0.01, b4 = 1.0))
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB4 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b3 = -0.02, b4 = 1.0, b4e = 0.0001, b4as = -0.035, b4ac = -0.04))
  psmeHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a4 * relativeDiameter)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4), psme2022, start = list(a1 = 25, a4 = 0.4, b1 = 0.3, b1as = -0.04, b1ac = 0.01, b2 = -0.02, b3 = -0.07, b4 = 1.0))
  psmeHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1 + a4 * relativeDiameter)*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3 * dbh))^b4, psme2022, start = list(a1 = 22, a4 = 0.4, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0)) # a4p, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rdp not significant, b1rd, b4rd not reliably significant
  # A5, B1, B2, B4 forms near interchangeable, B1 slightly preferred
  psmeHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), psme2022, start = list(a1 = 21, b1 = 0.3, b1as = 0.006, b1ac = 0.01, b2 = -0.02, b3 = 0, b4 = 1.0)) # b1e not significant # { a5, b1 } x { e, s, tr, tw }, a5ac, b1as, { b2, b3, b4 } x { s, tr, tw }, b3 not significant, { b2e, b3e }-b4e, { a4a, b1a, b2a, b3a }-b4a not mutually significant
  #psmeHeightFromDiameter$sharmaPartonPhysioA5 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a5ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4e * elevation), psme2022, start = list(a1 = 21, a5ac = 0.6, b1 = 0.3, b2 = -0.02, b3 = -0.01, b4 = 1.0, b4e = 0.0002)) # a5as not significant
  #psmeHeightFromDiameter$sharmaPartonPhysioB2 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2 + b2as * sin(3.14159/180 * aspect) + b2ac * cos(3.14159/180 * aspect))*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), psme2022, start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b2as = -0.001, b2ac = -0.001, b3 = 0, b4 = 1.0)) # b2e not significant
  #psmeHeightFromDiameter$sharmaPartonPhysioB3 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3as * sin(3.14159/180 * aspect) + b3ac * cos(3.14159/180 * aspect))*dbh))^(b4), psme2022, start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3e = -0.0001, b3as = 0.01, b3ac = 0.01, b4 = 1.0)) # b3e not significant, job singular gradient
  #psmeHeightFromDiameter$sharmaPartonPhysioB4 = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04))
  psmeHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psme2022, start = list(a1 = 13, b1 = 0.3, b1rd = 0, b2 = -0.02, b3 = 0.02, b4 = 1.0), significant = FALSE) # a4, a4p, b1rd, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rd, b4rdp not significant
  psmeHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1rd * relativeDiameter) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 18, b1 = 0.3, b1rd = 0.004, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04), significant = FALSE) # b1rd, b3 not significant
  psmeHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, psme2022, start = list(a1 = 45, a1p = -12, b1 = 0.08, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)) # b2p, b3p, b4p not significant, 10x10 AIC 3803
  psmeHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a2 * standBasalAreaPerHectare + a3 * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp((b2)*standTreesPerHectare^b3 * dbh))^(b4), psme2022, start = list(a1 = 50, a1p = -10, a2 = -0.1, a3 = 0.2, b1 = 0.1, b1p = 0.1, b2 = -0.03, b3 = -0.02, b4 = 1.1)) # a2p, a3p, { b1, b3 } x { bal, balp }, b1ba, b3, b4bap, b4balp not significant, b4ba, b4bal also significant, 10x10 AIC 3779
  #psmeHeightFromDiameter$sharmaZhangBalB4 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp((b2)*standTreesPerHectare^b3 * dbh))^(b4 + b4ba * standBasalAreaPerHectare + b4bal * basalAreaLarger), psme2022, start = list(a1 = 30, a1p = -10, b1 = 0.1, b1p = 0.1, b2 = -0.03, b3 = -0.02, b4 = 1.1, b4ba = 0.01, b4bal = -0.01))
  psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 0.2, b1 = 2.5, b1p = -0.2, b2 = -0.14, b2p = 0.02)) # a1p not significant
  psmeHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psme2022, start = list(a1 = 56, b1 = -0.01, b1p = 0.001, b2 = 1.1)) # a1p, b2p not significant
  # A3 and B2 forms MAE and AIC interchangeable
  psmeHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger))), psme2022, start = list(a1 = 60, b1 = -0.01, b2 = 1.0, b2bal = 0.002)) # a3p, b1ba, b1bap, b1bal, b2ba, b2balp, b2bal, b2balp not significant, a3 also significant
  #psmeHeightFromDiameter$weibullBalA3 = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a3 * basalAreaLarger) * (1 - exp(b1*dbh^(b2))), psme2022, start = list(a1 = 60, a3 = 0, b1 = -0.01, b2 = 1.0))
  #print(to_parameter_confidence_intervals(psmeHeightFromDiameter$weibullBal), n = 40)
  #to_fixed_coeffficients(psmeHeightFromDiameter$weibullBal)
  
  # all GAM forms NaN with Γ link
  psmeHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, bs = "ts", k = 9), data = psme2022, threads = 4) # 10x10 AIC 3918, MAE ~3.42 with power base, plantation not significant, specification of 1.37 m offset not supported (invalid model formula in ExtractVars, mgcv 1.9-3)
  #psmeHeightFromDiameter$gamHossfeld = fit_gam("REML GAM", height ~ I(71.61941/(1 + 158.72979 * dbh^(-1.26072 + 0.02926 * isPlantation))) + s(dbh, bs = "ts", k = 7), data = psme2022, threads = 4) # 10x10 AIC ~3924, MAE ~3.40 with Hossfeld IV
  #psmeHeightFromDiameter$gamMichaelis = fit_gam("REML GAM", height ~ I(72.299 * dbh^1.231/(143.174 + 18.141 * isPlantation + dbh^1.231)) + s(dbh, bs = "ts", k = 7), data = psme2022, threads = 4) # 10x10 AIC 3902, MAE ~3.44 with Michaelis-Menten
  #psmeHeightFromDiameter$gamS0 = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = psme2022, threads = 4) # 10x10 AIC ~3897, MAE ~3.46
  #psmeHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, threads = 4) # 10x10 AIC 3881
  psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, threads = 4) # 10x10 AIC 3825
  #psmeHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 4) # 10x10 AIC 3836
  psmeHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 14), data = psme2022, threads = 4) # 10x10 AIC 3818, plantation not significant
  psmeHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, threads = 4) # 10x10 AIC 3898
  
  if (psmeOptions$fitPhysioGams)
  {
    psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, threads = 4) # 10x10 AIC 3789
    psmeHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 19), data = psme2022, threads = 4) # 10x10 AIC 3848
    # modest reductions in MAE from slope, aspect, roughness, and wetness, including elevation is detrimental
    #psmeHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, elevation, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, threads = 8, significant = FALSE) # 10x10 AIC 4125, support at k ~= 33 but clearly overfit, k = 10 minimum for curvature
    #psmeHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 15), data = psme2022, threads = 8) # 10x10 AIC 3856
    #psmeHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 8) # 10x10 AIC 3937
    #psmeHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 8) # 10x10 AIC3948
    #psmeHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 8) # 10x10 AIC 3940
    #psmeHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 8) # 10x10 AIC 3941
    #psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 85), data = psme2022, threads = 4) # 10x10 AIC 4220 + higher MAE than slope or wetness alone @ k ~= 61 -> arguably overfit, fails on dplyr errors with threads = 8, warnings at threads = 4
    #psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), data = psme2022, threads = 4) # 10x10 AIC 3945 + higher MAE than slope or wetness alone
    psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 10), data = psme2022, threads = 8) # 10x10 AIC 3941
    psmeHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, topographicWetnessFD8f, bs = "ts", k = 12, by = as.factor(isPlantation)), data = psme2022, threads = 4) # 10x10 AIC 3925
    #lapply(psmeHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
    #lapply(psmeHeightFromDiameter$gamRelDbhPhysio$fit, summary)
    #psmeHeightFromDiameter$gam$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
    
    saveRDS(list(gamBalPhysio = psmeHeightFromDiameter$gamBalPhysio,
                 gamBalPhysioRelDbh = psmeHeightFromDiameter$gamBalPhysioRelDbh,
                 gamPhysio = psmeHeightFromDiameter$gamPhysio,
                 gamRelDbhPhysio = psmeHeightFromDiameter$gamRelDbhPhysio), 
            paste0("trees/height-diameter/data/PSME height physio GAMs ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  } else {
    psmePrimaryHeightGams = readRDS(paste0("trees/height-diameter/data/PSME height physio GAMs ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
    psmeHeightFromDiameter$gamBalPhysio = psmePrimaryHeightGams$gamBalPhysio
    psmeHeightFromDiameter$gamBalPhysioRelDbh = psmePrimaryHeightGams$psmeHeightFromDiameter$gamBalPhysioRelDbh
    psmeHeightFromDiameter$gamPhysio = psmePrimaryHeightGams$gamPhysio
    psmeHeightFromDiameter$gamRelDbhPhysio = psmePrimaryHeightGams$gamRelDbhPhysio
    rm(psmePrimaryHeightGams)
  }
  
  psmeHeightFromDiameter$randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"),
                                                   psme2022, mtry = 2, minNodeSize = 48, sampleFraction = 0.230) # from setup.R tuneRanger @ 500 trees
  psmeHeightFromDiameter$randomForestBalPhysioRelDbh = fit_ranger("random forest BA+L RelDbh physio", c("height", "dbh", "relativeDiameter", "standQmd", "topHeight", "standBasalAreaPerHectare", "basalAreaLarger", "terrainRoughness", "slope", "elevation", "topographicWetnessFD8f"), # from setup.R VSURF
                                                                  psme2022, mtry = 2, minNodeSize = 2, sampleFraction = 0.785) # from setup.R tuneRanger @ 500 trees
  
  # saving primary Rds at the end of this block ensures physio GAMs are included
  saveRDS(psmeHeightFromDiameter, paste0("trees/height-diameter/data/PSME height ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## Douglas-fir height mixed regressions
if (psmeOptions$fitHeightMixed)
{
  psmeHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1r) *(1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, # b1r, b2r singularity in backsolve
                                                                fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 45, b1 = -0.02, b2 = 1.2, b2p = -0.1))))
  psmeHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1r + a2 * standBasalAreaPerHectare) * (1 - exp(b1*dbh))^(b2 + b2bal * basalAreaLarger), psme2022, # b1r, b2r singularity in backsolve, b2balp not significant
                                                            fixedFormula = a1 + a2 + b1 + b2 + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 38, a2 = 0.2, b1 = -0.02, b2 = 1.2, b2bal = -0.003)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalRelDbh = fit_nlme("Chapman-Richards BA+L RelDbh", height ~ 1.37 + (a1 + a1r + a4 * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2bal * basalAreaLarger), psme2022,
                                                                  fixedFormula = a1 + a4 + b1 + b2 + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 70, a4 = -1.7, b1 = -0.015, b2 = 1.3, b2bal = -0.005)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + a2 * standBasalAreaPerHectare) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation + b2bal * basalAreaLarger), psme2022, # b1r, b2r singularity in backsolve, b2balp not significant
                                                                  fixedFormula = a1 + a2 + b1 + b1ac + b2 + b2p + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 40, a2 = 0.3, b1 = -0.03, b1ac = -0.002, b2 = 1.4, b2p = 0.1, b2bal = -0.005)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = fit_nlme("Chapman-Richards BA+L RelDbh physio", height ~ 1.37 + (a1 + a1r) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger), psme2022, # b1r, b2r singularity in backsolve
                                                                        fixedFormula = a1 + b1 + b1ac + b2 + b2rd + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 48, b1 = -0.03, b1ac = -0.001, b2 = 1.3, b2rd = 0.3, b2bal = -0.003)))
  psmeHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), psme2022, 
                                                               fixedFormula = a1 + b1 + b1ac + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 45, b1 = -0.02, b1ac = -0.002, b2 = 1.2, b2p = 0.1)))
  psmeHeightFromDiameterMixed$chapmanRichardsRelDbh = fit_nlme("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1r) *(1 - exp(b1*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psme2022, # b1r, b2r singularity in backsolve
                                                                fixedFormula = a1 + b1 + b2 + b2rd + b2rdp ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 46, b1 = -0.03, b2 = 1.2, b2rd = 0.4, b2rdp = -0.1)))
  psmeHeightFromDiameterMixed$chapmanRichardsRelDbhPhysio = fit_nlme("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1r) *(1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter), psme2022, # b1r singularity in backsolve, b2r max iterations, b2rdp not significant
                                                                     fixedFormula = a1 + b1 + b1ac + b2 + b2rd ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 46, b1 = -0.03, b1ac = -0.002, b2 = 1.2, b2rd = 0.4)))
  psmeHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation + b1r), psme2022, # a1r also viable but slower to converge
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 3, a1p = -0.8, b1 = 0.4, b1p = -0.1)))
  psmeHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + (a2) * dbh^(b1 + b1p * isPlantation + b1r)), psme2022, # b1r clearly preferable to a1r
                                                  fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 75, a2 = 80, b1 = -1.1, b1p = 0.05)))
  #psmeHeightFromDiameterMixed$hossfeldA1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1r) / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psme2022,
  #                                                fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 50, a2 = 90, b1 = -1.3, b1p = 0.01)))
  #psmeHeightFromDiameterMixed$hossfeldA2 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + (a2 + a2r) * dbh^(b1 + b1p * isPlantation)), psme2022, # job false convergence
  #                                                fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a2r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 75, a2 = 75, b1 = -1.0, b1p = 0.06)))
  psmeHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psme2022, # b1r max iterations, b2r step halving
                                              fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                              start = list(fixed = c(a1 = 90, b1 = -8, b1p = 2, b2 = -0.5, b2p = 0.1)))
  psmeHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1) / (a2 + a2p * isPlantation + a2r + dbh^(b1)), psme2022, # a2r clearly preferable to a1r
                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a2r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 50, a2 = 100, a2p = -20, b1 = 1.3)))
  #psmeHeightFromDiameterMixed$michaelisMentenA1 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1r)*dbh^(b1) / (a2 + a2p * isPlantation + dbh^(b1)), psme2022,
  #                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 50, a2 = 100, a2p = -17, b1 = 1.3)))
  #psmeHeightFromDiameterMixed$michaelisMentenB1 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + b1r) / (a2 + a2p * isPlantation + dbh^(b1 + b1r)), psme2022, # job false convergence
  #                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 75, a2 = 65, a2p = 11, b1 = 1.0)))
  psmeHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1)*dbh^2 + (a2 + a2r)*dbh + a3), psme2022, # a3r max iterations, a3, a3p not significant
                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.014, a2 = 0.9, a3 = 1)))
  #psmeHeightFromDiameterMixed$prodanA1 = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1r)*dbh^2 + (a2)*dbh + a3), psme2022, # a2r clearly preferable
  #                                              fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 0.018, a2 = 0.7, a3 = 2.7)))
  psmeHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation + a1r)*dbh^(b1 + b1p * isPlantation), psme2022, # a1r, b1r not distinguishable
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 1.9, a1p = -0.5, b1 = 0.7, b1p = 0.07)))
  #psmeHeightFromDiameterMixed$powerB1 = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation + b1r), psme2022,
  #                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2.2, a1p = -0.7, b1 = 0.67, b1p = 0.08)))
  psmeHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1 + b1p * isPlantation + b1r)/(dbh + b2)), psme2022, # b1r clearly preferable
                                                   fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 63, b1 = -35, b1p = -4, b2 = 8)))
  #psmeHeightFromDiameterMixed$ratkowskyA1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1r)*exp((b1 + b1p * isPlantation)/(dbh + b2)), psme2022,
  #                                                   fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 53, b1 = -29, b1p = 3, b2 = 5)))
  psmeHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), psme2022, # kUr max iterations, dr step halving, dp not significant
                                                   fixedFormula = Ha + Hap + d + kU ~ 1, randomFormula = Har ~ 1|stand/plot,
                                                   start = list(fixed = c(Ha = 47, Hap = -5, d = 0.4, kU = 0.018)))
  psmeHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1r)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, psme2022, # a1r, b2r, b4r singularity in backsolve, b3 not significant
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonB3 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3r)*dbh))^b4, psme2022, # false convergence
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^b4, psme2022, # a1r, b2r, b4r singularity in backsolve
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 40, b1 = 0.1, b2 = -0.02, b3 = -0.02, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalB3 = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3r)*dbh))^b4, psme2022, # b1r clearly preferable
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 40, b1 = 0.1, b2 = -0.02, b3 = -0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a5ac * cos(3.14159/180 * aspect))*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation), psme2022, # a1r, b24, b4r singularity in backsolve
                                                               fixedFormula = a1 + a5ac + b1 + b2 + b3 + b4 + b4e ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 40, a5ac = 1.5, b1 = 0.1, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4e = 0.0003)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysioB3 = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a5ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3r)*dbh))^(b4 + b4e * elevation), psme2022, # inclination to a1-b1 evaporation
  #                                                             fixedFormula = a1 + a5ac + b1 + b2 + b3 + b4 + b4e ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 25, a5ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0, b4e = 0.0002)))
  psmeHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1r + b1rd * relativeDiameter + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4), psme2022, # a1r, b2r, b4r singularity in backsolve, b3r max iterations, a4 not significant, b2rd, b3rd, b4rd all slow to converge compared to b1rd
                                                                     fixedFormula = a1 + b1 + b1rd + b1as + b1ac + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 35, b1 = 0.1, b1rd = 0.002, b1as = 0.005, b1ac = 0.01, b2 = -0.02, b3 = -0.15, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton BA+L RelDbh", height ~ 1.37 + (a1)*topHeight^(b1 + b1r) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4rd * relativeDiameter), psme2022, # b2r max iterations, b3r, b4r singularity in backsolve, a4 not significant
                                                               fixedFormula = a1 + b1 + b2 + b3 + b4 + b4rd ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 40, b1 = 0.1, b2 = -0.02, b3 = -0.2, b4 = 1.0, b4rd = -0.03)))
  psmeHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, # a1r, b4r singularity in backsolve, b2r step halving, b3 not significant
                                                            fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac + b4as ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  #psmeHeightFromDiameterMixed$sharmaPartonPhysioB3 = fit_nlme("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3r)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, # job false convergence
  #                                                          fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac + b4as ~ 1, randomFormula = b3r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  #psmeHeightFromDiameterMixed$sharmaPartonRelDbh = fit_nlme("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1r + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, psme2022, # fixed effects not significant
  #                                                          fixedFormula = a1 + b1 + b1rd + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 13, b1 = 0.3, b1rd = 0, b2 = -0.02, b3 = 0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonRelDbh = create_fit_statistics("Sharma-Parton RelDbh", fitting = "nlme", significant = FALSE)
  #psmeHeightFromDiameterMixed$sharmaPartonRelDbhPhysio = fit_nlme("Sharma-Parton RelDbh physio", height ~ 1.37 + a1*topHeight^(b1 + b1r + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, # fixed effects not significant
  #                                                                fixedFormula = a1 + b1 + b1rd + b2 + b3 + b4 + b4as + b4ac ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 18, b1 = 0.3, b1rd = 0.004, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  psmeHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*standTreesPerHectare^(b3 + b3r)*dbh))^(b4), psme2022, # a1r, b4r singularity in backsolve, b2r step halving, a1p not significant
                                                     fixedFormula = a1 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b3r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 35, b1 = 0.1, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaZhangB1 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1r)*(1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4), psme2022, # b3r preferable
  #                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 35, a1p = -10, b1 = 0.1, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a2 * standBasalAreaPerHectare + a3 * basalAreaLarger) * standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1r) * (1 - exp(b2 * standTreesPerHectare^(b3)*dbh))^(b4), psme2022, # a1r, b4r singularity in backsolve, b2r step halving, b2ba, b2bap not significant
                                                        fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 50, a1p = -70, a2 = 40, a3 = 11, b1 = -0.9, b1p = 0.01, b2 = -0.02, b3 = -0.03, b4 = 1.04)))
  psmeHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + (a1)*dbh^((b1 + b1p * isPlantation + b1r)*dbh^b2), psme2022, # a1r step halving, b2p not significant
                                                  fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.5, b1 = 1.8, b1p = -0.05, b2 = -0.13)))
  #psmeHeightFromDiameterMixed$sibbesenB2 = fit_nlme("Sibbesen", height ~ 1.37 + (a1)*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2r)), psme2022, # b1r preferable
  #                                                fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.5, b1 = 1.8, b1p = -0.05, b2 = -0.13)))
  psmeHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1r)*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psme2022, # b1r, b2r max iterations
                                                 fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 45, b1 = -0.01, b1p = 0.001, b2 = 1.1)))
  psmeHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger + b2r))), psme2022, # a1r max iterations, b1r step halving
                                                    fixedFormula = a1 + b1 + b2 + b2bal ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 100, b1 = -0.01, b2 = 0.9, b2bal = 0.003)))
  #to_fixed_coeffficients(psmeHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(psmeHeightFromDiameterMixed$sharmaZhang), n = 40)

  # fit with gamm() because gam() with random effect basis functions is unhappy, gamma links not used as convergence fails
  # basis  fit time, 9900X  threads  warnings
  # fs     ~80 s            4        yes
  # re     ~80 s            4        yes
  # sz     ~100 s           4        yes
  # mrf is impractical due to requiring lists of coordinates forming plot polygons
  psmeHeightFromDiameterMixed$gam = fit_gam("REML GAM", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, bs = "ts", k = 8), random = list(uniquePlotID = ~1), data = psme2022) # 10x10 AIC 4144, median MAE 3.60 m vs ~3.41 m for fixed effects
  psmeHeightFromDiameterMixed$gamBal = fit_gam("REML GAM BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13), random = list(uniquePlotID = ~1), data = psme2022)
  psmeHeightFromDiameterMixed$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 13), random = list(uniquePlotID = ~1), data = psme2022)
  psmeHeightFromDiameterMixed$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(uniquePlotID = ~1), data = psme2022)
  psmeHeightFromDiameterMixed$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 13), random = list(uniquePlotID = ~1), data = psme2022)
  psmeHeightFromDiameterMixed$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 18), random = list(uniquePlotID = ~1), data = psme2022)
  psmeHeightFromDiameterMixed$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(uniquePlotID = ~1), data = psme2022)
  psmeHeightFromDiameterMixed$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, topographicWetnessFD8f, bs = "ts", k = 13, by = as.factor(isPlantation)), random = list(uniquePlotID = ~1), data = psme2022)
  #lapply(psmeHeightFromDiameterMixed$gam$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(psmeHeightFromDiameterMixed$gam$fit, function(fit) { return(summary(fit$gam)) })
  #psmeHeightFromDiameterMixed$gam$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(psmeHeightFromDiameterMixed, paste0("trees/height-diameter/data/PSME height mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## Douglas-fir diameter regressions
if (psmeOptions$fitDbhPrimary)
{
  # linear regressions
  psmeDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37) + I(isPlantation*(height - 1.37)), psme2022))
  psmeDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I((height - 1.37)^2) + I(isPlantation*(height - 1.37)^2), psme2022) # plantation not significant
  
  # nonlinear regressions
  psmeDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation)*(exp(b1*(height - 1.37)) - 1)^b2, psme2022, start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8)) # b1p, b2p not significant
  #psmeDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 75, a1p = -35, a2 = -0.07, a3 = 0.5, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07))
  psmeDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1)*(exp((b1)*(height - 1.37)^(b2 + b2rh * relativeHeight)) - 1), psme2022, start = list(a1 = 45, b1 = 0.04, b2 = 0.8, b2rh = 0.02)) # a1p-b1rh, b2rh not mutually significant
  psmeDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -100, b1 = 0.012, b2 = 0.9)) # a1p, b1p, b2p not significant
  #psmeDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, start = list(a1 = -175, a1p = 100, a2 = 0.14, a3 = 1.2, b1 = 0.007, b1p = 0.0057, b2 = 0.79))
  psmeDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), psme2022, start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8, b2tw = 0), significant = FALSE) # { a5, b1, b2 } x { e, s, a, tr, tw } not significant
  # no clear differentiation between A4, B1, B2 forms, B2 maybe somewhat more accurate
  psmeDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), psme2022, start = list(a1 = -90, b1 = 0.012, b2 = 0.9, b2rh = -0.07)) # a4, b1rh also significant
  #psmeDiameterFromHeight$chapmanRichardsRelHtA4 = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a4 * relativeHeight)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), psme2022, start = list(a1 = -90, a4 = -11, b1 = 0.012, b2 = 0.9))
  #psmeDiameterFromHeight$chapmanRichardsRelHtB1 = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1)*log(1 - pmin(((b1 + b1rh * relativeHeight)*(height - 1.37))^(b2), 0.9999)), psme2022, start = list(a1 = -100, b1 = 0.01, b1rh = 0.001, b2 = 0.9))
  psmeDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 130, a2 = 70, b1 = 0.9, b1p = -0.008)) # a1p, a2p not significant
  psmeDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psme2022, start = list(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006))
  psmeDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2))
  #psmeDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + (a3 + a3p * isPlantation) * basalAreaTaller)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, start = list(a1 = 2.14, a1p = -0.051, a2 = 0.00085, a3 = -0.0065, a3p = -0.0038, b1 = 0.963, b1p = -0.053))
  psmeDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1e * elevation), psme2022, start = list(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2, b1e = 0), significant = FALSE) # { a5, b1 } x { e, s, a, tr, tw } not significant
  psmeDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight), psme2022, start = list(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)) # a1p-a4 not mutually significant
  psmeDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2p = -0.001)) # a1p not significant, b1p-b2p not mutually significant
  #psmeDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a3 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013))
  #psmeDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a2 = -0.004, a3 = -0.017, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012))
  #psmeDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", dbh ~ (a1 + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.5, a3 = -0.02, a4 = 0.2, b1 = 0.91, b1p = -0.16, b2 = 0, b2p = 0.009))
  #psmeDiameterFromHeight$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark ABA+T RelHt physio", dbh ~ (a1 + a2 * basalAreaTaller + a3 * bootstrapStandBasalAreaPerHectare + a6 * cos(3.14159/180 * aspect) + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.4, a2 = -0.013, a3 = -0.0037, a6 = -0.05, a4 = 0.1, b1 = 0.94, b1p = -0.16, b2 = 0.001, b2p = 0.01))
  # A5 modestly preferable to B1 and B2
  psmeDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ (a1 + a5ac * cos(pi/180*aspect))*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psme2022, start = list(a1 = 2, a5ac = -0.1, b1 = 0.8, b2 = 0.01)) # { a6, b1, b2 } x { e, s, as, tr, tw }, b2p not significant, b1ac, b2ac also  significant
  #psmeDiameterFromHeight$ruarkPhysioB1 = fit_gsl_nls("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1ac * cos(pi/180*aspect)) * exp((b2) * (height - 1.37)), psme2022, start = list(a1 = 2, b1 = 0.8, b1ac = -0.015, b2 = 0.01))
  #psmeDiameterFromHeight$ruarkPhysioB2 = fit_gsl_nls("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2ac * cos(pi/180*aspect)) * (height - 1.37)), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2ac = -0.001))
  # A4, B1, B2 not entirely separable but B2 preferred
  psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)) # a4, b1rh also significant, b2rh-b4p not mutually significant
  #psmeDiameterFromHeight$ruarkRelHtA4 = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a4 * relativeHeight)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psme2022, start = list(a1 = 2, a4 = 0.3, b1 = 0.8, b2 = 0.008))
  #psmeDiameterFromHeight$ruarkRelHtB1 = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2) * (height - 1.37)), psme2022, start = list(a1 = 2, b1 = 0.8, b1rh = 0.03, b2 = 0.008))
  psmeDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a5ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psme2022, start = list(a1 = 2, a5ac = -0.09, b1 = 0.8, b2 = 0.009, b2rh = 0.003))
  #psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1), 0.9999)), psme2022, start = list(a1 = 0.0001, a2 = 0.55, b1 = 1.05, Ha = 90)) # always nan-inf or singular gradient
  psmeDiameterFromHeight$schnute = create_fit_statistics(name = "Schnute inverse", fitting = "gsl_nls")
  psmeDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^b4, psme2022, start = list(a1 = 3, b1 = 0.7, b2 = 0.1, b3 = 0.03, b4 = 0.3)) # a1p, b1p, b2p, b3p, b4p not significant, some a1-b3 evaporation
  psmeDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)) # b1p not significant
  #psmeDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 4.4, a1p = -1.6, a2 = -0.006, a3 = -0.03, b1 = 0.61, b2 = 0.071, b2p = 0.021))
  #psmeDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077))
  #psmeDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.006, a3 = -0.01, a4 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021))
  #psmeDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022, start = list(a1 = 3.1, a1p = -0.5, a2 = -0.005, a3 = -0.01, a6 = -0.06, a4 = 0.1, b1 = 0.65, b2 = 0.07))
  # A5 and B2 similar, B2 probably slightly preferred
  psmeDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2ac * cos(pi/180*aspect))), psme2022, start = list(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2ac = -0.003)) # { a6, b1, b2 } x { e, s, as, tr, tw } not significant, a5ac also significant
  #psmeDiameterFromHeight$sibbesenReplacePhysioA5 = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation + a5ac * cos(pi/180*aspect))*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 1.6, a1p = 0.4, a5ac = -0.07, b1 = 0.7, b2 = 0.1, b2p = -0.02)) # { a6, b1, b2 } x { e, s, as, tr, tw } not significant, a5ac also significant
  # A4, B1, B2 not clearly separable
  psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight)), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01)) # a1p not significant, a4, b1rh also significant
  #psmeDiameterFromHeight$sibbesenReplaceRelHtA4 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a4 * relativeHeight)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 2, a4 = 0.4, b1 = 0.7, b2 = 0.07, b2p = -0.006))
  #psmeDiameterFromHeight$sibbesenReplaceRelHtB1 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1 + b1rh * relativeHeight)*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, start = list(a1 = 2, b1 = 0.8, b1rh = 0.03, b2 = 0.06, b2p = -0.006))
  psmeDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2ac * cos(pi/180*aspect))), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01, b2ac = -0.003)) # a1p not significant
  psmeDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, psme2022, start = list(a1 = -140, b1 = 0.01, b2 = 0.9)) # a1p, b1p, b2p not significant
  #print(to_parameter_confidence_intervals(psmeDiameterFromHeight$weibull), n = 40)
  #to_fixed_coeffficients(psmeDiameterFromHeight$weibull)
  
  # all GAM forms NaN with Γ link
  psmeDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), data = psme2022, threads = 2) # 10x10 MAE 14.9, AIC 5258
  #psmeDiameterFromHeight$gamRuark = fit_gam("REML GAM", dbh ~ I(1.92785 * (height - 1.37)^(0.84484) * exp((0.01328 - 0.00111 * isPlantation) * (height - 1.37))) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), data = psme2022, threads = 2) # 10x10 MAE 15.1, AIC 5297
  #psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 28), data = psme2022, threads = 2)
  psmeDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 8), data = psme2022, threads = 2) # 10x10 AIC 5221

  if (psmeOptions$fitPhysioGams)
  {
    #psmeDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), data = psme2022, threads = 2)
    #psmeDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ I(height^) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496), data = psme2022, threads = 2)
    # plausible but modest improvement with any physiographic predictor but all predictor min k = 331 > ~150 -> drop elevation + cos(aspect) on highest RMSE min k = 57 > ~42 -> drop slope
    #psmeDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, elevation, bs = "ts", by = as.factor(isPlantation), k = 9), data = psme2022, threads = 4) # 10x10 AIC 5280
    #psmeDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, slope, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 13), data = psme2022, threads = 4) # 10x10 AIC 5304
    #psmeDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 8), data = psme2022, threads = 4) # 10x10 AIC 5261
    #psmeDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 8), data = psme2022, threads = 4) # 10x10 AIC 5242
    #psmeDiameterFromHeight$gamRoughness = fit_gam("REML GAM roughness", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 8), data = psme2022, threads = 4) # 10x10 AIC 5200
    #psmeDiameterFromHeight$gamWetness = fit_gam("REML GAM wetness", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 9), data = psme2022, threads = 4) # 10x10 AIC 5119
    psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 18), data = psme2022, threads = 4) # 10x10 AIC 5398
    # minimum k = 57 > ~45 -> drop wetness
    psmeDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), data = psme2022, threads = 4) # minimum k = 16 acceptable, 10x10 AIC 
    #lapply(psmeDiameterFromHeight$gamRelHtPhysio$fit, k.check)
    #lapply(psmeDiameterFromHeight$gamRelHtPhysio$fit, summary)
    #psmeDiameterFromHeight$gamRelHtPhysio$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
    
    saveRDS(list(gamAbatPhysio = psmeDiameterFromHeight$gamAbatPhysio,
                 gamAbatPhysioRelHt = psmeDiameterFromHeight$gamAbatPhysioRelHt,
                 gamPhysio = psmeDiameterFromHeight$gamPhysio,
                 gamRelHtPhysio = psmeDiameterFromHeight$gamRelHtPhysio), 
            paste0("trees/height-diameter/data/PSME DBH physio GAMs ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  } else {
    psmePrimaryDbhGams = readRDS(paste0("trees/height-diameter/data/PSME DBH physio GAMs ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
    psmeDiameterFromHeight$gamAbatPhysio = psmePrimaryDbhGams$gamAbatPhysio
    psmeDiameterFromHeight$gamAbatPhysioRelHt = psmePrimaryDbhGams$gamAbatPhysioRelHt
    psmeDiameterFromHeight$gamPhysio = psmePrimaryDbhGams$gamPhysio
    psmeDiameterFromHeight$gamRelHtPhysio = psmePrimaryDbhGams$gamRelHtPhysio
    rm(psmePrimaryDbhGams)
  }
  
  psmeDiameterFromHeight$randomForest = fit_ranger("random forest", c("dbh", "height"),
                                                   psme2022, mtry = 1, minNodeSize = 97, sampleFraction = 0.307) # from setup.R tuneRanger @ 500 trees
  #psmeDiameterFromHeight$randomForestAbatRelHt = fit_ranger("random forest ABA+T RelHt", c("dbh", "height", "bootstrapStandBasalAreaPerHectare", "relativeHeight"), # from setup.R VSURF
  #                                                          psme2022, mtry = 3, minNodeSize = 2, sampleFraction = 0.328) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(psmeDiameterFromHeight, paste0("trees/height-diameter/data/PSME DBH ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}

if (psmeOptions$fitDbhMixed)
{
  #psmeDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation) * (exp((b1)*(height - 1.37)) - 1)^(b2 + b2r), psme2022, # a1r, b2r max iterations, b1r singularity in backsolve
  #                                                             fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot, # |stand, |uniquePlotID max iterations
  #                                                             start = list(fixed = c(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8))))
  psmeDiameterFromHeightMixed = list(chapmanReplace = create_fit_statistics(name = "Chapman-Richards replace", fitting = "nlme"))
  #psmeDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller)*(exp((b1 + b1p * isPlantation)*(height - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2022, 
  #                                                          fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 75, a1p = -35, a2 = -0.07, a3 = 0.5, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07)))
  #psmeDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ a1*(exp(b1*(height - 1.37)^(b2 + b2r + b2rh * relativeHeight)) - 1), psme2022, # a1r, b1r, b2r step halving
  #                                                           fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 45, b1 = 0.04, b2 = 0.8, b2rh = 0.02)))
  psmeDiameterFromHeightMixed$chapmanReplaceRelHt = create_fit_statistics(name = "Chapman-Richards replace RelHt", fitting = "nlme")
  #psmeDiameterFromHeightMixed = list(chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation + a1r)*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), psme2022, # a1r job step halving, b1r singular precision, b2r step halving
  #                                                              fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = -10, a1p = -1, b1 = 0.03, b2 = 0.3))))
  psmeDiameterFromHeightMixed$chapmanRichards = create_fit_statistics(name = "Chapman-Richards inverse", fitting = "nlme")
  #psmeDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller)*log(1 - pmin(((b1 + b1p * isPlantation)*(height - 1.37))^b2, 0.9999)), psme2022, 
  #                                                           fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = -175, a1p = 100, a2 = 0.14, a3 = 1.2, b1 = 0.007, b1p = 0.0057, b2 = 0.79)))
  #psmeDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1r + a1p * isPlantation)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), psme2022, # fixed effects not significant
  #                                                             fixedFormula = a1 + a1p + a5 + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8, b2tw = 0)), significant = FALSE)
  psmeDiameterFromHeightMixed$chapmanRichardsPhysio = create_fit_statistics(name = "Chapman-Richards inverse physio", fitting = "nlme", significant = FALSE)
  psmeDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r)*log(1 - pmin((b1*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), psme2022, # b1r step halving, b2r singularity in backsolve
                                                              fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = -20, b1 = 0.018, b2 = 0.5, b2rh = -0.3)))
  psmeDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation + b1r) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1r)), psme2022, # a1r step halving, a2r singularity in backsolve
                                                                fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 40, a2 = 18, b1 = 0.6, b1p = 0.008)))
  psmeDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1r + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psme2022, # a2r also tractable but false convergence with both a1r and a2r
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006)))
  #psmeDiameterFromHeightMixed$naslundA2 = fit_nlme("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation + a2r) * sqrt(height - 1.37)), psme2022, # 30x higher RMSE than a1r
  #                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a2r ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006)))
  psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation + b1r), psme2022, # a1r also tractable
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 0.7, a1p = 0.4, b1 = 1.2, b1p = -0.2)))
  #psmeDiameterFromHeightMixed$powerA1 = fit_nlme("power", dbh ~ (a1 + a1p*isPlantation + a1r)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, # false convergence but possibly lower mean or median error
  #                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.7, a1p = 0.4, b1 = 1.2, b1p = -0.2)))
  #psmeDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + (a3 + a3p * isPlantation) * basalAreaTaller)*(height - 1.37)^(b1 + b1p*isPlantation), psme2022, 
  #                                                 fixedFormula = a1 + a1p + a2 + a2p + a3 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 2.14, a1p = -0.051, a2 = 0.00085, a3 = -0.0065, a3p = -0.0038, b1 = 0.963, b1p = -0.053)))
  #psmeDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation), psme2022, # fixed effects not significant
  #                                                   fixedFormula = a1 + a1p + a6 + b1 + b1p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2, a6 = 0)), significant = FALSE)
  psmeDiameterFromHeightMixed$powerPhysio = create_fit_statistics(name = "power physio", fitting = "nlme", significant = FALSE)
  psmeDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight + b1r), psme2022, # a1r also tractable
                                                    fixedFormula = a1 + a1p + b1 + b1p + b1rh ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)))
  #psmeDiameterFromHeightMixed$powerRelHtA1 = fit_nlme("power RelHt", dbh ~ (a1 + a1p * isPlantation + a1r)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight), psme2022, # false convergence
  #                                                    fixedFormula = a1 + a1p + b1 + b1p + b1rh ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)))
  psmeDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1 + b1r) * exp((b2) * (height - 1.37)), psme2022, # a1r, b2r also tractable, b2p not significant
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkA1 = fit_nlme("Ruark", dbh ~ (a1 + a1r) * (height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psme2022, # false convergence
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkB2 = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1) * exp((b2 + b2r) * (height - 1.37)), psme2022, # not discriminable from b1r
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1r + a2 * basalAreaTaller)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
  #                                                 fixedFormula = a1 + a2 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 2.5, a3 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013)))
  #psmeDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1r + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
  #                                                       fixedFormula = a1 + a2 + a3 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 2.5, a2 = -0.004, a3 = -0.017, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012)))
  #psmeDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", dbh ~ (a1 + a1r + a2 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
  #                                                      fixedFormula = a1 + a2 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 2.5, a3 = -0.027, a4 = 0.3, b1 = 0.91, b1p = -0.2, b2 = 0, b2p = 0.013)))
  #psmeDiameterFromHeightMixed$ruarkAbatRelHtPhysio = fit_nlme("Ruark ABA+T RelHt physio", dbh ~ (a1 + a1r + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect) + a4*relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022,
  #                                                            fixedFormula = a1 + a2 + a3 + a6 + a4 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 2.4, a2 = -0.004, a3 = -0.017, a6 = -0.05, a4 = 0.3, b1 = 0.92, b1p = -0.19, b2 = 0.001, b2p = 0.012)))
  psmeDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a5ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2r) * (height - 1.37)), psme2022, # a1r singular precision, b1r also tractable
                                                     fixedFormula = a1 + a5ac + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 2, a5ac = -0.1, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkPhysioB1 = fit_nlme("Ruark physio", dbh ~ (a1 + a5ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1r) * exp((b2) * (height - 1.37)), psme2022, # near identical to b2r
  #                                                     fixedFormula = a1 + a5ac + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 2, a5ac = -0.1, b1 = 0.8, b2 = 0.01)))
  psmeDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight + b2r) * (height - 1.37)), psme2022, # a1r, b1r also tractable
                                                    fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtA1 = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psme2022, # false convergence
  #                                                    fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtB1 = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1r) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psme2022, # likely slightly less accurate than b1r
  #                                                    fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  psmeDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a5ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight + b2r) * (height - 1.37)), psme2022, # a1r, b1r also tractable
                                                          fixedFormula = a1 + a5ac + b1 + b2 + b2rh ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 2, a5ac = -0.09, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtPhysioA1 = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a5ac * cos(3.14159/180 * aspect) + a1r)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psme2022, # false convergence
  #                                                          fixedFormula = a1 + a5ac + b1 + b2 + b2rh ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 2, a5ac = -0.09, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtPhysioB1 = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a5ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1r) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psme2022, # less accurate than b2r
  #                                                          fixedFormula = a1 + a5ac + b1 + b2 + b2rh ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 2, a5ac = -0.09, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.37^b1), 0.9999)), psme2022, # fixed form failed to fit
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90)))
  psmeDiameterFromHeightMixed$schnute = create_fit_statistics(name = "Schnute inverse", fitting = "nlme")
  #psmeDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1 * (height - 1.37)^(b1 + b1r)*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^(b4), psme2022, # a1r, b2r, b3r singularity in backsolve, b1r, b4r step halving
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve with
  #                                                    start = list(fixed = c(a1 = 3, b1 = 0.7, b2 = 0.1, b3 = 0.03, b4 = 0.3)))
  psmeDiameterFromHeightMixed$sharmaParton = create_fit_statistics(name = "modified Sharma-Parton", fitting = "nlme")
  
  psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2r)), psme2022, # a1r step halving, b1r also tractable
                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)))
  psmeDiameterFromHeightMixed$sibbesenReplaceB1 = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1r)*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, # a1r step halving, b1r also tractable
                                                           fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                           start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)))
  
  #psmeDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
  #                                                           fixedFormula = a1 + a1p + a2 + a3 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.006, a3 = -0.03, b1 = 0.61, b2 = 0.071, b2p = 0.021)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect))*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
  #                                                                 fixedFormula = a1 + a1p + a2 + a3 + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psme2022, 
  #                                                                fixedFormula = a1 + a1p + a2 + a3 + a4 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.006, a3 = -0.03, a4 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysio = fit_nlme("Sibbesen replace ABA+T RelHt physio", dbh ~ (a1 + a1r + a1p * isPlantation + a2 * bootstrapStandBasalAreaPerHectare + a3 * basalAreaTaller + a6 * cos(3.14159/180 * aspect) + a4 * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), psme2022,
  #                                                                      fixedFormula = a1 + a1p + a2 + a3 + a6 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                      start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = 0, a3 = 0, a6 = -0.06, a4 = 0.5, b1 = 0.60, b2 = 0.09)))
  
  psmeDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2ac * cos(3.14159/180 * aspect) + b2r)), psme2022, # a1r false convergence, false convergence but more accurate than b1r
                                                               fixedFormula = a1 + a1p + b1 + b2 + b2p + b2ac ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 1.2, a1p = 0.3, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2ac = -0.003)))
  #psmeDiameterFromHeightMixed$sibbesenReplacePhysioB1 = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1r)*(height - 1.37)^(b2 + b2p * isPlantation + b2ac * cos(3.14159/180 * aspect))), psme2022,
  #                                                               fixedFormula = a1 + a1p + b1 + b2 + b2p + b2ac ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = 1.2, a1p = 0.3, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2ac = -0.003)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2r)), psme2022, # a1r step halving, b1r singularity in backsolve
                                                              fixedFormula = a1 + b1 + b2 + b2p + b2rh ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2ac * cos(3.14159/180 * aspect) + b2r)), psme2022, # a1r, b1r false convergence
                                                                    fixedFormula = a1 + b1 + b2 + b2p + b2rh + b2ac ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01, b2ac = -0.003)))
  #psmeDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2 + b2r), psme2022, # a1r step halving, b1r singularity in backsolve, b2r job max iterations
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = -140, b1 = 0.01, b2 = 0.9)))
  psmeDiameterFromHeightMixed$weibull = create_fit_statistics(name = "Weibull inverse", fitting = "nlme")
  #to_fixed_coeffficients(psmeDiameterFromHeightMixed$weibull)
  #print(to_parameter_confidence_intervals(psmeDiameterFromHeightMixed$weibull), n = 40)
  
  psmeDiameterFromHeightMixed$gam = fit_gam("REML GAM", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), random = list(uniquePlotID = ~1), data = psme2022)
  #psmeDiameterFromHeightMixed$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 28), random = list(uniquePlotID = ~1), data = psme2022)
  psmeDiameterFromHeightMixed$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9), random = list(uniquePlotID = ~1), data = psme2022)
  #psmeDiameterFromHeightMixed$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 331), random = list(uniquePlotID = ~1), data = psme2022)
  #psmeDiameterFromHeightMixed$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", dbh ~ I(height^) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496), random = list(uniquePlotID = ~1), data = psme2022)
  psmeDiameterFromHeightMixed$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 18), random = list(uniquePlotID = ~1), data = psme2022)
  psmeDiameterFromHeightMixed$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 17), random = list(uniquePlotID = ~1), data = psme2022)
  #lapply(psmeDiameterFromHeightMixed$gamRelHtPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(psmeDiameterFromHeightMixed$gamRelHtPhysio$fit, function(fit) { return(summary(fit$gam)) })
  #psmeDiameterFromHeightMixed$gam$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(psmeDiameterFromHeightMixed, paste0("trees/height-diameter/data/PSME DBH mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}

if (psmeOptions$fitHeightPrimary & psmeOptions$fitHeightMixed & psmeOptions$fitDbhPrimary & psmeOptions$fitDbhMixed)
{
  # assemble parameter and results tibbles using incremental load and reduce to fit in memory
  # height primary + nlrob + gsl_nls + dbh primary + nlrob = 120 GB DDR at 10x10 cross validation with models dropped
  # Therefore, batch formation height and diameter results.
  if (exists("psmeHeightFromDiameter") == FALSE) 
  { 
    psmeHeightFromDiameter = readRDS(paste0("trees/height-diameter/data/PSME height ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  if (exists("psmeHeightFromDiameterMixed") == FALSE) 
  { 
    psmeHeightFromDiameterMixed = readRDS(paste0("trees/height-diameter/data/PSME height mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  if (exists("psmeDiameterFromHeight") == FALSE)
  { 
    psmeDiameterFromHeight = readRDS(paste0("trees/height-diameter/data/PSME DBH ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  if (exists("psmeDiameterFromHeightMixed") == FALSE) 
  { 
    psmeDiameterFromHeightMixed = readRDS(paste0("trees/height-diameter/data/PSME DBH mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }

  psmeFixedEffects = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, unnest_fixed_effects)),
                                         bind_rows(lapply(psmeHeightFromDiameterMixed, unnest_fixed_effects))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(psmeDiameterFromHeight, unnest_fixed_effects)),
                                         bind_rows(lapply(psmeDiameterFromHeightMixed, unnest_fixed_effects))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")
  psmeFitStats = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, unnest_validation_statistics)),
                                    bind_rows(lapply(psmeHeightFromDiameterMixed, unnest_validation_statistics))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(psmeDiameterFromHeight, unnest_validation_statistics)),
                                    bind_rows(lapply(psmeDiameterFromHeightMixed, unnest_validation_statistics))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")

  check_plot_fit_stats(psmeFitStats)
  saveRDS(psmeFixedEffects, paste0("trees/height-diameter/data/PSME coefficients ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  saveRDS(psmeFitStats, paste0("trees/height-diameter/data/PSME fit stats ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## preferred forms identified (results.R, Figure 8 and trees.R)
if (psmeOptions$fitHeightPrimary & psmeOptions$fitDbhPrimary)
{
  psmeHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = psme2022, threads = 4, folds = 1, repetitions = 1))
  psmeHeightFromDiameterPreferred$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 73, a2 = 150, b1 = -1.3, b1p = 0.03), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"), psme2022, mtry = 2, minNodeSize = 48, sampleFraction = 0.230, folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", height ~ 1.37 + (a1 + a4 * relativeDiameter)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psme2022, start = list(a1 = 20, a4 = 0.4, b1 = 0.3, b2 = -0.02, b3 = -0.05, b4 = 1.0, b4e = 0.0001, b4as = -0.04, b4ac = -0.04), folds = 1, repetitions = 1)
  saveRDS(psmeHeightFromDiameterPreferred, "trees/height-diameter/data/PSME preferred height models.Rds")
  
  psmeDiameterFromHeightPreferred = list(michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 130, a2 = 70, b1 = 0.9, b1p = -0.008), folds = 1, repetitions = 1))
  psmeDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2ac * cos(pi/180*aspect))), psme2022, start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01, b2ac = -0.003), folds = 1, repetitions = 1)
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
  summary(psmeHeightGam) # all significant
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
