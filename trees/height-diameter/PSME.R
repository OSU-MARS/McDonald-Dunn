# load libraries, functions, and abgr2022 from setup.R
# psmeOptions from PSME job.R
psmeCrossValidation = get_cross_validation_folds(psme2022)


## Douglas-fir height regressions
if (psmeOptions$fitHeightPrimary)
{
  # linear regressions
  psmeHeightFromDiameter = list(linear = fit_lm("linear", height ~ dbh + I(isPlantation*dbh), psmeCrossValidation))
  psmeHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ dbh + I(dbh^2) + I(isPlantation*dbh^2), psmeCrossValidation) # plantation, plantation^2 not mutually significant
  
  # nonlinear regressions
  psmeHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1 * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psmeCrossValidation, 
                                                       start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2p = 0.1)) # a1p not significant, 10x10 MAE 3.45, AIC 3897
  #psmeHeightFromDiameter$chapmanRichardsB2a1 = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1)*dbh^(b2))), psmeCrossValidation, 
  #                                                         start = list(a1 = 58, a1p = -3, b1 = -0.01, b2 = 1.1)) # 10x10 MAE 3.44, AIC 3938
  #psmeHeightFromDiameter$chapmanRichardsB2b1 = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1) * (1 - exp((b1 + b1p * isPlantation)*dbh^(b2))), psmeCrossValidation, 
  #                                                         start = list(a1 = 55, b1 = -0.01, b1p = 0.0001, b2 = 1.1)) # 10x10 MAE 3.48, AIC 3851
  #psmeHeightFromDiameter$chapmanRichardsB2b2 = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh^(b2 + b2p * isPlantation))), psmeCrossValidation, 
  #                                                         start = list(a1 = 55, b1 = -0.01, b2 = 1.1, b2p = -0.03)) # 10x10 MAE 3.47, AIC 3882
  #psmeHeightFromDiameter$chapmanRichardsB3 = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1)*dbh^(b2)))^(b3), psmeCrossValidation, 
  #                                                       start = list(a1 = 58, a1p = -3, b1 = -0.001, b2 = 1.2, b3 = 0.8), distinct = FALSE) # b3 not significant, 10x10 MAE 3.47, AIC 3938
  # linear BA+L
  #psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp((b1)*dbh))^(b2 + b2bal * basalAreaLarger), psmeCrossValidation, 
  #                                                        start = list(a1 = 57, a1ba = 0.1, b1 = -0.02, b2 = 1.3, b2bal = -0.005)) # a1bap, a1balp, b1ba, b1bap, b1balp, b2p, b2ba, b2bap, b2balp not significant, a1bal, b1bal also significant, 10x10 MAE 3.19, AIC 3829
  #psmeHeightFromDiameter$chapmanRichardsBalA1 = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare + a1bal * basalAreaLarger) * (1 - exp((b1)*dbh))^(b2), psmeCrossValidation, 
  #                                                          start = list(a1 = 57, a1ba = 0.1, a1bal = 0, a1balp = 0, b1 = -0.02, b2 = 1.3)) # job singular gradient
  #psmeHeightFromDiameter$chapmanRichardsBalB1 = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp((b1 + b1bal * basalAreaLarger)*dbh))^(b2), psmeCrossValidation, 
  #                                                          start = list(a1 = 57, a1ba = 0.1, b1bal = 0, b1balp = 0, b1 = -0.02, b2 = 1.3)) # job singular gradient
  #psmeHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp((b1)*dbh))^(b2 + b2bal * basalAreaLarger), psmeCrossValidation, 
  #                                                              start = list(a1 = 70, a1rd = -1.7, b1 = -0.015, b2 = 1.3, b2bal = -0.005)) # a1ba, a1rdp, b1rdp, b2rdp not significant, a1rd, b1rd also significant, 10x10 MAE 3.13, AIC 3790
  #psmeHeightFromDiameter$chapmanRichardsBalRelDbhB2 = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp((b1)*dbh))^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger), psmeCrossValidation, 
  #                                                                start = list(a1 = 57, a1ba = 0.05, b1 = -0.02, b2 = 1.25, b2rd = 0.1, b2bal = -0.005))
  #psmeHeightFromDiameter$chapmanRichardsBalRelDbhB1 = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + a1 * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^(b2 + b2bal * basalAreaLarger), psmeCrossValidation, 
  #                                                                start = list(a1 = 65, b1 = -0.02, b1rd = 0.001, b2 = 1.3, b2bal = -0.005))
  #psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation + b2bal * basalAreaLarger), psmeCrossValidation, 
  #                                                              start = list(a1 = 55, a1ba = 0.1, b1 = -0.01, b1ac = -0.001, b2 = 1.3, b2p = 0.1, b2bal = -0.008)) # 10x10 MAE 3.06, AIC 3820
  #psmeHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger), psmeCrossValidation, 
  #                                                                    start = list(a1 = 61, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2rd = 0.1, b2bal = -0.005)) # a1, b1ba, b2ba not significant, 10x10 MAE 3.07, AIC 3789
  # reciprocal BAL
  psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, 
                                                          start = list(a1 = 62, b1 = -0.02, b2 = 1.0, b2bal = 4, b2balr = 7)) # a1ba+r, b1ba+r, b2ba+r, b2balp not significant, 10x10 MAE 3.07, AIC 3745
  #psmeHeightFromDiameter$chapmanRichardsBalA1 = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger)) * (1 - exp((b1)*dbh))^(b2), psmeCrossValidation, 
  #                                                          start = list(a1 = 90, a1bal = -100, a1balr = 30, b1 = -0.01, b2 = 1.1)) # a1ba not significant, 10x10 MAE 3.13, AIC 3751
  #psmeHeightFromDiameter$chapmanRichardsBalB1 = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1bal / (b1balr + basalAreaLarger))*dbh))^(b2), psmeCrossValidation, 
  #                                                          start = list(a1 = 62, b1 = -0.02, b1bal = 0.1, b1balr = 20, b2 = 1.1)) # a1ba not significant, 10x10 MAE 3.08, AIC 3766
  psmeHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp((b1)*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, 
                                                                start = list(a1 = 70, a1rd = -1.5, b1 = -0.02, b2 = 1.0, b2bal = 4, b2balr = 8)) # a1ba not significant, 10x10 MAE 3.04, AIC 3814
  psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, 
                                                                start = list(a1 = 65, b1 = -0.01, b1ac = -0.001, b2 = 0.9, b2bal = 5, b2balr = 8)) # a1ba not significant, 10x10 MAE 2.96, AIC 3749
  psmeHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, 
                                                                      start = list(a1 = 61, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2rd = 0.1, b2bal = 4, b2balr = 8), distinct = FALSE) # b2rd not significant, 10x10 ME 3.01, AIC 3694
  psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), psmeCrossValidation, 
                                                             start = list(a1 = 60, b1 = -0.02, b1ac = -0.001, b2 = 1.1, b2p = 0.1)) # { a1, b1, b2 } x { e, s, tr, tw }, { a1, b2 } x a not significant, b1as not reliably significant
  psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psmeCrossValidation, 
                                                             start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2rd = 0.1, b2rdp = 0.05)) # a1rdp, b2p not significant, b1rd NaN-inf, a1rd also significant, 10x10 MAE 3.28, AIC 3802
  #psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1 + b1p * isPlantation)*relativeDiameter*dbh))^(b2), psmeCrossValidation, 
  #                                                           start = list(a1 = 50, b1 = -0.008, b1p = 0.004, b2 = 0.6)) # a1p not significant, b1p-b2p not mutually significant, 10x10 MAE 4.80, AIC 4402
  #psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1)*relativeDiameter*dbh))^(b2 + b2p * isPlantation), psmeCrossValidation, 
  #                                                           start = list(a1 = 50, b1 = -0.008, b2 = 0.4, b2p = 0.2)) # 10x10 MAE 4.77, AIC 4362
  #psmeHeightFromDiameter$chapmanRichardsRelDbhA1 = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp((b1)*dbh))^(b2), psmeCrossValidation, 
  #                                                             start = list(a1 = 72, a1rd = -2.5, b1 = -0.02, b2 = 1.1))
  psmeHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psmeCrossValidation, 
                                                                   start = list(a1 = 58, b1 = -0.02, b1ac = -0.001, b2 = 1.1, b2rd = 0.1, b2rdp = 0.06))
  psmeHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation), psmeCrossValidation, 
                                              start = list(a1 = 3, a1p = -1.4, b1 = 0.4, b1p = -0.16)) # 10x10 MAE 3.95, AIC 4085
  psmeHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psmeCrossValidation, 
                                                start = list(a1 = 73, a2 = 150, b1 = -1.3, b1p = 0.03)) # a1p not significant, a2p-b1p not mutually significant, 10x10 MAE 3.47, AIC 3877
  #psmeHeightFromDiameter$hossfeldBal = fit_gsl_nls("Hossfeld IV BA+L", height ~ 1.37 + (a1) / (1 + (a2 + a2bal / (a2balr + basalAreaLarger)) * dbh^(b1)), psmeCrossValidation, 
  #                                                 start = list(a1 = 90, a2 = 120, a2bal = 1000, a2balr = 12, b1 = -1.2)) # a2balp, a2balrp not significant, 10x10 MAE 3.07, AIC 3725
  #psmeHeightFromDiameter$hossfeldBal = fit_gsl_nls("Hossfeld IV BA+L", height ~ 1.37 + (a1) / (1 + (a2 + a2bal / (a2balr + basalAreaLarger)) * dbh^(b1 + b1ba / standBasalAreaPerHectare)), psmeCrossValidation, 
  #                                                 start = list(a1 = 73, a2 = 140, a2bal = 800, a2balr = 10, b1 = -1.25, b1ba = 0.0), distinct = FALSE) # b1ba not significant, 10x10 MAE 3.08, AIC 3738
  #psmeHeightFromDiameter$hossfeldBaA1 = fit_gsl_nls("Hossfeld IV BA", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) / (1 + (a2) * dbh^(b1)), psmeCrossValidation, 
  #                                                  start = list(a1 = 60, a1ba = 0.3, a2 = 130, b1 = -1.2)) # 10x10 MAE 3.37, AIC 3833
  #psmeHeightFromDiameter$hossfeldBaA2 = fit_gsl_nls("Hossfeld IV BA", height ~ 1.37 + (a1) / (1 + (a2 + a2ba * standBasalAreaPerHectare) * dbh^(b1)), psmeCrossValidation, 
  #                                                  start = list(a1 = 70, a2 = 160, a2ba = -0.9, b1 = -1.2)) # 10x10 MAE 3.31, AIC 3828
  #psmeHeightFromDiameter$hossfeldBaB1 = fit_gsl_nls("Hossfeld IV BA", height ~ 1.37 + (a1) / (1 + (a2) * dbh^(b1 + b1ba / standBasalAreaPerHectare)), psmeCrossValidation, 
  #                                                  start = list(a1 = 75, a2 = 110, b1 = -1.2, b1ba = 3)) # 10x10 MAE 3.29, AIC 3818
  #psmeHeightFromDiameter$hossfeldBalA1 = fit_gsl_nls("Hossfeld IV BAL", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger)) / (1 + (a2) * dbh^(b1)), psmeCrossValidation, 
  #                                                   start = list(a1 = 125, a1bal = -1000, a1balr = 25, a2 = 170, b1 = -1.1)) # 10x10 MAE 3.14, AIC 3722
  #psmeHeightFromDiameter$hossfeldBalB1 = fit_gsl_nls("Hossfeld IV BAL", height ~ 1.37 + (a1) / (1 + (a2) * dbh^(b1 + b1bal / (b1balr + basalAreaLarger))), psmeCrossValidation, 
  #                                                   start = list(a1 = 90, a2 = 170, b1 = -1.3, b1bal = 3, b1balr = 20)) # 10x10 MAE 3.11, AIC 3788
  #psmeHeightFromDiameter$hossfeldBalRelDbh = fit_gsl_nls("Hossfeld IV RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter) / (1 + (a2 + a2bal / (a2balr + basalAreaLarger)) * dbh^b1), psmeCrossValidation, 
  #                                                       start = list(a1 = 73, a1rd = 20, a2 = 140, a2bal = 800, a2balr = 10, b1 = -1.25), distinct = FALSE) # a1rd not significant
  #psmeHeightFromDiameter$hossfeldRelDbh = fit_gsl_nls("Hossfeld IV RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter) / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                    start = list(a1 = 90, a1rd = -3, a2 = 170, b1 = -1.2, b1p = 0.02)) # 10x10 MAE 3.37, AIC 3946
  #psmeHeightFromDiameter$hossfeldRelDbhA2 = fit_gsl_nls("Hossfeld IV RelDbh", height ~ 1.37 + (a1) / (1 + (a2 + a2rd * relativeDiameter) * dbh^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                      start = list(a1 = 80, a2 = 170, a2rd = 30, b1 = -1.3, b1p = 0.02)) # 10x10 MAE 3.38, AIC 3882
  #psmeHeightFromDiameter$hossfeldRelDbhB1 = fit_gsl_nls("Hossfeld IV RelDbh", height ~ 1.37 + (a1) / (1 + (a2) * dbh^(b1 + b1rd * relativeDiameter)), psmeCrossValidation, 
  #                                                      start = list(a1 = 90, a2 = 200, b1 = -1.3, b1rd = 0.02)) # 10x10 MAE, 3.40, AIC 3890
  psmeHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psmeCrossValidation, 
                                            start = list(a1 = 200, b1 = -6.8, b1p = -1.2, b2 = -0.3, b2p = -0.04)) # a1p-b1p not mutually significant, 10x10 MAE 3.53, AIC 3850
  psmeHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psmeCrossValidation, 
                                                       start = list(a1 = 75, a2 = 140, a2p = 20, b1 = 1.2)) # a1p, b1p not significant, 10x10 MAE 3.44, AIC 3841
  psmeHeightFromDiameter$michaelisMentenBal = fit_gsl_nls("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1ba / standBasalAreaPerHectare)), psmeCrossValidation, 
                                                          start = list(a1 = 85, a2 = 115, a2bal = 900, a2balr = 10, b1 = 1.2, b1ba = -1)) # a2balp not significant, 10x10 MAE 3.06, AIC 3715
  psmeHeightFromDiameter$michaelisMentenBalRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + a2rd * relativeDiameter + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1ba / standBasalAreaPerHectare)), psmeCrossValidation, 
                                                                start = list(a1 = 85, a2 = 115, a2rd = 25, a2bal = 900, a2balr = 10, b1 = 1.2, b1ba = -1), distinct = FALSE) # a1rd, a2rd, b1rd not significant
  #psmeHeightFromDiameter$michaelisMentenBaB1 = fit_gsl_nls("Michaelis-Menten BA", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1ba / standBasalAreaPerHectare)), psmeCrossValidation, 
  #                                                         start = list(a1 = 75, a2 = 110, b1 = 1.2, b1ba = -3)) # a1ba+r, b1bar not significant, a2ba+r numerically unstable, 10x10 MAE 3.28, AIC 3833
  #psmeHeightFromDiameter$michaelisMentenBalA2 = fit_gsl_nls("Michaelis-Menten BAL", height ~ 1.37 + (a1)*dbh^(b1) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1)), psmeCrossValidation, 
  #                                                          start = list(a1 = 90, a2 = 125, a2bal = 1000, a2balr = 12, b1 = 1.2)) # a1bal+r not significant, 10x10 MAE 3.08, AIC 3718
  #psmeHeightFromDiameter$michaelisMentenBalB1 = fit_gsl_nls("Michaelis-Menten BAL", height ~ 1.37 + (a1)*dbh^(b1 + b1bal / (b1balr + basalAreaLarger)) / (a2 + dbh^(b1 + b1bal / (b1balr + standBasalAreaPerHectare))), psmeCrossValidation, 
  #                                                          start = list(a1 = 120, a2 = 190, b1 = 1.2, b1bal = -3, b1balr = 25)) # 10x10 MAE 3.14, AIC 3743
  #psmeHeightFromDiameter$michaelisMentenBaBalB1 = fit_gsl_nls("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare + b1bal / (b1balr + basalAreaLarger)) / (a2 + dbh^(b1 + b1ba / standBasalAreaPerHectare + b1bal / (b1balr + standBasalAreaPerHectare))), psmeCrossValidation, 
  #                                                            start = list(a1 = 110, a2 = 150, b1 = 1.2, b1ba = -1.5, b1bal = -3, b1balr = 30)) # 10x10 MAE 3.11, AIC 3850
  psmeHeightFromDiameter$michaelisMentenRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^b1 / (a2 + a2rd * relativeDiameter + a2p * isPlantation + dbh^b1), psmeCrossValidation, 
                                                             start = list(a1 = 80, a2 = 150, a2rd = 25, a2p = 15, b1 = 1.2)) # a2rdp not significant, 10x10 MAE 3.35, AIC 3898
  #psmeHeightFromDiameter$michaelisMentenRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^b1 / (a2 + (relativeDiameter * dbh)^b1), psmeCrossValidation, 
  #                                                           start = list(a1 = 420, a2 = 480, b1 = 0.9)) # 10x10 MAE 3.57, AIC 3944
  #psmeHeightFromDiameter$michaelisMentenRelDbhA1 = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psmeCrossValidation, 
  #                                                             start = list(a1 = 90, a1rd = -3, a2 = 140, a2p = 15, b1 = 1.2)) # 10x10 MAE 3.36, AIC 3851
  #psmeHeightFromDiameter$michaelisMentenRelDbhB1 = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + dbh^(b1 + b1rd * relativeDiameter)), psmeCrossValidation, 
  #                                                             start = list(a1 = 90, a2 = 170, b1 = 1.2, b1rd = -0.02)) # 10x10 MAE 3.41, AIC 3918
  psmeHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation), psmeCrossValidation, 
                                             start = list(a1 = 2.5, a1p = -1, b1 = 0.65, b1p = 0.14)) # 10x10 MAE 4.01, AIC 4047
  psmeHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3 + a3p* isPlantation), psmeCrossValidation, 
                                              start = list(a1 = 0.012, a2 = 0.9, a3 = 1.4, a3p = 3.5)) # a1p, a2p not significant, 10x10 MAE 3.49, AIC 3989
  psmeHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(dbh + b2)), psmeCrossValidation, 
                                                 start = list(a1 = 73, b1 = -46, b1p = -3, b2 = 10)) # a1p, b2p not significant, 10x10 MAE 3.44, AIC 3943
  psmeHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-kU * dbh)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psmeCrossValidation, 
                                                 start = list(Ha = 55, Hap = -3, d = 0.3, dp = 0.2, kU = 0.013)) # kUp not significant, 10x10 MAE 3.45, AIC 4166, 2x50, 2x250, 2x500 NaN-inf
  #psmeHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psmeCrossValidation, 
  #                                                  start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0)) # a1p, b1p, b2p, b3, b3p, b4p not significant, 10x10 MAE 3.13, AIC 3677
  psmeHeightFromDiameter$sharmaPartonB4 = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh^(b4))), psmeCrossValidation, 
                                                      start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0)) # a1p, b1p, b2p, b3, b3p, b4, b4p not significant, 10x10 MAE 3.04, AIC 3695
  psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh))^b4, psmeCrossValidation, 
                                                       start = list(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3bal = -2, b3balr = 14, b4 = 1.0)) # 10x10 MAE 2.91, AIC 3705
  #psmeHeightFromDiameter$sharmaPartonBaBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + b2balr * basalAreaLarger))^(b3) * dbh))^b4, psmeCrossValidation, 
  #                                                       start = list(a1 = 30, b1 = 0.2, b2 = -0.01, b2balr = 35, b3 = -0.05, b4 = 1.0)) # 10x10 MAE 2.95, AIC 3766 but potentially fitting issues with NaN-inf due to b2balr runaway and b2balr not significant, Levenberg scale more stable but lower accuracy, moving b2 into or splitting it within b3 scope leads to loss of significance and numerical issues, large b2balr acts to remove influence of stand basal area but modulation by stand density provides an accuracy edge as expected based on GAM smooths
  #psmeHeightFromDiameter$sharmaParton0.1Ba10Bal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(0.1 * standBasalAreaPerHectare + 10 * basalAreaLarger))^b3 * dbh))^b4, psmeCrossValidation, 
  #                                                            start = list(a1 = 30, b1 = 0.2, b2 = -0.01, b3 = -0.05, b4 = 1.0)) # 10x10 MAE 2.95, AIC 3704
  #psmeHeightFromDiameter$sharmaParton1Bal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(1 + basalAreaLarger))^b3 * dbh))^b4, psmeCrossValidation, 
  #                                                      start = list(a1 = 30, b1 = 0.2, b2 = -0.02, b3 = -0.05, b4 = 1.0)) # 10x10 MAE 2.98, AIC 3766
  #psmeHeightFromDiameter$sharmaPartonBalA1 = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger))*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psmeCrossValidation, 
  #                                                       start = list(a1 = 35, a1bal = -100, a1balr = 15, b1 = 0.2, b2 = -0.02, b3 = 0, b4 = 1.0)) # b3 not significant, 10x10 MAE 2.93, AIC 3734
  #psmeHeightFromDiameter$sharmaPartonBalB1 = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1bal / (b1balr + basalAreaLarger)) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psmeCrossValidation, 
  #                                                       start = list(a1 = 24, b1 = 0.3, b1bal = -1, b1balr = 12, b2 = -0.02, b3 = -0.02, b4 = 1.0)) # 10x10 MAE 2.96, AIC 3747
  #psmeHeightFromDiameter$sharmaPartonBalB2 = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2 + b2bal / (b2balr + basalAreaLarger))*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psmeCrossValidation, 
  #                                                       start = list(a1 = 24, b1 = 0.3, b2 = -0.02, b2bal = 0.1, b2balr = 12, b3 = -0.02, b4 = 1.0)) # 10x10 MAE 2.93, AIC 3688
  #psmeHeightFromDiameter$sharmaPartonBalB3 = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal * basalAreaLarger) * dbh))^b4, psmeCrossValidation, 
  #                                                       start = list(a1 = 30, b1 = 0.2, b2 = -0.02, b3 = 0, b3bal = 0.002, b4 = 1.1)) # b3 not significant, 10x10 MAE 2.97, AIC 3713  
  #psmeHeightFromDiameter$sharmaPartonBalB4 = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), psmeCrossValidation, 
  #                                                       start = list(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4bal = 2, b4balr = 5)) # b3 not significant, 10x10 MAE 2.96, AIC 3704
  #psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation), psmeCrossValidation, 
  #                                                           start = list(a1 = 25, a1ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0, b4e = 0.0002)) # { a1, b1 } x { e, s, tr, tw }, a1 x { as, ac }, b1as, { b2, b3, b4 } x { s, tr, tw }, b3 not significant, { b2e, b3e }-b4e, { a1a, b1a, b2a, b3a }-b4a not mutually significant
  psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4 + b4e * elevation), psmeCrossValidation, 
                                                             start = list(a1 = 30, a1ac = 1.0, b1 = 0.2, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 14, b4 = 1.0, b4e = 0.0001)) # 10x10 MAE 2.85, AIC 3721
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB1 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4 + b4e * elevation), psmeCrossValidation, 
  #                                                             start = list(a1 = 25, b1 = 0.2, b1as = 0.004, b1ac = 0.01, b2 = -0.02, b3 = -0.04, b3bal = -2, b3balr = 14, b4 = 1.0, b4e = 0.0002))
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB2 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2 + b2e * elevation + b2as * sin(3.14159/180 * aspect) + b2ac * cos(3.14159/180 * aspect))*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4), psmeCrossValidation, 
  #                                                             start = list(a1 = 25, b1 = 0.3, b2 = -0.02, b2e = 0.000001, b2as = -0.001, b2ac = -0.001, b3 = -0.0, b3bal = -2, b3balr = 145, b4 = 1.0))
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB3 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3e * elevation + b3as * sin(3.14159/180 * aspect) + b3ac * cos(3.14159/180 * aspect))*dbh))^(b4), psmeCrossValidation, 
  #                                                             start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3bal = -2, b3balr = 14, b3e = -0.0001, b3as = 0.01, b3ac = 0.01, b4 = 1.0))
  #psmeHeightFromDiameter$sharmaPartonBalPhysioB4 = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4 + b4e * elevation + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, 
  #                                                             start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3bal = -2, b3balr = 14, b4 = 1.0, b4e = 0.0001, b4as = -0.035, b4ac = -0.04))
  psmeHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4), psmeCrossValidation, 
                                                                   start = list(a1 = 30, a1rd = 0, b1 = 0.2, b1as = -0.04, b1ac = 0.01, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 14, b4 = 1.0), distinct = FALSE) # a1rd not significant
  psmeHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh))^(b4), psmeCrossValidation, 
                                                             start = list(a1 = 30, a1rd = 0, b1 = 0.2, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 15, b4 = 1.0), distinct = FALSE) # a1rd, a1rdp, b1rd, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rd, b4rdp not significant, a1rd closest to significance
  # a1ac, b1ac, b2ac, b4ac forms near interchangeable, b1ac slightly preferred
  psmeHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), psmeCrossValidation, 
                                                          start = list(a1 = 21, b1 = 0.3, b1as = 0.006, b1ac = 0.01, b2 = -0.02, b3 = 0, b4 = 1.0)) # { a1, b1 } x { e, s, tr, tw }, a1ac, b1as, { b2, b3, b4 } x { s, tr, tw }, b3 not significant, { b2e, b3e }-b4e, { a1a, b1a, b2a, b3a }-b4a not mutually significant, 10x10 MAE 3.06, AIC 3727
  #psmeHeightFromDiameter$sharmaPartonPhysioA1 = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4e * elevation), psmeCrossValidation, 
  #                                                          start = list(a1 = 21, a1ac = 0.6, b1 = 0.3, b2 = -0.02, b3 = -0.01, b4 = 1.0, b4e = 0.0002)) # a1as not significant
  #psmeHeightFromDiameter$sharmaPartonPhysioB2 = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2 + b2as * sin(3.14159/180 * aspect) + b2ac * cos(3.14159/180 * aspect))*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), psmeCrossValidation, 
  #                                                          start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b2as = -0.001, b2ac = -0.001, b3 = 0, b4 = 1.0)) # b2e not significant
  #psmeHeightFromDiameter$sharmaPartonPhysioB3 = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3as * sin(3.14159/180 * aspect) + b3ac * cos(3.14159/180 * aspect))*dbh))^(b4), psmeCrossValidation, 
  #                                                          start = list(a1 = 23, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3e = -0.0001, b3as = 0.01, b3ac = 0.01, b4 = 1.0)) # b3e not significant, job singular gradient
  #psmeHeightFromDiameter$sharmaPartonPhysioB4 = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, 
  #                                                          start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04))
  psmeHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^b4, psmeCrossValidation, 
                                                          start = list(a1 = 13, b1 = 0.3, b1rd = 0, b2 = -0.02, b3 = 0.02, b4 = 1.0), distinct = FALSE) # a1rd, a1rdp, b1rd, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rd, b4rdp not significant
  psmeHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1rd * relativeDiameter) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, 
                                                                start = list(a1 = 18, b1 = 0.3, b1rd = 0.004, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04), distinct = FALSE) # b1rd, b3 not significant
  #psmeHeightFromDiameter$sharmaPartonB4bal = fit_gsl_nls("Sharma-Parton b₄ BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh^b4)), psmeCrossValidation, 
  #                                                       start = list(a1 = 30, b1 = 0.2, b2 = -0.01, b3 = 0.1, b3bal = -3, b3balr = 15, b4 = 1.1)) # b3, b4 significant, 10x10 MAE 2.93, AIC 3696
  #psmeHeightFromDiameter$sharmaPartonB4balPhysio = fit_gsl_nls("Sharma-Parton b₄ BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh^(b4 + b4as * sin(pi/180*aspect)))), psmeCrossValidation, 
  #                                                             start = list(a1 = 33, a1ac = 1.3, b1 = 0.2, b2 = -0.01, b3 = 0.1, b3bal = -3, b3balr = 20, b4 = 1.1, b4as = 0.01)) # b4 x { e, s, ac, tw } not significant, 10x10 MAE 2.86, AIC 3652
  #psmeHeightFromDiameter$sharmaPartonB4balPhysioTR = fit_gsl_nls("Sharma-Parton b₄ BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh^(b4 + b4tr * terrainRoughness))), psmeCrossValidation, 
  #                                                               start = list(a1 = 30, a1ac = 1.1, b1 = 0.2, b2 = -0.01, b3 = 0.1, b3bal = -3, b3balr = 18, b4 = 1.1, b4tr = 0.05)) # 10x10 MAE 2.87, AIC 3704
  #psmeHeightFromDiameter$sharmaPartonB4balPhysioRelDbh = fit_gsl_nls("Sharma-Parton b₄ RelDbh BA+L physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh^(b4 + b4as * sin(pi/180*aspect)))), psmeCrossValidation, 
  #                                                                   start = list(a1 = 30, a1rd = 0.2, a1ac = 0, b1 = 0.2, b2 = -0.01, b3 = 0.1, b3bal = -3, b3balr = 16, b4 = 1.1, b4as = 0.008), distinct = FALSE) # a1rd not significant
  #psmeHeightFromDiameter$sharmaPartonB4balRelDbh = fit_gsl_nls("Sharma-Parton b₄ RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh^(b4))), psmeCrossValidation, 
  #                                                             start = list(a1 = 30, a1rd = 0, b1 = 0.2, b2 = -0.01, b3 = 0.1, b3bal = -2, b3balr = 15, b4 = 1.0), distinct = FALSE) # a1rd, a1rdp, b1rd, b1rdp, b2rd, b2rdp, b3rd, b3rdp, b4rd, b4rdp not significant, a1rd closest to significance
  #psmeHeightFromDiameter$sharmaPartonB4physio = fit_gsl_nls("Sharma-Parton b₄ physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^(b4))), psmeCrossValidation, 
  #                                                          start = list(a1 = 21, b1 = 0.3, b1as = 0.005, b1ac = 0.01, b2 = -0.02, b3 = 0, b4 = 1.0)) # b3, b4 not significant, 10x10 MAE 3.14, AIC 3708
  #psmeHeightFromDiameter$sharmaPartonB4relDbh = fit_gsl_nls("Sharma-Parton b₄ RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh^b4)), psmeCrossValidation, 
  #                                                          start = list(a1 = 19, b1 = 0.3, b1rd = 0, b2 = -0.02, b3 = 0, b4 = 1.0), distinct = FALSE) # a1rd, b3, b4 not significant
  #psmeHeightFromDiameter$sharmaPartonB4relDbhPhysio = fit_gsl_nls("Sharma-Parton b₄ RelDbh physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1rd * relativeDiameter) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, 
  #                                                                start = list(a1 = 20, b1 = 0.3, b1rd = 0.004, b2 = -0.02, b3 = -0.02, b4 = 1.0, b4ac = -0.05, b4as = -0.04), distinct = FALSE) # b1rd, b3, b4 not significant
  #psmeHeightFromDiameter$sharmaPartonB5 = fit_gsl_nls("Sharma-Parton b₅", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh^(b4)))^(b5), psmeCrossValidation, 
  #                                                    start = list(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b5 = 1.0), distinct = FALSE) # b2, b3, b5 not significant
  #psmeHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, psmeCrossValidation, 
  #                                                 start = list(a1 = 45, a1p = -12, b1 = 0.08, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)) # b2p, b3p, b4p not significant, 10x10 MAE 3.29, AIC 3840
  psmeHeightFromDiameter$sharmaZhangB4 = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh^(b4))), psmeCrossValidation, 
                                                     start = list(a1 = 30, b1 = 0.16, b2 = -0.03, b3 = -0.09, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant, 10x10 MAE 3.24, AIC 3807
  psmeHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3 * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), psmeCrossValidation, 
                                                      start = list(a1 = 55, b1 = 0.1, b2 = -0.02, b3 = -0.03, b4 = 1.0, b4bal = 4, b4balr = 6)) # a1ba+r, b1, b1ba+r, b3, b3ba+r, b4ba+r not significant, b2ba+r NaN-inf, 10x10 MAE 3.05, AIC 3749
  psmeHeightFromDiameter$sharmaZhangBalRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3 * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger) + b4rd * relativeDiameter), psmeCrossValidation, 
                                                            start = list(a1 = 55, b1 = 0.1, b2 = -0.02, b3 = -0.03, b4 = 1.0, b4bal = 4, b4balr = 6, b4rd = 0), distinct = FALSE) # b4rd not significant
  #psmeHeightFromDiameter$sharmaZhangBalLinear = fit_gsl_nls("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a1ba * standBasalAreaPerHectare + a1bal * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp((b2)*standTreesPerHectare^b3 * dbh))^(b4), psmeCrossValidation, 
  #                                                          start = list(a1 = 50, a1p = -10, a1ba = -0.1, a1bal = 0.2, b1 = 0.1, b1p = 0.1, b2 = -0.03, b3 = -0.02, b4 = 1.1)) # a1bap, a1balp, { b1, b3 } x { bal, balp }, b1ba, b3, b4bap, b4balp not significant, b4ba, b4bal also significant, 10x10 MAE 3.05, AIC 3789
  #psmeHeightFromDiameter$sharmaZhangBalA1 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger))*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^b3 * dbh))^(b4), psmeCrossValidation, 
  #                                                      start = list(a1 = 50, a1bal = -300, a1balr = 25, b1 = 0.1, b2 = -0.03, b3 = -0.02, b4 = 1.1)) # a1p, a1ba, b1, b3 not significant, 10x10 MAE 3.05, AIC 3754 before elimination -> MAE 3.11, AIC 3768 after
  #psmeHeightFromDiameter$sharmaZhangBalB1 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a1ba * standBasalAreaPerHectare)*standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1bal / (b1balr + basalAreaLarger)) * (1 - exp((b2)*standTreesPerHectare^b3 * dbh))^(b4), psmeCrossValidation, 
  #                                                      start = list(a1 = 50, a1p = -10, a1ba = -0.1, b1 = 0.1, b1p = 0.1, b1bal = -2, b1balr = 20, b2 = -0.03, b3 = -0.02, b4 = 1.1)) # 10x10 MAE 3.10, AIC 3714
  #psmeHeightFromDiameter$sharmaZhangBalB2 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare)*standBasalAreaPerHectare^(b1) * (1 - exp((b2 + b2bal / (b2balr + basalAreaLarger))*standTreesPerHectare^b3 * dbh))^(b4), psmeCrossValidation, 
  #                                                      start = list(a1 = 50, a1ba = -0.1, b1 = 0.1, b2 = -0.02, b2bal = 0.2, b2balr = 19, b3 = -0.03, b4 = 1.1)) # 10x10 MAE 3.06, AIC 3766
  #psmeHeightFromDiameter$sharmaZhangBalB3 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh))^(b4), psmeCrossValidation, 
  #                                                      start = list(a1 = 50, b1 = 0.02, b2 = -0.02, b3 = 0, b3bal = -1, b3balr = 15, b4 = 1.1)) # 10x10 MAE 3.06, AIC 3707
  #psmeHeightFromDiameter$sharmaZhangBalB4 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp((b2)*standTreesPerHectare^b3 * dbh))^(b4 + b4ba * standBasalAreaPerHectare + b4bal * basalAreaLarger), psmeCrossValidation, 
  #                                                      start = list(a1 = 30, a1p = -10, b1 = 0.1, b1p = 0.1, b2 = -0.03, b3 = -0.02, b4 = 1.1, b4ba = 0.01, b4bal = -0.01))
  psmeHeightFromDiameter$sharmaZhangRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^(b4 + b4rd * relativeDiameter), psmeCrossValidation, 
                                                         start = list(a1 = 45, a1p = -12, b1 = 0.08, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0, b4rd = 0), distinct = FALSE) # a1rd, b1rd, b2rd, b3rd, b4rd not significant
  #psmeHeightFromDiameter$sharmaZhangB4bal = fit_gsl_nls("Sharma-Zhang b₄ BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3 * dbh^(b4d + b4bal / (b4balr + basalAreaLarger)))), psmeCrossValidation, 
  #                                                      start = list(a1 = 55, b1 = 0.05, b2 = -0.02, b3 = -0.04, b4d = 1.1, b4bal = -2, b4balr = 16)) # 10x10 MAE 3.06, AIC 3757
  #psmeHeightFromDiameter$sharmaZhangB4balRelDbh = fit_gsl_nls("Sharma-Zhang b₄ RelDbh BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3 * dbh^(b4d + b4bal / (b4balr + basalAreaLarger) + b4rd * relativeDiameter))), psmeCrossValidation, 
  #                                                            start = list(a1 = 55, b1 = 0.05, b2 = -0.02, b3 = -0.04, b4d = 1.1, b4bal = -2, b4balr = 16, b4rd = 0), distinct = FALSE) # b4rd not significant
  #psmeHeightFromDiameter$sharmaZhangB4relDbh = fit_gsl_nls("Sharma-Zhang b₄ RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^(b4d + b4rd * relativeDiameter), psmeCrossValidation, 
  #                                                         start = list(a1 = 42, b1 = 0.08, b2 = -0.02, b3 = 0.06, b4d = 1.1, b4rd = 0), distinct = FALSE) # a1p, b1p, b3, b4rd not significant
  #psmeHeightFromDiameter$sharmaZhangB4relDbh = fit_gsl_nls("Sharma-Zhang b₄ RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp(b2*standTreesPerHectare^b3*(relativeDiameter*dbh)^(b4d))), psmeCrossValidation, 
  #                                                         start = list(a1 = 16, b1 = 0.3, b2 = -0.8, b3 = -0.4, b4d = 0.5), distinct = FALSE) # 10x10 MAE 3.34, AIC 3873 (3.43, 3876 without relative diameter within b4d, accuracy quite poor if dividing by relativeDiameter)
  #psmeHeightFromDiameter$sharmaZhangB4relDbh = fit_gsl_nls("Sharma-Zhang b₄ RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp(b2*standTreesPerHectare^b3*dbh^(b4d)*relativeDiameter^b5)), psmeCrossValidation, 
  #                                                         start = list(a1 = 50, b1 = 0.06, b2 = 0.003, b3 = 0.1, b4d = 1.3, b5 = -0.3), distinct = FALSE) # 10x10 MAE 3.95, AIC 4079
  #psmeHeightFromDiameter$sharmaZhangB5 = fit_gsl_nls("Sharma-Zhang b₅", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh^(b4d)))^(b5), psmeCrossValidation, 
  #                                                   start = list(a1 = 30, b1 = 0.16, b2 = -0.03, b3 = -0.09, b4d = 1.0, b5 = 1.0)) # b2, b4d, b5 not significant
  psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psmeCrossValidation, 
                                                start = list(a1 = 0.2, b1 = 2.5, b1p = -0.2, b2 = -0.14, b2p = 0.02)) # a1p not significant, 10x10 MAE 3.50, AIC 4012
  #psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1p * isPlantation)*(relativeDiameter*dbh)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                              start = list(a1 = 0.5, b1 = 1.3, b1p = -0.1, b2 = -0.05, b2p = 0.01)) # 10x10 MAE 3.54, AIC 3948
  #psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*(relativeDiameter*dbh)^((b1 + b1p * isPlantation)*(dbh)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                              start = list(a1 = 4.7, b1 = 0.7, b1p = -0.4, b2 = -0.1, b2p = 0.1)) # 10x10 MAE 5.68, AIC 4584
  psmeHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psmeCrossValidation, 
                                               start = list(a1 = 56, b1 = -0.01, b1p = 0.001, b2 = 1.1)) # a1p, b2p not significant, 10x10 MAE 3.44, AIC 3926
  # linear BA+L: a1bal and b2bal forms MAE and AIC interchangeable
  #psmeHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger))), psmeCrossValidation, 
  #                                                start = list(a1 = 60, b1 = -0.01, b2 = 1.0, b2bal = 0.002)) # a1balp, b1ba, b1bap, b1bal, b2ba, b2balp, b2bal, b2balp not significant, a1bal also significant, 10x10 MAE 3.13, AIC 3704
  # reciprocal BA+L: only slight increase in accuracy from linear
  #psmeHeightFromDiameter$weibullBaA1 = fit_gsl_nls("Weibull", height ~ 1.37 + (a1 + a1ba / standBasalAreaPerHectare)*(1 - exp((b1)*dbh^(b2))), psmeCrossValidation, 
  #                                                 start = list(a1 = 65, a1ba = -300, b1 = -0.01, b2 = 1.0)) # a1ba+r, b1ba+r, b2ba+r not significant, { a1ba, b1ba, b2ba }-b1p not mutually significant, 10x10 MAE 3.35, AIC 3821
  #psmeHeightFromDiameter$weibullBaB1 = fit_gsl_nls("Weibull", height ~ 1.37 + (a1)*(1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2))), psmeCrossValidation, 
  #                                                 start = list(a1 = 60, b1 = -0.02, b1ba = 0.1, b2 = 1.0)) # 10x10 MAE 3.31, AIC 3882
  #psmeHeightFromDiameter$weibullBaB2 = fit_gsl_nls("Weibull", height ~ 1.37 + (a1)*(1 - exp((b1)*dbh^(b2 + b2ba / standBasalAreaPerHectare))), psmeCrossValidation, 
  #                                                 start = list(a1 = 58, b1 = -0.01, b2 = 1.1, b2ba = -2)) # 10x10 MAE 3.32 AIC 3800
  #psmeHeightFromDiameter$weibullBalA1 = fit_gsl_nls("Weibull", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger))*(1 - exp((b1)*dbh^(b2))), psmeCrossValidation, 
  #                                                  start = list(a1 = 90, a1bal = -1000, a1balr = 30, b1 = -0.01, b2 = 1.1)) # b1bal+r not significant, { a1ba, b1ba, b2ba }-b1p not mutually significant, 10x10 MAE 3.14, AIC 3778
  #psmeHeightFromDiameter$weibullBalB2 = fit_gsl_nls("Weibull", height ~ 1.37 + (a1)*(1 - exp((b1)*dbh^(b2 + b2bal / (b2balr + basalAreaLarger)))), psmeCrossValidation, 
  #                                                  start = list(a1 = 65, b1 = -0.01, b2 = 1.1, b2bal = -3, b2balr = 20)) # 10x10 MAE 3.32, AIC 3845
  psmeHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2 + b2bal * basalAreaLarger))), psmeCrossValidation, 
                                                  start = list(a1 = 64, b1 = -0.01, b1ba = 0.03, b2 = 1.0, b2bal = 0.002)) # 10x10 MAE 3.12, AIC 3799
  #psmeHeightFromDiameter$weibullBalA1 = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp(b1*dbh^(b2))), psmeCrossValidation, 
  #                                                  start = list(a1 = 60, a1bal = 0, b1 = -0.01, b2 = 1.0))
  psmeHeightFromDiameter$weibullBalRelDbh = fit_gsl_nls("Weibull RelDbh BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger))), psmeCrossValidation, 
                                                        start = list(a1 = 65, b1 = -0.01, b1ba = 0.03, b2 = 1.0, b2rd = -0.005, b2bal = 0.002)) # 10x10 MAE 3.11, AIC 3791
  #psmeHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psmeCrossValidation, 
  #                                                   start = list(a1 = 67, a1rd = -2, b1 = -0.01, b1p = 0.001, b2 = 1.1)) # a1rdp not significant, 10x10 MAE 3.39, AIC 3863
  #psmeHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1)*(1 - exp((b1 + b1p * isPlantation + b1rd * relativeDiameter)*dbh^b2)), psmeCrossValidation, 
  #                                                   start = list(a1 = 58, b1 = -0.01, b1p = 0.001, b1rd = 0.0005, b2 = 1.1)) # b1rdp not significant, 10x10 MAE 3.39, AIC 3956
  psmeHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1)*(1 - exp((b1 + b1p * isPlantation)*dbh^(b2 + b2rd * relativeDiameter))), psmeCrossValidation, 
                                                     start = list(a1 = 60, b1 = -0.01, b1p = 0.001, b2 = 1.1, b2rd = -0.02)) # b2rdp not significant, 10x10 MAE 3.38, AIC 3852
  #print(to_parameter_confidence_intervals(psmeHeightFromDiameter$sharmaZhangBalRelDbh), n = 40)
  #to_fixed_coeffficients(psmeHeightFromDiameter$sharmaZhangBalRelDbh)
  #psmeHeightFromDiameter$weibullBal$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  # all GAM forms NaN with Γ link
  psmeHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, bs = "ts", k = 9), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3918, MAE ~3.42 with power base, plantation not significant, specification of 1.37 m offset not supported (invalid model formula in ExtractVars, mgcv 1.9-3)
  #psmeHeightFromDiameter$gamHossfeld = fit_gam("REML GAM", height ~ I(71.61941/(1 + 158.72979 * dbh^(-1.26072 + 0.02926 * isPlantation))) + s(dbh, bs = "ts", k = 7), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC ~3924, MAE ~3.40 with Hossfeld IV
  #psmeHeightFromDiameter$gamMichaelis = fit_gam("REML GAM", height ~ I(72.299 * dbh^1.231/(143.174 + 18.141 * isPlantation + dbh^1.231)) + s(dbh, bs = "ts", k = 7), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3902, MAE ~3.44 with Michaelis-Menten
  #psmeHeightFromDiameter$gamS0 = fit_gam("REML GAM", height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC ~3897, MAE ~3.46
  #psmeHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 15), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3881
  psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3825
  #psmeHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3836
  psmeHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 14), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3818, plantation not significant
  psmeHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 10), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3898
  
  if (psmeOptions$fitPhysioGams)
  {
    psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 15), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3789
    psmeHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 19), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3848
    # modest reductions in MAE from slope, aspect, roughness, and wetness, including elevation is detrimental
    #psmeHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, elevation, bs = "ts", by = as.factor(isPlantation), k = 10), speciesData = psmeCrossValidation, threads = 8, distinct = FALSE) # 10x10 AIC 4125, support at k ~= 33 but clearly overfit, k = 10 minimum for curvature
    #psmeHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 15), speciesData = psmeCrossValidation, threads = 8) # 10x10 AIC 3856
    #psmeHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 13), speciesData = psmeCrossValidation, threads = 8) # 10x10 AIC 3937
    #psmeHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 14), speciesData = psmeCrossValidation, threads = 8) # 10x10 AIC3948
    #psmeHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 13), speciesData = psmeCrossValidation, threads = 8) # 10x10 AIC 3940
    #psmeHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), speciesData = psmeCrossValidation, threads = 8) # 10x10 AIC 3941
    #psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 85), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 4220 + higher MAE than slope or wetness alone @ k ~= 61 -> arguably overfit, fails on dplyr errors with threads = 8, warnings at threads = 4
    #psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, slope, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 14), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3945 + higher MAE than slope or wetness alone
    psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 10), speciesData = psmeCrossValidation, threads = 8) # 10x10 AIC 3941
    psmeHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, topographicWetnessFD8f, bs = "ts", k = 12, by = as.factor(isPlantation)), speciesData = psmeCrossValidation, threads = 4) # 10x10 AIC 3925
    #lapply(psmeHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
    #lapply(psmeHeightFromDiameter$gamRelDbhPhysio$fit, summary)
    #psmeHeightFromDiameter$gam$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
    
    saveRDS(list(gamBalPhysio = psmeHeightFromDiameter$gamBalPhysio,
                 gamBalPhysioRelDbh = psmeHeightFromDiameter$gamBalPhysioRelDbh,
                 gamPhysio = psmeHeightFromDiameter$gamPhysio,
                 gamRelDbhPhysio = psmeHeightFromDiameter$gamRelDbhPhysio), 
            paste0("trees/height-diameter/data/PSME height physio GAMs ", get_cross_validation_data_suffix(), ".Rds"))
  } else {
    psmePrimaryHeightGams = readRDS(paste0("trees/height-diameter/data/PSME height physio GAMs ", get_cross_validation_data_suffix(), ".Rds"))
    psmeHeightFromDiameter$gamBalPhysio = psmePrimaryHeightGams$gamBalPhysio
    psmeHeightFromDiameter$gamBalPhysioRelDbh = psmePrimaryHeightGams$gamBalPhysioRelDbh
    psmeHeightFromDiameter$gamPhysio = psmePrimaryHeightGams$gamPhysio
    psmeHeightFromDiameter$gamRelDbhPhysio = psmePrimaryHeightGams$gamRelDbhPhysio
    rm(psmePrimaryHeightGams)
  }
  
  psmeHeightFromDiameter$randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"),
                                                   psmeCrossValidation, mtry = 2, minNodeSize = 48, sampleFraction = 0.230) # from setup.R tuneRanger @ 500 trees
  psmeHeightFromDiameter$randomForestBalPhysioRelDbh = fit_ranger("random forest RelDbh BA+L physio", c("height", "dbh", "relativeDiameter", "standQmd", "topHeight", "standBasalAreaPerHectare", "basalAreaLarger", "terrainRoughness", "slope", "elevation", "topographicWetnessFD8f"), # from setup.R VSURF
                                                                  psmeCrossValidation, mtry = 2, minNodeSize = 2, sampleFraction = 0.785) # from setup.R tuneRanger @ 500 trees
  
  # saving primary Rds at the end of this block ensures physio GAMs are included
  saveRDS(psmeHeightFromDiameter, paste0("trees/height-diameter/data/PSME height ", get_cross_validation_data_suffix(), ".Rds"))
}


## Douglas-fir height mixed regressions
if (psmeOptions$fitHeightMixed)
{
  psmeHeightFromDiameterMixed = list(linear = fit_lme("linear", fixedFormula = height ~ dbh + I(isPlantation*dbh), psmeCrossValidation, randomFormula = ~1|stand/plot))
  psmeHeightFromDiameterMixed$parabolic = fit_lme("parabolic", fixedFormula = height ~ dbh + I(dbh^2) + I(isPlantation*dbh), psmeCrossValidation, randomFormula = ~1|stand/plot) # I(isPlantation*dbh^2) not significant, 5x200, 5x200 random, 10x100 false convergence

  psmeHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1ra) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psmeCrossValidation, # b1ra, b2ra singularity in backsolve, 10x10 MAE 6.13, AIC 4178
                                                         fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 45, b1 = -0.02, b2 = 1.2, b2p = -0.1)))
  #psmeHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1) * (1 + a1rm) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psmeCrossValidation, # a1rm job leading minor, b1rm, b2rm singularity in backsolve
  #                                                       fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 45, b1 = -0.02, b2 = 1.2, b2p = -0.1)))
  # linear BA+L
  #psmeHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ra + a1ba * standBasalAreaPerHectare) * (1 - exp(b1*dbh))^(b2 + b2bal * basalAreaLarger), psmeCrossValidation, # b1ra, b2ra singularity in backsolve, b2balp not significant
  #                                                          fixedFormula = a1 + a1ba + b1 + b2 + b2bal ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 38, a1ba = 0.2, b1 = -0.02, b2 = 1.2, b2bal = -0.003)))
  #psmeHeightFromDiameterMixed$chapmanRichardsBalRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1ra + a1rd * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2bal * basalAreaLarger), psmeCrossValidation,
  #                                                                fixedFormula = a1 + a1rd + b1 + b2 + b2bal ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 70, a1rd = -1.7, b1 = -0.015, b2 = 1.3, b2bal = -0.005)))
  #psmeHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ra + a1ba * standBasalAreaPerHectare) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation + b2bal * basalAreaLarger), psmeCrossValidation, # b1ra, b2ra singularity in backsolve, b2balp not significant
  #                                                                fixedFormula = a1 + a1ba + b1 + b1ac + b2 + b2p + b2bal ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 40, a1ba = 0.3, b1 = -0.03, b1ac = -0.002, b2 = 1.4, b2p = 0.1, b2bal = -0.005)))
  #psmeHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger), psmeCrossValidation, # b1ra, b2ra singularity in backsolve
  #                                                                      fixedFormula = a1 + b1 + b1ac + b2 + b2rd + b2bal ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                      start = list(fixed = c(a1 = 48, b1 = -0.03, b1ac = -0.001, b2 = 1.3, b2rd = 0.3, b2bal = -0.003)))
  # reciprocal BA+L
  psmeHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1)*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, # b1ra, b2ra singularity in backsolve
                                                            fixedFormula = a1 + b1 + b2 + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 62, b1 = -0.02, b2 = 1.0, b2bal = 4, b2balr = 7)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1ra) * (1 - exp((b1)*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, # b1ra, b2ra singularity in backsolve
                                                                  fixedFormula = a1 + a1rd + b1 + b2 + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 70, a1rd = -1.5, b1 = -0.02, b2 = 1.0, b2bal = 4, b2balr = 8))) # 2x50 max iterations, 2x250 singularity in backsolve, 10x100 step halving
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation, # b1ra, b2ra singularity in backsolve
                                                                  fixedFormula = a1 + b1 + b1ac + b2 + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 65, b1 = -0.01, b1ac = -0.001, b2 = 0.9, b2bal = 5, b2balr = 8)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter + b2bal / (b2balr + basalAreaLarger)), psmeCrossValidation,
                                                                        fixedFormula = a1 + b1 + b1ac + b2 + b2rd + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 61, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2rd = 0.1, b2bal = 4, b2balr = 8)), distinct = FALSE)
  psmeHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), psmeCrossValidation, 
                                                               fixedFormula = a1 + b1 + b1ac + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 45, b1 = -0.02, b1ac = -0.002, b2 = 1.2, b2p = 0.1)))
  psmeHeightFromDiameterMixed$chapmanRichardsRelDbh = fit_nlme("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1ra) *(1 - exp(b1*dbh))^(b2 + (b2rd + b2rdp * isPlantation) * relativeDiameter), psmeCrossValidation, # b1ra, b2ra singularity in backsolve
                                                                fixedFormula = a1 + b1 + b2 + b2rd + b2rdp ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 46, b1 = -0.03, b2 = 1.2, b2rd = 0.4, b2rdp = -0.1)))
  psmeHeightFromDiameterMixed$chapmanRichardsRelDbhPhysio = fit_nlme("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1ra) *(1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2rd * relativeDiameter), psmeCrossValidation, # b1ra singularity in backsolve, b2ra max iterations, b2rdp not significant
                                                                     fixedFormula = a1 + b1 + b1ac + b2 + b2rd ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 46, b1 = -0.03, b1ac = -0.002, b2 = 1.2, b2rd = 0.4)))
  psmeHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^(b1 + b1p*isPlantation + b1ra), psmeCrossValidation, # a1ra job max iterations, 10x10 RMSE 5.36, AIC 4142
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 3, a1p = -0.8, b1 = 0.4, b1p = -0.1)))
  #psmeHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1p*isPlantation) * dbh / (1 + dbh)^((b1 + b1p*isPlantation)* (1 + b1rm)), psmeCrossValidation, # a1rm step halving, 10x10 RMSE 6.28, AIC 4286
  #                                              fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 3, a1p = -0.6, b1 = 0.4, b1p = -0.04)))
  psmeHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + (a2) * dbh^((b1 + b1p * isPlantation) * (1 + a1rm))), psmeCrossValidation, # a1rm computationally singular, a2rm step halving, 10x10 RMSE 4.73, AIC 4006
                                                  fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 75, a2 = 80, b1 = -1.1, b1p = 0.05))) # 10x40 step halving
  #psmeHeightFromDiameterMixed$hossfeldA1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1ra) / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psmeCrossValidation, # 10x10 RMSE 5.89, AIC 4144
  #                                                  fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 75, a2 = 80, b1 = -1.1, b1p = 0.05)))
  #psmeHeightFromDiameterMixed$hossfeldB1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + (a2) * dbh^(b1 + b1p * isPlantation + b1ra)), psmeCrossValidation, # 10x10 RMSE 4.96, AIC 4023
  #                                                  fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 75, a2 = 80, b1 = -1.1, b1p = 0.05)))
  #psmeHeightFromDiameterMixed$hossfeldA1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1ra) / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psmeCrossValidation,
  #                                                fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 50, a2 = 90, b1 = -1.3, b1p = 0.01)))
  #psmeHeightFromDiameterMixed$hossfeldA2 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + (a2 + a2ra) * dbh^(b1 + b1p * isPlantation)), psmeCrossValidation, # job false convergence
  #                                                fixedFormula = a1 + a2 + b1 + b1p ~ 1, randomFormula = a2ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 75, a2 = 75, b1 = -1.0, b1p = 0.06)))
  psmeHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1)*exp((b1 + b1p * isPlantation)*dbh^((b2 + b2p * isPlantation)*(1 + b2rm))), psmeCrossValidation, # a1rm job step halving, b1rm singularity in backsolve, 10x10 RMSE 4.98, AIC 4076
                                              fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                              start = list(fixed = c(a1 = 140, b1 = -7, b1p = 0.5, b2 = -0.4, b2p = 0.04))) # 2x250 step halving
  #psmeHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1ra)*exp((b1 + b1p * isPlantation)*dbh^(b2 + b2p * isPlantation)), psmeCrossValidation, # b1ra max iterations, b2ra step halving, 10x10 RMSE 5.48, AIC 4087
  #                                            fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                            start = list(fixed = c(a1 = 90, b1 = -8, b1p = 2, b2 = -0.5, b2p = 0.1)))
  psmeHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1) / (a2 + a2p * isPlantation + a2ra + dbh^(b1)), psmeCrossValidation, # 10x10 RMSE 4.89, AIC 3932
                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 50, a2 = 100, a2p = -20, b1 = 1.3))) # 2x50, 2x250 leading minor
  #psmeHeightFromDiameterMixed$michaelisMentenA1 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1ra)*dbh^(b1) / (a2 + a2p * isPlantation + dbh^(b1)), psmeCrossValidation, # 10x10 RMSE 6.49, AIC 4286
  #                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 50, a2 = 100, a2p = -20, b1 = 1.3)))
  #psmeHeightFromDiameterMixed$michaelisMentenB1 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + b1ra) / (a2 + a2p * isPlantation + dbh^(b1 + b1ra)), psmeCrossValidation, # job false convergence, 10x10 RMSE 4.95, AIC 3988
  #                                                         fixedFormula = a1 + a2 + a2p + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 75, a2 = 75, a2p = -8, b1 = 1.0)))
  #psmeHeightFromDiameterMixed$michaelisMentenB1 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1*(1 + a2rm)) / (a2 + dbh^(b1*(1 + a2rm))), psmeCrossValidation, # a1rm job step halving, a2p not significant, a2rm singularity in backsolve, 10x10 RMSE 4.90, AIC 4021
  #                                                         fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a2rm ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 70, a2 = 90, b1 = 1.1)))
  psmeHeightFromDiameterMixed$michaelisMentenBal = fit_nlme("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare + b1ra) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1ba / standBasalAreaPerHectare + b1ra)), psmeCrossValidation, # a1ra step halving, a2ra job step halving, b1ra job false convergence
                                                            fixedFormula = a1 + a2 + a2bal + a2balr + b1 + b1ba ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 85, a2 = 115, a2bal = 900, a2balr = 10, b1 = 1.2, b1ba = -1))) # 2x250 max iterations, 10x25 observations step halving
  psmeHeightFromDiameterMixed$michaelisMentenBalRelDbh = fit_nlme("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + a2rd * relativeDiameter + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1ba / standBasalAreaPerHectare)), psmeCrossValidation, 
                                                                  fixedFormula = a1 + a2 + a2rd + a2bal + a2balr + b1 + b1ba ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 85, a2 = 115, a2rd = 25, a2bal = 900, a2balr = 10, b1 = 1.2, b1ba = -1)), distinct = FALSE)
  #psmeHeightFromDiameterMixed$michaelisMentenRelDbh = fit_nlme("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a2ra)*dbh^b1 / (a2 + a2p * isPlantation + a2rd * relativeDiameter + dbh^b1), psmeCrossValidation, # a2ra singularity in backsolve, b1ra and b1rm computationally singular
  #                                                             fixedFormula = a1 + a2 + a2p + a2rd + b1 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 80, a2 = 150, a2rd = 25, a2p = 15, b1 = 1.2))) # 10x10 and 2x50 step halving
  psmeHeightFromDiameterMixed$michaelisMentenRelDbh = create_fit_statistics("Michaelis-Menten RelDbh", fitting = "nlme")
  psmeHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^((b1 + b1p * isPlantation)*(1 + b1rm)), psmeCrossValidation, # 10x10 MAE 4.08, AIC 4117
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2.0, a1p = -0.6, b1 = 0.7, b1p = 0.08)))
  #psmeHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation + a1ra)*dbh^(b1 + b1p * isPlantation), psmeCrossValidation, # 10x10 MAE 4.18, AIC 4152 vs 4.01, 4047 in fixed effects
  #                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.9, a1p = -0.5, b1 = 0.7, b1p = 0.07)))
  #psmeHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation)*(1 + a1ra)*dbh^(b1 + b1p * isPlantation), psmeCrossValidation, # step halving on (a1) * ((1 + a1ra) * dbh)^(b1), 10x10 MAE 4.17, AIC 4239
  #                                             fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 2.6, a1p = -0.9, b1 = 0.6, b1p = 0.07)))
  #psmeHeightFromDiameterMixed$powerB1 = fit_nlme("power", height ~ 1.37 + (a1 + a1p * isPlantation)*dbh^(b1 + b1p * isPlantation + b1ra), psmeCrossValidation, # 10x10 MAE 4.10, AIC 4135
  #                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2.2, a1p = -0.7, b1 = 0.67, b1p = 0.08)))
  psmeHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1)*dbh^2 + (a2 + a2ra)*dbh + a3), psmeCrossValidation, # a3ra max iterations, a3, a3p not significant, 10x10 RMSE 5.31, AIC 4074
                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.014, a2 = 0.9, a3 = 1)))
  #psmeHeightFromDiameterMixed$prodanA1 = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1ra)*dbh^2 + (a2)*dbh + a3), psmeCrossValidation, # 10x10 RMSE 6.67, AIC 4235
  #                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.018, a2 = 0.7, a3 = 2.6)))
  #psmeHeightFromDiameterMixed$prodanA2 = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1)*(1 + a2rm)*dbh^2 + (a2)*dbh + a3), psmeCrossValidation, # a2rm step halving, 10x10 RMSE 5.92, AIC 4412
  #                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.02, a2 = 0.7, a3 = 3)))
  #psmeHeightFromDiameterMixed$prodanA1 = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1ra)*dbh^2 + (a2)*dbh + a3), psmeCrossValidation, # a2ra clearly preferable
  #                                              fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 0.018, a2 = 0.7, a3 = 2.7)))
  psmeHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1 + b1p * isPlantation + b1ra)/(dbh + b2)), psmeCrossValidation, # 10x10 MAE 3.67, AIC 3972
                                                   fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 63, b1 = -35, b1p = -4, b2 = 8))) # 10x100 sep halving
  #psmeHeightFromDiameterMixed$ratkowskyA1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1ra)*exp((b1 + b1p * isPlantation)/(dbh + b2)), psmeCrossValidation, # 10x10 MAE 5.18, AIC 4393
  #                                                   fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 53, b1 = -29, b1p = 3, b2 = 5)))
  #psmeHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1)*(1 + b1rm)/(dbh + b2)), psmeCrossValidation, # a1rm, b1rm, b2rm step halving, dbh * (1 + b2rm) singularity in backsolve, 10x10 MAE 
  #                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 63, b1 = -35, b2 = 8)))
  psmeHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation) * (1 + Harm) * (1 + ((1.37/((Ha + Hap*isPlantation) * (1 + Harm)))^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), psmeCrossValidation, # drm step halving, kUrm job singularity in backsolve, 10x10 MAE 3.76, AIC 4175
                                                   fixedFormula = Ha + Hap + d + kU ~ 1, randomFormula = Harm ~ 1|stand/plot,
                                                   start = list(fixed = c(Ha = 47, Hap = -5, d = 0.4, kU = 0.018)))
  #psmeHeightFromDiameterMixed$richardsWha = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Hap*isPlantation + Hara) * (1 + ((1.37/(Ha + Hap*isPlantation + Hara))^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), psmeCrossValidation, # kUra max iterations, dra step halving, dp not significant, 10x10 MAE 4.02, AIC 4186
  #                                                   fixedFormula = Ha + Hap + d + kU ~ 1, randomFormula = Hara ~ 1|stand/plot,
  #                                                   start = list(fixed = c(Ha = 47, Hap = -5, d = 0.4, kU = 0.018)))
  psmeHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh^(b4d*(1 + b4rm)))), psmeCrossValidation, # b3 not significant, 10x10 MAE 3.10, AIC 3775
                                                        fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4rm ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = 0, b4d = 1.0))) # 2x50, 2x250 max iterations
  #psmeHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1 * (1 + a1rm) * topHeight^((b1))*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, psmeCrossValidation, # b1rm, b3rm, b4m step halving, 10x10 MAE 3.14, AIC 3819
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*(1 + a1rm)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh^(b4d))), psmeCrossValidation, # b1rm step halving, b3 not significant, b3rm leading minor, 10x10 MAE 3.16, AIC 3753
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 27, b1 = 0.3, b2 = -0.02, b3 = 0, b4d = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(1 + b2rm)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh^(b4d))), psmeCrossValidation, # b3 not significant, 10x10 MAE 3.18, AIC 3779
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b2rm ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4d = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonB1 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1ra)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4d, psmeCrossValidation, # a1ra, b2ra, b4ra singularity in backsolve, b3 not significant, 10x10 MAE 3.24, AIC 3814
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4d = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonB2 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1 * topHeight^((b1))*(1 - exp(b2 * (1 + b2rm) * (standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4d, psmeCrossValidation, # 10x10 MAE 3.22, AIC 3809
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b2rm ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 22, b1 = 0.3, b2 = -0.02, b3 = -0.01, b4d = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonB3 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3ra)*dbh))^b4, psmeCrossValidation, # false convergence
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1) * (1 + a1rm) * topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh))^b4, psmeCrossValidation, # 10x10 MAE 2.98, AIC 3709
                                                         fixedFormula = a1 + b1 + b2 + b3 + b3bal + b3balr + b4 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3bal = -2, b3balr = 14, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalA1 = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh))^b4, psmeCrossValidation, # a1ra, b4ra singularity in backsolve, b2ra step halving, 10x10 MAE 3.01, AIC 3669
  #                                                         fixedFormula = a1 + b1 + b2 + b3 + b3bal + b3balr + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3bal = -2, b3balr = 14, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalB3 = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger) + b3ra) * dbh))^b4, psmeCrossValidation, # 10x10 MAE 3.05, AIC 3727
  #                                                         fixedFormula = a1 + b1 + b2 + b3 + b3bal + b3balr + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 24, b1 = 0.3, b2 = -0.02, b3 = -0.02, b3bal = -2, b3balr = 14, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalB1 = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^b4, psmeCrossValidation, # a1ra, b2ra, b4ra singularity in backsolve
  #                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 40, b1 = 0.1, b2 = -0.02, b3 = -0.02, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalB3 = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3ra)*dbh))^b4, psmeCrossValidation, # b1ra clearly preferable
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 40, b1 = 0.1, b2 = -0.02, b3 = -0.02, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect)) * (1 + a1ra) * topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4 + b4e * elevation), # 10x10 MAE 2.86, AIC 3688
                                                               fixedFormula = a1 + a1ac + b1 + b2 + b3 + b3bal + b3balr + b4 + b4e ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               psmeCrossValidation, start = list(fixed = c(a1 = 38, a1ac = 1.7, b1 = 0.2, b2 = -0.02, b3 = 0.2, b3bal = -7, b3balr = 30, b4 = 1.1, b4e = 0.0002)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1 + b1ra) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4 + b4e * elevation), # a1ra step halving, b2ra, b4ra singularity in backsolve, 10x10 MAE 2.89, AIC 3679
  #                                                             fixedFormula = a1 + a1ac + b1 + b2 + b3 + b3bal + b3balr + b4 + b4e ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                             psmeCrossValidation, start = list(fixed = c(a1 = 30, a1ac = 1.0, b1 = 0.2, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 14, b4 = 1.0, b4e = 0.0001)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysioB3 = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh) + b3ra)^(b4 + b4e * elevation), # job singularity in backsolve
  #                                                               fixedFormula = a1 + a1ac + b1 + b2 + b3 + b3bal + b3balr + b4 + b4e ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                               psmeCrossValidation, start = list(fixed = c(a1 = 30, a1ac = 1.0, b1 = 0.2, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 14, b4 = 1.0, b4e = 0.0001)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation), psmeCrossValidation, # a1ra, b24, b4ra singularity in backsolve
  #                                                             fixedFormula = a1 + a1ac + b1 + b2 + b3 + b4 + b4e ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 40, a1ac = 1.5, b1 = 0.1, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4e = 0.0003)))
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysioB3 = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3ra)*dbh))^(b4 + b4e * elevation), psmeCrossValidation, # inclination to a1-b1 evaporation
  #                                                             fixedFormula = a1 + a1ac + b1 + b2 + b3 + b4 + b4e ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 25, a1ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0, b4e = 0.0002)))
  psmeHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*topHeight^(b1 + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4), psmeCrossValidation,
                                                                     fixedFormula = a1 + a1ac + b1 + b2 + b3 + b3bal + b3balr + b4 + b4e ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                     start = list(a1 = 30, a1rd = 0, b1 = 0.2, b1as = -0.04, b1ac = 0.01, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 14, b4 = 1.0), distinct = FALSE)
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra + b1rd * relativeDiameter + b1as * sin(3.14159/180 * aspect) + b1ac * cos(3.14159/180 * aspect)) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4), psmeCrossValidation, # a1ra, b2ra, b4ra singularity in backsolve, b3ra max iterations, a1rd not significant, b2rd, b3rd, b4rd all slow to converge compared to b1rd
  #                                                                   fixedFormula = a1 + b1 + b1rd + b1as + b1ac + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                                   start = list(fixed = c(a1 = 35, b1 = 0.1, b1rd = 0.002, b1as = 0.005, b1ac = 0.01, b2 = -0.02, b3 = -0.15, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*topHeight^(b1 + b1ra) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger)) * dbh))^(b4), psmeCrossValidation, 
                                                               fixedFormula = a1 + a1rd + b1 + b2 + b3 + b3bal + b3balr + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 30, a1rd = 0, b1 = 0.2, b2 = -0.02, b3 = 0.1, b3bal = -2, b3balr = 15, b4 = 1.0)), distinct = FALSE)
  #psmeHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4rd * relativeDiameter), psmeCrossValidation, # b2ra max iterations, b3ra, b4ra singularity in backsolve, a1rd not significant
  #                                                             fixedFormula = a1 + b1 + b2 + b3 + b4 + b4rd ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 40, b1 = 0.1, b2 = -0.02, b3 = -0.2, b4 = 1.0, b4rd = -0.03)))
  psmeHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, # a1ra, b4ra singularity in backsolve, b2ra step halving, b3 not significant, 10x10 MAE 3.19, AIC 3856
                                                            fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac + b4as ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  #psmeHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, # 10x10 MAE 3.22, AIC 3750
  #                                                          fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac + b4as ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 18, b1 = 0.3, b2 = -0.02, b3 = -0.03, b4 = 1.0, b4ac = -0.05, b4as = -0.05)))
  #psmeHeightFromDiameterMixed$sharmaPartonPhysioB3 = fit_nlme("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3ra)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, # job false convergence
  #                                                          fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac + b4as ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  psmeHeightFromDiameterMixed$sharmaPartonRelDbh = fit_nlme("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1ra + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, psmeCrossValidation, # fixed effects not significant
                                                            fixedFormula = a1 + b1 + b1rd + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 13, b1 = 0.3, b1rd = 0, b2 = -0.02, b3 = 0.02, b4 = 1.0)), distinct = FALSE)
  psmeHeightFromDiameterMixed$sharmaPartonRelDbhPhysio = fit_nlme("Sharma-Parton RelDbh physio", height ~ 1.37 + a1*topHeight^(b1 + b1ra + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4as * sin(3.14159/180 * aspect) + b4ac * cos(3.14159/180 * aspect)), psmeCrossValidation, # fixed effects not significant
                                                                  fixedFormula = a1 + b1 + b1rd + b2 + b3 + b4 + b4as + b4ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 18, b1 = 0.3, b1rd = 0.004, b2 = -0.02, b3 = 0, b4 = 1.0, b4ac = -0.05, b4as = -0.04)))
  psmeHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*standTreesPerHectare^(b3 + b3ra)*dbh))^(b4), psmeCrossValidation, # a1ra, b4ra singularity in backsolve, b2ra step halving, a1p not significant, 10x10 MAE 3.26, AIC 3874
                                                     fixedFormula = a1 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 35, b1 = 0.1, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaZhangB4 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1 + b1ra)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh^(b4d))), psmeCrossValidation, # a1ra job step halving, b2ra step halving, 10x10 MAE 3.43, AIC 3801
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 23, b1 = 0.2, b2 = -0.04, b3 = -0.1, b4d = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaZhangB4 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3 + b3ra)*dbh^(b4d))), psmeCrossValidation, # 10x10 MAE 3.38, AIC 3803
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 37, b1 = 0.14, b2 = -0.04, b3 = -0.1, b4d = 0.9)))
  #psmeHeightFromDiameterMixed$sharmaZhangB4 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh^(b4d + b4ra))), psmeCrossValidation, # 10x10 MAE 3.27, AIC 3847
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 33, b1 = 0.15, b2 = -0.04, b3 = -0.1, b4d = 0.9)))
  #psmeHeightFromDiameterMixed$sharmaZhangB4 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1)*(1 + a1rm)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh^(b4d))), psmeCrossValidation, # b1rm, b2rm step halving, b3rm singularity in backsolve, 10x10 MAE 3.40, AIC 3874
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 23, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4d = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaZhangB4 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh^(b4d*(1 + b4rm)))), psmeCrossValidation, # 10x10 MAE 3.30, AIC 3822
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 33, b1 = 0.15, b2 = -0.03, b3 = -0.1, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*(1 + a1rm)*standBasalAreaPerHectare^(b1)*(1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4), psmeCrossValidation, # b1p not significant, b1rm, b2rm, b3rm, b4rm singularity in backsolve, 10x10 MAE 3.40, AIC 3921
  #                                                   fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 24, b1 = 0.2, b2 = -0.03, b3 = -0.09, b4 = 1.0)))
  #psmeHeightFromDiameterMixed$sharmaZhangB1 = fit_nlme("Sharma-Zhang", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1ra)*(1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4), psmeCrossValidation, # b3ra preferable
  #                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 35, a1p = -10, b1 = 0.1, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0)))
  psmeHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1) * (1 - exp(b2*standTreesPerHectare^(b3 + b3ra) * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), psmeCrossValidation, # a1ra, b2ra, b4ra singularity in backsolve, 10x10 MAE 3.12, AIC 3760
                                                        fixedFormula = a1 + b1 + b2 + b3 + b4 + b4bal + b4balr ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 55, b1 = 0.1, b2 = -0.02, b3 = -0.03, b4 = 1.0, b4bal = 4, b4balr = 6)))
  #psmeHeightFromDiameterMixed$sharmaZhangBalB1 = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1ra) * (1 - exp(b2*standTreesPerHectare^(b3) * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), psmeCrossValidation, # 10x10 MAE 3.21, AIC 3757
  #                                                        fixedFormula = a1 + b1 + b2 + b3 + b4 + b4bal + b4balr ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 55, b1 = 0.1, b2 = -0.02, b3 = -0.03, b4 = 1.0, b4bal = 4, b4balr = 6)))
  #psmeHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1p * isPlantation + a1ba * standBasalAreaPerHectare + a1bal * basalAreaLarger) * standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1ra) * (1 - exp(b2 * standTreesPerHectare^(b3)*dbh))^(b4), psmeCrossValidation, # a1ra, b4ra singularity in backsolve, b2ra step halving, b2ba, b2bap not significant
  #                                                      fixedFormula = a1 + a1p + a1ba + a1bal + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 50, a1p = -70, a1ba = 40, a1bal = 11, b1 = -0.9, b1p = 0.01, b2 = -0.02, b3 = -0.03, b4 = 1.04)))
  psmeHeightFromDiameterMixed$sharmaZhangBalRelDbh = fit_nlme("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1ra) * (1 - exp(b2*standTreesPerHectare^b3 * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger) + b4rd * relativeDiameter), psmeCrossValidation, 
                                                              fixedFormula = a1 + b1 + b2 + b3 + b4 + b4bal + b4balr + b4rd ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 55, b1 = 0.1, b2 = -0.02, b3 = -0.03, b4 = 1.0, b4bal = 4, b4balr = 6, b4rd = 0)), distinct = FALSE)
  psmeHeightFromDiameterMixed$sharmaZhangRelDbh = fit_nlme("Sharma-Zhang RelDbh", height ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation + b1ra)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^(b4 + b4rd * relativeDiameter), psmeCrossValidation, 
                                                           fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 + b4rd ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                           start = list(fixed = c(a1 = 45, a1p = -12, b1 = 0.08, b1p = 0.08, b2 = -0.03, b3 = -0.06, b4 = 1.0, b4rd = 0)), distinct = FALSE)
  psmeHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + (a1)*dbh^((b1 + b1p * isPlantation)*(1 + b1rm)*dbh^b2), psmeCrossValidation, # a1rm job step halving, b2rm max iterations, 10x10 MAE 4.02, AIC 4016
                                                  fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.09, b1 = 3.0, b1p = -0.06, b2 = -0.17))) # 5x200 random step halving
  #psmeHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + (a1)*dbh^((b1 + b1p * isPlantation + b1ra)*dbh^b2), psmeCrossValidation, # a1ra step halving, b2p not significant, 10x10 MAE 4.09, AIC 4082
  #                                                fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.5, b1 = 1.8, b1p = -0.05, b2 = -0.13)))
  #psmeHeightFromDiameterMixed$sibbesenB2 = fit_nlme("Sibbesen", height ~ 1.37 + (a1)*dbh^((b1 + b1p * isPlantation)*dbh^(b2 + b2ra)), psmeCrossValidation, # b1ra preferable
  #                                                  fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 0.5, b1 = 1.8, b1p = -0.05, b2 = -0.13)))
  psmeHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1ra)*(1 - exp((b1 + b1p * isPlantation)*dbh^b2)), psmeCrossValidation, # b1ra, b2ra max iterations, 10x10 MAE 4.67, AIC 4262
                                                 fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 45, b1 = -0.01, b1p = 0.001, b2 = 1.1)))
  #psmeHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1)*(1 - exp((b1 + b1p * isPlantation)*(1 + b1rm)*dbh^b2)), psmeCrossValidation, # a1rm job step halving, b1p not significant, b1rm job step halving, b2rm max iterations, 10x10 MAE 
  #                                               fixedFormula = a1 + b1 + b1p + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 60, b1 = -0.02, b1p = 0.002, b2 = 1.0)))
  psmeHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2 + b2bal * basalAreaLarger + b2ra))), psmeCrossValidation,
                                                    fixedFormula = a1 + b1 + b1ba + b2 + b2bal ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 64, b1 = -0.01, b1ba = 0.03, b2 = 1.0, b2bal = 0.002))) # 5x200 random step halving
  #psmeHeightFromDiameterMixed$weibullBalA1 = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2 + b2bal * basalAreaLarger))), psmeCrossValidation, # b1ra step halving, job false convergence
  #                                                    fixedFormula = a1 + b1 + b1ba + b2 + b2bal ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 64, b1 = -0.01, b1ba = 0.03, b2 = 1.0, b2bal = 0.002)))
  #psmeHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger + b2ra))), psmeCrossValidation, # a1ra max iterations, b1ra step halving
  #                                                  fixedFormula = a1 + b1 + b2 + b2bal ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 100, b1 = -0.01, b2 = 0.9, b2bal = 0.003)))
  #psmeHeightFromDiameterMixed$weibullBalRelDbh = fit_nlme("Weibull RelDbh BA+L", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger))), psmeCrossValidation, # a1ra max iterations, b1ra, b2ra step halving
  #                                                        fixedFormula = a1 + b1 + b1ba + b2 + b2rd + b2bal ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 65, b1 = -0.01, b1ba = 0.03, b2 = 1.0, b2rd = -0.005, b2bal = 0.002)))
  #psmeHeightFromDiameterMixed$weibullBalRelDbh = fit_nlme("Weibull RelDbh BA+L", height ~ 1.37 + (a1) * (1 + a1rm) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^(b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger))), psmeCrossValidation, # b1rm step halving
  #                                                        fixedFormula = a1 + b1 + b1ba + b2 + b2rd + b2bal ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 67, b1 = -0.01, b1ba = 0, b2 = 1.0, b2rd = 0, b2bal = 0.002)), distinct = FALSE) # b1ba, b2rd not significant
  psmeHeightFromDiameterMixed$weibullBalRelDbh = fit_nlme("Weibull RelDbh BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ba / standBasalAreaPerHectare)*dbh^((b2 + b2rd * relativeDiameter + b2bal * basalAreaLarger) * (1 + b2rm)))), psmeCrossValidation,
                                                          fixedFormula = a1 + b1 + b1ba + b2 + b2rd + b2bal ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 65, b1 = -0.01, b1ba = 0.03, b2 = 1.0, b2rd = -0.005, b2bal = 0.002)), distinct = FALSE) # b1ba not significant
  psmeHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1)*(1 - exp((b1 + b1p * isPlantation)*dbh^(b2 + b2rd * relativeDiameter + b2ra))), psmeCrossValidation, # a1ra max iterations, b1ra step halving, 10x10 MAE 3.57, AIC 3931
                                                       fixedFormula = a1 + b1 + b1p + b2 + b2rd ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 70, b1 = -0.01, b1p = 0.001, b2 = 1.0, b2rd = -0.02))) # 2x50 max iterations, 2x250 step halving
  #psmeHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1)*(1 + a1rm)*(1 - exp((b1 + b1p * isPlantation)*dbh^((b2 + b2rd * relativeDiameter)))), psmeCrossValidation, # a1rm job max iterations, b1rm step halving, b2rm max iterations
  #                                                     fixedFormula = a1 + b1 + b1p + b2 + b2rd ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 54, b1 = -0.01, b1p = 0.001, b2 = 1.1, b2rd = -0.01)))
  #to_fixed_coeffficients(psmeHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(psmeHeightFromDiameterMixed$sharmaZhang), n = 40)
  
  # fit with gamm() because gam() with random effect basis functions is unhappy, gamma links not used as convergence fails
  # basis  fit time, 9900X  threads  warnings
  # fs     ~80 s            4        yes
  # re     ~80 s            4        yes
  # sz     ~100 s           4        yes
  # mrf is impractical due to requiring lists of coordinates forming plot polygons
  psmeHeightFromDiameterMixed$gam = fit_gam("REML GAM", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, bs = "ts", k = 8), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation) # 10x10 AIC 4144, median MAE 3.60 m vs ~3.41 m for fixed effects, niterPQL convergence with Γ link
  psmeHeightFromDiameterMixed$gamBal = fit_gam("REML GAM BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeHeightFromDiameterMixed$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 13), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeHeightFromDiameterMixed$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeHeightFromDiameterMixed$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 13), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeHeightFromDiameterMixed$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, basalAreaLarger, relativeDiameter, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 18), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeHeightFromDiameterMixed$gamPhysio = fit_gam("REML GAM physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeHeightFromDiameterMixed$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + s(dbh, relativeDiameter, topographicWetnessFD8f, bs = "ts", k = 13, by = as.factor(isPlantation)), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  #lapply(psmeHeightFromDiameterMixed$gam$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(psmeHeightFromDiameterMixed$gam$fit, function(fit) { return(summary(fit$gam)) })
  #psmeHeightFromDiameterMixed$gam$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(psmeHeightFromDiameterMixed, paste0("trees/height-diameter/data/PSME height mixed ", get_cross_validation_data_suffix(), ".Rds"))
}


## Douglas-fir diameter regressions
if (psmeOptions$fitDbhPrimary)
{
  # linear regressions
  psmeDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)), psmeCrossValidation))
  psmeDiameterFromHeight$linearAbat = fit_lm("linear ABA+T", dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)) + I(1/(1 + basalAreaTaller)), psmeCrossValidation) # a1aba, a1abar, a1aat, a1aatrp not significant
  psmeDiameterFromHeight$linearAbatRelHt = fit_lm("linear RelHt ABA+T", dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)) + relativeHeight + I(1/(1 + basalAreaTaller)), psmeCrossValidation) # a1rhp not significant
  psmeDiameterFromHeight$linearRelHt = fit_lm("linear RelHt", dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)) + relativeHeight + I(isPlantation*relativeHeight), psmeCrossValidation)
  psmeDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ I(height - 1.37) + I((height - 1.37)^2) + I(isPlantation*(height - 1.37)^2), psmeCrossValidation) # plantation not significant
  #lapply(psmeDiameterFromHeight$linearRelHt$fit, summary)
  
  # nonlinear regressions
  psmeDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation)*(exp(b1*(height - 1.37)) - 1)^b2, psmeCrossValidation, 
                                                      start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8)) # b1p, b2p not significant
  psmeDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation)*(exp((b1 + b1aba * bootstrapStandBasalAreaPerHectare + b1aat * basalAreaTaller)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
                                                          start = list(a1 = 50, a1p = -3, b1 = 0.03, b1aba = -0.0001, b1aat = -0.0001, b2 = 0.9)) # a1aat+r, b1aatp, b1aat+r, b1abap, b2aat+r not significant, 10x10 RMSE 14.5, AIC 5301
  #psmeDiameterFromHeight$chapmanReplaceAbatA1 = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a1aat * basalAreaTaller)*(exp((b1 + b1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                          start = list(a1 = 50, a1p = -3, a1aat = -0.1, b1 = 0.03, b1aba = -0.0001, b2 = 0.9)) # 10x10 RMSE 14.5, AIC 5300
  #psmeDiameterFromHeight$chapmanReplaceAba = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation)*(exp((b1 + b1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                       start = list(a1 = 40, a1p = -3, b1 = 0.02, b1aba = -0.0001, b2 = 0.9)) # 10x0 RMSE 14.8, AIC 5223
  #psmeDiameterFromHeight$chapmanReplaceAbatA1aba = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a1aba * bootstrapStandBasalAreaPerHectare)*(exp((b1)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                             start = list(a1 = 50, a1p = -3, a1aba = -0.1, b1 = 0.02, b2 = 0.9)) # b2aat not significant, 10x0 RMSE 14.7, AIC 5161
  #psmeDiameterFromHeight$chapmanReplaceAbatB2aba = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation)*(exp((b1)*(height - 1.37)) - 1)^(b2 + b2aba * bootstrapStandBasalAreaPerHectare), psmeCrossValidation, 
  #                                                             start = list(a1 = 150, a1p = -10, b1 = 0.01, b2 = 0.9, b2aba = 0.003)) # 10x0 RMSE 14.8, AIC 5189
  #psmeDiameterFromHeight$chapmanReplaceAbatA1bat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a1aat * basalAreaTaller)*(exp((b1)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                             start = list(a1 = 55, a1p = -3, a1aat = -0.3, b1 = 0.02, b2 = 0.9)) # 10x0 RMSE 14.6, AIC 5203
  #psmeDiameterFromHeight$chapmanReplaceAbatB1bat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation)*(exp((b1 + b1aat * basalAreaTaller)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                             start = list(a1 = 55, a1p = -3, b1 = 0.02, b1aat = -0.0001, b2 = 0.9)) # 10x0 RMSE 14.8, AIC 5213
  #psmeDiameterFromHeight$chapmanReplaceAbatA1batB1aba = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a1aat * basalAreaTaller)*(exp((b1 + b1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                                  start = list(a1 = 50, a1p = -3, a1aat = -0.1, b1 = 0.02, b1aba = -0.00005, b2 = 0.9)) # 10x0 RMSE 14.7, AIC 5168
  psmeDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1)*(exp((b1 + b1rh * relativeHeight)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
                                                           start = list(a1 = 45, b1 = 0.03, b1rh = 0.002, b2 = 0.8)) # b2rh not significant, a1p-b1rh not mutually significant, 10x10 RMSE 14.7, AIC 5227
  #psmeDiameterFromHeight$chapmanReplaceRelHtA1 = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(exp((b1)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, 
  #                                                           start = list(a1 = 35, a1rh = 5, b1 = 0.03, b2 = 0.8)) # 10x10 RMSE 15.2, AIC 5318
  psmeDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), psmeCrossValidation, 
                                                       start = list(a1 = -100, b1 = 0.012, b2 = 0.9)) # a1p, b1p, b2p not significant
  psmeDiameterFromHeight$chapmanRichardsAba = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1)*log(1 - pmin(((b1 + (b1aba + b1abap * isPlantation) * bootstrapStandBasalAreaPerHectare)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
                                                          start = list(a1 = -100, b1 = 0.01, b1aba = -0.00002, b1abap = -0.00001, b2 = 1.0)) # a1aat+r, b1aat+r, b2aat, b2aat+r not significant, 10x10 RMSE 14.7, AIC 5255
  #psmeDiameterFromHeight$chapmanRichardsAbaA1 = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + (a1aba + a1abap * isPlantation) * bootstrapStandBasalAreaPerHectare)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                          start = list(a1 = -120, a1aba = 0.01, a1abap = 0.1, b1 = 0.01, b2 = 1.0)) # 10x10 RMSE 14.7, AIC 5250
  #psmeDiameterFromHeight$chapmanRichardsAbaB2 = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + (b2aba + b2abap * isPlantation) * bootstrapStandBasalAreaPerHectare), 0.9999)), psmeCrossValidation, 
  #                                                          start = list(a1 = -120, b1 = 0.01, b2 = 1.0, b2aba = 0.001, b2abap = 0.001)) # 10x10 RMSE 14.8, AIC 5209
  #psmeDiameterFromHeight$chapmanRichardsAbatA1 = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + (a1aba + a1abap * isPlantation) * bootstrapStandBasalAreaPerHectare + a1aat * basalAreaTaller)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                           start = list(a1 = -120, a1aba = 0.2, a1abap = 0.2, a1aat = 0.3, b1 = 0.01, b2 = 1.0)) # 10x10 RMSE 14.7, AIC 5124
  #psmeDiameterFromHeight$chapmanRichardsAbatA1 = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*log(1 - pmin(((b1 + (b1aba + b1abap * isPlantation) * bootstrapStandBasalAreaPerHectare)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                           start = list(a1 = -115, a1aat = 0.3, b1 = 0.01, b1aba = -0.00002, b1abap = -0.00001, b2 = 1.0)) # 10x10 RMSE 14.8, AIC 5236
  #psmeDiameterFromHeight$chapmanRichardsAbatB1 = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1)*log(1 - pmin(((b1 + (b1aba + b1abap * isPlantation) * bootstrapStandBasalAreaPerHectare + b1aat * basalAreaTaller)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                           start = list(a1 = -100, b1 = 0.01, b1aba = -0.00002, b1abap = -0.00001, b1aat = -0.00002, b2 = 1.0)) # 10x10 RMSE 14.6, AIC 5193
  #psmeDiameterFromHeight$chapmanRichardsAbatA1aba = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                              start = list(a1 = -115, a1aba = 0.2, b1 = 0.01, b2 = 1.0)) # 10x10 RMSE 14.9, AIC 5239
  #psmeDiameterFromHeight$chapmanRichardsAbatB2aba = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2aba * bootstrapStandBasalAreaPerHectare), 0.9999)), psmeCrossValidation, 
  #                                                              start = list(a1 = -115, b1 = 0.01, b2 = 0.9, b2aba = 0.002)) # 10x10 RMSE 14.8, AIC 5179
  #psmeDiameterFromHeight$chapmanRichardsAbatA1bat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                              start = list(a1 = -110, a1aat = 0.4, b1 = 0.01, b2 = 0.9)) # 10x10 RMSE 150, AIC 5283
  #psmeDiameterFromHeight$chapmanRichardsAbatB1bat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1)*log(1 - pmin(((b1 + b1aat * basalAreaTaller)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                              start = list(a1 = -110, b1 = 0.01, b1aat = -0.00004, b2 = 1.0)) # 10x10 RMSE 14.9, AIC 5259
  #psmeDiameterFromHeight$chapmanRichardsAbatB2bat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2aat * basalAreaTaller), 0.9999)), psmeCrossValidation, 
  #                                                              start = list(a1 = -100, b1 = 0.01, b2 = 0.9, b2aat = 0.003)) # 10x10 RMSE 15.0, AIC 5291
  #psmeDiameterFromHeight$chapmanRichardsAbatB1aba = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*log(1 - pmin(((b1 + b1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                              start = list(a1 = -100, a1aat = 0.3, b1 = 0.01, b1aba = -0.00002, b2 = 1.0)) # 10x10 RMSE 14.9, AIC 5121
  psmeDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), psmeCrossValidation, 
                                                             start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8, b2tw = 0), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  # no clear differentiation between a1rh, b1rh, b2rh forms, b2rh maybe somewhat more accurate
  psmeDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), psmeCrossValidation, 
                                                            start = list(a1 = -90, b1 = 0.012, b2 = 0.9, b2rh = -0.07)) # a1rh, b1rh also significant
  #psmeDiameterFromHeight$chapmanRichardsRelHtA1 = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1rh * relativeHeight)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                            start = list(a1 = -90, a1rh = -11, b1 = 0.012, b2 = 0.9))
  #psmeDiameterFromHeight$chapmanRichardsRelHtB1 = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1)*log(1 - pmin(((b1 + b1rh * relativeHeight)*(height - 1.37))^(b2), 0.9999)), psmeCrossValidation, 
  #                                                            start = list(a1 = -100, b1 = 0.01, b1rh = 0.001, b2 = 0.9))
  psmeDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
                                                              start = list(a1 = 130, a2 = 70, b1 = 0.9, b1p = -0.008)) # a1p, a2p not significant, 10x10 RMSE 14.9, AIC 5224, 2x500 step halving
  psmeDiameterFromHeight$michaelisMentenReplaceAbat = fit_gsl_nls("Michaelis-Menten replace ABA+T", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, 
                                                                  start = list(a1 = 120, a2 = 70, b1 = 0.9, b1aba = -0.0005, b1p = -0.01)) # 10x10 RMSE 14.4, AIC 5241
  psmeDiameterFromHeight$michaelisMentenReplaceAbatRelHt = fit_gsl_nls("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight) * (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, 
                                                                       start = list(a1 = 110, a1rh = 10, a2 = 60, b1 = 0.9, b1aba = -0.0005, b1p = -0.01))
  #psmeDiameterFromHeight$michaelisMentenReplaceAbaA1 = fit_gsl_nls("Michaelis-Menten replace ABA", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                 start = list(a1 = 170, a1aba = -0.4, a2 = 90, b1 = 0.9, b1p = -0.01)) # 10x10 RMSE 14.6, AIC 5128
  #psmeDiameterFromHeight$michaelisMentenReplaceAbaA2 = fit_gsl_nls("Michaelis-Menten replace ABA", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2aba * bootstrapStandBasalAreaPerHectare - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                 start = list(a1 = 130, a2 = 70, a2aba = 0.2, b1 = 0.9, b1p = -0.01)) # 10x10 RMSE 14.5, AIC 5206
  #psmeDiameterFromHeight$michaelisMentenReplaceAatA1 = fit_gsl_nls("Michaelis-Menten replace AAT", dbh ~ (a1 + a1aat * basalAreaTaller) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                 start = list(a1 = 170, a1aat = -0.8, a2 = 90, b1 = 0.9, b1p = -0.01)) # 10x10 RMSE 14.6, AIC 5265
  #psmeDiameterFromHeight$michaelisMentenReplaceAatA2 = fit_gsl_nls("Michaelis-Menten replace AAT", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2aat * basalAreaTaller - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                 start = list(a1 = 170, a2 = 100, a2aat = 0.4, b1 = 0.9, b1p = -0.01)) # 10x10 RMSE 14.5, AIC 5189
  #psmeDiameterFromHeight$michaelisMentenReplaceAatB1 = fit_gsl_nls("Michaelis-Menten replace AAT", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation + b1aat * basalAreaTaller) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1aat * basalAreaTaller)), psmeCrossValidation, 
  #                                                                 start = list(a1 = 180, a2 = 100, b1 = 0.9, b1aat = -0.001, b1p = -0.01)) # 10x10 RMSE 14.9, AIC 5228
  #psmeDiameterFromHeight$michaelisMentenReplaceAbat = fit_gsl_nls("Michaelis-Menten replace ABA+T", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 + a2aat * basalAreaTaller - (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, 
  #                                                                start = list(a1 = 120, a2 = 70, a2aat = 0.4, b1 = 0.9, b1aba = -0.0005, b1p = -0.01), distinct = FALSE) # a2aat not reliably significant
  psmeDiameterFromHeight$michaelisMentenReplaceRelHt = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
                                                                   start = list(a1 = 130, a1rh = 20, a2 = 70, b1 = 0.9, b1p = -0.01)) # 10x10 RMSE 14.8, AIC 5214
  #psmeDiameterFromHeight$michaelisMentenReplaceRelHtA2 = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2rh * relativeHeight - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                   start = list(a1 = 170, a2 = 85, a2rh = -8, b1 = 0.9, b1p = -0.01)) # 10x10 RMSE 15.1, AIC 5181
  #psmeDiameterFromHeight$michaelisMentenReplaceRelHtB1 = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight)), psmeCrossValidation, 
  #                                                                   start = list(a1 = 190, a2 = 100, b1 = 0.9, b1p = -0.01, b1rh = 0.02)) # 10x10 RMSE 14.9, AIC 5239
  #psmeDiameterFromHeight$michaelisMentenReplaceRelHtA1 = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh / relativeHeight) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                   start = list(a1 = 120, a1rh = -2.5, a2 = 55, b1 = 0.8, b1p = -0.007)) # 10x10 RMSE 14.9, AIC 5291
  #psmeDiameterFromHeight$michaelisMentenReplaceRelHtA2 = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2rh / relativeHeight - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, 
  #                                                                   start = list(a1 = 110, a2 = 45, a2rh = 3, b1 = 0.8, b1p = -0.01)) # 10x10 RMSE 14.9, AIC 5206
  #psmeDiameterFromHeight$michaelisMentenReplaceRelHtB1 = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1) * (height - 1.37)^(b1 + b1p * isPlantation + b1rh / relativeHeight) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1rh / relativeHeight)), psmeCrossValidation, 
  #                                                                   start = list(a1 = 150, a2 = 80, b1 = 0.9, b1p = -0.007, b1rh = -0.02)) # 10x10 RMSE 15.0, AIC 5230
  psmeDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psmeCrossValidation, 
                                               start = list(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006)) # 10x10 RMSE 15.6, AIC 5514
  # a1*(height - 1.37)^b1 * a2*topHeight^b2 singular gradient
  # a1*(height - 1.37)^(b1*topHeight^b2) not significant since b2 not significant
  # a1*(height - 1.37)^(b1*((height - 1.37)/relativeHeight)^b2) not significant since b2 not significant, seems logical since greater relative height implies larger rather than smaller DBH
  psmeDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, 
                                             start = list(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2)) # 10x10 RMSE 15.1, AIC 5311
  psmeDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, 
                                                 start = list(a1 = 1.2, a1aat = -0.01, b1 = 1.1, b1p = -0.2)) # a1p, a1aatp not significant, 10x10 RMSE 14.8, AIC 5241
  #psmeDiameterFromHeight$powerAbatA1aba = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1p * isPlantation + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, 
  #                                                    start = list(a1 = 0.7, a1p = 0.6, a1aba = 0, b1 = 1.2, b1p = -0.2)) # 10x10 RMSE 14.9, AIC 5291
  #psmeDiameterFromHeight$powerAbatB1aba = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1p*isPlantation + b1aba * bootstrapStandBasalAreaPerHectare), psmeCrossValidation, 
  #                                                    start = list(a1 = 1.0, b1 = 1.2, b1p = -0.02, b1aba = -0.001)) # a1p not significant, 10x10 RMSE 15.2, AIC 5273
  #psmeDiameterFromHeight$powerAbatB1bat = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1p*isPlantation + b1aat * basalAreaTaller), psmeCrossValidation, 
  #                                                    start = list(a1 = 1.2, b1 = 1.1, b1p = -0.03, b1aat = -0.002)) # a1p not significant, 10x10 RMSE 14.9, AIC 5205
  #psmeDiameterFromHeight$powerAbatA1abaBat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1aat * basalAreaTaller)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, 
  #                                                       start = list(a1 = 0.7, a1aba = 0, a1aat = -0.01, b1 = 1.2, b1p = -0.2), distinct = FALSE) # a1p not significant, a1aba-a1aat not mutually significant
  psmeDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1e * elevation), psmeCrossValidation, 
                                                   start = list(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2, b1e = 0), distinct = FALSE) # { a1, b1 } x { e, s, a, tr, tw } not significant
  psmeDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight), psmeCrossValidation, 
                                                  start = list(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)) # a1p-a1rh not mutually significant, 10x10 RMSE 14.9, AIC 5277 (17.0, 5498 if multiplying height by relative height)
  psmeDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
                                             start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2p = -0.001)) # a1p not significant, b1p-b2p not mutually significant, 10x10 RMSE 14.9, AIC 5247
  psmeDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
                                                 start = list(a1 = 2.0, a1aba = -0.006, b1 = 0.9, b2 = 0.01, b2p = -0.002)) # a1abap, a1aat+r, b2aat+r not significant, b1aat+r NaN-Inf, 10x10 RMSE 14.6, AIC 5148
  #psmeDiameterFromHeight$ruarkAbaB1 = fit_gsl_nls("Ruark ABA", dbh ~ (a1)*(height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
  #                                                start = list(a1 = 1.9, b1 = 0.9, b1aba = -0.001, b2 = 0.01, b2p = -0.002)) # 10x10 RMSE 14.6, AIC 5272
  #psmeDiameterFromHeight$ruarkAbaB2 = fit_gsl_nls("Ruark ABA", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + b2aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)), psmeCrossValidation, 
  #                                                start = list(a1 = 2.0, b1 = 0.8, b2 = 0.02, b2p = -0.002, b2aba = -0.0001)) # 10x10 RMSE 15.6, AIC 5343
  #psmeDiameterFromHeight$ruarkAbatA1 = fit_gsl_nls("Ruark ABAT", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
  #                                                 start = list(a1 = 2.0, a1aat = -0.01, b1 = 0.9, b2 = 0.01, b2p = -0.002)) # 10x10 RMSE 14.7, AIC 5178
  #psmeDiameterFromHeight$ruarkAbatB1 = fit_gsl_nls("Ruark ABAT", dbh ~ (a1)*(height - 1.37)^(b1 + b1aat * basalAreaTaller) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
  #                                                 start = list(a1 = 2.0, b1 = 0.9, b1aat = -0.001, b2 = 0.01, b2p = -0.002)) # 10x10 RMSE 15.1, AIC 5238
  #psmeDiameterFromHeight$ruarkAbatB2 = fit_gsl_nls("Ruark ABAT", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + b2aat * basalAreaTaller) * (height - 1.37)), psmeCrossValidation, 
  #                                                 start = list(a1 = 2.0, b1 = 0.9, b2 = 0.01, b2p = -0.002, b2aat = -0.0001)) # 10x10 RMSE 14.5, AIC 5359, b2aatp not significant, but not consistently discriminable from AbaA1
  #psmeDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + (b2aat + b2aatp * isPlantation) * basalAreaTaller) * (height - 1.37)), psmeCrossValidation, 
  #                                               start = list(a1 = 2.0, b1 = 0.9, b2 = 0.01, b2p = -0.002, b2aat = -0.0001, b2aatp = 0))
  psmeDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
                                                       start = list(a1 = 2.0, a1aba = -0.006, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002))
  psmeDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
                                                      start = list(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, b1 = 0.9, b2 = 0.01, b2p = -0.002))
  psmeDiameterFromHeight$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1rh*relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, 
                                                            start = list(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002))
  # a1ac modestly preferable to b1ac and b2ac
  psmeDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ (a1 + a1ac * cos(pi/180*aspect))*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psmeCrossValidation, 
                                                   start = list(a1 = 2, a1ac = -0.1, b1 = 0.8, b2 = 0.01)) # { a1, b1, b2 } x { e, s, as, tr, tw }, b2p not significant, b1ac, b2ac also  significant
  #psmeDiameterFromHeight$ruarkPhysioB1 = fit_gsl_nls("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1ac * cos(pi/180*aspect)) * exp((b2) * (height - 1.37)), psmeCrossValidation, 
  #                                                   start = list(a1 = 2, b1 = 0.8, b1ac = -0.015, b2 = 0.01))
  #psmeDiameterFromHeight$ruarkPhysioB2 = fit_gsl_nls("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2ac * cos(pi/180*aspect)) * (height - 1.37)), psmeCrossValidation, 
  #                                                   start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2ac = -0.001))
  # a1rh, b1rh, b2rh not entirely separable but b2rh preferred
  psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psmeCrossValidation, 
                                                  start = list(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)) # a1rh, b1rh also significant, b2rh-b4p not mutually significant, exp(b2 * rh), b2rhe not significant, b2rhb singular gradient, 10x10 RMSE 14.9, AIC 5287
  #psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rhr / relativeHeight) * (height - 1.37)), psmeCrossValidation, 
  #                                                start = list(a1 = 2, b1 = 0.8, b2 = 0.008, b2rhr = -0.006)) # 10x10 RMSE 14.8, AIC 5279
  #psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rh * if_else(relativeHeight < 1.2, 0, relativeHeight - 1.2)) * (height - 1.37)), psmeCrossValidation, 
  #                                                start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2rh = 0.003)) # 10x10 RMSE 14.5, AIC 5300
  #psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight^2) * (height - 1.37)), psmeCrossValidation, 
  #                                                start = list(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.0001)) # 10x10 RMSE 14.8, AIC 5282
  #psmeDiameterFromHeight$ruarkRelHtA1 = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psmeCrossValidation, 
  #                                                  start = list(a1 = 2, a1rh = 0.3, b1 = 0.8, b2 = 0.008))
  #psmeDiameterFromHeight$ruarkRelHtB1 = fit_gsl_nls("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2) * (height - 1.37)), psmeCrossValidation, 
  #                                                  start = list(a1 = 2, b1 = 0.8, b1rh = 0.03, b2 = 0.008))
  psmeDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psmeCrossValidation, 
                                                        start = list(a1 = 2, a1ac = -0.09, b1 = 0.8, b2 = 0.009, b2rh = 0.003))
  psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1), 0.9999)), psmeCrossValidation, 
                                               start = list(a1 = 0.000001, a2 = 0.0003, b1 = 1.2, Ha = 130))
  psmeDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^b4, psmeCrossValidation, 
                                                    start = list(a1 = 3, b1 = 0.7, b2 = 0.1, b3 = 0.03, b4 = 0.3)) # a1, a1p, b1, b1p, b2, b2p, b3, b3p, b4, b4p not significant, some a1-b3 evaporation, 10x10 RMSE 15.1, AIC 5291
  psmeDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
                                                       start = list(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)) # b1p not significant, 10x10 RMSE 15.1, AIC 5242, 2x500 max iterations
  psmeDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, 
                                                           start = list(a1 = 2.0, b1 = 0.7, b2 = 0.1, b2p = -0.007, b2aba = -0.0002)) # a1p, a1aat+r, b1aat+r, b2abap, b2aat+r not significant, 10x10 RMSE 14.6, AIC 5136, 2x500 max iterations
  #psmeDiameterFromHeight$sibbesenReplaceAbatA1aba = fit_gsl_nls("Sibbesen replace ABA", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                              start = list(a1 = 2.1, a1aba = -0.006, b1 = 0.7, b2 = 0.1, b2p = -0.005)) # a1p not significant, 10x10 RMSE 14.6, AIC 5241
  #psmeDiameterFromHeight$sibbesenReplaceAbatB1aba = fit_gsl_nls("Sibbesen replace ABA", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                              start = list(a1 = 1.6, a1p = 0.3, b1 = 0.7, b1aba = -0.001, b2 = 0.1, b2p = 0.02)) # 10x10 RMSE 14.9, AIC 5276
  #psmeDiameterFromHeight$sibbesenReplaceAbatA1bat = fit_gsl_nls("Sibbesen replace ABAT", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                              start = list(a1 = 2.0, a1aat = -0.01, b1 = 0.7, b2 = 0.08, b2p = -0.005)) # a1p not significant, 10x10 RMSE 14.7, AIC 5264
  #psmeDiameterFromHeight$sibbesenReplaceAbatB1bat = fit_gsl_nls("Sibbesen replace ABAT", dbh ~ (a1)*(height - 1.37)^((b1 + b1aat * basalAreaTaller)*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                              start = list(a1 = 2.0, b1 = 0.8, b1aat = -0.001, b2 = 0.07, b2p = 0.006)) # a1p not significant, 10x10 RMSE 14.8, AIC 5180
  #psmeDiameterFromHeight$sibbesenReplaceAbatB2bat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2aat * basalAreaTaller)), psmeCrossValidation, 
  #                                                              start = list(a1 = 1.9, b1 = 0.7, b2 = 0.07, b2p = -0.005, b2aat = -0.0005)) # a1p not significant, 10x10 RMSE 15.0, AIC 5268
  #psmeDiameterFromHeight$sibbesenReplaceAbatB2abaA1bat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, 
  #                                                                   start = list(a1 = 2.0, a1aat = -0.01, b1 = 0.7, b2 = 0.09, b2p = -0.007, b2aba = -0.0002)) # a1p not significant, 10x10 RMSE 14.7, AIC 5274
  psmeDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2aba * bootstrapStandBasalAreaPerHectare + b2ac * cos(3.14159/180 * aspect))), psmeCrossValidation, 
                                                                 start = list(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2aba = -0.0002, b2ac = -0.003))
  psmeDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, 
                                                                start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.006, b2rh = 0.01, b2aba = -0.0002)) # 2x500 max iterations
  psmeDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare + b2ac * cos(3.14159/180 * aspect))), psmeCrossValidation, 
                                                                      start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.006, b2rh = 0.01, b2aba = -0.0002, b2ac = -0.003)) # 2x500 max iterations
  # a1ac and b2ac similar, b2ac probably slightly preferred
  psmeDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2ac * cos(pi/180*aspect))), psmeCrossValidation, 
                                                             start = list(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2ac = -0.003)) # { a1, b1, b2 } x { e, s, as, tr, tw } not significant, a1ac also significant
  #psmeDiameterFromHeight$sibbesenReplacePhysioA1 = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation + a1ac * cos(pi/180*aspect))*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                             start = list(a1 = 1.6, a1p = 0.4, a1ac = -0.07, b1 = 0.7, b2 = 0.1, b2p = -0.02)) # { a1, b1, b2 } x { e, s, as, tr, tw } not significant, a1ac also significant
  # a1rh, b1rh, b2rh not clearly separable
  #psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight)), psmeCrossValidation, 
  #                                                          start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01)) # a1p not significant, a1rh, b1rh also significant, 10x10 RMSE 14.8, AIC 5243
  #psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * if_else(relativeHeight < 1.2, 0, relativeHeight - 1.2))), psmeCrossValidation, 
  #                                                          start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01)) # 10x10 RMSE 14.8, AIC 5267
  #psmeDiameterFromHeight$sibbesenReplaceRelHtA1 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                            start = list(a1 = 2, a1rh = 0.4, b1 = 0.7, b2 = 0.07, b2p = -0.006))
  #psmeDiameterFromHeight$sibbesenReplaceRelHtB1 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1 + b1rh * relativeHeight)*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                            start = list(a1 = 2, b1 = 0.8, b1rh = 0.03, b2 = 0.06, b2p = -0.006))
  psmeDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2ac * cos(pi/180*aspect))), psmeCrossValidation, 
                                                                  start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01, b2ac = -0.003)) # a1p not significant
  # (a1 + a1p*isPlantation)^(b1*topHeight^b2) numerically tractable but dramatically less accurate than power, Ruark, and Sibbesen replace
  #psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight)), psmeCrossValidation, 
  #                                                          start = list(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01)) # 10x10 MAE 14.8, AIC 5254
  #psmeDiameterFromHeight$sibbesenReplaceRelHtA1 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1p * isPlantation)*(relativeHeight * (height - 1.37))^((b1)*((height - 1.37))^(b2)), psmeCrossValidation, 
  #                                                            start = list(a1 = 5.2, a1p = -0.7, b1 = 0.4, b2 = 0.1)) # 10x10 MAE 15.6, AIC 5346
  #psmeDiameterFromHeight$sibbesenReplaceRelHtB1 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(relativeHeight * (height - 1.37))^((b1 + b1p * isPlantation)*((height - 1.37))^(b2)), psmeCrossValidation, 
  #                                                            start = list(a1 = 4.7, b1 = 0.5, b1p = -0.2, b2 = 0.12)) # 10x10 MAE 15.7, AIC 5423
  #psmeDiameterFromHeight$sibbesenReplaceRelHtB2 = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(relativeHeight * (height - 1.37))^((b1)*((height - 1.37))^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                            start = list(a1 = 4.7, b1 = 0.4, b2 = 0.13, b2p = -0.01)) # 10x10 MAE 15.8, AIC 344
  #psmeDiameterFromHeight$sibbesenReplaceRelHtA1p = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(relativeHeight * (height - 1.37))^(b2)), psmeCrossValidation, 
  #                                                             start = list(a1 = 2.1, a1p = -0.1, b1 = 0.8, b2 = 0.06)) # 10x10 RMSE 14.7, AIC 5233
  psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1 + b1p * isPlantation)*(relativeHeight * (height - 1.37))^(b2)), psmeCrossValidation, 
                                                            start = list(a1 = 2.0, b1 = 0.8, b1p = -0.01, b2 = 0.05)) # a1p-b1p, b1p-b2p not mutually significant, 10x10 RMSE 14.6, AIC 5312
  #psmeDiameterFromHeight$sibbesenReplaceRelHtB2p = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(relativeHeight * (height - 1.37))^(b2 + b2p * isPlantation)), psmeCrossValidation, 
  #                                                             start = list(a1 = 2.0, b1 = 0.8, b2 = 0.06, b2p = -0.004)) # 10x10 RMSE 14.8, AIC 5313
  #psmeDiameterFromHeight$sibbesenReplaceRelHtA1p = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1p * isPlantation)*(relativeHeight*(height - 1.37))^((b1)*(height - 1.37)^(b2)), psmeCrossValidation, 
  #                                                             start = list(a1 = 5.2, a1p = -0.6, b1 = 0.4, b2 = 0.1)) # b1p, b2p not significant, 10x10 RMSE 15.8, AIC 5439
  psmeDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, psmeCrossValidation, 
                                               start = list(a1 = -140, b1 = 0.01, b2 = 0.9)) # a1p, b1p, b2p not significant
  #print(to_parameter_confidence_intervals(psmeDiameterFromHeight$weibull), n = 40)
  #to_fixed_coeffficients(psmeDiameterFromHeight$weibull)
  #psmeDiameterFromHeight$power$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  # all GAM forms NaN with Γ link
  psmeDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.0, AIC 5233
  #psmeDiameterFromHeight$gamRuark = fit_gam("REML GAM", dbh ~ I(1.92785 * (height - 1.37)^(0.84484) * exp((0.01328 - 0.00111 * isPlantation) * (height - 1.37))) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), speciesData = psmeCrossValidation, threads = 2) # 10x10 RMSE 15.1, AIC 5297
  psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 9), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.7, AIC 5263
  #psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 12), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.7, AIC 5287
  #psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 11), speciesData = psmeCrossValidation, threads = 4, distinct = FALSE) # min k = 11, 10x10 RMSE 14.9, AIC 5266
  psmeDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bs = "ts", k = 4), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.4, AIC 5243, relative height marginally significant, plantation not significant
  
  if (psmeOptions$fitPhysioGams)
  {
    psmeDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 17), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.5, AIC 5260
    psmeDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bootstrapStandBasalAreaPerHectare, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), speciesData = psmeCrossValidation, threads = 4) # min k = 16, 10x10 RMSE 14.4, AIC 5173
    #psmeDiameterFromHeight$gamAbaPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, basalAreaTaller, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), speciesData = psmeCrossValidation, threads = 4) # min k = 16, 10x10 RMSE 14.6, AIC 5295
    #psmeDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, basalAreaTaller, bootstrapStandBasalAreaPerHectare, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 57), speciesData = psmeCrossValidation, threads = 4, distinct = FALSE) # min k 57 > 40
    # plausible but modest improvement with any physiographic predictor but all predictor min k = 331 > ~150 -> drop elevation + cos(aspect) on highest RMSE min k = 57 > ~42 -> drop slope
    #psmeDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, elevation, bs = "ts", by = as.factor(isPlantation), k = 4), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.4, AIC 5269
    #psmeDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, slope, bs = "ts", by = as.factor(isPlantation), k = 9), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.0, AIC 5309
    #psmeDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 8), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.8, AIC 5275
    #psmeDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 8), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.8, AIC 5259
    #psmeDiameterFromHeight$gamRoughness = fit_gam("REML GAM roughness", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 8), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.2, AIC 5308
    #psmeDiameterFromHeight$gamWetness = fit_gam("REML GAM wetness", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 9), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.1, AIC 5305
    #psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.8, AIC 5341, k = 16 minimum only weakly supported
    #psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, slope, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 11), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 15.3, AIC 5376
    #psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, slope, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 11), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.9, AIC 5291, weak support for significance
    #psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 11), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.9, AIC 5375, weak support for significance
    psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 8), speciesData = psmeCrossValidation, threads = 4)
    psmeDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 11), speciesData = psmeCrossValidation, threads = 4) # 10x10 RMSE 14.7, AIC 5335
    #lapply(psmeDiameterFromHeight$gamRelHtPhysio$fit, k.check)
    #lapply(psmeDiameterFromHeight$gamRelHtPhysio$fit, summary)
    #psmeDiameterFromHeight$gamRelHtPhysio$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
    
    saveRDS(list(gamAbatPhysio = psmeDiameterFromHeight$gamAbatPhysio,
                 gamAbatPhysioRelHt = psmeDiameterFromHeight$gamAbatPhysioRelHt,
                 gamPhysio = psmeDiameterFromHeight$gamPhysio,
                 gamRelHtPhysio = psmeDiameterFromHeight$gamRelHtPhysio), 
            paste0("trees/height-diameter/data/PSME DBH physio GAMs ", get_cross_validation_data_suffix(), ".Rds"))
  } else {
    psmePrimaryDbhGams = readRDS(paste0("trees/height-diameter/data/PSME DBH physio GAMs ", get_cross_validation_data_suffix(), ".Rds"))
    psmeDiameterFromHeight$gamAbatPhysio = psmePrimaryDbhGams$gamAbatPhysio
    psmeDiameterFromHeight$gamAbatPhysioRelHt = psmePrimaryDbhGams$gamAbatPhysioRelHt
    psmeDiameterFromHeight$gamPhysio = psmePrimaryDbhGams$gamPhysio
    psmeDiameterFromHeight$gamRelHtPhysio = psmePrimaryDbhGams$gamRelHtPhysio
    rm(psmePrimaryDbhGams)
  }
  
  psmeDiameterFromHeight$randomForest = fit_ranger("random forest", c("dbh", "height"),
                                                   psmeCrossValidation, mtry = 1, minNodeSize = 97, sampleFraction = 0.307) # from setup.R tuneRanger @ 500 trees
  psmeDiameterFromHeight$randomForestAbatRelHtPhysio = fit_ranger("random forest RelHt ABA+T physio", c("dbh", "height", "topHeight", "basalAreaTaller", "bootstrapStandBasalAreaPerHectare", "elevation", "terrainRoughness", "slope", "aspect", "topographicWetnessFD8f"), # from setup.R VSURF
                                                                  psmeCrossValidation, mtry = 3, minNodeSize = 4, sampleFraction = 0.670) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(psmeDiameterFromHeight, paste0("trees/height-diameter/data/PSME DBH ", get_cross_validation_data_suffix(), ".Rds"))
}

if (psmeOptions$fitDbhMixed)
{
  psmeDiameterFromHeightMixed = list(linear = fit_lme("linear", fixedFormula = dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)), psmeCrossValidation, randomFormula = ~1|stand/plot)) # 10x100 false convergence
  psmeDiameterFromHeightMixed$linearAbat = fit_lme("linear ABA+T", fixedFormula = dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)) + I(1/(1 + basalAreaTaller)), psmeCrossValidation, randomFormula = ~1|stand/plot)
  psmeDiameterFromHeightMixed$linearAbatRelHt = fit_lme("linear RelHt ABA+T", fixedFormula = dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)) + relativeHeight + I(1/(1 + basalAreaTaller)), psmeCrossValidation, randomFormula = ~1|stand/plot)
  psmeDiameterFromHeightMixed$linearRelHt = fit_lme("linear RelHt", fixedFormula = dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)), psmeCrossValidation, randomFormula = ~1|stand/plot) # a1rhp not significant
  psmeDiameterFromHeightMixed$parabolic = fit_lme("parabolic", fixedFormula = dbh ~ I(height - 1.37) + I(isPlantation*(height - 1.37)) + I((height - 1.37)^2), psmeCrossValidation, randomFormula = ~1|stand/plot) # I(isPlantation*(height - 1.37)^2) not significant
  #lapply(psmeDiameterFromHeightMixed$linearRelHt$fit, summary)

  #psmeDiameterFromHeightMixed$chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation) * (exp((b1)*(height - 1.37)) - 1)^(b2 + b2ra), psmeCrossValidation, # a1ra, b2ra max iterations, b1ra singularity in backsolve
  #                                                      fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot, # |stand, |uniquePlotID max iterations
  #                                                      start = list(fixed = c(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8))))
  psmeDiameterFromHeightMixed$chapmanReplace = create_fit_statistics(name = "Chapman-Richards replace", fitting = "nlme")
  #psmeDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1p * isPlantation + a1ra)*(exp((b1 + b1aba * bootstrapStandBasalAreaPerHectare + b1aat * basalAreaTaller)*(height - 1.37)) - 1)^(b2), psmeCrossValidation, # a1ra, b1ra singularity in backsolve, b2ra max iterations
  #                                                          fixedFormula = a1 + a1p + b1 + b1aba + b1aat + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 50, a1p = -3, b1 = 0.03, b1aba = -0.0001, b1aat = -0.0001, b2 = 0.9)))
  psmeDiameterFromHeightMixed$chapmanReplaceAbat = create_fit_statistics(name = "Chapman-Richards replace ABA+T", fitting = "nlme")
  psmeDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1)*(exp((b1 + b1rh * relativeHeight)*(height - 1.37)) - 1)^(b2), psmeCrossValidation,
                                                             fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 45, b1 = 0.03, b1rh = 0.002, b2 = 0.8)), distinct = FALSE)
  #psmeDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation + a1ra)*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), psmeCrossValidation, # a1ra job step halving, b1ra singular precision, b2ra step halving
  #                                                       fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = -10, a1p = -1, b1 = 0.03, b2 = 0.3)))
  psmeDiameterFromHeightMixed$chapmanRichards = create_fit_statistics(name = "Chapman-Richards inverse", fitting = "nlme")
  #psmeDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1ra)*log(1 - pmin(((b1 + (b1aba + b1abap * isPlantation) * bootstrapStandBasalAreaPerHectare)*(height - 1.37))^b2, 0.9999)), psmeCrossValidation, # a1ra false convergence, b1ra, b2ra singularity in backsolve
  #                                                           fixedFormula = a1 + b1 + b1aba + b1abap + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = -100, b1 = 0.01, b1aba = -0.00002, b1abap = -0.00001, b2 = 1.0)))
  psmeDiameterFromHeightMixed$chapmanRichardsAbat = create_fit_statistics(name = "Chapman-Richards inverse ABA+T", fitting = "nlme")
  psmeDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1ra + a1p * isPlantation)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), psmeCrossValidation, # fixed effects not significant
                                                               fixedFormula = a1 + a1p + b1 + b2 + b2tw ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8, b2tw = 0)), distinct = FALSE)
  psmeDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1ra)*log(1 - pmin((b1*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), psmeCrossValidation, # a1ra 2x50 job, b1ra step halving, b2ra singularity in backsolve
                                                              fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = -20, b1 = 0.018, b2 = 0.5, b2rh = -0.3))) # 2x250, 5x100, 5x200, 10x100 step halving
  psmeDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1) * (1 + a1rm) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), psmeCrossValidation, # b1p not significant, lower height step halving, lower height power max iterations, upper height power singular precision, both heights max iterations, 10x10 RMSE 16.4, AIC 5404, 10x50, 10x100 step halving
                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 150, a2 = 90, b1 = 0.9))) # 2x250, 5x200, 10x100 step halving
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceB1 = fit_nlme("Michaelis-Menten replace", dbh ~ (a1) * (height - 1.37)^(b1 + b1ra) / (a2 - (height - 1.37)^(b1 + b1ra)), psmeCrossValidation, # a1ra step halving, a2ra singularity in backsolve, b1p not significant, 10x10 RMSE 26.0, AIC 5539 vs 14.9, 5224 fixed effect
  #                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 40, a2 = 22, b1 = 0.7)))
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceB1 = fit_nlme("Michaelis-Menten replace", dbh ~ (a1) * (height - 1.37)^(b1 * (1 + b1rm)) / (a2 - (height - 1.37)^(b1 * (1 + b1rm))), psmeCrossValidation, # b1p not significant, 10x10 RMSE 18.2, AIC 5416
  #                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 80, a2 = 45, b1 = 0.8)))
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceAbat = fit_nlme("Michaelis-Menten replace ABA+T", dbh ~ (a1 + a1ra) * (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, # a1ra, a1rm, b1ra step halving, a2ra singularity in backsolve
  #                                                                  fixedFormula = a1 + a2 + b1 + b1p + b1aba ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 120, a2 = 70, b1 = 0.9, b1aba = -0.0005, b1p = -0.01)))
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceAbat = fit_nlme("Michaelis-Menten replace ABA+T", dbh ~ (a1) * (1 + a1rm) * (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, # a1rm, b1rm step halving
  #                                                                  fixedFormula = a1 + a2 + b1 + b1p + b1aba ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 120, a2 = 70, b1 = 0.9, b1aba = -0.0005, b1p = -0.01)))
  psmeDiameterFromHeightMixed$michaelisMentenReplaceAbat = create_fit_statistics(name = "Michaelis-Menten replace ABA+T", fitting = "nlme")
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceAbatRelHt = fit_nlme("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight) * (1 + a1rm) * (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, # 
  #                                                                     fixedFormula = a1 + a1rh + a2 + b1 + b1aba ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                                     start = list(fixed = c(a1 = 110, a1rh = 10, a2 = 60, b1 = 0.9, b1aba = -0.0005)), distinct = FALSE) # b1aba not significant
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceAbatRelHt = fit_nlme("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight + a1ra) * (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation + b1aba * bootstrapStandBasalAreaPerHectare)), psmeCrossValidation, # a1ra maximum iterations, a2ra, b1ra step halving
  #                                                                       fixedFormula = a1 + a1rh + a2 + b1 + b1p + b1aba ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                       start = list(fixed = c(a1 = 110, a1rh = 10, a2 = 60, b1 = 0.9, b1aba = -0.0005, b1p = -0.01)))
  psmeDiameterFromHeightMixed$michaelisMentenReplaceAbatRelHt = create_fit_statistics(name = "Michaelis-Menten replace RelHt ABA+T", fitting = "nlme")
  #psmeDiameterFromHeightMixed$michaelisMentenReplaceRelHt = fit_nlme("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight + a1ra) * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psmeCrossValidation, # a1ra, a2ra singularity in backsolve
  #                                                                   fixedFormula = a1 + a1rh + a2 + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                   start = list(fixed = c(a1 = 130, a1rh = 20, a2 = 70, b1 = 0.9, b1p = -0.01)))
  psmeDiameterFromHeightMixed$michaelisMentenReplaceRelHt = fit_nlme("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight) * (1 + a1rm) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), psmeCrossValidation, # b1p not significant, b1rm step halving, 10x10 RMSE 15.8, AIC 5422
                                                                     fixedFormula = a1 + a1rh + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 130, a1rh = 25, a2 = 70, b1 = 0.9)))
  psmeDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1ra + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psmeCrossValidation, # a2ra also tractable but false convergence with both a1ra and a2ra, 10x10 RMSE 17.6, AIC 5518 vs 15.6, 5514 fixed effect
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006)))
  #psmeDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * (1 + a1rm) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), psmeCrossValidation, # a2rm max iterations, 10x10 RMSE 18.0, AIC 5606
  #                                               fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 4.3, a1p = -0.6, a2 = -0.1, a2p = -0.01)))
  #psmeDiameterFromHeightMixed$naslundA2 = fit_nlme("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation + a2ra) * sqrt(height - 1.37)), psmeCrossValidation, # 30x higher RMSE than a1ra
  #                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a2ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 4.0, a1p = -0.6, a2 = -0.1, a2p = -0.006)))
  psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^((b1 + b1p*isPlantation)*(1 + b1rm)), psmeCrossValidation, # 10x10 RMSE 15.1, AIC 5346
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 0.6, a1p = 0.5, b1 = 1.3, b1p = -0.2)))
  #psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1p*isPlantation)*(height - 1.37)^(b1 + b1p*isPlantation + b1ra), psmeCrossValidation, # a1ra also tractable, 10x10 RMSE 15.3, AIC 5321 vs 15.1, 5311 fixed
  #                                             fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 0.7, a1p = 0.4, b1 = 1.2, b1p = -0.2)))
  #psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1 + a1p*isPlantation)*(1 + a1ra)*(height - 1.37)^(b1), psmeCrossValidation, # a1p-b1p not mutually significant, step halving or singularity in backsolve with (a1)*((1 + a1ra)*(height - 1.37))^(b1) regardess of a1p or b1p, 10x10 RMSE 16.4, AIC 5424
  #                                             fixedFormula = a1 + a1p + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 1.2, a1p = -0.1, b1 = 1.1)))
  #psmeDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1)*(1 + a1ra)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, # 10x10 RMSE 16.4, AIC 5352
  #                                             fixedFormula = a1 + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 1.1, b1 = 1.2, b1p = -0.02)))
  #psmeDiameterFromHeightMixed$powerA1 = fit_nlme("power", dbh ~ (a1 + a1p*isPlantation + a1ra)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, # false convergence but possibly lower mean or median error
  #                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.7, a1p = 0.4, b1 = 1.2, b1p = -0.2)))
  psmeDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^((b1 + b1p*isPlantation)*(1 + b1rm)), psmeCrossValidation, # 10x10 RMSE 15.0, AIC 5266, maybe lower AIC than a1ra
                                                   fixedFormula = a1 + a1aat + b1 + b1p ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 1.2, a1aat = -0.01, b1 = 1.1, b1p = -0.2)))
  #psmeDiameterFromHeightMixed$powerAbatA1 = fit_nlme("power ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller + a1ra)*(height - 1.37)^(b1 + b1p*isPlantation), psmeCrossValidation, # 10x10 RMSE 14.8, AIC 5338
  #                                                   fixedFormula = a1 + a1aat + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1.2, a1aat = -0.01, b1 = 1.1, b1p = -0.2)))
  #psmeDiameterFromHeightMixed$powerAbatB1 = fit_nlme("power ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1 + b1p*isPlantation + b1ra), psmeCrossValidation, # 10x10 RMSE 15.1, AIC 5403
  #                                                   fixedFormula = a1 + a1aat + b1 + b1p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1.2, a1aat = -0.01, b1 = 1.1, b1p = -0.2)))
  psmeDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1ra + a1p * isPlantation + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation), psmeCrossValidation, # fixed effects not significant
                                                     fixedFormula = a1 + a1p + a1ac + b1 + b1p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 0.7, a1p = 0.6, b1 = 1.2, b1p = -0.2, a1ac = 0)), distinct = FALSE)
  psmeDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1p * isPlantation + b1rh * relativeHeight)*(1 + b1rm)), psmeCrossValidation, # 10x10 RMSE 15.0, AIC 5305
                                                    fixedFormula = a1 + a1p + b1 + b1p + b1rh ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)))
  #psmeDiameterFromHeightMixed$powerRelHtA1 = fit_nlme("power RelHt", dbh ~ (a1 + a1p * isPlantation + a1ra)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight), psmeCrossValidation, # false convergence
  #                                                    fixedFormula = a1 + a1p + b1 + b1p + b1rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)))
  #psmeDiameterFromHeightMixed$powerRelHtB1 = fit_nlme("power RelHt", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^(b1 + b1p * isPlantation + b1rh * relativeHeight + b1ra), psmeCrossValidation, # 10x10 RMSE 15.2, AIC 5271
  #                                                    fixedFormula = a1 + a1p + b1 + b1p + b1rh ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 1.1, a1p = 0.4, b1 = 1.1, b1p = -0.1, b1rh = 0.05)))
  psmeDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1 + b1ra) * exp((b2) * (height - 1.37)), psmeCrossValidation, # b2p not significant, 10x10 RMSE 15.0, AIC 5288
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkA1 = fit_nlme("Ruark", dbh ~ (a1)* (1 + a1rm) * (height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psmeCrossValidation, # b2rm leading minor, 10x10 RMSE 16.6, AIC 5429
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.8, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkB1 = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^((b1)*(1 + b1rm)) * exp((b2) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 15.7, AIC 5290
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.1, b1 = 1.0, b2 = 0)))
  #psmeDiameterFromHeightMixed$ruarkA1 = fit_nlme("Ruark", dbh ~ (a1 + a1ra) * (height - 1.37)^(b1) * exp((b2) * (height - 1.37)), psmeCrossValidation, # false convergence
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkB2 = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1) * exp((b2 + b2ra) * (height - 1.37)), psmeCrossValidation, # not discriminable from b1ra
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.01)))
  psmeDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.6, AIC 5267
                                                   fixedFormula = a1 + a1aba + b1 + b2 + b2p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 2.0, a1aba = -0.006, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  #psmeDiameterFromHeightMixed$ruarkAbatA1 = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.7, AIC 5279
  #                                                   fixedFormula = a1 + a1aba + b1 + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 2.0, a1aba = -0.006, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  #psmeDiameterFromHeightMixed$ruarkAbatB2 = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + b2ra) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.7, AIC 5177
  #                                                   fixedFormula = a1 + a1aba + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 2.0, a1aba = -0.006, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  psmeDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, # a1ra false convergence, 10x10 RMSE 14.3, AIC 5240
                                                         fixedFormula = a1 + a1aba + a1ac + b1 + b2 + b2p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  #psmeDiameterFromHeightMixed$ruarkAbatPhysioB2 = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + b2ra) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.4, AIC 5254
  #                                                       fixedFormula = a1 + a1aba + a1ac + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  psmeDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.4, AIC, 5262
                                                        fixedFormula = a1 + a1rh + a1aba + b1 + b2 + b2p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  #psmeDiameterFromHeightMixed$ruarkAbatRelHtA1 = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1 + a1ra + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, # job false convergence, 10x10 RMSE 14.7, AIC 5281
  #                                                        fixedFormula = a1 + a1rh + a1aba + b1 + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  #psmeDiameterFromHeightMixed$ruarkAbatRelHtB2 = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + b2ra) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.7, AIC 5265
  #                                                      fixedFormula = a1 + a1rh + a1aba + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  psmeDiameterFromHeightMixed$ruarkAbatRelHtPhysio = fit_nlme("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psmeCrossValidation, # a1ra singular precision, 10x10 RMSE 14.3, AIC 5273
                                                              fixedFormula = a1 + a1aba + a1ac + a1rh + b1 + b2 + b2p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  #psmeDiameterFromHeightMixed$ruarkAbatRelHtPhysioB2 = fit_nlme("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation + b2ra) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.6, AIC 5176
  #                                                            fixedFormula = a1 + a1aba + a1ac + a1rh + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002)))
  psmeDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2ra) * (height - 1.37)), psmeCrossValidation, # a1ra singular precision, b1ra also tractable
                                                     fixedFormula = a1 + a1ac + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 2, a1ac = -0.1, b1 = 0.8, b2 = 0.01)))
  #psmeDiameterFromHeightMixed$ruarkPhysioB1 = fit_nlme("Ruark physio", dbh ~ (a1 + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1ra) * exp((b2) * (height - 1.37)), psmeCrossValidation, # near identical to b2ra
  #                                                     fixedFormula = a1 + a1ac + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 2, a1ac = -0.1, b1 = 0.8, b2 = 0.01)))
  
  psmeDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight + b2ra) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 14.8, AIC 5197
                                                    fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtA1 = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1ra)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psmeCrossValidation, # false convergence, 10x10 RMSE 14.9, AIC 5341
  #                                                      fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtB1 = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psmeCrossValidation, # 10x10 RMSE 15.1, AIC 5282
  #                                                    fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  psmeDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight + b2ra) * (height - 1.37)), psmeCrossValidation, # a1ra, b1ra also tractable
                                                          fixedFormula = a1 + a1ac + b1 + b2 + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 2, a1ac = -0.09, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtPhysioA1 = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1ac * cos(3.14159/180 * aspect) + a1ra)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psmeCrossValidation, # false convergence
  #                                                          fixedFormula = a1 + a1ac + b1 + b2 + b2rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 2, a1ac = -0.09, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$ruarkRelHtPhysioB1 = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), psmeCrossValidation, # less accurate than b2ra
  #                                                          fixedFormula = a1 + a1ac + b1 + b2 + b2rh ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 2, a1ac = -0.09, b1 = 0.8, b2 = 0.008, b2rh = 0.003)))
  #psmeDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/((Ha + Hara)^b1 - 1.37^b1), 0.9999)), psmeCrossValidation, # a1ra, a1rm, a2ra, a2rm, b1ra, b1rm, Hara, Harm step halving
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Hara ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.000001, a2 = 0.0003, b1 = 1.2, Ha = 130)))
  psmeDiameterFromHeightMixed$schnute = create_fit_statistics(name = "Schnute inverse", fitting = "nlme")
  #psmeDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1 * (height - 1.37)^(b1 + b1ra)*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^(b4), psmeCrossValidation, # a1ra, b2ra, b3ra, b3rm singularity in backsolve, a1rm, b1ra, b1rm, b3ra, b3rm, b4rm step halving
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve
  #                                                    start = list(fixed = c(a1 = 3, b1 = 0.7, b2 = 0.1, b3 = 0.03, b4 = 0.3)))
  psmeDiameterFromHeightMixed$sharmaParton = create_fit_statistics(name = "modified Sharma-Parton", fitting = "nlme")
  psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2ra)), psmeCrossValidation, # a1ra step halving, b2rm max iterations, 10x10 RMSE 15.0, AIC 5310
                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02))) # 2x250 step halving
  #psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), psmeCrossValidation, # a1p, b2p not significant, 10x10 RMSE 16.5, AIC 5403
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1.9, b1 = 0.7, b2 = 0.09)))
  #psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1)*(height - 1.37)^((b1)*(1 + b1rm)*(height - 1.37)^(b2)), psmeCrossValidation, # a1p, b2p not significant, 10x10 RMSE 15.3, AIC 5329
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1.4, b1 = 0.9, b2 = 0.05)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceB1 = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1ra)*(height - 1.37)^(b2 + b2p * isPlantation)), psmeCrossValidation, # a1ra step halving, b1ra also tractable
  #                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2aba * bootstrapStandBasalAreaPerHectare + b2ra)), psmeCrossValidation, # a1ra step halving, b1ra false convergence
                                                             fixedFormula = a1 + b1 + b2 + b2p + b2aba ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 2.0, b1 = 0.7, b2 = 0.1, b2p = -0.007, b2aba = -0.0002)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1p * isPlantation + a1ra)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2aba * bootstrapStandBasalAreaPerHectare + b2ac * cos(3.14159/180 * aspect))), psmeCrossValidation,
  #                                                                 fixedFormula = a1 + a1p + b1 + b2 + b2p + b2aba + b2ac + ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2aba = -0.0002, b2ac = -0.003)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = create_fit_statistics(name = "Sibbesen replace ABA+T physio", fitting = "nlme") # parsing failure: object `a1` not found
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace RelHt ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare + b2ra)), psmeCrossValidation, # a1ra singular precision, b1ra job max iterations
                                                                  fixedFormula = a1 + b1 + b2 + b2p + b2rh + b2aba ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.006, b2rh = 0.01, b2aba = -0.0002)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysio = fit_nlme("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare + b2ac * cos(3.14159/180 * aspect) + b2ra)), psmeCrossValidation, # a1ra singular precision, job false convergence, 10x0 RMSE 14.6, AIC 5267
                                                                        fixedFormula = a1 + b1 + b2 + b2p + b2rh + b2aba + b2ac ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.006, b2rh = 0.01, b2aba = -0.0002, b2ac = -0.003))) # 2x250, 10x10 observation step halving
  #psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysioB1 = fit_nlme("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^((b1 + b1ra)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare + b2ac * cos(3.14159/180 * aspect))), psmeCrossValidation, # job false convergence, 10X10 RMSE 16.0, AIC 5295
  #                                                                        fixedFormula = a1 + b1 + b2 + b2p + b2rh + b2aba + b2ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                                        start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.006, b2rh = 0.01, b2aba = -0.0002, b2ac = -0.003)))
  psmeDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2ac * cos(3.14159/180 * aspect) + b2ra)), psmeCrossValidation, # a1ra false convergence, false convergence but more accurate than b1ra
                                                               fixedFormula = a1 + a1p + b1 + b2 + b2p + b2ac ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 1.2, a1p = 0.3, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2ac = -0.003)))
  #psmeDiameterFromHeightMixed$sibbesenReplacePhysioB1 = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1ra)*(height - 1.37)^(b2 + b2p * isPlantation + b2ac * cos(3.14159/180 * aspect))), psmeCrossValidation,
  #                                                               fixedFormula = a1 + a1p + b1 + b2 + b2p + b2ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = 1.2, a1p = 0.3, b1 = 0.7, b2 = 0.1, b2p = -0.02, b2ac = -0.003)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2ra)), psmeCrossValidation, # a1ra step halving, b1ra singularity in backsolve, 10x10 RMSE 1.50, AIC 5775
  #                                                            fixedFormula = a1 + b1 + b2 + b2p + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1 + b1ra)*(relativeHeight * (height - 1.37))^(b2)), psmeCrossValidation, # a1ra job step halving, 10x10 RMSE 14.9, AIC 5267
                                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2.0, b1 = 0.7, b2 = 0.08))) # 10x100 step halving
  #psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(relativeHeight * (height - 1.37))^(b2 + b2ra)), psmeCrossValidation, # b1p not significant, 10x10 RMSE 15.3, AIC 5331
  #                                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 1.6, b1 = 0.9, b2 = 0.04)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^((b1)*(relativeHeight * (height - 1.37))^(b2)), psmeCrossValidation, # 10x10 RMSE 15.7, AIC 5425
  #                                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 1.9, b1 = 0.8, b2 = 0.05)))
  #psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1*(1 + b1rm))*(relativeHeight * (height - 1.37))^(b2)), psmeCrossValidation, # b2rm max iterations, 10x10 RMSE 15.1, AIC 5283
  #                                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 1.5, b1 = 1.0, b2 = 0.04)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2rh * relativeHeight + b2ac * cos(3.14159/180 * aspect) + b2ra)), psmeCrossValidation, # a1ra, b1ra false convergence
                                                                    fixedFormula = a1 + b1 + b2 + b2p + b2rh + b2ac ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 2, b1 = 0.8, b2 = 0.05, b2p = -0.002, b2rh = 0.01, b2ac = -0.003)))
  #psmeDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2 + b2ra), psmeCrossValidation, # a1ra step halving, b1ra singularity in backsolve, b2ra job max iterations
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = -140, b1 = 0.01, b2 = 0.9)))
  psmeDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2*(1 + b2rm)), psmeCrossValidation, # a1rm singularity in backsolve
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = -135, b1 = 0.01, b2 = 0.97)))
  #to_fixed_coeffficients(psmeDiameterFromHeightMixed$weibull)
  #print(to_parameter_confidence_intervals(psmeDiameterFromHeightMixed$weibull), n = 40)
  
  psmeDiameterFromHeightMixed$gam = fit_gam("REML GAM", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation, distinct = FALSE) # collapses to linear, plantation not significant
  psmeDiameterFromHeightMixed$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 8), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation)
  psmeDiameterFromHeightMixed$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation) # min k = 16
  psmeDiameterFromHeightMixed$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bootstrapStandBasalAreaPerHectare, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation) # min k = 16
  psmeDiameterFromHeightMixed$gamPhysio = fit_gam("REML GAM physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, cos(3.14159/180 * aspect), bs = "ts", k = 8), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation) # plantation not significant
  psmeDiameterFromHeightMixed$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation) # plantation not significant
  psmeDiameterFromHeightMixed$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, relativeHeight, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), speciesData = psmeCrossValidation) # k = 16 minimum -> only weakly significant
  #lapply(psmeDiameterFromHeightMixed$gamRelHtPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(psmeDiameterFromHeightMixed$gamRelHtPhysio$fit, function(fit) { return(summary(fit$gam)) })
  #psmeDiameterFromHeightMixed$gam$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(psmeDiameterFromHeightMixed, paste0("trees/height-diameter/data/PSME DBH mixed ", get_cross_validation_data_suffix(), ".Rds"))
}

if (psmeOptions$recalcResultCollections)
{
  # assemble parameter and results tibbles using incremental load and reduce to fit in memory
  # height primary + nlrob + gsl_nls + dbh primary + nlrob = 120 GB DDR at 10x10 cross validation with models dropped
  # Therefore, batch formation height and diameter results.
  if (exists("psmeHeightFromDiameter") == FALSE) 
  { 
    psmeHeightFromDiameter = readRDS(paste0("trees/height-diameter/data/PSME height ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("psmeHeightFromDiameterMixed") == FALSE) 
  { 
    psmeHeightFromDiameterMixed = readRDS(paste0("trees/height-diameter/data/PSME height mixed ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("psmeDiameterFromHeight") == FALSE)
  { 
    psmeDiameterFromHeight = readRDS(paste0("trees/height-diameter/data/PSME DBH ", get_cross_validation_data_suffix(), ".Rds")) 
    #psmeDiameterFromHeight = list()
  }
  if (exists("psmeDiameterFromHeightMixed") == FALSE) 
  { 
    psmeDiameterFromHeightMixed = readRDS(paste0("trees/height-diameter/data/PSME DBH mixed ", get_cross_validation_data_suffix(), ".Rds")) 
    #psmeDiameterFromHeightMixed = list()
  }
  
  psmeCoefficients = unnest_and_bind("PSME", psmeHeightFromDiameter, psmeHeightFromDiameterMixed, psmeDiameterFromHeight, psmeDiameterFromHeightMixed, unnest_coefficients)
  psmeFitStats = unnest_and_bind("PSME", psmeHeightFromDiameter, psmeHeightFromDiameterMixed, psmeDiameterFromHeight, psmeDiameterFromHeightMixed, unnest_validation_statistics)
  psmePredictors = unnest_and_bind("PSME", psmeHeightFromDiameter, psmeHeightFromDiameterMixed, psmeDiameterFromHeight, psmeDiameterFromHeightMixed, unnest_predictors) %>%
    mutate(across(where(is.logical), ~replace_na(.x, FALSE)))
  
  check_plot_fit_stats(psmeFitStats)
  saveRDS(psmeCoefficients, paste0("trees/height-diameter/data/PSME coefficients ", get_cross_validation_data_suffix(), ".Rds"))
  saveRDS(psmeFitStats, paste0("trees/height-diameter/data/PSME fit stats ", get_cross_validation_data_suffix(), ".Rds"))
  saveRDS(psmePredictors, paste0("trees/height-diameter/data/PSME predictors ", get_cross_validation_data_suffix(), ".Rds"))
}


## preferred forms identified (results.R, Figure 8 and trees.R)
if (psmeOptions$recalcPreferredModels)
{
  # model selections
  # responseVariable species       cross  maeName                  rmseName                 aicName          nseName                  isGam isRandomForest maeFitting rmseFitting aicFitting nseFitting
  # height           Douglas-fir   10x50  random forest            Hossfeld IV              Michaelis-Menten Hossfeld IV                  0              1 ranger     gsl_nls     gsl_nls    gsl_nls  
  #                                2x250  Ratkowsky                Ratkowsky                Hossfeld IV      Ratkowsky                    0              0 gsl_nls    gsl_nls     gsl_nls    gsl_nls 
  # DBH              Douglas-fir   10x50  Chapman-Richards inverse Chapman-Richards replace Ruark            Chapman-Richards replace     0              0 gsl_nls    gsl_nls     gsl_nls    gsl_nls
  #                                2x250  Michaelis-Menten replace Ruark                    Ruark            Ruark                        0              0 gsl_nls    gsl_nls     gsl_nls    gsl_nls     
  #psmeHeightFromDiameterPreferred = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1 * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), psme2022, start = list(a1 = 60, b1 = -0.02, b2 = 1.1, b2p = 0.1)))
  #psmeHeightFromDiameterPreferred$gam = fit_gam("REML GAM", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), data = psme2022, threads = 4)
  psmeHeightFromDiameterPreferred = list(hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + (a2) * dbh^(b1 + b1p * isPlantation)), psme2022, 
                                                                start = list(a1 = 73, a2 = 150, b1 = -1.3, b1p = 0.03)))
  psmeHeightFromDiameterPreferred$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + a2p * isPlantation + dbh^b1), psme2022, 
                                                                start = list(a1 = 75, a2 = 140, a2p = 20, b1 = 1.2))
  psmeHeightFromDiameterPreferred$randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"), psme2022, 
                                                            mtry = 2, minNodeSize = 48, sampleFraction = 0.230)
  psmeHeightFromDiameterPreferred$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(dbh + b2)), psme2022, 
                                                          start = list(a1 = 73, b1 = -46, b1p = -3, b2 = 10))
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4e * elevation), psme2022, 
                                                                      start = list(a1 = 25, a1ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = -0.04, b4 = 1.0, b4e = 0.0002))
  #psmeHeightFromDiameterPreferred$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, psme2022, start = list(a1 = -140, b1 = 0.01, b2 = 0.9))
  saveRDS(psmeHeightFromDiameterPreferred, "trees/height-diameter/data/PSME preferred height models.Rds")
  
  psmeDiameterFromHeightPreferred = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1 + a1p * isPlantation)*(exp(b1*(height - 1.37)) - 1)^b2, psme2022, 
                                                                      start = list(a1 = 40, a1p = -1, b1 = 0.03, b2 = 0.8)))
  psmeDiameterFromHeightPreferred$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), psme2022, 
                                                                start = list(a1 = -100, b1 = 0.012, b2 = 0.9))
  psmeDiameterFromHeightPreferred$gam = fit_gam("REML GAM", dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 6), data = psme2022, threads = 4)
  psmeDiameterFromHeightPreferred$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, 
                                                                       start = list(a1 = 130, a2 = 70, b1 = 0.9, b1p = -0.008))
  psmeDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                      start = list(a1 = 2, b1 = 0.8, b2 = 0.01, b2p = -0.001))
  psmeDiameterFromHeightPreferred$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1rh*relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare + a1ac * cos(3.14159/180 * aspect))*(height - 1.37)^(b1) * exp((b2 + b2p * isPlantation) * (height - 1.37)), psme2022, 
                                                                     start = list(a1 = 2.0, a1aba = -0.006, a1rh = 0.2, a1ac = -0.1, b1 = 0.9, b2 = 0.01, b2p = -0.002))
  saveRDS(psmeDiameterFromHeightPreferred, "trees/height-diameter/data/PSME preferred diameter models.Rds")
}


## error residual distributions
if (htDiaOptions$includeInvestigatory)
{
  psmeDbhMichaelisMentenReplace = gsl_nls(dbh ~ a1 * (height - 1.37)^(b1 + b1p * isPlantation) / (a2 - (height - 1.37)^(b1 + b1p * isPlantation)), psme2022, start = list(a1 = 130, a2 = 70, b1 = 0.9, b1p = -0.008))
  ggplot() +
    geom_point(aes(x = relativeHeight, y = (dbh - residuals(psmeDbhMichaelisMentenReplace)) / dbh), psme2022, alpha = 0.1, shape = 16) +
    geom_smooth(aes(x = relativeHeight, y = (dbh - residuals(psmeDbhMichaelisMentenReplace)) / dbh), psme2022, formula = y ~ x, method = "loess") + # GAM linear, loess linear/parabolic
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "relative height", y = "scale error")
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  psmeHeightGam = gam(height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + 
                               s(dbh, bs = "ts", by = as.factor(isPlantation), k = 8) + 
                               s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 3) + # linear
                               s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 8) + # broadly linear-ish
                               s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 7) + # curved, effects to 30-50 m²/ha
                               s(elevation, bs = "ts", k = 3) + # linear
                               s(slope, bs = "ts", k = 3) + # linear
                               s(aspect, bs = "ts", k = 3) +  # parabolic, discontinuous
                               # s(terrainRoughness, bs = "ts", k = 6) + # not significant
                               s(topographicWetnessFD8f, bs = "ts", k = 6), 
                      data = psme2022, method = "REML", select = TRUE, weights = dbhWeight, threads = 8)
  k.check(psmeHeightGam)
  summary(psmeHeightGam) # 0.954 R²
  plot_gam_smooths(psmeHeightGam)

  psmeHeightGamSimplified = gam(height ~ I((2.4201 - 1.0816 * isPlantation) * dbh^(0.6516 + 0.1398 * isPlantation)) + 
                                           s(dbh, bs = "ts", by = as.factor(isPlantation), k = 8) + 
                                           s(relativeDiameter, bs = "ts", k = 3) + 
                                           s(standBasalAreaPerHectare, bs = "ts", k = 8) +
                                           s(basalAreaLarger, bs = "ts", k = 7) + 
                                           s(elevation, bs = "ts", k = 3) +
                                           s(slope, bs = "ts", k = 3) + 
                                           s(aspect, bs = "ts", k = 3) +
                                           #s(terrainRoughness, bs = "ts", k = 6) + # not significant
                                           s(topographicWetnessFD8f, bs = "ts", k = 3), 
                                  data = psme2022, method = "REML", select = TRUE, weights = dbhWeight, threads = 8)
  k.check(psmeHeightGamSimplified)
  summary(psmeHeightGamSimplified) # 0.953 R²
  plot_gam_smooths(psmeHeightGamSimplified)
  
  psmeDbhGam = gam(dbh ~ I((0.668 + 0.607 * isPlantation) * height^(1.283 - 0.200 * isPlantation)) + 
                         s(height, bs = "ts", by = as.factor(isPlantation), k = 6) + # piecewise linear, upwards from ~1.2
                         s(relativeHeight, bs = "ts", k = 5) + # plantation not significant
                         s(bootstrapStandBasalAreaPerHectare, bs = "ts", k = 6) + # plantation not significant, concave down, arguably piecewise linear with break ~65 m²/ha
                         s(basalAreaTaller, bs = "ts", k = 6) + # plantation not significant, linear-ish overall, arguably piecewise linear with breakpoint ~15 m²/ha
                         #s(elevation, bs = "ts", k = 4) + # not significant
                         # s(slope, bs = "ts", k = 3) + # not significant
                         s(aspect, bs = "ts", k = 4), # +4 mm southsouthwest, -4 mm north, not quite continuous
                         #s(terrainRoughness, bs = "ts", k = 3) + # not significant
                         # s(topographicWetnessFD8f, bs = "ts", k = 3), # not significant 
                   data = psme2022, method = "REML", select = TRUE, weights = heightWeight, threads = 8)
  k.check(psmeDbhGam)
  summary(psmeDbhGam)
  plot_gam_smooths(psmeDbhGam)
}
