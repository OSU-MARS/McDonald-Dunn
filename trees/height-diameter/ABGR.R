# load libraries, functions, and abgr2022 from setup.R
# abgrOptions from ABGR job.R

## grand fir height regressions
if (abgrOptions$fitHeightPrimary)
{
  # linear regressions
  abgrHeightFromDiameter = list(linear = fit_lm("linear", height ~ dbh, abgr2022)) # plantation mostly not significant
  abgrHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ dbh + I(dbh^2), abgr2022) # plantation, plantation^2 only sometimes significant

  # nonlinear regressions
  # reciprocal BAL
  abgrHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), abgr2022, 
                                                       start = list(a1 = 64, a1p = -6, b1 = -0.02, b2 = 1.3, b2p = -0.1)) # 10x10 MAE 1.86, RMSE 609
  abgrHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
                                                          start = list(a1 = 70, b1 = -0.02, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6)) # a1ba+r, a1bal+rp, b1ba+r, b2ba+r, b2bap, b2balp+rp not significant, b1bal+r NaN-inf, 10x10 MAE 1.69, AIC 575
  #abgrHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BAL", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
  #                                                        start = list(a1 = 70, b1 = -0.02, b2 = 1.2, b2bal = 1.2, b2balr = 2)) # 10x10 MAE 1.77, AIC 584
  abgrHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
                                                                start = list(a1 = 65, b1 = -0.02, b1ac = -0.001, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6), distinct = FALSE) # b1ac not significant
  abgrHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
                                                                      start = list(a1 = 65, b1 = -0.02, b1rd = 0.001, b1ac = -0.001, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6), distinct = FALSE) # b1ac, b1rd not significant
  abgrHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
                                                                start = list(a1 = 80, b1 = -0.02, b1rd = 0.001, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6), distinct = FALSE) # { a1rd, b1rd, b2rd }-b2bal+r not mutually significant
  # linear BAL: a1ba, a1bal, and b1bal interchangeable
  #abgrHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + (a1bal + a1balp * isPlantation) * basalAreaLarger) * (1 - exp((b1)*dbh))^(b2), abgr2022, 
  #                                                        start = list(a1 = 70, a1bal = 0.2, a1balp = 0.2, b1 = -0.02, b2 = 1.3)) # a1p, a1bap, b1ba, b2p, b2ba, b2bal not significant, a1ba-{ a1bal, b1bal } not mutually significant, 10x10 MAE 1.77, AIC 598
  #abgrHeightFromDiameter$chapmanRichardsBalA1ba = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp((b1)*dbh))^(b2), abgr2022, 
  #                                                            start = list(a1 = 55, a1ba = 0.1, b1 = -0.02, b2 = 1.3)) # 10x1) MAE 1.77, AIC 601
  #abgrHeightFromDiameter$chapmanRichardsBalB1bal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + (b1bal + b1balp * isPlantation) * basalAreaLarger)*dbh))^(b2), abgr2022, 
  #                                                             start = list(a1 = 75, b1 = -0.02, b1bal = -0.00005, b1balp = -0.00005, b2 = 1.3)) # 10x10 MAE 1.83, AIC 609
  #abgrHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter) * (1 - exp((b1)*dbh))^(b2), abgr2022, 
  #                                                              start = list(a1 = 80, a1bal = 0.2, a1rd = -4, b1 = -0.02, b2 = 1.2), distinct = FALSE) # a1balp, a1rdp, b1rd not significant, a1bal-a1rd not reliably mutually significant
  #abgrHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare + a1bal * basalAreaLarger) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^b2, abgr2022, 
  #                                                              start = list(a1 = 65, a1ba = -0.1, a1bal = 0.3, b1 = -0.01, b1ac = 0, b2 = 1.3), distinct = FALSE) # b2p, b1bal, b1ac, b2bal not significant
  #abgrHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp((b1 + b1rd * relativeDiameter + b1ac * cos(3.14159/180 * aspect))*dbh))^b2, abgr2022, 
  #                                                                    start = list(a1 = 61, a1bal = 0.2, b1 = -0.02, b1rd = 0.001, b1ac = -0.001, b2 = 1.3), distinct = FALSE) # a1ba, b1ac not significant
  abgrHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), abgr2022, 
                                                             start = list(a1 = 60, a1p = -6, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2p = -0.1)) # { a1, b1, b2 } x { e, s, as, tr, tw }, { a1, b2 } x ac not significant
  # no clear differentiation between a1rd and b1rd
  abgrHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^b2, abgr2022, 
                                                             start = list(a1 = 75, b1 = -0.02, b1rd = 0.001, b2 = 1.3)) # a1p, a1rdp, b1rdp, b2p, b2rd, b2rdp not significant
  #abgrHeightFromDiameter$chapmanRichardsRelDbhA1 = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp((b1)*dbh))^b2, abgr2022, 
  #                                                             start = list(a1 = 85, a1rd = -6, b1 = -0.02, b2 = 1.3))
  abgrHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect) + b1rd * relativeDiameter)*dbh))^(b2 + b2p * isPlantation), abgr2022, 
                                                                   start = list(a1 = 66, a1p = -6, b1 = -0.02, b1ac = 0, b1rd = 0.001, b2 = 1.3, b2p = -0.1), distinct = FALSE) # b1ac not significant
  abgrHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1 * dbh / (1 + dbh)^b1, abgr2022, 
                                              start = list(a1 = 0.6, b1 = 0.03), distinct = FALSE) # a1p, b1, b1p not significant -> collapses to linear
  abgrHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), abgr2022, 
                                                start = list(a1 = 80, a2 = 200, b1 = -1.3)) # a1p, a2p, b1p not significant
  abgrHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1*dbh^b2), abgr2022, 
                                            start = list(a1 = 200, b1 = -8, b2 = -0.2)) # a1, a1p, b1p, b2p not significant
  abgrHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), abgr2022, 
                                                       start = list(a1 = 85, a2 = 200, b1 = 1.3)) # a1p, a2p, b1p not significant
  abgrHeightFromDiameter$michaelisMentenBal = fit_gsl_nls("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1)), abgr2022, 
                                                          start = list(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2)) # a1bal+r, a2balp, b1bal+r not significant, 10x10 MAE 1.69, AIC 574
  #abgrHeightFromDiameter$michaelisMentenBaA1 = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare)*dbh^(b1) / (a2 + dbh^(b1)), abgr2022, 
  #                                                         start = list(a1 = 85, a1ba = 0, a2 = 200, b1 = 1.3)) # a2ba not significant, 10x10 MAE 1.91, AIC 609
  #abgrHeightFromDiameter$michaelisMentenBaB1 = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + b1ba * standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1ba * standBasalAreaPerHectare)), abgr2022, 
  #                                                         start = list(a1 = 85, a2 = 200, b1 = 1.3, b1ba = 0)) # 10x10 MAE 1.82, AIC 605
  #abgrHeightFromDiameter$michaelisMentenBarA1 = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1 + a1ba / standBasalAreaPerHectare)*dbh^(b1) / (a2 + dbh^(b1)), abgr2022, 
  #                                                          start = list(a1 = 95, a1ba = -200, a2 = 200, b1 = 1.3)) # 10x10 MAE 1.88, AIC 618
  #abgrHeightFromDiameter$michaelisMentenBarA2 = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1) / (a2 + a2ba / standBasalAreaPerHectare + dbh^(b1)), abgr2022, 
  #                                                          start = list(a1 = 85, a2 = 190, a2ba = 700, b1 = 1.2)) # 10x10 MAE 1.85, AIC 605
  #abgrHeightFromDiameter$michaelisMentenBarB1 = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1ba / standBasalAreaPerHectare)), abgr2022, 
  #                                                          start = list(a1 = 80, a2 = 200, b1 = 1.3, b1ba = -2)) # 10x10 MAE 1.82, AIC 589
  #abgrHeightFromDiameter$michaelisMentenBal = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1ba / standBasalAreaPerHectare)), abgr2022, 
  #                                                        start = list(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2, b1ba = 0), distinct = FALSE) # b1ba not significant
  abgrHeightFromDiameter$michaelisMentenBalRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1rd * relativeDiameter)), abgr2022, 
                                                                start = list(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2, b1rd = 0), distinct = FALSE) # b1rd not significant
  abgrHeightFromDiameter$michaelisMentenRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + dbh^(b1 + b1rd * relativeDiameter)), abgr2022, 
                                                             start = list(a1 = 110, a2 = 300, b1 = 1.3, b1rd = -0.03)) # b1rdp not significant, 10x10 MAE 1.80, AIC 583
  #abgrHeightFromDiameter$michaelisMentenRelDbhA1$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #abgrHeightFromDiameter$michaelisMentenRelDbhA2$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #abgrHeightFromDiameter$michaelisMentenRelDbhB1$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  abgrHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, abgr2022, 
                                             start = list(a1 = 0.6, b1 = 1.0)) # a1p, b1p not significant
  abgrHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3), abgr2022, 
                                              start = list(a1 = 0.01, a2 = 1.1, a3 = 3)) # a1p, a2p, a3p not significant
  abgrHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2)), abgr2022, 
                                                 start = list(a1 = 70, b1 = -40, b2 = 7)) # a1p, b1p, b2p not significant
  abgrHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), abgr2022, 
                                                 start = list(Ha = 45, d = 1.3, kU = 0.02)) # dp, Hap, kUp not significant
  abgrHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh))^(b4 + b4p * isPlantation), abgr2022, 
                                                    start = list(a1 = 25, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4p = -0.1)) # a1p, b1p, b2p, b3, b3p not significant, 10x10 RMSE 3.10, AIC 622
  #abgrHeightFromDiameter$sharmaPartonB4 = fit_gsl_nls("Sharma-Parton b₄", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3 * dbh^(b4))), abgr2022, 
  #                                                    start = list(a1 = 30, b1 = 0.3, b2 = -0.008, b3 = 0, b4 = 1.3)) # b3, b4p not significant, 10x10 RMSE 3.22, AIC 586
  abgrHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh))^(b4 + b4p * isPlantation + b4bal / (b4balr + basalAreaLarger)), abgr2022, 
                                                       start = list(a1 = 35, b1 = 0.2, b2 = -0.02, b3 = -0.04, b4 = 1.3, b4p = -0.1, b4bal = 0, b4balr = 2), distinct = FALSE) # a1bal+r, b1bal+r, b2bal+r, b3, b3bal+r, b4bal+r not significant with or without b4p, 5x200, 5x100, 10x100 NaN-inf
  #abgrHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh))^(b4 + b4p * isPlantation + b4bal * basalAreaLarger), abgr2022, 
  #                                                     start = list(a1 = 35, b1 = 0.2, b2 = -0.02, b3 = -0.04, b4 = 1.3, b4p = -0.1, b4bal = 0)) # b1bal, b2bal, b3bal, b4bal not significant
  abgrHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4bal / (b4balr + basalAreaLarger) + b4ac * cos(pi/18 *aspect)), abgr2022, 
                                                             start = list(a1 = 40, b1 = 0.2, b2 = -0.02, b3 = 0, b4 = 1.2, b4ac = 0, b4bal = 0, b4balr = 2), distinct = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw }, b4bal+r not significant, job NaN-inf
  abgrHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4bal / (b4balr + basalAreaLarger) + b4ac * cos(pi/18 *aspect)), abgr2022, 
                                                                   start = list(a1 = 50, b1 = 0.1, b2 = -0.02, b3 = 0.1, b3rd = -0.04, b4 = 1.3, b4ac = 0, b4bal = 0, b4balr = 2), distinct = FALSE) # b1, b4ac, b4bal+r not significant
  abgrHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter) * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), abgr2022, 
                                                             start = list(a1 = 50, b1 = 0.1, b2 = -0.02, b3 = 0.1, b3rd = -0.03, b4 = 1.3, b4bal = 0, b4balr = 3), distinct = FALSE) # b4bal+r not significant, 2x50 NaN-Inf
  # a1rd and b4rd quite similar with conflicting MAE and AIC signals
  abgrHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a1ac * cos(pi/180*aspect))*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), abgr2022, 
                                                          start = list(a1 = 25, a1ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3)) # { a1, b1, b2, b3, b4 } x { e, s, as, tr, tw }, { b1, b2, b3 } x ac, b3, b4p not significant
  #abgrHeightFromDiameter$sharmaPartonPhysioB4 = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4ac * cos(pi/180*aspect)), abgr2022, 
  #                                                          start = list(a1 = 25, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4ac = 0.8))
  # b1 and b3 interchangeable
  abgrHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh))^b4, abgr2022, 
                                                          start = list(a1 = 25, b1 = 0.3, b1rd = -0.03, b2 = -0.02, b3 = 0.1, b4 = 1.3)) # a1rd, b2rd NaN-inf, b4p, b4rd not significant
  #abgrHeightFromDiameter$sharmaPartonRelDbhB3 = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter) * dbh))^b4, abgr2022, 
  #                                                          start = list(a1 = 30, b1 = 0.3, b2 = -0.02, b3 = 0.1, b3rd = -0.03, b4 = 1.3))
  abgrHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a1ac * cos(pi/180*aspect))*topHeight^(b1 + b1rd * relativeDiameter) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4), abgr2022, 
                                                                start = list(a1 = 25, a1ac = 0, b1 = 0.3, b1rd = -0.03, b2 = -0.02, b3 = 0.1, b4 = 1.3), distinct = FALSE) # a1ac not significant
  abgrHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*standTreesPerHectare^b3*dbh))^(b4 + b4p * isPlantation), abgr2022, 
                                                   start = list(a1 = 33, b1 = 0.15, b2 = -0.03, b3 = -0.06, b4 = 1.2, b4p = -0.1)) # a1p, b1p, b2p, b3p not significant, 10x10 MAE 1.73, AIC 579
  #abgrHeightFromDiameter$sharmaZhangB4 = fit_gsl_nls("Sharma-Zhang b₄", height ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*standTreesPerHectare^b3*dbh))^(b4 + b4p * isPlantation), abgr2022, 
  #                                                   start = list(a1 = 33, b1 = 0.15, b2 = -0.03, b3 = -0.08, b4 = 1.3, b4p = -0.1)) # b3 not significant, 10x10 MAE 3.17, AIC 597
  abgrHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3) * dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, 
                                                      start = list(a1 = 20, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4ba = 0.007)) # 10x10 MAE 1.73, AIC 579
  #abgrHeightFromDiameter$sharmaZhangBaB2 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^(b3) * dbh))^(b4), abgr2022, 
  #                                                     start = list(a1 = 15, b1 = 0.4, b2 = -0.04, b2ba = 0.0003, b3 = -0.05, b4 = 1.3))
  #abgrHeightFromDiameter$sharmaZhangBaB2r = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2 / standBasalAreaPerHectare)*standTreesPerHectare^(b3) * dbh))^(b4), abgr2022, 
  #                                                      start = list(a1 = 3, b1 = 0.9, b2 = -1.4, b3 = -0.1, b4 = 1.2)) # a1ba, a1ba+r, b1ba, b1ba+r, b2ba+r, b3ba, b3ba+r, b4ba+r not significant, b2-b2ba not mutually significant, 10x10 MAE 1.85, AIC 612
  #abgrHeightFromDiameter$sharmaZhangBaB3 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3 + b3ba * standBasalAreaPerHectare) * dbh))^(b4), abgr2022, 
  #                                                     start = list(a1 = 15, b1 = 0.4, b2 = -0.03, b3 = -0.1, b3ba = -0.001, b4 = 1.3)) # b3 not significant
  #abgrHeightFromDiameter$sharmaZhangBaB4 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3) * dbh))^(b4 + b4ba / standBasalAreaPerHectare), abgr2022, 
  #                                                     start = list(a1 = 30, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4 = 1.4, b4ba = -4)) # 10x10 MAE 1.78, AOC 593
  #abgrHeightFromDiameter$sharmaZhangBalB4 = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3) * dbh))^(b4 + b4bal * basalAreaLarger), abgr2022, 
  #                                                      start = list(a1 = 25, b1 = 0.2, b2 = -0.04, b2ba = 0, b3 = -0.1, b4 = 1.0, b4bal = 0.005)) # a1ba, a1bal, b1ba, b1bal, b2bal, b3bal not significant, b4bal x { b2ba, b3ba, b4ba } not mutually significant
  #abgrHeightFromDiameter$sharmaZhangBalB4r = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), abgr2022, 
  #                                                       start = list(a1 = 33, b1 = 0.15, b2 = -0.03, b3 = -0.06, b4 = 1.2, b4bal = 0, b4balr = 1), distinct = FALSE) # a1bal+r, b1bal+r, b2bal, b3bal, b4bal+r not significant
  abgrHeightFromDiameter$sharmaZhangBalRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3 + b3rd * relativeDiameter) * dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, 
                                                            start = list(a1 = 25, b1 = 0.2, b2 = -0.03, b3 = 0, b3rd = 0, b4 = 1.1, b4ba = 0.007), distinct = FALSE) # b3rd not significant
  abgrHeightFromDiameter$sharmaZhangRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3 + b3rd * relativeDiameter)*dbh))^(b4), abgr2022, 
                                                         start = list(a1 = 65, b1 = 0.05, b2 = -0.01, b3 = 0, b3rd = -0.02, b4 = 1.2)) # a1rd, b1rd, b3, b4rd, b4p not significant, b2rd NaN-inf
  abgrHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1*dbh^b2), abgr2022, 
                                                start = list(a1 = 0.2, b1 = 1.8, b2 = -0.1)) # a1p, b1p, b2p not significant
  abgrHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1*(1 - exp(b1*dbh^b2)), abgr2022, 
                                               start = list(a1 = 58, b1 = -0.01, b2 = 1.2)) # a1p, b1p, b2p not significant, 10x10 MAE 1.83, AIC 598
  abgrHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1bal / (b1balr + basalAreaLarger))*dbh^(b2))), abgr2022, 
                                                  start = list(a1 = 70, b1 = -0.006, b1bal = 0.02, b1balr = 5, b2 = 1.2)) # a1bal+r, b2bal+r not significant, 10x10 MAE 1.73, AIC 593
  #abgrHeightFromDiameter$weibullBaA1 = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a1ba / standBasalAreaPerHectare) * (1 - exp((b1)*dbh^(b2))), abgr2022, 
  #                                                 start = list(a1 = 70, a1ba = -200, b1 = -0.006, b2 = 1.2)) # a1ba+r, b1ba, b1ba+r, b2ba+r not significant, 10x10 MAE 1.85, AIC 575
  #abgrHeightFromDiameter$weibullBaB2 = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh^(b2 + b2ba / standBasalAreaPerHectare))), abgr2022, 
  #                                                 start = list(a1 = 60, b1 = -0.008, b2 = 1.2, b2ba = -2)) # 10x10 MAE 1.80, AIC 602
  # linear BA+L: a1bal, b1ba, b2ba, b2bal conflicting indications but slightly against b2ba and maybe weakly for a1bal
  #abgrHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp((b1)*dbh^(b2))), abgr2022, 
  #                                                start = list(a1 = 58, a1bal = 0, b1 = -0.01, b2 = 1.2)) # a1ba, b1ba not significant, { a1bal, b1bal, b2bal }-b2ba not mutually significant, 10x10 MAE 1.77, AIC 595
  #abgrHeightFromDiameter$weibullBalB1bal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1bal * basalAreaLarger)*dbh^(b2))), abgr2022, 
  #                                                     start = list(a1 = 58, b1 = -0.01, b1bal = 0, b2 = 1.2))
  #abgrHeightFromDiameter$weibullBalB2ba = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh^(b2 + b2ba * standBasalAreaPerHectare))), abgr2022, 
  #                                                    start = list(a1 = 57, b1 = -0.01, b2 = 1.2, b2ba = 0.001))
  #abgrHeightFromDiameter$weibullBalB2bal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh^(b2 + b2bal * basalAreaLarger))), abgr2022, 
  #                                                     start = list(a1 = 58, b1 = -0.01, b2 = 1.2, b2bal = 0))
  abgrHeightFromDiameter$weibullBalRelDbh = fit_gsl_nls("Weibull RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp((b1 + b1bal / (b1balr + basalAreaLarger))*dbh^(b2))), abgr2022, 
                                                        start = list(a1 = 70, a1rd = 0, b1 = -0.006, b1bal = 0.02, b1balr = 5, b2 = 1.2), distinct = FALSE) # a1rd-b1bal+r not mutually significant
  abgrHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp(b1*dbh^b2)), abgr2022, 
                                                     start = list(a1 = 78, a1rd = -5, b1 = -0.01, b2 = 1.2)) # 10x10 MAE 1.81, AIC 578
  #abgrHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1)*(1 - exp((b1 + b1rd * relativeDiameter)*dbh^b2)), abgr2022, 
  #                                                   start = list(a1 = 65, b1 = -0.01, b1rd = 0.001, b2 = 1.2)) # 10x10 MAE 1.83, AIC 578
  #abgrHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1)*(1 - exp((b1)*dbh^(b2 + b1rd * relativeDiameter))), abgr2022, 
  #                                                   start = list(a1 = 70, b1 = -0.01, b1rd = -0.03, b2 = 1.3)) # 10x10 MAE 1.83, AIC 590
  #print(to_parameter_confidence_intervals(abgrHeightFromDiameter$sharmaZhang), n = 40)
  #to_fixed_coeffficients(abgrHeightFromDiameter$sibbesen)
  #abgrHeightFromDiameter$weibullBalA1$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  # GAM NaNs with Γ link
  abgrHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I(0.6539 * dbh^0.9757) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = abgr2022, threads = 2) # 10x10 MAE 1.87, AIC 610
  #abgrHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(0.6539 * dbh^0.9757) + s(dbh, standBasalAreaPerHectare, bs = "ts", k = 12), data = abgr2022, threads = 2) # 10x10 MAE 1.84, AIC 599
  abgrHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13), data = abgr2022, threads = 2) # 10x10 MAE 1.71, AIC 602
  #abgrHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(0.6539 * dbh^0.9757) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13), data = abgr2022, threads = 2, distinct = FALSE) # 10x10 MAE 1.73, AIC 637
  abgrHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 14), data = abgr2022, threads = 2) # 10x10 MAE 1.71, AIC 595
  abgrHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 15), data = abgr2022, threads = 2)
  abgrHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, relativeDiameter, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), data = abgr2022, threads = 2, distinct = FALSE) # k = 16 minimum > 0-10.6
  # indication against elevation, slope, cos(aspect), roughness, no clear signal between sin(aspect) and wetness but indication against both -> select sin(aspect) on tendency to lower AIC
  #abgrHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I(0.6539 * dbh^0.9757) + s(dbh, elevation, bs = "ts", by = as.factor(isPlantation), k = 4), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I(0.6539 * dbh^0.9757) + s(dbh, slope, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 4), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I(0.6539 * dbh^0.9757) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", k = 9), data = abgr2022, threads = 2) # 10x10 MAE 1.81, AIC 607
  #abgrHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I(0.6539 * dbh^0.9757) + s(dbh, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 8), data = abgr2022, threads = 2, distinct = FALSE) # 10x10 MAE 1.94, AIC 627
  #abgrHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I(0.6539 * dbh^0.9757) + s(dbh, terrainRoughness, bs = "ts", k = 8), data = abgr2022, threads = 2, distinct = FALSE) # 10x10 MAE 1.94, AIC 607
  #abgrHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I(0.6539 * dbh^0.9757) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 8), data = abgr2022, threads = 2) # 10x10 MAE 1.81, AIC 619
  abgrHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", k = 9), data = abgr2022, threads = 2)
  #abgrHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, sin(3.14159/180 * aspect), topographicWetnessFD8f, bs = "ts", k = 11), data = abgr2022, threads = 2, distinct = FALSE) # k = 11 minimum, collapses to linear
  abgrHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(0.6539 * dbh^0.9757) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 8), data = abgr2022, threads = 2) # 10x10 MAE 1.72, AIC 595
  abgrHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, relativeDiameter, sin(3.14159/180 * aspect), bs = "ts", k = 12), data = abgr2022, threads = 2) # 10x10 MAE 1.79, AIC 588, debatably significant
  #lapply(abgrHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(abgrHeightFromDiameter$gamRelDbhPhysio$fit, summary)
  #abgrHeightFromDiameter$gamRelDbhPhysio$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))

  abgrHeightFromDiameter$randomForest = fit_ranger("random forest", c("height", "dbh"), # 10x10 MAE 1.86
                                                   abgr2022, mtry = 1, minNodeSize = 15, sampleFraction = 0.220) # from setup.R tuneRanger @ 500 trees
  abgrHeightFromDiameter$randomForestBalPhysioRelDbh = fit_ranger("random forest RelDbh BA+L physio", c("height", "dbh", "relativeDiameter", "basalAreaLarger", "standQmd", "topHeight", "standBasalAreaPerHectare", "terrainRoughness", "slope"), # from setup.R VSURF
                                                                  abgr2022, mtry = 4, minNodeSize = 2, sampleFraction = 0.513) # from setup.R tuneRanger @ 500 trees
  
  # saving primary Rds at the end of this block ensures physio GAMs are included
  saveRDS(abgrHeightFromDiameter, paste0("trees/height-diameter/data/ABGR height ", get_cross_validation_data_suffix(), ".Rds"))
}


## grand fir height mixed regressions
if (abgrOptions$fitHeightMixed)
{
  abgrHeightFromDiameterMixed = list(linear = fit_lme("linear", fixedFormula = height ~ dbh, abgr2022, randomFormula = ~1|stand/plot))
  abgrHeightFromDiameterMixed$parabolic = fit_lme("parabolic", fixedFormula = height ~ dbh + I(dbh^2), abgr2022, randomFormula = ~1|stand/plot)

  abgrHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 + a1rm) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), abgr2022, # b2rm singularity in backsolve, 10x10 MAE 1.91, AIC 633
                                                                fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 54, a1p = -8, b1 = -0.02, b2 = 1.4, b2p = -0.1)))
  #abgrHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation + a1ra)*(1 - exp(b1*dbh))^(b2 + b2p * isPlantation), abgr2022, # b1ra, b2ra singularity in backsolve, 10x10 MAE 2.08, AIC 657
  #                                                       fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 64, a1p = -6, b1 = -0.02, b2 = 1.3, b2p = -0.1)))
  #abgrHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*(1 + b1rm)*dbh))^(b2), abgr2022, # 10x10 MAE 1.95, AIC 614
  #                                                       fixedFormula = a1 + a1p + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 64, a1p = -6, b1 = -0.02, b2 = 1.3)))
  abgrHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 + a1rm) * (1 - exp((b1)*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, # 10x10 MAE 1.72, AIC 597
                                                            fixedFormula = a1 + b1 + b2 + b2ba + b2bal + b2balr ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 70, b1 = -0.02, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6)))
  #abgrHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1)*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, # b1ra, b2ra singularity in backsolve, 10x10 MAE 2.04, AIC 597
  #                                                          fixedFormula = a1 + b1 + b2 + b2ba + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 70, b1 = -0.02, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6)))
  #abgrHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + (a1bal + a1balp * isPlantation) * basalAreaLarger + a1ra) * (1 - exp(b1*dbh))^b2, abgr2022, # a1ra 2x50 job max iterations, b1ra, b2ra singularity in backsolve
  #                                                          fixedFormula = a1 + a1bal + a1balp + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 70, a1bal = 0.2, a1balp = 0.2, b1 = -0.02, b2 = 1.3)))
  #abgrHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
  #                                                                     fixedFormula = a1 + b1 + b1ac + b2 + b2ba + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 65, b1 = -0.02, b1ac = -0.001, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6)), distinct = FALSE)
  abgrHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ra + a1ba * standBasalAreaPerHectare + a1bal * basalAreaLarger) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^b2, abgr2022,
                                                                  fixedFormula = a1 + a1ba + a1bal + b1 + b1ac + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 65, a1ba = -0.1, a1bal = 0.3, b1 = -0.01, b1ac = 0, b2 = 1.3)))
  #abgrHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1rd * relativeDiameter + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
  #                                                                     fixedFormula = a1 + a1bal + a1balp + b1 + b1rd + b1ac + b2 + b1bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                      start = list(fixed = c(a1 = 65, b1 = -0.02, b1rd = 0.001, b1ac = -0.001, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6)), distinct = FALSE)
  abgrHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1ra + a1bal * basalAreaLarger) * (1 - exp((b1 + b1rd * relativeDiameter + b1ac * cos(3.14159/180 * aspect))*dbh))^b2, abgr2022,
                                                                        fixedFormula = a1 + a1bal + b1 + b1rd + b1ac + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 61, a1bal = 0.2, b1 = -0.02, b1rd = 0.001, b1ac = -0.001, b2 = 1.3)), distinct = FALSE)
  abgrHeightFromDiameterMixed$chapmanRichardsBalRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
                                                                  fixedFormula = a1 + a1bal + b1 + b1rd + b2 + b2ba + b2bal + b2balr ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 80, b1 = -0.02, b1rd = 0.001, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6)), distinct = FALSE)
  abgrHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1p * isPlantation + a1ra) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), abgr2022, # b1ra, b2ra singularity in backsolve, 10x10 MAE 1.72, AIC 597
                                                               fixedFormula = a1 + a1p  + b1 + b1ac + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 60, a1p = -6, b1 = -0.02, b1ac = -0.001, b2 = 1.3, b2p = -0.1)))
  #abgrHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 + a1rm) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect))*dbh))^(b2 + b2p * isPlantation), abgr2022, # 10x10 MAE 1.86, AIC 606
  #                                                             fixedFormula = a1 + a1p  + b1 + b1ac + b2 + b2p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 55, a1p = -6, b1 = -0.02, b1ac = -0.002, b2 = 1.3, b2p = -0.1)))
  abgrHeightFromDiameterMixed$chapmanRichardsRelDbh = fit_nlme("Chapman-Richards RelDbh", height ~ 1.37 + (a1) * (1 + a1rm) * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^b2, abgr2022, # 10x10 MAE 1.81, AIC 599
                                                               fixedFormula = a1 + b1 + b1rd + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 65, b1 = -0.02, b1rd = 0.002, b2 = 1.3))) # 2x250 max iterations
  #abgrHeightFromDiameterMixed$chapmanRichardsRelDbh = fit_nlme("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1rd * relativeDiameter)*dbh))^b2, abgr2022, # a1ra 2x50 step halving, b1ra, b2ra singularity in backsolve, 10x01 MAE 2.65, AIC 732
  #                                                             fixedFormula = a1 + b1 + b1rd + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 75, b1 = -0.02, b1rd = 0.001, b2 = 1.3)))
  abgrHeightFromDiameterMixed$chapmanRichardsRelDbhPhysio = fit_nlme("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1ac * cos(3.14159/180 * aspect) + b1rd * relativeDiameter)*dbh))^(b2 + b2p * isPlantation), abgr2022,
                                                                     fixedFormula = a1 + a1p + b1 + b1ac + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 66, a1p = -6, b1 = -0.02, b1ac = 0, b1rd = 0.001, b2 = 1.3, b2p = -0.1)), distinct = FALSE)
  abgrHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1) * dbh / (1 + dbh)^(b1 + b1ra), abgr2022,
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.6, b1 = 0.03)), distinct = FALSE)
  abgrHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + a2 * dbh^(b1 + b1ra)), abgr2022, # a2ra max iterations, 10x10 MAE 1.85, AIC 609
                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 80, a2 = 200, b1 = -1.3)))
  #abgrHeightFromDiameterMixed$hossfeldA1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) * (1 + a1rm) / (1 + a2 * dbh^b1), abgr2022, # 10x10 MAE 1.96, AIC 620
  #                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 70, a2 = 180, b1 = -1.3)))
  #abgrHeightFromDiameterMixed$hossfeldA2 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + a2 * (1 + a2rm) * dbh^b1), abgr2022, # 10x10 MAE 2.01, AIC 611
  #                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a2rm ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 100, a2 = 210, b1 = -1.1)))
  #abgrHeightFromDiameterMixed$hossfeldB1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1) / (1 + a2 * dbh^(b1*(1 + b1rm))), abgr2022, # 10x10 MAE 1.91, AIC 588
  #                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 80, a2 = 200, b1 = -1.3)))
  #abgrHeightFromDiameterMixed$hossfeldA1 = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1ra) / (1 + a2 * dbh^(b1)), abgr2022, # b1ra more accurate
  #                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 80, a2 = 200, b1 = -1.3)))
  #abgrHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + a1*exp(b1*dbh^(b2 + b2ra)), abgr2022, # a1ra max iterations, b1ra, b2ra step halving
  #                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot, # |stand, |uniquePlotID step halving
  #                                            start = list(fixed = c(a1 = 200, b1 = -8, b2 = -0.2)))
  #abgrHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + a1*exp(b1*dbh^(b2 + b2ra)), abgr2022, # a1ra max iterations, b1ra, b2ra step halving
  #                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot, # |stand, |uniquePlotID step halving
  #                                            start = list(fixed = c(a1 = 200, b1 = -8, b2 = -0.2)))
  abgrHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + a1*exp(b1*dbh^(b2 * (1 + b2rm))), abgr2022, # b1rm step halving, 10x10 MAE 1.93, AIC 632
                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot, 
                                              start = list(fixed = c(a1 = 200, b1 = -8, b2 = -0.2))) # 5x200, 10x100 step halving
  #abgrHeightFromDiameterMixed$korfA1 = fit_nlme("Korf", height ~ 1.37 + a1 * (1 + a1rm) * exp(b1 * dbh^(b2)), abgr2022, # 10x10 MAE 1.97, AIC 639
  #                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot, 
  #                                              start = list(fixed = c(a1 = 200, b1 = -8, b2 = -0.2)))
  abgrHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*(1 + a1rm)*dbh^b1 / (a2 + dbh^b1), abgr2022, # b1rm max iterations, 10x10 MAE 1.85, AIC 627
                                                         fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 70, a2 = 180, b1 = 1.3)))
  #abgrHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^(b1 + b1ra) / (a2 + dbh^(b1 + b1ra)), abgr2022, # a1ra 2x50 job max iterations, a2ra max iterations, 10x10 MAE 1.91, AIC 611
  #                                                       fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 85, a2 = 200, b1 = 1.3)))
  #abgrHeightFromDiameterMixed$michaelisMentenA2 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1)*dbh^b1 / (a2 + (1 + a2rd)*dbh^b1), abgr2022, # 10x10 MAE 2.15, AIC 649
  #                                                         fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a2rd ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 50, a2 = 170, b1 = 1.4)))
  #abgrHeightFromDiameterMixed$michaelisMentenA1 = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1ra)*dbh^(b1) / (a2 + dbh^(b1)), abgr2022, # b1ra more accurate
  #                                                         fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 85, a2 = 200, b1 = 1.3)))
  abgrHeightFromDiameterMixed$michaelisMentenBal = fit_nlme("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*(1 + a1rm)*dbh^(b1) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1)), abgr2022, # 10x10 MAE 1.68, AIC 591
                                                            fixedFormula = a1 + a2 + a2bal + a2balr + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2)), distinct = FALSE) # a2bal, a2balr not significant
  #abgrHeightFromDiameterMixed$michaelisMentenBalA2 = fit_nlme("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1) / (a2 + a2ra + a2bal / (a2balr + basalAreaLarger) + dbh^(b1)), abgr2022, # a1ra max iterations, a2ra job step halving
  #                                                          fixedFormula = a1 + a2 + a2bal + a2balr + b1 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2)))
  #abgrHeightFromDiameterMixed$michaelisMentenBalB1 = fit_nlme("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ra) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1ra)), abgr2022, # 10x10 MAE 1.80, AIC 591
  #                                                          fixedFormula = a1 + a2 + a2bal + a2balr + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2)))
  abgrHeightFromDiameterMixed$michaelisMentenBalRelDbh = fit_nlme("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + a2bal / (a2balr + basalAreaLarger) + dbh^(b1 + b1rd * relativeDiameter)), abgr2022,
                                                                  fixedFormula = a1 + a2 + a2bal + a2balr + b1 + b1rd ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 100, a2 = 230, a2bal = 600, a2balr = 3, b1 = 1.2, b1rd = 0)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$michaelisMentenRelDbh = fit_nlme("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a1ra)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + dbh^(b1 + b1rd * relativeDiameter)), abgr2022, # a1ra job step halving + singularity in backsolve, a2ra max iterations, b1ra step halving
  #                                                             fixedFormula = a1 + a2 + b1 + b1rd ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 110, a2 = 300, b1 = 1.3, b1rd = -0.03)))
  abgrHeightFromDiameterMixed$michaelisMentenRelDbh = fit_nlme("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*(1 + a1ra)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + dbh^(b1 + b1rd * relativeDiameter)), abgr2022, # 10x10 MAE 1.80, AIC 604
                                                               fixedFormula = a1 + a2 + b1 + b1rd ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 90, a2 = 240, b1 = 1.3, b1rd = -0.03)))
  #abgrHeightFromDiameterMixed$michaelisMentenRelDbh = fit_nlme("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^((b1 + b1rd * relativeDiameter)*(1 + a1ra)) / (a2 + dbh^((b1 + b1rd * relativeDiameter)*(1 + a1ra))), abgr2022, # 10x10 MAE, AIC 594
  #                                                             fixedFormula = a1 + a2 + b1 + b1rd ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 200, a2 = 400, b1 = 1.3, b1rd = -0.05)))
  abgrHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1)*dbh^(b1*(1 + b1ra)), abgr2022, # 10x10 MAE 2.20, AIC 629
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 0.6, b1 = 1.0)))
  #abgrHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1 + a1ra)*dbh^(b1), abgr2022, # 10x10 MAE 2.29, AIC 649
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.6, b1 = 1.0)))
  #abgrHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1)*(1 + a1ra)*dbh^(b1), abgr2022, # 10x10 MAE 2.19, AIC 648
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.6, b1 = 1.0)))
  #abgrHeightFromDiameterMixed$powerB1 = fit_nlme("power", height ~ 1.37 + (a1)*dbh^(b1 + b1ra), abgr2022, # same MAE but a1ra slightly lower AIC
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.6, b1 = 1.0)))
  abgrHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1)*dbh^2 + (a2 + a2ra)*dbh + a3), abgr2022, # 10x10 MAE 1.93, AIC 622
                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 0.01, a2 = 1.1, a3 = 3))) # 2x250 max iterations
  #abgrHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1)*dbh^2 + (a2)*(1 + a2rm)*dbh + a3), abgr2022, # a1rm singularity in backsolve, a3rm singular precision, 10x10 MAE 2.12, AIC 649
  #                                              fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2rm ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 0.01, a2 = 1.4, a3 = 31)))
  #abgrHeightFromDiameterMixed$prodanA1 = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1ra)*dbh^2 + (a2)*dbh + a3), abgr2022, # a2ra slighly more accurate
  #                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.01, a2 = 1.1, a3 = 3)))
  #abgrHeightFromDiameterMixed$prodanA3 = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1)*dbh^2 + (a2)*dbh + a3 + a3ra), abgr2022, # job max iterations
  #                                                fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a3ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.01, a2 = 1.1, a3 = 3)))
  abgrHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1)/(dbh + b2 + b2ra)), abgr2022, # 10x10 MAE 1.84, AIC 590
                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 70, b1 = -40, b2 = 7)))
  #abgrHeightFromDiameterMixed$ratkowskyA1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*(1 + a1rm)*exp((b1)/(dbh + b2)), abgr2022, # 10x10 MAE 2.00, AIC 605
  #                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 63, b1 = -33, b2 = 6)))
  #abgrHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1*(1 + b1rm))/(dbh + b2)), abgr2022, # 10x10 MAE 1.94, AIC 622
  #                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 72, b1 = -46, b2 = 8)))
  #abgrHeightFromDiameterMixed$ratkowskyA1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1 + a1ra)*exp((b1)/(dbh + b2)), abgr2022, # b1ra and b2ra more accurate
  #                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 70, b1 = -40, b2 = 7)))
  #abgrHeightFromDiameterMixed$ratkowskyB1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1 + b1ra)/(dbh + b2)), abgr2022, # b2ra modestly more accurate
  #                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 70, b1 = -40, b2 = 7)))
  abgrHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha + Hara) * (1 + ((1.37/(Ha + Hara))^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), abgr2022, # conflicting MAE and AIC signals between Har and kUr, 10x10 MAE 2.44, AIC 558
                                                   fixedFormula = Ha + d + kU ~ 1, randomFormula = Hara ~ 1|stand/plot,
                                                   start = list(fixed = c(Ha = 45, d = 1.3, kU = 0.02))) # 2x500, 5x100, 5x200 step halving
  #abgrHeightFromDiameterMixed$richardsWkU = fit_nlme("unified Richards", height ~ 1.37 + (Ha) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * (1 + kUrm) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), abgr2022, # Harm, drm singular precision, 10x10 MAE 2.40, AIC 1160
  #                                                   fixedFormula = Ha + d + kU ~ 1, randomFormula = kUrm ~ 1|stand/plot,
  #                                                   start = list(fixed = c(Ha = 45, d = 1.3, kU = 0.02)))
  #abgrHeightFromDiameterMixed$richardsWkU = fit_nlme("unified Richards", height ~ 1.37 + (Ha) * (1 + ((1.37/(Ha))^(1 - d) - 1) * exp((-(kU + kUr) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), abgr2022,
  #                                                   fixedFormula = Ha + d + kU ~ 1, randomFormula = kUr ~ 1|stand/plot,
  #                                                   start = list(fixed = c(Ha = 45, d = 1.3, kU = 0.02)))
  #abgrHeightFromDiameterMixed$richardsWd = fit_nlme("unified Richards", height ~ 1.37 + (Ha) * (1 + ((1.37/(Ha))^(1 - (d + dr)) - 1) * exp((-kU * dbh)/(d + dr)^((d + dr)/(1 - (d + dr)))))^(1/(1 - (d + dr))), abgr2022, # job step halving
  #                                                  fixedFormula = Ha + d + kU ~ 1, randomFormula = dr ~ 1|stand/plot,
  #                                                  start = list(fixed = c(Ha = 45, d = 1.3, kU = 0.02)))
  abgrHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2 * (1 + b2rm) * (standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022, # b1rm, b4rm singularity in backsolve, b3 not significant, b3rm leading minor, 10x10 MAE 1.77, AIC 610
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 30, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4p = -0.1))) # 2x250, 2x500 max iterations
  #abgrHeightFromDiameterMixed$sharmaPartonA1 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1 * (1 + a1rm) * topHeight^(b1)*(1 - exp(b2 * (standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022, # 10x10 MAE 1.81, AIC 630
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 20, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4p = -0.1)))
  #abgrHeightFromDiameterMixed$sharmaPartonB3 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3ra)*dbh))^(b4 + b4p * isPlantation), abgr2022, # a1ra singularity in backsolve, b2ra, b3ra 2x50 job max iterations, b4ra step halving, 10x10 MAE 1.78, AIC 599
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 25, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4p = -0.1)))
  #abgrHeightFromDiameterMixed$sharmaPartonB1 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1ra)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022, # b3ra more accurate
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 25, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4p = -0.1)))
  abgrHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3) * dbh))^(b4 + b4p * isPlantation + b4bal / (b4balr + basalAreaLarger)), abgr2022, 
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p + b4bal + b4balr ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 35, b1 = 0.2, b2 = -0.02, b3 = -0.04, b4 = 1.3, b4p = -0.1, b4bal = 0, b4balr = 2)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022,
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 35, b1 = 0.2, b2 = -0.02, b3 = -0.04, b4 = 1.3, b4p = -0.1)))
  abgrHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4 + b4bal / (b4balr + basalAreaLarger) + b4ac * cos(pi/18 *aspect)), abgr2022, 
                                                               fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 40, b1 = 0.2, b2 = -0.02, b3 = 0, b4 = 1.2, b4ac = 0, b4bal = 0, b4balr = 2)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1ac * cos(3.14159/180 * aspect))*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4ac * cos(pi/18 *aspect)), abgr2022,
  #                                                             fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 25, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3, b4ac = 0)))
  abgrHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4bal / (b4balr + basalAreaLarger) + b4ac * cos(pi/18 *aspect)), abgr2022, 
                                                                     fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 + b4ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 50, b1 = 0.1, b2 = -0.02, b3 = 0.1, b3rd = -0.04, b4 = 1.3, b4ac = 0, b4bal = 0, b4balr = 2)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b1rd * relativeDiameter)*dbh))^(b4 + b4ac * cos(pi/18 *aspect)), abgr2022,
  #                                                                   fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 + b4ac ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                                   start = list(fixed = c(a1 = 25, b1 = 0.3, b2 = -0.02, b3 = 0, b3rd = -0.04, b4 = 1.3, b4ac = 0)))
  abgrHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter) * dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), abgr2022, 
                                                               fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 50, b1 = 0.1, b2 = -0.02, b3 = 0.1, b3rd = -0.03, b4 = 1.3, b4bal = 0, b4balr = 3)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4), abgr2022, # a1ra, b2ra, b4ra singularity in backsolve
  #                                                             fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 33, b1 = 0.2, b2 = -0.02, b3 = -0.08, b3rd = -0.04, b4 = 1.3)))
  #abgrHeightFromDiameterMixed$sharmaPartonBalRelDbhB3 = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter + b3ra)*dbh))^(b4), abgr2022, # b1ra more accurate
  #                                                               fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = 33, b1 = 0.2, b2 = -0.02, b3 = -0.08, b3rd = -0.04, b4 = 1.3)))
  abgrHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1ac * cos(pi/180*aspect))*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3ra)*dbh))^(b4), abgr2022, # a1ra, b4ra singularity in backsolve, b1ra step halving, 10x10 MAE 1.78, AIC 622
                                                            fixedFormula = a1 + a1ac + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 25, a1ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3))) # 2x500 max iterations
  #abgrHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1ac * cos(pi/180*aspect))*topHeight^(b1) * (1 - exp(b2*(1 + b2rm)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), abgr2022, # a1ac not significant, 10x10 MAE 1.81, AIC 590
  #                                                          fixedFormula = a1 + a1ac + b1 + b2 + b3 + b4 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 35, a1ac = 0, b1 = 0.2, b2 = -0.02, b3 = 0.05, b4 = 1.2)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$sharmaPartonPhysioB1 = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1ac * cos(pi/180*aspect))*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), abgr2022, # b3ra more accurate
  #                                                            fixedFormula = a1 + a1ac + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 25, a1ac = 0.8, b1 = 0.3, b2 = -0.02, b3 = 0, b4 = 1.3)))
  abgrHeightFromDiameterMixed$sharmaPartonRelDbh = fit_nlme("Sharma-Parton RelDbh", height ~ 1.37 + (a1)*(1 + a1rm)*topHeight^(b1 + b1rd * relativeDiameter)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4), abgr2022, # b1rm, b2rm, b3rm, b4rm singularity in backsolve, 10x10 MAE 1.78, AIC 592
                                                            fixedFormula = a1 + b1 + b1rd + b2 + b3 + b4 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 35, b1 = 0.3, b1rd = -0.03, b2 = -0.02, b3 = 0.1, b4 = 1.3)))
  #abgrHeightFromDiameterMixed$sharmaPartonRelDbh = fit_nlme("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1rd * relativeDiameter + b1ra)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, abgr2022, # a1ra step halving, b2ra, b2rm, b3ra, b4ra singularity in backsolve, 10x10 MAE 1.85, AIC 582
  #                                                          fixedFormula = a1 + b1 + b1rd + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 25, b1 = 0.3, b1rd = -0.03, b2 = -0.02, b3 = 0.1, b4 = 1.3)))
  abgrHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_nlme("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a1ac * cos(pi/180*aspect))*topHeight^(b1 + b1rd * relativeDiameter + b1ra) * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^(b4), abgr2022, 
                                                             fixedFormula = a1 + a1ac + b1 + b1rd + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                             start = list(a1 = 25, a1ac = 0, b1 = 0.3, b1rd = -0.03, b2 = -0.02, b3 = 0.1, b4 = 1.3), distinct = FALSE)
  abgrHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1)*(1 - exp(b2*standTreesPerHectare^(b3 + b3ra)*dbh))^(b4 + b4p * isPlantation), abgr2022, # a1ra, b2ra, b4ra singularity in backsolve, 10x10 MAE 1.78, AIC 611
                                                     fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 33, b1 = 0.15, b2 = -0.03, b3 = -0.06, b4 = 1.2, b4p = -0.1)))
  #abgrHeightFromDiameterMixed$sharmaZhangA1 = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*(1 + a1rm)*standBasalAreaPerHectare^(b1)*(1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022, # b3 not significant, 10x10 MAE 1.88, AIC 622
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 33, b1 = 0.1, b2 = -0.03, b3 = 0, b4 = 1.3, b4p = -0.1)))
  #abgrHeightFromDiameterMixed$sharmaZhangB2 = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1)*(1 - exp(b2*(1 + b2rm)*standTreesPerHectare^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022, # b2rm, b4rm singularity in backsolve, b3 not signficant, b3rm max iterations, 10x10 MAE 1.81, AIC 616
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b2rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 38, b1 = 0.1, b2 = -0.03, b3 = 0, b4 = 1.3, b4p = -0.1)))
  abgrHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3 + b3ra) * dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, # a1ra, b2ra, b4ra singularity in backsolve, 10x10 MAE 1.80, AIC 583
                                                        fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ba ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 20, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4ba = 0.007))) # 2x500, 5x200 singularity in backsolve
  #abgrHeightFromDiameterMixed$sharmaZhangBalB1 = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1 + b1ra) * (1 - exp((b2)*standTreesPerHectare^(b3) * dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, # 10x10 MAE 1.91, AIC 591
  #                                                        fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ba ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 20, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4ba = 0.007)))
  #abgrHeightFromDiameterMixed$sharmaZhangB1 = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1ra)*(1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4 + b4p * isPlantation), abgr2022, # job false convergence
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 33, b1 = 0.15, b2 = -0.03, b3 = -0.06, b4 = 1.2, b4p = -0.1)))
  #abgrHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1 * standBasalAreaPerHectare^(b1) * (1 - exp(b2 * standTreesPerHectare^(b3 + b3ra)*dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, # a1ra, b2ra, b4ra singularity in backsolve
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ba ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 20, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4ba = 0.007)))
  #abgrHeightFromDiameterMixed$sharmaZhangBalB1 = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1 * standBasalAreaPerHectare^(b1 + b1ra) * (1 - exp(b2 * standTreesPerHectare^(b3)*dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, # b3ra more accurate
  #                                                        fixedFormula = a1 + b1 + b2 + b3 + b4 + b4ba ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 20, b1 = 0.2, b2 = -0.03, b3 = -0.1, b4 = 1.0, b4ba = 0.007)))
  abgrHeightFromDiameterMixed$sharmaZhangBalRelDbh = fit_nlme("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3 + b3rd * relativeDiameter) * dbh))^(b4 + b4ba * standBasalAreaPerHectare), abgr2022, 
                                                              fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 + b4ba ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 25, b1 = 0.2, b2 = -0.03, b3 = 0, b3rd = 0, b4 = 1.1, b4ba = 0.007)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$sharmaZhangRelDbh = fit_nlme("Sharma-Zhang RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1 + b1ra)*(1 - exp((b2)*standTreesPerHectare^(b3 + b3rd * relativeDiameter)*dbh))^(b4), abgr2022, # a1ra, b2ra, b3ra, b4ra singularity in backsolve
  #                                                         fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 65, b1 = 0.05, b2 = -0.01, b3 = 0, b3rd = -0.02, b4 = 1.2)), distinct = FALSE) # b3, b3rd not significant
  abgrHeightFromDiameterMixed$sharmaZhangRelDbh = fit_nlme("Sharma-Zhang RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*(1 + b2rm)*standTreesPerHectare^(b3 + b3rd * relativeDiameter)*dbh))^(b4), abgr2022, # a1rm-b3rd not mutually significant, b1rm, b4rm singularity in backsolve, b3rm step halving
                                                           fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                           start = list(fixed = c(a1 = 65, b1 = 0, b2 = -0.01, b3 = 0.2, b3rd = -0.03, b4 = 1.4))) # b1, b2 not significant, 2x50, 2x250, 2x500 singularity in backsolve
  abgrHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1)*(1 + b1rm)*dbh^(b2)), abgr2022, # b2rm max iterations, 10x10 MAE 1.90, AIC 630
                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.2, b1 = 2.3, b2 = -0.1)))
  #abgrHeightFromDiameterMixed$sibbesenA1 = fit_nlme("Sibbesen", height ~ 1.37 + a1*(1 + a1rm)*dbh^((b1)*dbh^(b2)), abgr2022, # 10x10 MAE 1.98, AIC 630
  #                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 0.2, b1 = 2.0, b2 = -0.1)))
  #abgrHeightFromDiameterMixed$sibbesenB2 = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1)*dbh^(b2 + b2ra)), abgr2022, # a1ra singular precision, b2ra 2x50 job step halving, 10x10 MAE 1.94, AIC 631
  #                                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 0.2, b1 = 1.8, b2 = -0.1)))
  #abgrHeightFromDiameterMixed$sibbesenB1 = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1ra)*dbh^(b2)), abgr2022, # b2ra more accurate but higher AIC
  #                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 0.2, b1 = 1.8, b2 = -0.1)))
  abgrHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1)*(1 - exp(b1*dbh^(b2*(1 + b2rm)))), abgr2022, # 10x10 MAE 1.88, AIC 612
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 61, b1 = -0.01, b2 = 1.2))) # 2x50, 2x250, 2x500, 5x100, 5x200, 10x25, 10x100 max iterations, 5x200 max iterations
  #abgrHeightFromDiameterMixed$weibullB1 = fit_nlme("Weibull", height ~ 1.37 + (a1)*(1 - exp(b1*(1 + b1rm)*dbh^(b2))), abgr2022, # a1rm max iterations, 10x10 MAE 1.94, AIC604
  #                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 64, b1 = -0.01, b2 = 1.2)))
  #abgrHeightFromDiameterMixed$weibullB2 = fit_nlme("Weibull", height ~ 1.37 + (a1)*(1 - exp(b1*dbh^(b2 + b2ra))), abgr2022, # b1ra step halving, b2ra 2x50 job max iterations, 10x10 MAE 1.90, AIC 613
  #                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 58, b1 = -0.01, b2 = 1.2)))
  #abgrHeightFromDiameterMixed$weibullA1 = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1ra)*(1 - exp(b1*dbh^b2)), abgr2022, # b2ra more accurate
  #                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 58, b1 = -0.01, b2 = 1.2)))
  abgrHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1bal / (b1balr + basalAreaLarger))*dbh^(b2*(1 + b2rm)))), abgr2022, # 10x10 MAE 1.75, AIC 598
                                                    fixedFormula = a1 + b1 + b1bal + b1balr + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 75, b1 = -0.006, b1bal = 0.02, b1balr = 10, b2 = 1.2))) # 2x50, 2x500 max iterations, 5x200 step halving
  #abgrHeightFromDiameterMixed$weibullBalA1 = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1bal / (b1balr + basalAreaLarger))*dbh^(b2))), abgr2022, # b1ra max iterations, 10x10 MAE 1.94, AIC 601
  #                                                    fixedFormula = a1 + b1 + b1bal + b1balr + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 70, b1 = -0.006, b1bal = 0.02, b1balr = 5, b2 = 1.2)))
  #abgrHeightFromDiameterMixed$weibullBalB2 = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger + b2ra))), abgr2022,
  #                                                    fixedFormula = a1 + b1 + b2 + b2bal ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 100, b1 = -0.01, b2 = 0.9, b2bal = 0.003)))
  abgrHeightFromDiameterMixed$weibullBalRelDbh = fit_nlme("Weibull RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp((b1 + b1bal / (b1balr + basalAreaLarger))*dbh^(b2))), abgr2022, 
                                                          fixedFormula = a1 + a1rd + b1 + b1bal + b1balr + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 70, a1rd = 0, b1 = -0.006, b1bal = 0.02, b1balr = 5, b2 = 1.2)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp((b1)*dbh^(b2 + b2ra))), abgr2022, # a1ra max iterations, b1ra step halving, b2ra job max iterations, a1rd+b1+b2ra not mutually significant, signs of parameter evaporation, 10x10 MAE 1.83, AIC 584
  #                                                     fixedFormula = a1 + a1rd + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 80, a1rd = -6, b1 = -0.01, b2 = 1.2)))
  #abgrHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 + a1rm)*(1 - exp((b1)*dbh^(b2))), abgr2022, # a1rd-a1rm not mutually significant
  #                                                     fixedFormula = a1 + a1rd + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 100, a1rd = -10, b1 = -0.01, b2 = 1.2)), distinct = FALSE)
  #abgrHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp((b1)*(1 + b1rm)*dbh^(b2))), abgr2022, # 10x10 MAE 1.92, AIC 608
  #                                                     fixedFormula = a1 + a1rd + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 100, a1rd = -10, b1 = -0.01, b2 = 1.2)))
  abgrHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp((b1)*dbh^(b2*(1 + b2rm)))), abgr2022, # 10x10 MAE 1.82, AIC 578
                                                       fixedFormula = a1 + a1rd + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 90, a1rd = -6, b1 = -0.01, b2 = 1.2))) # 2x50, 2x250, 2x500, 5x100, 5x200 max iterations
  #to_fixed_coeffficients(abgrHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(abgrHeightFromDiameterMixed$sharmaZhang), n = 40)
  #abgrHeightFromDiameterMixed$weibullBal$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  abgrHeightFromDiameterMixed$gam = fit_gam("REML GAM", height ~ I(0.6539 * dbh^0.9757) + s(dbh, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # collapses to linear
  abgrHeightFromDiameterMixed$gamBal = fit_gam("REML GAM BA+L", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 4), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # collapses to linear
  abgrHeightFromDiameterMixed$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 14), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # by propagation
  abgrHeightFromDiameterMixed$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 15), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # by propagation
  abgrHeightFromDiameterMixed$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, basalAreaLarger, relativeDiameter, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # by propagation
  abgrHeightFromDiameterMixed$gamPhysio = fit_gam("REML GAM physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # collapses to linear
  abgrHeightFromDiameterMixed$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(0.6539 * dbh^0.9757) + s(dbh, relativeDiameter, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # collapses to linear
  abgrHeightFromDiameterMixed$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(0.6539 * dbh^0.9757) + s(dbh, relativeDiameter, sin(3.14159/180 * aspect), bs = "ts", k = 12), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # by propagation
  #lapply(abgrHeightFromDiameterMixed$gamPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(abgrHeightFromDiameterMixed$gamPhysio$fit, function(fit) { return(summary(fit$gam)) })
  #abgrHeightFromDiameterMixed$gam$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(abgrHeightFromDiameterMixed, paste0("trees/height-diameter/data/ABGR height mixed ", get_cross_validation_data_suffix(), ".Rds"))
}


## grand fir diameter regressions
if (abgrOptions$fitDbhPrimary)
{
  # linear regressions
  abgrDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ I(height - 1.37), abgr2022)) # plantation not significant
  #abgrDiameterFromHeight$linearRelHtTopHt = fit_lm("linear", dbh ~ I(height - 1.37) + relativeHeight + topHeight, abgr2022, distinct = FALSE) # relativeHeight and topHeight not significant
  abgrDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ I(height - 1.37) + I((height - 1.37)^2), abgr2022, distinct = FALSE) # dbh^2, plantation^2 not significant
  #lapply(abgrDiameterFromHeight$linearRelHt$fit, summary)
  
  # nonlinear regressions
  abgrDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ (a1)*(exp(b1*(height - 1.37)) - 1)^b2, abgr2022, 
                                                      start = list(a1 = 100, b1 = 0.01, b2 = 0.9)) # a1p, b1p, b2p not significant
  abgrDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(exp((b1)*(height - 1.37)) - 1)^(b2), abgr2022, 
                                                          start = list(a1 = 100, a1aba = 0, b1 = 0.01, b2 = 0.9), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  abgrDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ (a1)*(exp((b1 + b1rh * relativeHeight)*(height - 1.37)) - 1)^(b2), abgr2022, 
                                                           start = list(a1 = 100, b1 = 0.01, b1rh = 0, b2 = 0.9), distinct = FALSE) # a1rh, b1rh, b2rh not significant but b1rh closest
  abgrDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), abgr2022, 
                                                       start = list(a1 = -100, b1 = 0.01, b2 = 0.9)) # a1p, b1p, b2p not significant
  abgrDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), abgr2022, 
                                                           start = list(a1 = -100, a1aba = 0, b1 = 0.01, b2 = 0.9), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  abgrDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ (a1)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), abgr2022, 
                                                             start = list(a1 = 100, b1 = 0.01, b2 = 0.9, b2tw = 0), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  abgrDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ a1*log(1 - pmin((b1*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), abgr2022, 
                                                            start = list(a1 = -90, b1 = 0.012, b2 = 0.9, b2rh = -0.0), distinct = FALSE) # a1rh, b1rh, b2rh not significant
  abgrDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), abgr2022, 
                                                              start = list(a1 = 200, a2 = 100, b1 = 0.9)) # a1p, a2p, b2p not significant
  abgrDiameterFromHeight$michaelisMentenReplaceAbat = fit_gsl_nls("Michaelis-Menten replace ABA+T", dbh ~ (a1) * (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare)), abgr2022, 
                                                                  start = list(a1 = 200, a2 = 100, b1 = 0.9, b1aba = -0.0005)) # a1aba, a1abar, a1aat, a1aat+r, a2aba, a2abar, a2aat, a2aat+r, b1abap, b1abar, b1aat, b1aat+r not significant
  abgrDiameterFromHeight$michaelisMentenReplaceAbatRelHt = fit_gsl_nls("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight) * (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare)), abgr2022, 
                                                                       start = list(a1 = 200, a1rh = 0, a2 = 100, b1 = 0.9, b1aba = -0.0005), distinct = FALSE) # a1rh, a2rh, b1rh not significant
  abgrDiameterFromHeight$michaelisMentenReplaceRelHt = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), abgr2022, 
                                                                   start = list(a1 = 200, a1rh = 0, a2 = 100, b1 = 0.9), distinct = FALSE) # a1rh, a2rh, b1rh not significant
  abgrDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), abgr2022, 
                                               start = list(a1 = 3.2, a2 = -0.1, a2p = -0.005)) # a1p not significant
  abgrDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, abgr2022, 
                                             start = list(a1 = 2.0, b1 = 0.9)) # a1p, b1p not significant, 10x10 RMSE 6.75, AIC 921
  abgrDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aat * basalAreaTaller), abgr2022, 
                                                 start = list(a1 = 2.1, b1aat = -0.001, b1 = 0.9)) # a1aba, b1aba not significant, 10x10 RMSE 6.51, AIC 919
  #abgrDiameterFromHeight$powerAbatA1 = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1), abgr2022, 
  #                                                 start = list(a1 = 2.2, a1aat = -0.01, b1 = 0.9)) # 10x10 RMSE 6.51, AIC 933
  abgrDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ a1*(height - 1.37)^(b1 + b1e * topographicWetnessFD8f), abgr2022, 
                                                   start = list(a1 = 2.0, b1 = 0.9, b1e = 0.01)) # { a1, b1 } x { e, s, a, tr }, a1e not significant, 10x10 RMSE 6.46, AIC 957
  abgrDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight), abgr2022, 
                                                  start = list(a1 = 2.0, b1 = 0.9, b1rh = 0.0)) # 10x10 RMSE 6.65, AIC 902
  #abgrDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^(b1), abgr2022, 
  #                                                start = list(a1 = 2.0, a1rh = 0.4, b1 = 0.9)) # 10x10 RMSE 6.68, AIC 942
  #abgrDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ a1*(height - 1.37)^b1*relativeHeight^b1rh, abgr2022, 
  #                                                start = list(a1 = 2.0, b1 = 0.9, b1rh = 0), distinct = FALSE) # b1rh not significant, 10x10 RMSE 6.72, AIC 956
  #abgrDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ a1*(height - 1.37)^b1*(1 + relativeHeight)^b1rh, abgr2022, 
  #                                                start = list(a1 = 2.0, b1 = 0.9, b1rh = 0.4)) # 10x10 RMSE 6.85, AIC 930
  #abgrDiameterFromHeight$powerTopHt = fit_gsl_nls("power", dbh ~ a1*(topHeight - 1.37)^b1th*(height - 1.37)^b1, abgr2022, 
  #                                                start = list(a1 = 2.0, b1th = 0, b1 = 0.9), distinct = FALSE) # b1th not significant, 10x10 RMSE 6.72, AIC 910
  #abgrDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh / (b1rhr + relativeHeight)), abgr2022, 
  #                                                start = list(a1 = 2.0, b1 = 0.9, b1rh = 0, b1rhr = 1), distinct = FALSE) # a1rh+r, b1rh+r not significant
  abgrDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), abgr2022, 
                                             start = list(a1 = 2.2, b1 = 0.8, b2 = 0.007)) # a1p, b1p, b2p not significant
  abgrDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aat * basalAreaTaller) * (height - 1.37)), abgr2022, 
                                                 start = list(a1 = 2.1, b1 = 0.8, b2 = 0.004, b2aat = -0.0001)) # 10x10 RMSE 6.48, AIC 927
  #abgrDiameterFromHeight$ruarkAbatA1aba = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), abgr2022, 
  #                                                    start = list(a1 = 2.3, a1aba = -0.003, b1 = 0.8, b2 = 0.005)) # 10x10 RMSE 6.85, AIC 940
  #abgrDiameterFromHeight$ruarkAbatB1aba = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) * exp((b2) * (height - 1.37)), abgr2022, 
  #                                                    start = list(a1 = 2.2, b1 = 0.9, b2 = 0.007, b1aba = -0.001)) # 10x10 RMSE 6.75, AIC 876
  #abgrDiameterFromHeight$ruarkAbatB2aba = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)), abgr2022, 
  #                                                    start = list(a1 = 2.1, b1 = 0.8, b2 = 0.01, b2aba = -0.0001)) # 10x10 RMSE 6.82, AIC 935
  #abgrDiameterFromHeight$ruarkAbatA1bat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a1aat * basalAreaTaller)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), abgr2022, 
  #                                                    start = list(a1 = 2.3, a1aat = -0.01, b1 = 0.8, b2 = 0.005)) # 10x10 RMSE 6.59, AIC 935
  #abgrDiameterFromHeight$ruarkAbatB1bat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aat * basalAreaTaller) * exp((b2) * (height - 1.37)), abgr2022, 
  #                                                    start = list(a1 = 2.2, b1 = 0.9, b2 = 0.003, b1aat = -0.001)) # 10x10 RMSE 6.65, AIC 929
  #abgrDiameterFromHeight$ruarkAbatB1abaB2bat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) * exp((b2 + b2aat * basalAreaTaller) * (height - 1.37)), abgr2022, 
  #                                                         start = list(a1 = 2.2, b1 = 0.9, b1aba = -0.001, b2 = 0.007, b2aat = -0.0001), distinct = FALSE) # b1aba-b2aat not mutually significant
  abgrDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aat * basalAreaTaller + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, 
                                                       start = list(a1 = 2.2, b1 = 0.8, b2 = 0.01, b2aat = -0.0001, b2tw = 0), distinct = FALSE) # b2tw not significant
  abgrDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark RelHt ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2aat * basalAreaTaller) * (height - 1.37)), abgr2022, 
                                                      start = list(a1 = 2.2, b1 = 0.8, b1rh = 0.03, b2 = 0.002, b2aat = 0), distinct = FALSE) # b1rh-b2aat not mutually significant
  abgrDiameterFromHeight$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2aat * basalAreaTaller + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, 
                                                            start = list(a1 = 2.2, b1 = 0.8, b1rh = 0, b2 = 0.01, b2aat = -0.0001, b2tw = 0), distinct = FALSE) # by propagation
  abgrDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^(b1) * exp((b2 + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, 
                                                   start = list(a1 = 2.2, b1 = 0.8, b2 = 0.01, b2tw = 0.001)) # { a1, b1, b2 } x { e, s, a, tr }, a1tw, b2 not significant
  #abgrDiameterFromHeight$ruarkPhysioB1 = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^(b1 + b1tw * topographicWetnessFD8f) * exp((b2) * (height - 1.37)), abgr2022, 
  #                                                   start = list(a1 = 2.2, b1 = 0.8, b1tw = 0.007, b2 = 0.007)) # b2tw possibly slightly lower error
  abgrDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ a1*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp(b2 * (height - 1.37)), abgr2022, 
                                                  start = list(a1 = 2.2, b1 = 0.8, b1rh = 0.03, b2 = 0.002)) # a1rh, b2rh not significant
  abgrDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ a1*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, 
                                                        start = list(a1 = 2.2, b1 = 0.8, b1rh = 0.03, b2 = 0.002, b2tw = 0.001), distinct = FALSE) # b1rh-b2tw not mutually significant
  abgrDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1), 0.9999)), abgr2022, 
                                               start = list(a1 = 0.000001, a2 = 0.0003, b1 = 0.9, Ha = 130))
  abgrDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^b4, abgr2022, 
                                                    start = list(a1 = 2.2, b1 = 0.9, b2 = 2.5, b3 = 0.2, b4 = 0.002)) # a1p, b1p, b2, b2p, b3, b3p, b4p not significant, 2x50, 2x500 NaN-inf
  abgrDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), abgr2022, 
                                                       start = list(a1 = 2.2, b1 = 0.8, b2 = 0.05)) # a1p, b1p, b2p not significant, 10x10 RMSE 6.73, AIC 940
  abgrDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), abgr2022, 
                                                           start = list(a1 = 2.2, a1aba = 0, b1 = 0.8, b2 = 0.05), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  abgrDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tw * topographicWetnessFD8f)*(height - 1.37)^(b1*(height - 1.37)^b2), abgr2022, 
                                                                 start = list(a1 = 1.8, a1aba = 0, a1tw = 0.06, b1 = 0.8, b2 = 0.04), distinct = FALSE) # by propagation
  abgrDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), abgr2022, 
                                                                start = list(a1 = 2.2, a1aba = 0, b1 = 0.8, b2 = 0.04, b2rh = 0), distinct = FALSE) # by propagation
  abgrDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tw * topographicWetnessFD8f)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), abgr2022, 
                                                                      start = list(a1 = 2.2, a1aba = 0, a1tw = 0, b1 = 0.8, b2 = 0.04, b2rh = 0), distinct = FALSE) # by propagation
  abgrDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(height - 1.37)^((b1)*(height - 1.37)^b2), abgr2022, 
                                                             start = list(a1 = 1.8, a1tw = 0.06, b1 = 0.8, b2 = 0.04)) # { b1, b2 } x { e, s, a, tr }, b2tw not significant
  #abgrDiameterFromHeight$sibbesenReplacePhysioB1 = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1)*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^b2), abgr2022, 
  #                                                             start = list(a1 = 2.2, b1 = 0.8, b1tw = 0.08, b2 = 0.04)) # a1tw and b1tw near interchangeable, a1tw maybe slightly lower AIC
  abgrDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^(b1*(relativeHeight*(height - 1.37))^b2), abgr2022, 
                                                            start = list(a1 = 2.2, b1 = 0.8, b2 = 0.04), distinct = FALSE) # 10x10 RMSE 6.68, AIC 905
  #abgrDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2rh * relativeHeight)), abgr2022, 
  #                                                          start = list(a1 = 2.2, b1 = 0.8, b2 = 0.04, b2rh = 0), distinct = FALSE) # a1rh, b1rh, b2rh not significant, 10x10 RMSE 6.81, AIC 918
  #abgrDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ a1*(relativeHeight*(height - 1.37))^(b1*(height - 1.37)^b2), abgr2022, 
  #                                                          start = list(a1 = 2.2, b1 = 0.8, b2 = 0.04), distinct = FALSE) # 10x10 RMSE 7.25, AIC 975
  abgrDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), abgr2022, 
                                                                  start = list(a1 = 2.2, a1tw = 0.06, b1 = 0.8, b2 = 0.04, b2rh = 0), distinct = FALSE) # by propagation, also a1tw-b2rh not mutually significant
  #abgrDiameterFromHeight$sibbesenTopHt = fit_gsl_nls("Sibbesen TopHt", dbh ~ a1*(height - 1.37)^(b1*(topHeight - 1.37)^b2), abgr2022, 
  #                                                   start = list(a1 = 1.9, b1 = 1.0, b2 = 0), distinct = FALSE) # b2 not significant, 10x10 RMSE 6.82, AIC 930
  #abgrDiameterFromHeight$sibbesenTopHt = fit_gsl_nls("Sibbesen TopHt", dbh ~ a1*(topHeight - 1.37)^(b1*(height - 1.37)^b2), abgr2022, 
  #                                                   start = list(a1 = 0.9, b1 = 0.4, b2 = 0.3)) # 10x10 RMSE 11.3, AIC 1183
  abgrDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, abgr2022, 
                                               start = list(a1 = -250, b1 = 0.01, b2 = 0.9)) # a1p, b1p, b2p not significant, signs of a1-b1 evaporation
  #abgrDiameterFromHeight$weibull$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #print(to_parameter_confidence_intervals(abgrDiameterFromHeight$weibull), n = 20)
  #to_fixed_coeffficients(abgrDiameterFromHeight$weibull)
  
  # GAM NaNs with Γ link
  abgrDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, bs = "ts", k = 3), data = abgr2022, threads = 2) # plantation not significant, 10x10 RMSE 6.76 AIC 919
  abgrDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 11), data = abgr2022, threads = 2, distinct = FALSE) # min k = 11 > 3-8
  #abgrDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bs = "ts", by = as.factor(isPlantation), k = 3), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 3), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  abgrDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), data = abgr2022, threads = 2, distinct = FALSE) # min k = 16 > 3-10
  abgrDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), data = abgr2022, threads = 2, distinct = FALSE) # min k = 16 > 12
  # plausible but modest improvement with any physiographic predictor but all predictor min k = 331 > ~150 -> drop elevation + cos(aspect) on highest RMSE min k = 57 > ~42 -> drop slope
  #abgrDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, elevation, bs = "ts", k = 4), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, slope, bs = "ts", k = 4), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, sin(3.14159/180 * aspect), bs = "ts", k = 4), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, cos(3.14159/180 * aspect), bs = "ts", k = 4), data = abgr2022, threads = 2, distinct = FALSE) # collapses to linear
  #abgrDiameterFromHeight$gamRoughness = fit_gam("REML GAM roughness", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, terrainRoughness, bs = "ts", k = 10), data = abgr2022, threads = 2) # 10x10 RMSE 6.70 AIC 919
  #abgrDiameterFromHeight$gamWetness = fit_gam("REML GAM wetness", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, topographicWetnessFD8f, bs = "ts", k = 9), data = abgr2022, threads = 2) # 10x10 RMSE 6.54 AIC 954
  #abgrDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, terrainRoughness, topographicWetnessFD8f, bs = "ts", k = 11), data = abgr2022, threads = 2, distinct = FALSE) # minimum k = 11 > ~6
  abgrDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, terrainRoughness, bs = "ts", k = 10), data = abgr2022, threads = 2)
  abgrDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, relativeHeight, bs = "ts", k = 4), data = abgr2022, threads = 2, distinct = FALSE) # k = 4 minimum, plantation not significant, 10x10 RMSE 6.82 AIC 948
  abgrDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, relativeHeight, terrainRoughness, bs = "ts", k = 11), data = abgr2022, threads = 2, distinct = FALSE) # minimum k = 11 > ~5
  #lapply(abgrDiameterFromHeight$gamRelHtPhysio$fit, k.check)
  #lapply(abgrDiameterFromHeight$gamRelHtPhysio$fit, summary)
  #abgrDiameterFromHeight$gamRelHtPhysio$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
    
  abgrDiameterFromHeight$randomForest = fit_ranger("random forest", c("dbh", "height"),
                                                   abgr2022, mtry = 1, minNodeSize = 36, sampleFraction = 0.526) # from setup.R tuneRanger @ 500 trees
  abgrDiameterFromHeight$randomForestAbatRelHtPhysio = fit_ranger("random forest RelHt ABA+T physio", c("dbh", "height", "relativeHeight", "topHeight", "bootstrapStandBasalAreaPerHectare", "slope"), # from setup.R VSURF
                                                                  abgr2022, mtry = 3, minNodeSize = 4, sampleFraction = 0.216) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(abgrDiameterFromHeight, paste0("trees/height-diameter/data/ABGR DBH ", get_cross_validation_data_suffix(), ".Rds"))
}

if (abgrOptions$fitDbhMixed)
{
  abgrDiameterFromHeightMixed = list(linear = fit_lme("linear", fixedFormula = dbh ~ I(height - 1.37), abgr2022, randomFormula = ~1|stand/plot)) # 2x500 false convergence, 5x200 singular convergence
  abgrDiameterFromHeightMixed$parabolic = fit_lme("parabolic", fixedFormula = dbh ~ I((height - 1.37)^2), abgr2022, randomFormula = ~1|stand/plot)
  
  #abgrDiameterFromHeightMixed$chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1) * (exp((b1)*(height - 1.37)) - 1)^(b2 + b2ra), abgr2022, # a1ra max iterations, b1ra, b2ra singularity in backsolve
  #                                                      fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 100, b1 = 0.01, b2 = 0.9)))
  abgrDiameterFromHeightMixed$chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1)*(1 + a1rm)*(exp((b1)*(height - 1.37)) - 1)^(b2), abgr2022, # b2rm singularity in backsolve, 10x10 RMSE 7.61, AIC 931
                                                        fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 40, b1 = 0.03, b2 = 0.8)))
  #abgrDiameterFromHeightMixed$chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1) * (exp((b1)*(1 + b1rm)*(height - 1.37)) - 1)^(b2), abgr2022, # 10x10 RMSE 8.27, AIC 953
  #                                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 30, b1 = 0.03, b2 = 0.8)))
  abgrDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(exp((b1)*(height - 1.37)) - 1)^(b2), abgr2022, 
                                                            fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 100, a1aba = 0, b1 = 0.01, b2 = 0.9)), distinct = FALSE)
  abgrDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1)*(exp((b1 + b1rh * relativeHeight)*(height - 1.37)) - 1)^(b2), abgr2022,
                                                             fixedFormula = a1 + b1 + b1rh + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 100, b1 = 0.01, b1rh = 0, b2 = 0.9)), distinct = FALSE)
  #abgrDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1ra)*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), abgr2022, # a1ra job max iterations, b1ra, b2ra singularity in backsolve, 10x10 RMSE 8.99, AIC 1016
  #                                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = -100, b1 = 0.01, b2 = 0.9)))
  abgrDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1)*(1 + a1rm)*log(1 - pmin((b1*(height - 1.37))^b2, 0.9999)), abgr2022, # 10x10 RMSE 7.18, AIC 957
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = -100, b1 = 0.01, b2 = 0.8)))
  #abgrDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1)*log(1 - pmin((b1*(1 + b1rm)*(height - 1.37))^b2, 0.9999)), abgr2022, # b2rm singularity in backsolve, 10x10 RMSE 41.5, AIC 1413
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = -45, b1 = 0.01, b2 = 0.8)))
  abgrDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*log(1 - pmin(((b1)*(height - 1.37))^(b2), 0.9999)), abgr2022, 
                                                             fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = -100, a1aba = 0, b1 = 0.01, b2 = 0.9)), distinct = FALSE)
  abgrDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1ra)*log(1 - pmin(((b1)*(height - 1.37))^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), abgr2022,
                                                               fixedFormula = a1+ b1 + b2 + b2tw ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 100, b1 = 0.01, b2 = 0.9, b2tw = 0)))
  abgrDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1ra)*log(1 - pmin((b1*(height - 1.37))^(b2 + b2rh * relativeHeight), 0.9999)), abgr2022,
                                                              fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = -90, b1 = 0.012, b2 = 0.9, b2rh = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), abgr2022, # 10x10 RMSE 7.40, AIC 932
                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 200, a2 = 100, b1 = 0.9)))
  #abgrDiameterFromHeightMixed$michaelisMentenReplaceB1 = fit_nlme("Michaelis-Menten replace", dbh ~ (a1) * (height - 1.37)^(b1 + b1ra) / (a2 - (height - 1.37)^(b1 + b1ra)), abgr2022, # a1ra, a2ra singularity in backsolve, 10x10 RMSE 9.41, AIC 1008
  #                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 140, a2 = 60, b1 = 0.8)))
  #abgrDiameterFromHeightMixed$michaelisMentenReplaceB1 = fit_nlme("Michaelis-Menten replace", dbh ~ (a1) * (height - 1.37)^(b1*(1 + b1rm)) / (a2 - (height - 1.37)^(b1*(1 + b1rm))), abgr2022, # 10x10 RMSE 8.34, AIC 963
  #                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 80, a2 = 35, b1 = 0.8)))
  abgrDiameterFromHeightMixed$michaelisMentenReplaceAbat = fit_nlme("Michaelis-Menten replace ABA+T", dbh ~ (a1)*(1 + a1rm)* (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare)), abgr2022, # 10x10 RMSE 7.03, AIC 951
                                                                    fixedFormula = a1 + a2 + b1 + b1aba ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 120, a2 = 60, b1 = 0.8, b1aba = -0.0005)))
  #abgrDiameterFromHeightMixed$michaelisMentenReplaceAbat = fit_nlme("Michaelis-Menten replace ABA+T", dbh ~ (a1) * (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare + b1ra) / (a2 - (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare + b1ra)), abgr2022, # a1ra step halving, a2ra max iterations, 10x10 RMSE 9.46, AIC 970
  #                                                                  fixedFormula = a1 + a2 + b1 + b1aba ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 200, a2 = 100, b1 = 0.9, b1aba = -0.0005)))
  abgrDiameterFromHeightMixed$michaelisMentenReplaceAbatRelHt = fit_nlme("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight) * (1 + a1rm) * (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare) / (a2 - (height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare)), abgr2022, 
                                                                         fixedFormula = a1 + a1rh + a2 + b1 + b1aba ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                         start = list(fixed = c(a1 = 200, a1rh = 0, a2 = 100, b1 = 0.9, b1aba = -0.0005)), distinct = FALSE)
  abgrDiameterFromHeightMixed$michaelisMentenReplaceRelHt = fit_nlme("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight) * (1 + a1rm) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), abgr2022, 
                                                                     fixedFormula = a1 + a1rh + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 200, a1rh = 0, a2 = 100, b1 = 0.9)), distinct = FALSE)
  abgrDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1)*(1 + a1rm)*sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), abgr2022, # 10x10 RMSE 10.4, AIC 1132
                                                 fixedFormula = a1 + a2 + a2p ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 3.2, a2 = -0.1, a2p = -0.01))) # 2x500 max iterations
  #abgrDiameterFromHeightMixed$naslundA1 = fit_nlme("Näslund inverse", dbh ~ (a1 + a1ra) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), abgr2022, # a1ra sometimes job false convergence, 10x10 RMSE 10.6, AIC 1136
  #                                                 fixedFormula = a1 + a2 + a2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 3.2, a2 = -0.1, a2p = -0.005)))
  #abgrDiameterFromHeightMixed$naslundA2 = fit_nlme("Näslund inverse", dbh ~ (a1) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation + a2ra) * sqrt(height - 1.37)), abgr2022, # 10x10 RMSE 26, AIC 1317
  #                                                 fixedFormula = a1 + a2 + a2p ~ 1, randomFormula = a2ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 3.2, a2 = -0.1, a2p = -0.005)))
  #abgrDiameterFromHeightMixed$naslundA2 = fit_nlme("Näslund inverse", dbh ~ (a1)*sqrt(height - 1.37) / (1 + a2 * (1 + a2rm) * sqrt(height - 1.37)), abgr2022, # a2p not significant, 10x10 RMSE 105, AIC 1692
  #                                                 fixedFormula = a1 + a2 ~ 1, randomFormula = a2rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 2.2, a2 = -0.1)))
  #abgrDiameterFromHeightMixed$naslundA2 = fit_nlme("Näslund inverse", dbh ~ (a1) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation + a2ra) * sqrt(height - 1.37)), abgr2022, # a1ra 2.8x lower RMSE
  #                                                 fixedFormula = a1 + a2 + a2p ~ 1, randomFormula = a2ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 3.2, a2 = -0.1, a2p = -0.005)))
  abgrDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^(b1), abgr2022, # 10x10 RMSE 6.63, AIC 928
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2.1, b1 = 0.9))) # 10x100 step halving
  #abgrDiameterFromHeightMixed$powerA1 = fit_nlme("power", dbh ~ (a1 + a1ra)*(height - 1.37)^(b1), abgr2022, # 10x10 RMSE 6.70, AIC 939
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2.0, b1 = 0.9)))
  #abgrDiameterFromHeightMixed$powerB1 = fit_nlme("power", dbh ~ (a1)*(height - 1.37)^(b1*(1 + b1rm)), abgr2022, # 10x10 RMSE 6.73, AIC 936
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2.0, b1 = 0.9)))
  #abgrDiameterFromHeightMixed$powerB1 = fit_nlme("power", dbh ~ (a1)*(height - 1.37)^(b1 + b1ra), abgr2022, 
  #                                             fixedFormula = a1 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot, # a1ra lower RMSE, same AIC
  #                                             start = list(fixed = c(a1 = 2.0, b1 = 0.9)))
  #abgrDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1ra)*(height - 1.37)^(b1 + b1aat * basalAreaTaller), abgr2022, # a1ra singular precision, b1ra singularity in backsolve
  #                                                 fixedFormula = a1 + b1 + b1aat ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 2.1, b1aat = -0.001, b1 = 0.9)))
  abgrDiameterFromHeightMixed$powerAbat = create_fit_statistics(name = "power ABA+T", fitting = "nlme")
  abgrDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^(b1 + b1e * topographicWetnessFD8f), abgr2022, # 10x10 RMSE 6.68, AIC 904
                                                     fixedFormula = a1 + b1 + b1e ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 2.1, b1 = 0.9, b1e = 0.01)))
  #abgrDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1e * topographicWetnessFD8f + b1ra), abgr2022, # 10x10 RMSE 6.70, AIC 947
  #                                                   fixedFormula = a1 + b1 + b1e ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 2.0, b1 = 0.9, b1e = 0.01)))
  #abgrDiameterFromHeightMixed$powerPhysioA1 = fit_nlme("power physio", dbh ~ (a1 + a1ra)*(height - 1.37)^(b1 + b1e * topographicWetnessFD8f), abgr2022, # b1ra more accurate
  #                                                     fixedFormula = a1 + b1 + b1e ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 2.0, b1 = 0.9, b1e = 0.01)))
  abgrDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^(b1 + b1ra), abgr2022, # a1ra max iterations, 10x10 RMSE 6.67, AIC 932
                                                    fixedFormula = a1 + a1rh + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2.0, a1rh = 0.4, b1 = 0.9)))
  #abgrDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(1 + a1rm)*(height - 1.37)^(b1), abgr2022, # 10x10 RMSE 6.75, AIC 946
  #                                                  fixedFormula = a1 + a1rh + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 2.3, a1rh = 1.0, b1 = 0.8)))
  abgrDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1*(1 + b1rm)) * exp((b2) * (height - 1.37)), abgr2022, # b2rm max iterations, 10x10 RMSE 6.69, AIC 926
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 2.1, b1 = 0.9, b2 = 0.007)))
  #abgrDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1 + b1ra) * exp((b2) * (height - 1.37)), abgr2022, # a1ra max iterations, 10x10 RMSE 6.93, AIC 923
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 2.3, b1 = 0.8, b2 = 0.007)))
  #abgrDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), abgr2022, # 10x10 RMSE 7.13, AIC 921
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 2.1, b1 = 0.9, b2 = 0.007)))
  #abgrDiameterFromHeightMixed$ruarkB2 = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1 + b2ra) * exp((b2) * (height - 1.37)), abgr2022, # b1ra more accurate
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.007)))
  abgrDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2aat * basalAreaTaller) * (height - 1.37)), abgr2022, # a1ra max iterations, 10x10 RMSE 6.87, AIC 933
                                                   fixedFormula = a1 + b1 + b2 + b2aat ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 2.1, b1 = 0.8, b2 = 0.004, b2aat = -0.0001)))
  #abgrDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1*(1 + b1rm)) * exp((b2 + b2aat * basalAreaTaller) * (height - 1.37)), abgr2022, # a1ra max iterations, 10x10 RMSE 6.93, AIC 932
  #                                                 fixedFormula = a1 + b1 + b2 + b2aat ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 2.1, b1 = 0.9, b2 = 0.004, b2aat = -0.0001)))
  #abgrDiameterFromHeightMixed$ruarkAbatB2 = fit_nlme("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aat * basalAreaTaller + b2ra) * (height - 1.37)), abgr2022, # 10x10 RMSE 7.00, AIC 923
  #                                                 fixedFormula = a1 + b1 + b2 + b2aat ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 2.1, b1 = 0.8, b2 = 0.004, b2aat = -0.0001)))
  abgrDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aat * basalAreaTaller + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022,
                                                         fixedFormula = a1 + b1 + b2 + b2aat + b2tw ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.01, b2aat = -0.0001, b2tw = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2aat * basalAreaTaller) * (height - 1.37)), abgr2022, 
                                                        fixedFormula = a1 + b1 + b1rh + b2 + b2aat ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 2.2, b1 = 0.8, b1rh = 0.03, b2 = 0.002, b2aat = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$ruarkAbatRelHtPhysio = fit_nlme("Ruark RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2aat * basalAreaTaller + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022,
                                                              fixedFormula = a1 + b1 + b2 + b2aat + b2tw ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2.2, b1 = 0.8, b1rh = 0, b2 = 0.01, b2aat = -0.0001, b2tw = 0)), distinct = FALSE)
  #abgrDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1*(1 + b1rm)) * exp((b2 + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, # b2, b2tw not significant, 10x10 RMSE 6.72, AIC 917
  #                                                   fixedFormula = a1 + b1 + b2 + b2tw ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 2.1, b1 = 0.9, b2 = 0, b2tw = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1ra) * exp((b2 + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, # 10x10 RMSE 6.93, AIC 932
                                                     fixedFormula = a1 + b1 + b2 + b2tw ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.01, b2tw = 0.001)))
  #abgrDiameterFromHeightMixed$ruarkPhysioA1 = fit_nlme("Ruark physio", dbh ~ (a1 + a1ra)*(height - 1.37)^(b1) * exp((b2 + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, # b1ra more accurate
  #                                                     fixedFormula = a1 + b1 + b2 + b2tw ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.01, b2tw = 0.001)))
  #abgrDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2ra + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022, # false convergence
  #                                                   fixedFormula = a1 + b1 + b2 + b2tw ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.01, b2tw = 0.001)))
  abgrDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2ra) * (height - 1.37)), abgr2022, # b1rh, b2 not significant, 10x10 RMSE 6.76, AIC 928
                                                    fixedFormula = a1 + b1 + b1rh + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 2.2, b1 = 0.8, b1rh = 0.05, b2 = 0.002)), distinct = FALSE)
  #abgrDiameterFromHeightMixed$ruarkRelHtB1 = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight + b1ra) * exp((b2) * (height - 1.37)), abgr2022, # a1ra max iterations, b1rh, b2 not significant, 10x10 RMSE 7.03, AIC 935
  #                                                    fixedFormula = a1 + b1 + b1rh + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2.2, b1 = 0.8, b1rh = 0.05, b2 = 0.002)), distinct = FALSE)
  #abgrDiameterFromHeightMixed$ruarkRelHtB1 = fit_nlme("Ruark RelHt", dbh ~ (a1)*(height - 1.37)^((b1 + b1rh * relativeHeight)*(1 + b1rm)) * exp((b2) * (height - 1.37)), abgr2022, # 10x10 RMSE 7.00, AIC 946
  #                                                    fixedFormula = a1 + b1 + b1rh + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 2.2, b1 = 0.9, b1rh = 0.1, b2 = -0.005)))
  abgrDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1)*(height - 1.37)^(b1 + b1rh * relativeHeight) * exp((b2 + b2rh * relativeHeight + b2tw * topographicWetnessFD8f) * (height - 1.37)), abgr2022,
                                                          fixedFormula = a1 + a1ac + b1 + b2 + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 2.2, b1 = 0.8, b1rh = 0.03, b2 = 0.002, b2tw = 0.001)), distinct = FALSE)
  #abgrDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/((Ha + Hara)^b1 - 1.37^b1), 0.9999)), abgr2022, # fixed form failed to fit
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Hara ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90)))
  abgrDiameterFromHeightMixed$schnute = create_fit_statistics(name = "Schnute inverse", fitting = "nlme")
  #abgrDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1 * (height - 1.37)^(b1 + b1ra)*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^(b4), abgr2022,
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve
  #                                                    start = list(fixed = c(a1 = 3, b1 = 0.7, b2 = 0.1, b3 = 0.03, b4 = 0.3)))
  #abgrDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1 * (height - 1.37)^(b1)*(exp(b2*(standTreesPerHectare/topHeight)^(b3)*(height - 1.37)) - 1)^(b4*(1 + b4rm)), abgr2022, # a1rm, b2rm, b3rm, b4rm singularity in backsolve, b1rm stephalving
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b4rm ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 3, b1 = 0.7, b2 = 0.1, b3 = 0.03, b4 = 0.3)))
  abgrDiameterFromHeightMixed$sharmaParton = create_fit_statistics(name = "modified Sharma-Parton", fitting = "nlme")
  abgrDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), abgr2022, # a1p, b1p not significant, 10x10 RMSE 6.73, AIC 935, sometimes job false convergence
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 2.3, b1 = 0.8, b2 = 0.05)))
  #abgrDiameterFromHeightMixed$sibbesenReplaceB1 = fit_nlme("Sibbesen replace", dbh ~ (a1)*(height - 1.37)^((b1)*(1 + b1rm)*(height - 1.37)^(b2)), abgr2022, # a1p, b1p not significant, b2rm max iterations, 10x10 RMSE 6.78, AIC 935
  #                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 2.1, b1 = 0.7, b2 = 0.05)))
  #abgrDiameterFromHeightMixed$sibbesenReplaceB2 = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2p * isPlantation + b2ra)), abgr2022, # 10x10 RMSE 6.85, AIC 924
  #                                                       fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)))
  #abgrDiameterFromHeightMixed$sibbesenReplaceB1 = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1p * isPlantation)*(height - 1.37)^((b1 + b1ra)*(height - 1.37)^(b2 + b2p * isPlantation)), abgr2022, # b2ra more accurate
  #                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                         start = list(fixed = c(a1 = 1.6, a1p = 0.4, b1 = 0.7, b2 = 0.1, b2p = -0.02)))
  abgrDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), abgr2022, 
                                                             fixedFormula = a1 + a1aba  + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 2.2, a1aba = 0, b1 = 0.8, b2 = 0.05)))
  abgrDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), abgr2022,
                                                                   fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                   start = list(fixed = c(a1 = 2.2, a1aba = 0, b1 = 0.8, b2 = 0.05)), distinct = FALSE)
  abgrDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), abgr2022, 
                                                                  fixedFormula = a1 + a1aba + b1 + b2 + b2rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                  start = list(fixed = c(a1 = 2.2, a1aba = 0, b1 = 0.8, b2 = 0.04, b2rh = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysio = fit_nlme("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tw * topographicWetnessFD8f + a1ra)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), abgr2022,
                                                                        fixedFormula = a1 + a1p + a1aba + a1aat + a1ac + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                        start = list(fixed = c(a1 = 2.2, a1aba = 0, a1tw = 0, b1 = 0.8, b2 = 0.04, b2rh = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2ra)), abgr2022, # 10x10 RMSE 6.67, AIC 922
                                                               fixedFormula = a1 + a1tw + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 1.8, a1tw = 0.06, b1 = 0.8, b2 = 0.04)))
  #abgrDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(1 + a1rm)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), abgr2022, # 10x10 RMSE 6.82, AIC 921
  #                                                             fixedFormula = a1 + a1tw + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 1.9, a1tw = 0.06, b1 = 0.8, b2 = 0.05)), distinct = FALSE) # a1tw not significant
  #abgrDiameterFromHeightMixed$sibbesenReplacePhysioA1 = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1ra + a1tw * topographicWetnessFD8f)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), abgr2022, # false convergence
  #                                                             fixedFormula = a1 + a1tw + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 1.8, a1tw = 0.06, b1 = 0.8, b2 = 0.04)))
  #abgrDiameterFromHeightMixed$sibbesenReplacePhysioB1 = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(height - 1.37)^((b1 + b1ra)*(height - 1.37)^(b2)), abgr2022, # b2ra more accurate
  #                                                             fixedFormula = a1 + a1tw + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 1.8, a1tw = 0.06, b1 = 0.8, b2 = 0.04)))
  abgrDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2rh * relativeHeight + b2ra)), abgr2022,
                                                              fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.04, b2rh = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2rh * relativeHeight + b2ra)), abgr2022,
                                                                    fixedFormula = a1 + a1tw + b1 + b2 + b2rh ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 2.2, a1tw = 0.06, b1 = 0.8, b2 = 0.04, b2rh = 0)), distinct = FALSE)
  abgrDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2*(1 + b2rm)), abgr2022, # a1rm job singularity in backsolve, b1rm job step halving, 10x10 RMSE 6.93, AIC 953
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = -210, b1 = 0.01, b2 = 0.9)))
  #abgrDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2 + b2ra), abgr2022, # a1ra max iterations, b1ra singularity in backsolve, 10x10 RMSE 7.12, AIC 950
  #                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = -250, b1 = 0.01, b2 = 0.9)))
  #to_fixed_coeffficients(abgrDiameterFromHeightMixed$weibull)
  #print(to_parameter_confidence_intervals(abgrDiameterFromHeightMixed$weibull), n = 40)
  #abgrDiameterFromHeightMixed$weibull$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  abgrDiameterFromHeightMixed$gam = fit_gam("REML GAM", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, bs = "ts", k = 3), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # collapses to linear
  abgrDiameterFromHeightMixed$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE)
  abgrDiameterFromHeightMixed$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE)
  abgrDiameterFromHeightMixed$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE)
  abgrDiameterFromHeightMixed$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, terrainRoughness, bs = "ts", k = 10), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE) # collapses to linear
  abgrDiameterFromHeightMixed$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, relativeHeight, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE)
  abgrDiameterFromHeightMixed$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, relativeHeight, terrainRoughness, bs = "ts", k = 11), random = list(stand = ~1, plot = ~1), data = abgr2022, distinct = FALSE)
  #lapply(abgrDiameterFromHeightMixed$gamRelHtPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(abgrDiameterFromHeightMixed$gamRelHtPhysio$fit, function(fit) { return(summary(fit$gam)) })
  #abgrDiameterFromHeightMixed$gam$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  saveRDS(abgrDiameterFromHeightMixed, paste0("trees/height-diameter/data/ABGR DBH mixed ", get_cross_validation_data_suffix(), ".Rds"))
}

if (abgrOptions$fitHeightPrimary & abgrOptions$fitHeightMixed & abgrOptions$fitDbhPrimary & abgrOptions$fitDbhMixed)
{
  # assemble parameter and results tibbles using incremental load and reduce to fit in memory
  # height primary + nlrob + gsl_nls + dbh primary + nlrob = 120 GB DDR at 10x10 cross validation with models dropped
  # Therefore, batch formation height and diameter results.
  if (exists("abgrHeightFromDiameter") == FALSE) 
  { 
    abgrHeightFromDiameter = readRDS(paste0("trees/height-diameter/data/ABGR height ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("abgrHeightFromDiameterMixed") == FALSE) 
  { 
    abgrHeightFromDiameterMixed = readRDS(paste0("trees/height-diameter/data/ABGR height mixed ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("abgrDiameterFromHeight") == FALSE)
  { 
    abgrDiameterFromHeight = readRDS(paste0("trees/height-diameter/data/ABGR DBH ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("abgrDiameterFromHeightMixed") == FALSE) 
  { 
    abgrDiameterFromHeightMixed = readRDS(paste0("trees/height-diameter/data/ABGR DBH mixed ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  
  abgrCoeffcients = bind_rows(bind_rows(bind_rows(lapply(abgrHeightFromDiameter, unnest_coefficients)),
                                         bind_rows(lapply(abgrHeightFromDiameterMixed, unnest_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(abgrDiameterFromHeight, unnest_coefficients)),
                                         bind_rows(lapply(abgrDiameterFromHeightMixed, unnest_coefficients))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "ABGR")
  abgrFitStats = bind_rows(bind_rows(bind_rows(lapply(abgrHeightFromDiameter, unnest_validation_statistics)),
                                    bind_rows(lapply(abgrHeightFromDiameterMixed, unnest_validation_statistics))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(abgrDiameterFromHeight, unnest_validation_statistics)),
                                    bind_rows(lapply(abgrDiameterFromHeightMixed, unnest_validation_statistics))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "ABGR")

  check_plot_fit_stats(abgrFitStats)
  saveRDS(abgrCoeffcients, paste0("trees/height-diameter/data/ABGR coefficients ", get_cross_validation_data_suffix(), ".Rds"))
  saveRDS(abgrFitStats, paste0("trees/height-diameter/data/ABGR fit stats ", get_cross_validation_data_suffix(), ".Rds"))
}


## preferred forms identified (results.R, Figure 8 and trees.R)
if (abgrOptions$recalcPreferredModels)
{
  # 10x10 model selections
  # responseVariable species       folds  maeName                  rmseName                 aicName          nseName                  isGam isRandomForest maeFitting rmseFitting aicFitting nseFitting
  # height           grand fir     10x10  random forest            Chapman-Richards         Hossfeld IV      Chapman-Richards             0              1 ranger     gsl_nls     gsl_nls    gsl_nls
  #                                2x250  Chapman-Richards         Chapman-Richards         Chapman-Richards Chapman-Richards             0              0 gsl_nls    gsl_nls     gsl_nls    gsl_nls  
  # DBH              grand fir     10x10  Schnute inverse          linear                   Weibull inverse  linear                       0              0 gsl_nls    lm          gsl_nls    lm      
  #                                2x250  Schnute inverse          linear                   Sibbesen replace linear                       0              0 gsl_nls    lm          gsl_nls    lm
  abgrHeightFromDiameterPreferred = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), abgr2022, 
                                                                       start = list(a1 = 64, a1p = -6, b1 = -0.02, b2 = 1.3, b2p = -0.1), folds = 1, repetitions = 1))
  abgrHeightFromDiameterPreferred$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + b2ba * standBasalAreaPerHectare + b2bal / (b2balr + basalAreaLarger)), abgr2022, 
                                                                   start = list(a1 = 70, b1 = -0.02, b2 = 1.2, b2ba = 0.004, b2bal = 3.6, b2balr = 3.6), folds = 1, repetitions = 1)
  abgrHeightFromDiameterPreferred$gam = fit_gam("REML GAM", height ~ I(0.6539 * dbh^0.9757) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = abgr2022, 
                                                threads = 2, folds = 1, repetitions = 1)
  abgrHeightFromDiameterPreferred$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), abgr2022, 
                                                         start = list(a1 = 80, a2 = 200, b1 = -1.3), folds = 1, repetitions = 1)
  abgrHeightFromDiameterPreferred$randomForest = fit_ranger("random forest", c("height", "dbh"), abgr2022, mtry = 1, minNodeSize = 15, sampleFraction = 0.220, folds = 1, repetitions = 1)
  #abgrHeightFromDiameterPreferred$randomForestBalPhysioRelDbh = fit_ranger("random forest RelDbh BA+L physio", c("height", "dbh", "relativeDiameter", "basalAreaLarger", "standQmd", "topHeight", "standBasalAreaPerHectare", "terrainRoughness", "slope"),
  #                                                                         abgr2022, mtry = 4, minNodeSize = 2, sampleFraction = 0.513, folds = 1, repetitions = 1)
  abgrHeightFromDiameterPreferred$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1*(1 - exp(b1*dbh^b2)), abgr2022, 
                                                        start = list(a1 = 58, b1 = -0.01, b2 = 1.2), folds = 1, repetitions = 1)
  saveRDS(abgrHeightFromDiameterPreferred, "trees/height-diameter/data/abgr preferred height models.Rds")
  
  abgrDiameterFromHeightPreferred = list(linear = fit_lm("linear", dbh ~ I(height - 1.37), abgr2022, folds = 1, repetitions = 1))
  abgrDiameterFromHeightPreferred$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aat * basalAreaTaller), abgr2022, start = list(a1 = 2.1, b1aat = -0.001, b1 = 0.9), folds = 1, repetitions = 1)
  #abgrDiameterFromHeightPreferred = list(gam = fit_gam("REML GAM", dbh ~ I(1.9505 * (height - 1.37)^0.9556) + s(height, bs = "ts", k = 3), data = abgr2022, threads = 2, folds = 1, repetitions = 1))
  #abgrDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), abgr2022, start = list(a1 = 2.2, b1 = 0.8, b2 = 0.007), folds = 1, repetitions = 1)
  #abgrDiameterFromHeightPreferred$ruarkMixed = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1 + b1ra) * exp((b2) * (height - 1.37)), abgr2022,
  #                                                      fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 2.2, b1 = 0.8, b2 = 0.007)), folds = 1, repetitions = 1)
  abgrDiameterFromHeightPreferred$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1), 0.9999)), abgr2022, 
                                                        start = list(a1 = 0.000001, a2 = 0.0003, b1 = 0.9, Ha = 130), folds = 1, repetitions = 1)
  #abgrDiameterFromHeightPreferred$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), abgr2022, start = list(a1 = 2.2, b1 = 0.8, b2 = 0.05), folds = 1, repetitions = 1)
  #abgrDiameterFromHeightPreferred$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ (a1 + a1tw * topographicWetnessFD8f)*(height - 1.37)^((b1)*(height - 1.37)^b2), abgr2022, start = list(a1 = 1.8, a1tw = 0.06, b1 = 0.8, b2 = 0.04), folds = 1, repetitions = 1)
  abgrDiameterFromHeightPreferred$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, abgr2022, 
                                                        start = list(a1 = -250, b1 = 0.01, b2 = 0.9), folds = 1, repetitions = 1)
  saveRDS(abgrDiameterFromHeightPreferred, "trees/height-diameter/data/abgr preferred diameter models.Rds")
}


## error residual distributions
if (htDiaOptions$includeInvestigatory)
{
  abgrHeightWeibull = gsl_nls(height ~ 1.37 + a1*(1 - exp(b1*dbh^b2)), abgr2022, start = list(a1 = 58, b1 = -0.01, b2 = 1.2))
  ggplot() +
    geom_point(aes(x = relativeDiameter, y = residuals(abgrHeightWeibull)), abgr2022, alpha = 0.1, shape = 16) +
    geom_smooth(aes(x = relativeDiameter, y = (height - residuals(abgrHeightWeibull)) / height), abgr2022, formula = y ~ x, method = "loess", span = 0.03) + # GAM and loess both linear
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "relative diameter", y = "scale error")

  abgrDbhLinear = lm(dbh ~ I(height - 1.37), abgr2022)
  ggplot() +
    geom_point(aes(x = relativeHeight, y = (dbh - residuals(abgrDbhLinear)) / dbh), abgr2022, alpha = 0.1, shape = 16) +
    geom_smooth(aes(x = relativeHeight, y = (dbh - residuals(abgrDbhLinear)) / dbh), abgr2022, formula = y ~ x, method = "loess") + # GAM linear, loess suggests 1/rh
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "relative height", y = "scale error")
  ggplot() +
    geom_point(aes(x = topHeight, y = relativeHeight, color = 100 * residuals(abgrDbhLinear) / dbh), abgr2022, alpha = 1, shape = 16) +
    labs(x = bquote(H[100]*", m"), y = "relative height", color = "% error") +
    scico::scale_color_scico(palette = "bam", limits = c(-100, 100))

  ggplot() +
    geom_point(aes(x = height, y = relativeHeight, color = dbh), abgr2022, alpha = 0.3, shape = 16) +
    labs(x = "height, m", y = bquote(H[100]*", m"), color = "DBH, cm") +
    scale_color_viridis_c(limits = c(0, NA))
  ggplot() +
    geom_point(aes(x = topHeight, y = relativeHeight, color = residuals(abgrDbhLinear)), abgr2022, alpha = 0.5, shape = 16) +
    labs(x = bquote(H[100]*", m"), y = "relative height", color = "DBH error, cm") +
    scico::scale_color_scico(palette = "bam", limits = c(-50, 50), trans = scales::pseudo_log_trans(), breaks = c(-50, -10, -5, 0, 5, 10, 50))
  
  abgrPca = princomp(abgr2022 %>% select(dbh, height, relativeHeight, topHeight))
  summary(abgrPca)
  abgrPca$loadings
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  abgrHeightGam = gam(height ~ s(dbh, bs = "ts", by = as.factor(isPlantation), k = 8) + 
                               s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 6) + # sort of linearish
                               #s(standBasalAreaPerHectare, bs = "ts", k = 4) + # not significant
                               s(basalAreaLarger, bs = "ts", k = 12), # plantation not significant, piecewise linear, effects 0-10 m²/ha
                               #s(elevation, bs = "ts", k = 6) + # not significant
                               #s(slope, bs = "ts", k = 6) + # not significant
                               #s(aspect, bs = "ts", k = 8) + # not significant
                               #s(terrainRoughness, bs = "ts", k = 6) + # not significant
                               #s(topographicWetnessFD8f, bs = "ts", k = 6) + # not significant
                      data = abgr2022, method = "REML", select = TRUE, weights = dbhWeight)
  k.check(abgrHeightGam)
  summary(abgrHeightGam)
  plot_gam_smooths(abgrHeightGam)
  
  #ggplot() + geom_line(aes(x = seq(0, 60), y = -1/(seq(0, 60) + 1)))
  
  abgrDbhGam = gam(dbh ~ I(1.9505 * (height - 1.37)^0.9556) + 
                         s(height, bs = "ts", by = as.factor(isPlantation), k = 7) +
                         s(relativeHeight, bs = "ts", k = 6) + # plantation not significant, roughly linear
                         #s(bootstrapStandBasalAreaPerHectare, bs = "ts", k = 3) + # not significant
                         #s(basalAreaTaller, bs = "ts", k = 3) + # not significant
                         #s(elevation, bs = "ts", k = 9) + # not significant
                         #s(slope, bs = "ts", k = 11) + # not significant
                         #s(aspect, bs = "ts", k = 25) + # not significant
                         #s(terrainRoughness, bs = "ts", k = 15) + # not significant
                         s(topographicWetnessFD8f, bs = "ts", k = 3), # significant but basically linear
                   data = abgr2022, method = "REML", select = TRUE, weights = heightWeight)
  k.check(abgrDbhGam)
  summary(abgrDbhGam)
  plot_gam_smooths(abgrDbhGam)

  #abgrHeightGamBase = gam(height ~ s(dbh, bs = "ts", k = 3), data = abgr2022, method = "REML", select = TRUE, weights = dbhWeight)
  #abgrHeightGamBaseGuide = gam(height ~ I(0.6539 * dbh^0.9757) + s(dbh, bs = "ts", k = 3), data = abgr2022, method = "REML", select = TRUE, weights = dbhWeight)
  #abgrHeightGamBaseGuideOffset = gam(height ~ offset(breastHeight) + I(0.6539 * dbh^0.9757) + s(dbh, bs = "ts", k = 3), data = abgr2022, method = "REML", select = TRUE, weights = dbhWeight)
  #max(abs(predict(abgrHeightGamBaseGuide) - predict(abgrHeightGamBaseGuideOffset)))
  #ggplot() +
  #  geom_point(aes(x = dbh, y = height), abgr2022, alpha = 0.1, shape = 16) +
  #  geom_line(aes(x = dbh, y = predict(abgrHeightGamBase), color = "base"), abgr2022) +
  #  geom_line(aes(x = dbh, y = predict(abgrHeightGamBaseGuide), color = "base guide"), abgr2022) +
  #  geom_line(aes(x = dbh, y = predict(abgrHeightGamBaseGuideOffset), color = "base offset"), abgr2022) +
  #  labs(x = "DBH, cm", y = "height, m", color = NULL)
}
