# load libraries, functions, and other2022 from setup.R
# otherOptions from other job.R

## other species height regressions
if (otherOptions$fitHeight)
{
  otherHeightFromDiameter = list(linear = fit_lm("linear", height ~ dbh, other2022)) # a1p not significant in 50% of fits
  otherHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ dbh + I(dbh^2), other2022) # a1p, a2p not significant in 50% of fits

  otherHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1 * dbh))^b2, other2022, 
                                                        start = list(a1 = 30, b1 = -0.02, b2 = 0.84)) # a1p, b1p, b2p (25%) not significant
  # linear BA+L
  #otherHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                         start = list(a1 = 30, a1bal = 0, b1 = -0.02, b2 = 0.84), distinct = FALSE) # a1p, a1ba, a2p, a1bal, a1balp, b1p, b2p, exponent BAL not significant
  #otherHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1tw * topographicWetnessFD8f) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                               start = list(a1 = 30, a1bal = 0, a1tw = 0, b1 = -0.02, b2 = 0.8), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  #otherHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                                     start = list(a1 = 80, a1bal = 0.3, a1rd = 0, a1e = 0, b1 = -0.01, b2 = 0.8), distinct = FALSE) # by propagation
  #otherHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter) * (1 - exp(b1 * dbh))^b2, other2022, 
  #                                                               start = list(a1 = 30, a1bal = 0, a1rd = 0, b1 = -0.02, b2 = 0.84), distinct = FALSE) # singular gradient with a1bal, a1rd, a1rdp not significant, exponents NaN-Inf or singular gradient
  # reciprocal BA+L
  otherHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger)) * (1 - exp((b1) * dbh))^(b2), other2022, 
                                                           start = list(a1 = 30, a1bal = 0, a1balr = 20, b1 = -0.02, b2 = 0.84), distinct = FALSE) # a1ba, a1ba+r, a1bal+r, b1ba, b1ba+r, b2ba, b2ba+r, b2bal+r not significant, b1bal+r NaN-Inf
  otherHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger) + a1tw * topographicWetnessFD8f) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                 start = list(a1 = 30, a1bal = 0, a1balr = 20, a1tw = 0, b1 = -0.02, b2 = 0.8), distinct = FALSE) # by propagation, a1tw+e not significant
  otherHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger) + a1rd * relativeDiameter + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                       start = list(a1 = 80, a1bal = 0, a1balr = 20, a1rd = 0, a1e = 0, b1 = -0.01, b2 = 0.8), distinct = FALSE) # by propagation
  otherHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger) + a1rd * relativeDiameter) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                 start = list(a1 = 30, a1bal = 0, a1balr = 20, a1rd = 0, b1 = -0.02, b2 = 0.84), distinct = FALSE) # by propagation
  otherHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp(b1 * dbh))^(b2 + b2ac * cos(pi/180 * aspect)), other2022, 
                                                              start = list(a1 = 30, b1 = -0.02, b2 = 0.81, b2ac = -0.08), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant, though b2as + b2ac intermittently test significant
  otherHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp(b1 * dbh))^b2, other2022, 
                                                              start = list(a1 = 30, a1rd = 0, b1 = -0.02, b2 = 0.8), distinct = FALSE) # a1rd, a1rdp, b1rd, b2rd not significant, signs of a1-b1 evaporation
  otherHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                    start = list(a1 = 30, a1rd = 0, a1e = 0, b1 = -0.02, b2 = 0.81), distinct = FALSE) # by propagation
  otherHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1*dbh / (1 + dbh)^b1, other2022, 
                                               start = list(a1 = 2, b1 = 0.4))
  otherHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), other2022, 
                                                 start = list(a1 = 44, a2 = 33, b1 = -0.9)) # prone to singular gradient
  otherHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1 * dbh^b2), other2022, 
                                             start = list(a1 = 300, b1 = -5, b2 = -0.1)) # a1p, b1p, b2p not significant, form not well posed (should probably be a1 * (exp() - 1))
  otherHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), other2022, 
                                                        start = list(a1 = 50, a2 = 40, b1 = 0.9)) # a1p, a2p, b1p not significant
  otherHeightFromDiameter$michaelisMentenBal = fit_gsl_nls("Michaelis-Menten BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger)*dbh^(b1) / (a2 + dbh^(b1)), other2022, 
                                                           start = list(a1 = 50, a1bal = 0, a2 = 40, b1 = 0.9), distinct = FALSE) # a1ba, a1bar, a1bal, a1bal+r, a2ba, a2bar, a2bal, a2bal+r, b1ba, b1bar, b1bal, b1bal+r not significant
  otherHeightFromDiameter$michaelisMentenBalRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1bal * basalAreaLarger)*dbh^(b1) / (a2 + dbh^(b1)), other2022, 
                                                                 start = list(a1 = 50, a1bal = 0, a1rd= 0, a2 = 40, b1 = 0.9), distinct = FALSE) # by propagation
  otherHeightFromDiameter$michaelisMentenRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*dbh^(b1) / (a2 + dbh^(b1)), other2022, 
                                                              start = list(a1 = 50, a1rd = 0, a2 = 40, b1 = 0.9), distinct = FALSE) # a1rd, a2rd, b1rd not significant
  otherHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, other2022, 
                                              start = list(a1 = 1.5, b1 = 0.7)) # a1p, b1p not significant
  otherHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + a2*dbh + a3), other2022, 
                                               start = list(a1 = -0.01, a2 = 2.5, a3 = -5)) # a1p, a2p, a3p not significant
  otherHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2)), other2022, 
                                                  start = list(a1 = 30, b1 = -11, b2 = 3)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), other2022, 
                                                  start = list(Ha = 28, d = 0.2, kU = 0.025)) # NaN-inf @ 2x50, 2x250, 2x500
  #otherHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^b4, other2022, 
  #                                                   start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant, 10x10MAE 2.18, AIC 762
  otherHeightFromDiameter$sharmaPartonB4 = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh^b4d)), other2022, 
                                                       start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4d = 1.0)) # b4p not significant, 10x10 MAE 2.13, AIC 719
  # fixed BAL
  #otherHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, other2022, 
  #                                                      start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant
  #otherHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, 
  #                                                            start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0.1)) # { a1, b1, b2, b3, b4 } x { e, s, ac, tr, tw }, { a1, b1, b2, b3 } x as not significant
  #otherHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, 
  #                                                                  start = list(a1 = 9, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.06, b4 = 1.0, b4as = -0.1)) # a1rd, b1rd, b4rd not significant, b2rd often significant
  #otherHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, 
  #                                                            start = list(a1 = 10, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.07, b4 = 1.0)) # a1rd, b1rd, b2rd, b4rd not significant
  # reciprocal BAL, not b4 as none are distinct
  otherHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4 + b4bal / (b4balr + basalAreaLarger)), other2022, 
                                                           start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4bal = 0, b4balr = 20), distinct = FALSE) # a1bal+r, b1bal+r, b3bal+r, b4bal+r not significant, b2bal+r NaN-Inf
  otherHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, 
                                                              start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0.1), distinct = FALSE) # by propagation
  otherHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, 
                                                                    start = list(a1 = 9, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.06, b4 = 1.0, b4as = -0.1), distinct = FALSE) # by propagation
  otherHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, 
                                                              start = list(a1 = 10, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.07, b4 = 1.0), distinct = FALSE) # by propagation
  otherHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + (a1 + a1tw * topographicWetnessFD8f)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/(standBasalAreaPerHectare))^(b3)*dbh))^(b4), other2022, 
                                                           start = list(a1 = 5, a1tw = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0), distinct = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw }, a1tw+e not significant
  otherHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + (a1)*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, 
                                                           start = list(a1 = 8, b1 = 0.6, b2 = -0.01, b3 = 0.4, b3rd = -0.05, b4 = 1.0), distinct = FALSE) # a1rd not significant
  otherHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1e * elevation)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, 
                                                                 start = list(a1 = 10, a1rd = 0, a1e = 0, b1 = 0.5, b2 = -0.007, b3 = 0.4, b4 = 1.1), distinct = FALSE) # by propagation
  otherHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, other2022, 
                                                    start = list(a1 = 25, b1 = 0, b2 = -0.01, b3 = 0.2, b4 = 0.9)) # a1p, b1, b1p, b2p, b3, b3p, b4p not significant, NaN-inf @ 2x50, 10x10 RMSE 3.69, AIC 767, 2x250, 2x500, 5x200 NaN-inf
  #otherHeightFromDiameter$sharmaZhangB4 = fit_gsl_nls("Sharma-Zhang b₄", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^b3*dbh^b4d)), other2022, 
  #                                                    start = list(a1 = 25, b1 = 0, b2 = -0.001, b3 = 0.2, b4d = 0.9), distinct = FALSE) # b1, b3 not significant -> collapses to Chapman-Richards, 10x10 RMSE 3.66, AIC 775
  otherHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2bal * basalAreaLarger)*standTreesPerHectare^b3*dbh))^b4, other2022, 
                                                       start = list(a1 = 25, b1 = 0.05, b2 = -0.01, b2bal = 0, b3 = 0.2, b4 = 0.9), distinct = FALSE) # a1ba, a1bar, a1bal+r, a2p, a1bal, b1ba, b1bar, b1bal, b1bal+r, b2bar, b2bal, b2bal+r, b3ba-a1, b3bal, b4ba, b4bar, b4bal, b4bal+r not significant, a1ba, a1bal+r, b1bal+r, b2ba, b3bar, b3bal+r NaN-inf, 2x50, 2x250 NaN-inf
  otherHeightFromDiameter$sharmaZhangBalRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2 + b2bal * basalAreaLarger)*standTreesPerHectare^(b3)*dbh))^(b4 + b4rd * relativeDiameter), other2022, 
                                                             start = list(a1 = 25, b1 = 0, b2 = -0.01, b2bal = 0, b3 = 0.2, b4 = 0.9, b4rd = 0), distinct = FALSE) # by propagation, 2x50, 2x250, 5x100 NaN-inf
  otherHeightFromDiameter$sharmaZhangRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3)*dbh))^(b4 + b4rd * relativeDiameter), other2022, 
                                                          start = list(a1 = 25, b1 = 0, b2 = -0.01, b3 = 0.2, b4 = 0.9, b4rd = 0), distinct = FALSE) # a1rd, b1rd, b2rd, b3rd NaN-inf, b4rd not significant, 2x50, 2x250 NaN-inf
  otherHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1 * dbh^b2), other2022, 
                                                 start = list(a1 = 1.2, b1 = 1.0, b2 = -0.08)) # a1p, b1p, b2p not significant
  otherHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + a1 * (1 - exp(b1 * dbh^b2)), other2022, 
                                                start = list(a1 = 30, b1 = -0.05, b2 = 0.9)) # a1p, b1p, b2p not significant
  #otherHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + a1 * (1 - exp(b1 * dbh^(b2 + b2bal * basalAreaLarger))), other2022, 
  #                                                 start = list(a1 = 40, b1 = -0.04, b2 = 0.9, b2bal = 0), distinct = FALSE) # a1ba, a1bal, b1ba, b1bal, b2ba, b2bal not significant
  otherHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1bal / standBasalAreaPerHectare) * dbh^(b2))), other2022, 
                                                   start = list(a1 = 40, b1 = -0.04, b1bal = 0, b2 = 0.9), distinct = FALSE) # { a1, b1, b2 } x { ba, ba+r, bal+r } not significant
  otherHeightFromDiameter$weibullBalRelDbh = fit_gsl_nls("Weibull RelDbh BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter + b1bal / standBasalAreaPerHectare) * dbh^(b2))), other2022, 
                                                         start = list(a1 = 40, b1 = -0.04, b1rd = 0, b1bal = 0, b2 = 0.9), distinct = FALSE) # b1rd-b1bal not mutually significant
  otherHeightFromDiameter$weibullRelDbh = fit_gsl_nls("Weibull RelDbh", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter) * dbh^b2)), other2022, 
                                                      start = list(a1 = 30, b1 = -0.05, b1rd = 0.002, b2 = 0.9)) # a1rd, b2rd not significant
  #print(to_parameter_confidence_intervals(otherHeightFromDiameter$weibullRelDbh), n = 24)
  #to_fixed_coeffficients(otherHeightFromDiameter$weibullRelDbh)
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
  #otherHeightFromDiameter$gamSibbesen = fit_gam("REML GAM", height ~ I(1.24408 * dbh^(1.01570 * dbh^-0.08286)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022, distinct = FALSE) # 10x10 MAE 2.25, AIC 673 but collapses to linear
  #otherHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, bs = "ts", k = 4), data = other2022, distinct = FALSE) # k = 4 minimum collapses to linear, plantation not significant
  #otherHeightFromDiameter$gamBal = fit_gam("REML GAM BAL", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, bs = "ts", k = 4), data = other2022, distinct = FALSE) # k = 4 minimum collapses to linear, plantation not significant
  otherHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", k = 11), data = other2022, distinct = FALSE) # k = 11 minimum > ~0, plantation not significant, 10x0 AIC 
  otherHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, bs = "ts", by = isPlantation, k = 11), data = other2022, distinct = FALSE) # by propagation
  otherHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), data = other2022, distinct = FALSE) # by propagation
  otherHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 11), data = other2022) # k = 16 minimum possible -> drop basal area -> k = 11 minimum, significant but not by much, plantation not reliably significant, 10x10 AIC 790
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
  otherHeightFromDiameter$randomForestBalPhysioRelDbh = fit_ranger("random forest RelDbh BA+L physio", c("height", "dbh", "relativeDiameter", "standQmd", "topHeight", "standBasalAreaPerHectare", "basalAreaLarger", "standTreesPerHectare", "slope", "elevation", "terrainRoughness", "topographicWetnessFD8f", "aspect", "isPlantation"), # from setup.R VSURF
                                                                   other2022, mtry = 3, minNodeSize = 2, sampleFraction = 0.839) # from setup.R tuneRanger @ 500 trees
 
  saveRDS(otherHeightFromDiameter, paste0("trees/height-diameter/data/other height ", get_cross_validation_data_suffix(), ".Rds"))
}

if (otherOptions$fitHeightMixed)
{
  otherHeightFromDiameterMixed = list(linear = fit_lme("linear", fixedFormula = height ~ dbh, other2022, randomFormula = ~1|stand/plot))
  otherHeightFromDiameterMixed$parabolic = fit_lme("parabolic", fixedFormula = height ~ dbh + I(dbh^2), other2022, randomFormula = ~1|stand/plot)

  #otherHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1 * dbh))^(b2 + b2ra), other2022, # a1ra, b1ra, b2ra singularity in backsolve
  #                                                        fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve
  #                                                        start = list(fixed = c(a1 = 32, b1 = -0.5, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 1E-4))
  #otherHeightFromDiameterMixed$chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + a1*(1 + a1rm)*(1 - exp(b1 * dbh))^(b2), other2022, # a1rm, b1rm, b2rm singularity in backsolve
  #                                                        fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 32, b1 = -0.5, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 1E-4))
  otherHeightFromDiameterMixed$chapmanRichards = create_fit_statistics(name = "Chapman-Richards", fitting = "nlme")
  otherHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp((b1) * dbh))^(b2 + b2ra), other2022, # a1ra, b1ra, b24 singularity in backsolve
                                                             fixedFormula = a1 + a1bal + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot, # |stand singularity in backsolve
                                                             start = list(fixed = c(a1 = 30, a1bal = 0, b1 = -0.02, b2 = 0.84)), distinct = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1ra + a1bal * basalAreaLarger + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                   fixedFormula = a1 + a1bal + a1e + b1 + b2 ~ 1, randomFormula = a1ra ~ 1,
                                                                   start = list(fixed = c(a1 = 80, a1bal = 0.3, a1e = 0.005, b1 = -0.005, b2 = 0.8)), distinct = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger) + a1rd * relativeDiameter + a1e * elevation + a1ra) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                         fixedFormula = a1 + a1bal + a1balr + a1rd + a1e + b1 + b2 ~ 1, randomFormula = a1ra ~ 1,
                                                                         start = list(a1 = 80, a1bal = 0, a1balr = 20, a1rd = 0, a1e = 0, b1 = -0.01, b2 = 0.8), distinct = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsBalRelDbh = fit_nlme("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1bal / (a1balr + basalAreaLarger) + a1rd * relativeDiameter + a1ra) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                   fixedFormula = a1 + a1bal + a1balr + a1rd + b1 + b2 ~ 1, randomFormula = a1ra ~ 1,
                                                                   start = list(a1 = 30, a1bal = 0, a1balr = 20, a1rd = 0, b1 = -0.02, b2 = 0.84), distinct = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1ra + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                fixedFormula = a1 + a1e + b1 + b2 ~ 1, randomFormula = a1ra ~ 1,
                                                                start = list(fixed = c(a1 = 70, a1e = -0.01, b1 = -0.006, b2 = 0.81)), distinct = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsRelDbh = fit_nlme("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1ra)*(1 - exp(b1 * dbh))^b2, other2022, 
                                                                fixedFormula = a1 + a1rd + b1 + b2 ~ 1, randomFormula = a1ra ~ 1,
                                                                start = list(a1 = 30, a1rd = 0, b1 = -0.02, b2 = 0.8), distinct = FALSE)
  otherHeightFromDiameterMixed$chapmanRichardsRelDbhPhysio = fit_nlme("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1e * elevation) * (1 - exp(b1 * dbh))^b2, other2022, 
                                                                      fixedFormula = a1 + a1rd + a1e + b1 + b2 ~ 1, randomFormula = a1ra ~ 1,
                                                                      start = list(a1 = 30, a1rd = 0, a1e = 0, b1 = -0.02, b2 = 0.81), distinct = FALSE)
  otherHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1)*dbh / (1 + dbh)^(b1 + b1ra), other2022, # 10x10 MAE 2.30, AIC 742
                                                 fixedFormula = a1 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 1.0, b1 = 0.18)))
  #otherHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1ra)*dbh / (1 + dbh)^(b1), other2022, # 10x10 MAE 2.75, AIC 683
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.0, b1 = 0.18)))
  #otherHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1)*(1 + a1rm)*dbh / (1 + dbh)^(b1), other2022, # b1rm job step halving, 10x10 MAE 2.68, AIC 777
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.9, b1 = 0.4)))
  #otherHeightFromDiameterMixed$curtisA1 = fit_nlme("Curtis", height ~ 1.37 + (a1 + a1ra)*dbh / (1 + dbh)^(b1), other2022,
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 1.0, b1 = 0.18)))
  otherHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^(b1 + b1ra)), other2022, # a1ra, a2ra max iterations, b1ra 2x50 job computationally singular, 10x10 MAE 2.31, AIC 784
                                                   fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 34, b1 = 25, b2 = -0.9))) # 2x50, 2x250, 2x500, 5x100 step halving, computational singularity
  #otherHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + a1*(1 + a1rm) / (1 + a2 * dbh^(b1)), other2022, # a2rm step halving, 10x10 MAE 2.59, AIC 754
  #                                                 fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 55, b1 = 40, b2 = -0.9)))
  #otherHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^(b1*(1 + b1rm))), other2022, # 10x10 MAE 2.43, AIC 744
  #                                                 fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 60, b1 = 60, b2 = -0.9)))
  #otherHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1)*exp((b1 + b1ra) * dbh^(b2)), other2022, # |stand/plot a1ra job singularity in backsolve, b1ra singularity in backsolve, b2ra step halving, |plot b2ra step halving, a1ra tractable but excluded
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot, # |stand, |uniquePlotID step halving
  #                                             start = list(fixed = c(a1 = 100, b1 = -6, b2 = -0.2)))
  otherHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1)*exp((b1) * dbh^(b2*(1 + b2rm))), other2022, # 10x10 MAE 2.33, AIC 689
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 400, b1 = -6, b2 = -0.2)))
  #otherHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1)*(1 + b1rm)*exp((b1) * dbh^(b2)), other2022, # b1rm step halving, 10x10 MAE 2.43, AIC 784
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 200, b1 = -4, b2 = -0.2)))
  #otherHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + a1*dbh^(b1 + b1ra) / (a2 + dbh^(b1 + b1ra)), other2022, # a1ra max iterations, a2ra singularity in backsolve, potential condition number and job step halving with b1ra|stand/plot -> reduced to stand -> job max iterations
  #                                                        fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot, # |uniquePlotID step halving
  #                                                        start = list(fixed = c(a1 = 100, a2 = 100, b1 = 0.8)))
  otherHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + a1*dbh^(b1*(1 + b1rm)) / (a2 + dbh^(b1*(1 + b1rm))), other2022, # 10x10 MAE 2.41, AIC 749
                                                          fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 60, a2 = 55, b1 = 1.0))) # 2x50, 2x500 max iterations
  #otherHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + a1*(1 + a1rm)*dbh^(b1) / (a2 + dbh^(b1)), other2022, # 10x10 MAE 2.53, AIC 751
  #                                                        fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 60, a2 = 55, b1 = 0.9)))
  otherHeightFromDiameterMixed$michaelisMentenBal = fit_nlme("Michaelis-Menten BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger)*dbh^(b1) / (a2 + dbh^(b1)), other2022, 
                                                             fixedFormula = a1 + a1bal + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                             start = list(a1 = 50, a1bal = 0, a2 = 40, b1 = 0.9), distinct = FALSE)
  otherHeightFromDiameterMixed$michaelisMentenBalRelDbh = fit_nlme("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1bal * basalAreaLarger)*dbh^(b1) / (a2 + dbh^(b1)), other2022, 
                                                                   fixedFormula = a1 + a1rd + a1bal + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                                   start = list(a1 = 50, a1bal = 0, a1rd = 0, a2 = 40, b1 = 0.9), distinct = FALSE)
  otherHeightFromDiameterMixed$michaelisMentenRelDbh = fit_nlme("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*dbh^(b1) / (a2 + dbh^(b1)), other2022, 
                                                                fixedFormula = a1 + a1rd + a2 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                                start = list(a1 = 50, a1rd = 0, a2 = 40, b1 = 0.9), distinct = FALSE)
  otherHeightFromDiameterMixed$power = fit_nlme("power", height ~ 1.37 + (a1)*dbh^(b1 + b1ra), other2022, # 10x10 MAE 2.39, AIC 790
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.0, b1 = 0.85)))
  #otherHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1 + a1ra)*dbh^(b1), other2022, # 10x10 MAE 2.67, AIC 712
  #                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 1.0, b1 = 0.85)))
  #otherHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1)*(1 + a1rm)*dbh^(b1), other2022, # 10x10 2.97, AIC 819
  #                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 1.5, b1 = 0.75)))
  #otherHeightFromDiameterMixed$powerB1 = fit_nlme("power", height ~ 1.37 + (a1)*dbh^(b1*(1 + b1rm)), other2022, # 10x10 2.72, AIC 779
  #                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 1.3, b1 = 0.8)))
  #otherHeightFromDiameterMixed$powerA1 = fit_nlme("power", height ~ 1.37 + (a1 + a1ra)*dbh^(b1), other2022,
  #                                              fixedFormula = a1 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 1.0, b1 = 0.85)))
  otherHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / (a1*dbh^2 + (a2 + a2ra)*dbh + a3), other2022, # 10x10 MAE 2.25, AIC 704
                                                 fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 0.025, a2 = 1.08, a3 = -0.01))) # 2x50, 2x250, 2x500 max iterations
  #otherHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / ((a1 + a1ra)*dbh^2 + (a2)*dbh + a3), other2022, # a3ra max iterations, 10x10 MAE 2.35, AIC 783
  #                                               fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.025, a2 = 1.08, a3 = -0.01)))
  #otherHeightFromDiameterMixed$prodan = fit_nlme("Prodan", height ~ 1.37 + dbh^2 / (a1*(1 + a1rm)*dbh^2 + (a2)*dbh + a3), other2022, # a1rm singular precision, a2rm step halving
  #                                               fixedFormula = a1 + a2 + a3 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.025, a2 = 1.08, a3 = -0.01)))
  otherHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1 + b1ra)/(dbh + b2)), other2022, # a1ra less accurate, b1ra 2x50 job step halving, b2ra step halving, 10x10 MAE 2.38, AIC 747
                                                    fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 20, b1 = -10, b2 = 2))) # 2x250 singularity in backsolve, 2x500 max iterations
  #otherHeightFromDiameterMixed$ratkowskyA1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*(1 + a1rm)*exp((b1)/(dbh + b2)), other2022, # 10x10 MAE 2.70, AIC 773
  #                                                    fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 30, b1 = -13, b2 = 3)))
  #otherHeightFromDiameterMixed$ratkowskyB1 = fit_nlme("Ratkowsky", height ~ 1.37 + (a1)*exp((b1*(1 + b1rm))/(dbh + b2)), other2022, # 10x10 MAE 2.75, AIC 884
  #                                                    fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 40, b1 = -35, b2 = 7)))
  #otherHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU + kUra) * dbh)/d^(d/(1 - d))))^(1/(1 - d)), other2022, # Hara, kUra, dra singularity in backsolve
  #                                                  fixedFormula = Ha + d + kU ~ 1, randomFormula = kUra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(Ha = 50, d = 0.2, kU = 0.01)))
  #otherHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", height ~ 1.37 + (Ha) * (1 + ((1.37/(Ha))^(1 - (d*(1 + drm))) - 1) * exp((-kU * dbh)/(d*(1 + drm))^((d*(1 + drm))/(1 - (d*(1 + drm))))))^(1/(1 - (d*(1 + drm)))), other2022, # drm, Harm singularity in backsolve, kUrm max iterations
  #                                                  fixedFormula = Ha + d + kU ~ 1, randomFormula = kUrm ~ 1|stand/plot,
  #                                                  start = list(fixed = c(Ha = 50, d = 0.2, kU = 0.01)))
  otherHeightFromDiameterMixed$richardsW = create_fit_statistics(name = "unified Richards", fitting = "nlme")
  otherHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1 + b1ra)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^b4, other2022, # a1ra, b2ra, b4ra singularity in backsolve, 10x10 MAE 2.13, AIC 717
                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0))) # 2x500 singularity in backsolve
  #otherHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^b4d)), other2022, # a1ra singularity in backsolve, b2ra, b3ra step halving, 10x10 MAE 2.18, AIC 714
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 9, b1 = 0.4, b2 = -0.01, b3 = 0.3, b4d = 0.9)))
  #otherHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^(b4d + b4dra))), other2022, # 10x10 MAE 2.26, AIC 753
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4dra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 30, b1 = 0.2, b2 = -0.01, b3 = 0.3, b4d = 0.9)))
  #otherHeightFromDiameterMixed$sharmaPartonB3 = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3ra)*dbh))^b4, other2022, # 10x10 MAE 2.25, AIC 683
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  #otherHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^(b4d + b4dra))), other2022, # a1, b1, b2, b3 not significant, 10x10 MAE 2.32, AIC 716
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4dra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 50, b1 = 0.4, b2 = -0.01, b3 = 0.3, b4d = 0.9)))
  #otherHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^(b4d*(1 + b4drm)))), other2022, # a1rm singularity in backsolve, b1rm, b3rm computationally singular, b2rm step halving, 10x10 MAE 2.36, AIC 655
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4drm ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 8, b1 = 0.4, b2 = -0.01, b3 = 0.3, b4d = 0.9)))
  #otherHeightFromDiameterMixed$sharmaPartonB4 = fit_nlme("Sharma-Parton", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^(b4d*(1 + b4rm)))), other2022, # a1rm, b1rm, b2rm, b3rm singularity in backsolve, b4rm job max iterations
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4d ~ 1, randomFormula = b4rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 10, b1 = 0.4, b2 = -0.01, b3 = 0.3, b4d = 1.0)))
  #otherHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", height ~ 1.37 + a1*topHeight^(b1)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh))^(b4*(1 + b4rm)), other2022, # a1rm, b1rm, b3rm, b4rm singularity in backsolve, b2rm step halving
  #                                                     fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b4rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)))
  otherHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, other2022, # a1ra, b2ra, b4ra singularity in backsolve, b3ra job max iterations
                                                          fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 12, b1 = 0.3, b2 = -0.02, b3 = 0.3, b4 = 0.9)), distinct = FALSE)
  otherHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, # a1ra, b2ra, b4ra singularity in backsolve, b3ra less accurate
                                                                fixedFormula = a1 + b1 + b2 + b3 + b4 + b4as ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0)), distinct = FALSE)
  #otherHeightFromDiameterMixed$sharmaPartonBalPhysioB3 = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3ra)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022,
  #                                                                fixedFormula = a1 + b1 + b2 + b3 + b4 + b4as ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0)))
  otherHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, # a1ra, b4ra singularity in backsolve, b2ra, b3ra step halving
                                                                      fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 + b4as ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                      start = list(fixed = c(a1 = 9, b1 = 0.4, b2 = -0.02, b3 = 0.3, b3rd = -0.06, b4 = 1.0, b4as = -0.1)), distinct = FALSE)
  otherHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1ra) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, # a1ra, b2ra, b4ra singularity in backsolve, b3ra less accurate
                                                                fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 30, b1 = 0.1, b2 = -0.02, b3 = 0.2, b3rd = -0.04, b4 = 1.0)), distinct = FALSE)
  #otherHeightFromDiameterMixed$sharmaPartonBalRelDbhB3 = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^(b1) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter + b3ra)*dbh))^b4, other2022,
  #                                                                fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 30, b1 = 0.1, b2 = -0.02, b3 = 0.2, b3rd = -0.04, b4 = 1.0)))
  otherHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1ra + a1tw * topographicWetnessFD8f)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, 
                                                             fixedFormula = a1 + a1tw + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 5, a1tw = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0)), distinct = FALSE)
  otherHeightFromDiameterMixed$sharmaPartonRelDbh = fit_nlme("Sharma-Parton RelDbh", height ~ 1.37 + (a1)*topHeight^(b1 + b1ra)*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^b4, other2022, 
                                                             fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                             start = list(a1 = 8, b1 = 0.6, b2 = -0.01, b3 = 0.4, b3rd = -0.05, b4 = 1.0), distinct = FALSE)
  otherHeightFromDiameterMixed$sharmaPartonRelDbhPhysio = fit_nlme("Sharma-Parton RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1e * elevation)*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, other2022, 
                                                                   fixedFormula = a1 + a1rd + a1e + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                                   start = list(a1 = 10, a1rd = 0, a1e = 0, b1 = 0.5, b2 = -0.007, b3 = 0.4, b4 = 1.1), distinct = FALSE)
  otherHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1) * (1 - exp(b2*standTreesPerHectare^(b3 + b3ra)*dbh))^b4, other2022, # a1ra step halving, b2ra, b4ra singularity in backsolve, b1ra of similar accuracy with notably conflicting signals between MAE and AIC
                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 20, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9))) # 2x50 max iterations, 2x250, 2x500 singularity in backsolve
  #otherHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*(1 + a1rm)*standBasalAreaPerHectare^(b1) * (1 - exp(b2*standTreesPerHectare^(b3)*dbh))^(b4), other2022, # a1rm, b1rm, b2rm, b3rm, b4rm singularity in backsolve
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 20, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  #otherHeightFromDiameterMixed$sharmaZhangB1 = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1ra) * (1 - exp(b2*standTreesPerHectare^(b3)*dbh))^b4, other2022,
  #                                                      fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 20, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  #otherHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*standTreesPerHectare^(b3 + a1ra)*dbh))^b4, other2022, # a1ra, b1ra, b2ra, b3ra, b4ra singularity in backsolve, presumably due to a1bal lacking significance
  #                                                       fixedFormula = a1 + a1bal + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 20, a1bal = 0, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 0.9)))
  otherHeightFromDiameterMixed$sharmaZhangBal = create_fit_statistics("Sharma-Zhang BA+L", fitting = "nlme", distinct = FALSE)
  otherHeightFromDiameterMixed$sharmaZhangBalRelDbh = fit_nlme("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + (a1 + a1ra)*standBasalAreaPerHectare^(b1) * (1 - exp((b2 + b2bal * basalAreaLarger)*standTreesPerHectare^(b3)*dbh))^(b4 + b4rd * relativeDiameter), other2022, 
                                                               fixedFormula = a1 + a1bal + b1 + b2 + b2bal + b3 + b4 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 25, b1 = 0, b2 = -0.01, b2bal = 0, b3 = 0.2, b4 = 0.9, b4rd = 0)), distinct = FALSE)
  otherHeightFromDiameter$sharmaZhangRelDbh = fit_nlme("Sharma-Zhang RelDbh", height ~ 1.37 + (a1 + a1ra)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3)*dbh))^(b4 + b4rd * relativeDiameter), other2022, 
                                                       fixedFormula = a1 + a1bal + b1 + b2 + b3 + b4 + b4rd ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 25, b1 = 0, b2 = -0.01, b3 = 0.2, b4 = 0.9, b4rd = 0)), distinct = FALSE)
  otherHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1ra) * dbh^b2), other2022, # a1ra, b2ra max iterations, 10x10 MAE 2.22, AIC 731, probably slightly more accurate than b1rm
                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 1.2, b1 = 0.8, b2 = -0.05)))
  #otherHeightFromDiameterMixed$sibbesenA1 = fit_nlme("Sibbesen", height ~ 1.37 + a1*(1 + a1rm)*dbh^((b1) * dbh^(b2)), other2022, # 10x10 MAE 2.45, AIC 812
  #                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1.3, b1 = 1.1, b2 = -0.1)))
  #otherHeightFromDiameterMixed$sibbesenB1 = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1*(1 + b1rm)) * dbh^b2), other2022, # b2rm max iterations, 10x10 MAE 2.23, AIC 732
  #                                                   fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1.0, b1 = 1.3, b2 = -0.12)))
  otherHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + a1*(1 - exp(b1 * dbh^(b2 + b2ra))), other2022, # a1ra max iterations, b1ra step halving, b2ra 2x50 job step halving, 10x10 MAE 2.40, AIC 795
                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 90, b1 = -0.015, b2 = 0.8))) # 2x50, 2x250, 2x500, 5x100 step halving
  #otherHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + a1*(1 + b2rm)*(1 - exp(b1 * dbh^(b2))), other2022, # b1rm step halving, b2rm max iterations, 10x10 MAE 2.64, AIC 775
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 40, b1 = -0.03, b2 = 0.9)))
  otherHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare) * (1 - exp(b1 * dbh^(b2 + b2ra))), other2022, # a1ra max iterations, b1ra singular, a1bal not significant
                                                     fixedFormula = a1 + a1ba + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 90, a1ba = 0, b1 = -0.015, b2 = 0.8)), distinct = FALSE)
  #to_fixed_coeffficients(otherHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(otherHeightFromDiameterMixed$weibullBal), n = 16)
  #otherHeightFromDiameterMixed$sharmaZhang$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  otherHeightFromDiameterMixed$weibullBalRelDbh = fit_nlme("Weibull RelDbh BA+L", height ~ 1.37 + (a1) * (1 - exp((b1 + b1rd * relativeDiameter + b1bal / standBasalAreaPerHectare) * (1 + b1rm) * dbh^(b2))), other2022, 
                                                           fixedFormula = a1 + b1 + b1rd + b1bal + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                           start = list(a1 = 40, b1 = -0.04, b1rd = 0, b1bal = 0, b2 = 0.9), distinct = FALSE)
  #otherHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1 + a1ra) * (1 - exp((b1 + b1rd * relativeDiameter) * dbh^(b2))), other2022, # a1ra max iterations, b1ra, b2ra step halving
  #                                                      fixedFormula = a1 + b1 + b1rd + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 30, b1 = -0.05, b1rd = 0.002, b2 = 0.9)))
  otherHeightFromDiameterMixed$weibullRelDbh = fit_nlme("Weibull RelDbh", height ~ 1.37 + (a1)*(1 - exp((b1 + b1rd * relativeDiameter)*(1 + b1rm)*dbh^(b2))), other2022, # a1rm-b1rd, b1rd-b2rm not mutually significant, b1rm step halving
                                                        fixedFormula = a1 + b1 + b1rd + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 30, b1 = -0.05, b1rd = 0.002, b2 = 0.9)), distinct = FALSE)
  #print(to_parameter_confidence_intervals(otherHeightFromDiameterMixed$weibullRelDbh), n = 40)
  #to_fixed_coeffficients(otherHeightFromDiameterMixed$weibullRelDbh)
  #otherHeightFromDiameterMixed$weibullRelDbh$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #htDiaOptions$folds = 2
  #htDiaOptions$repetitions = 2
  #htDiaOptions$folds = 10
  #htDiaOptions$repetitions = 10

  otherHeightFromDiameterMixed$gam = fit_gam("REML GAM", height ~ I(1.523 * dbh^0.697) + s(dbh, bs = "ts", k = 3), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE) # plantation not significant, collapses to linear
  otherHeightFromDiameterMixed$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", k = 11), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherHeightFromDiameterMixed$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, bs = "ts", by = isPlantation, k = 11), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherHeightFromDiameterMixed$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, terrainRoughness, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherHeightFromDiameterMixed$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I(1.523 * dbh^0.697) + s(dbh, basalAreaLarger, relativeDiameter, bs = "ts", k = 11), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE) # k = 11 minimum> ~6
  otherHeightFromDiameterMixed$gamPhysio = fit_gam("REML GAM physio", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE) # k = 4 minimum, collapses to linear
  
  if (htDiaOptions$folds == 2)
  {
    reset_future() # work around chronic future failures @ 2x50 on Future (<unnamed-112>) of class MultisessionFuture interrupted, while running on 'localhost' 
  }
  otherHeightFromDiameterMixed$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(1.523 * dbh^0.697) + s(dbh, relativeDiameter, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = other2022) # k = 4 minimum, only marginally significant
  otherHeightFromDiameterMixed$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(1.523 * dbh^0.697) + s(dbh, terrainRoughness, relativeDiameter, bs = "ts", k = 11, by = as.factor(isPlantation)), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE) # k = 11 minimum > ~5
  #lapply(otherHeightFromDiameterMixed$gamRelDbhPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(otherHeightFromDiameterMixed$gamRelDbhPhysio$fit, function(fit) { return(summary(fit$gam)) })
  
  saveRDS(otherHeightFromDiameterMixed, paste0("trees/height-diameter/data/other height mixed ", get_cross_validation_data_suffix(), ".Rds"))
}


## other species diameter regressions
if (otherOptions$fitDbh)
{
  otherDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ I(height - 1.37), other2022))
  otherDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ I(height - 1.37) + I((height - 1.37)^2), other2022)

  otherDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^b2, other2022, 
                                                       start = list(a1 = 4, b1 = 0.3, b2 = 1.0), control = gsl_nls_control(scale = "levenberg")) # a1-b1 evaporation, a1p, b1p, b2p not significant
  otherDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(exp((b1)*(height - 1.37)) - 1)^(b2), other2022, 
                                                           start = list(a1 = 4, a1aba = 0.0, b1 = 0.3, b2 = 1.0), distinct = FALSE) # a1-b1 evaporation, { a1, b1, b2 } x { aba, bat } not significant
  otherDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^(b2 + b2rh * relativeHeight), other2022, 
                                                            start = list(a1 = 4, b1 = 0.3, b2 = 1.0, b2rh = 0), distinct = FALSE) # a1rh, b1rh, b2rh not significant
  otherDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, 
                                                        start = list(a1 = -6.5, a1p = 2, b1 = 0.3, b2 = 0.4, b2p = 0.08)) # b1p not significant
  otherDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*log(1 - pmin((b1)*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, 
                                                            start = list(a1 = -7, a1aba = 0.0, b1 = 0.35, b2 = 0.37, b2p = -0.025), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  otherDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2as * sin(pi/180 * aspect) + b2e * elevation), 0.9999)), other2022, 
                                                              start = list(a1 = -5, b1 = 0.3, b2 = 0.4, b2as = 0.0, b2e = 0.0001)) # { a1, b1, b2 } x { s, ac, tr, tw }, a1as, b1as not significant, b1e often significant
  otherDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1p * isPlantation + a1rh * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, 
                                                             start = list(a1 = -6.5, a1p = 2, a1rh = 0, b1 = 0.3, b2 = 0.4, b2p = 0.08), distinct = FALSE) # a1rh, b1rh, b2rh not significant
  otherDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, 
                                                               start = list(a1 = -100, a2 = -200, b1 = 1.5)) # a1p, a2p, b1p not significant, form not well posed
  otherDiameterFromHeight$michaelisMentenReplaceAbat = fit_gsl_nls("Michaelis-Menten replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), other2022, 
                                                                   start = list(a1 = -100, a1aba = 0, a2 = -200, b1 = 1.5), distinct = FALSE) # a1aba, a1abar, a1aat, a1aat+r, a2aba, a2abar, a2aat, a2aat+r, b1aba, b1abar, b1aat, b1aat+r not significant
  otherDiameterFromHeight$michaelisMentenReplaceAbatRelHt = fit_gsl_nls("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), other2022, 
                                                                        start = list(a1 = -100, a1rh = 0, a1aba = 10, a2 = -200, b1 = 1.5), distinct = FALSE)
  otherDiameterFromHeight$michaelisMentenReplaceRelHt = fit_gsl_nls("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), other2022, 
                                                                    start = list(a1 = -100, a1rh = 0, a2 = -200, b1 = 1.5), distinct = FALSE) # a1rh, a2rh, b1rh not significant
  otherDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + a2*sqrt(height - 1.37)), other2022, 
                                                start = list(a1 = 2.3, a2 = -0.13)) # a1p, a2p not significant
  otherDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, other2022, 
                                              start = list(a1 = 0.8, b1 = 1.3)) # a1p, b1p not significant, 10x10 RMSE 9.35, AIC 928
  otherDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1), other2022, 
                                                  start = list(a1 = 0.8, a1aba = 0, b1 = 1.3), distinct = FALSE) # { a1, b1 } x { aba, bat } not significant
  otherDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ (a1 + a1s * slope)*(height - 1.37)^b1, other2022, 
                                                    start = list(a1 = 2.7, a1s = 0, b1 = 0.8), distinct = FALSE) # { a1, b1 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^b1, other2022, 
                                                   start = list(a1 = 0.8, a1rh = 0, b1 = 1.3), distinct = FALSE) # a1rh, b1rh not significant, 10x10 RMSE 9.49, AIC 931
  #otherDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1)*(height - 1.37)^b1*relativeHeight^b1rh, other2022, 
  #                                                 start = list(a1 = 0.8, b1 = 1.3, b1rh = 0), distinct = FALSE) # b1rh not significant
  #otherDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ (a1)*(height - 1.37)^b1*(1 + relativeHeight)^b1rh, other2022, 
  #                                                 start = list(a1 = 0.8, b1 = 1.3, b1rh = 0), distinct = FALSE) # b1rh not significant
  otherDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
                                              start = list(a1 = 1, b1 = 1.5, b2 = 0), distinct = FALSE) # a1p, b1p, b2, b2p not significant -> collapses to power
  otherDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), other2022, 
                                                  start = list(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # { a1, b1, b2 } x { aba, bat }, b2 not significant
  otherDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
                                                        start = list(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0, b2s = 0), distinct = FALSE) # by propagation
  otherDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
                                                             start = list(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0, b2s = 0), distinct = FALSE) # by propagation
  otherDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
                                                       start = list(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # by propagation
  otherDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
                                                    start = list(a1 = 1, b2s = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
                                                   start = list(a1 = 1, a1rh = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # a1rh, a1rhp, b1rh, b1rhp, b2rh, b2rhp not significant
  otherDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ (a1 + a1rh * relativeHeight + a1e * elevation)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
                                                         start = list(a1 = 1, a1rh = 0, a1e = 0, b1 = 1.5, b2 = 0), distinct = FALSE)
  otherDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), other2022, 
                                                start = list(a1 = 0.00003, a2 = 0.01, b1 = 1.4, Ha = 70)) # a1p, a2p, b1p not significant
  otherDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1), other2022, 
                                                     start = list(a1 = 10, b1 = 0.2, b2 = 0.1, b3 = -0.1), control = gsl_nls_control(scale = "levenberg")) # a1p, b1p, b2p, b3p not significant, a1-b2 evaporation, nan-inf with b4
  otherDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                        start = list(a1 = 1, b1 = 1.5, b2 = 0)) # a1p, b1p, b2p not significant, 10x10 RMSE 9.58, AIC 928, 2x250 max iterations
  otherDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), other2022, 
                                                            start = list(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  otherDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                                  start = list(a1 = 1, a1aba = 0, a1tr = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # by propagation
  otherDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                                       start = list(a1 = 1, a1aba = 0, a1tr = 0, a1rh = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # by propagation
  otherDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                                 start = list(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0), distinct = FALSE)
  otherDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2tw * topographicWetnessFD8f)), other2022, 
                                                              start = list(a1 = 1, b1 = 1.5, b2 = 0, b2tw = 0), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  otherDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                             start = list(a1 = 1, a1rh = 0, b1 = 1.5, b2 = 0), distinct = FALSE) # a1rh, b1rh, b2rh not significant, 10x10 RMSE 
  #otherDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(height - 1.37)^(b1*(relativeHeight*(height - 1.37))^b2), other2022, 
  #                                                           start = list(a1 = 1, b1 = 1.5, b2 = 0), distinct = FALSE) # b2 not significant
  #otherDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ (a1)*(relativeHeight*(height - 1.37))^(b1*(height - 1.37)^b2), other2022, 
  #                                                           start = list(a1 = 4.5, b1 = 0.5, b2 = 0.1)) # 10x10 RMSE 10.5, AIC 956
  otherDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ (a1 + a1rh * relativeHeight + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                                   start = list(a1 = 3.5, a1rh = 0, a1tr = 0, b1 = 0.3, b2 = 0.33), distinct = FALSE) # by propagation
  otherDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, other2022, 
                                                start = list(a1 = -10000, b1 = 0.0001, b2 = 1.3))
  #to_fixed_coeffficients(otherDiameterFromHeight$weibull)
  #print(to_parameter_confidence_intervals(otherDiameterFromHeight$weibull), n = 16)
  #ggplot() +
  #  geom_point(aes(x = dbh, y = height), other2022, alpha = 0.1, shape = 16) +
  #  geom_line(aes(x = 6*(exp(0.3*(height - 1.37)) - 1)^0.3, y = height), other2022) + # Chapman-Richards replace
  #  geom_line(aes(x = -160 * (height - 1.37)^1.2 / (-250 - (height - 1.37)^1.5), y = height), other2022) # Michaelis-Menten replace
  #  geom_line(aes(x = (-10000*log(1 - pmin(0.0001*(height - 1.37), 0.9999)))^1.2, y = height), other2022) # Weibull

  # all significant GAM forms NaN with Γ link
  otherDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, bs = "ts", k = 6), data = other2022) # plantation not significant, 10x10 RMSE 8.67, AIC 1005
  #otherDiameterFromHeight$gamMichaelis = fit_gam("REML GAM", dbh ~ I(-202.573 * (height - 1.37)^1.479/(-339.456 - (height - 1.37)^1.479)) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022) # 10x10 RMSE 9.27, AIC 963
  #otherDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, basalAreaTaller, bs = "ts", k = 12), data = other2022, distinct = FALSE) # 10x10 RMSE 8.85, AIC 1022
  #otherDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, bootstrapStandBasalAreaPerHectare, bs = "ts", k = 3), data = other2022, distinct = FALSE) # collapses to linear
  otherDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 11), data = other2022, distinct = FALSE) # min k = 11 > ~7
  otherDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, basalAreaTaller, cos(3.141593/180*aspect), bs = "ts", k = 11), data = other2022, distinct = FALSE) # by propagation
  otherDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, relativeHeight, basalAreaTaller, cos(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 16), data = other2022, distinct = FALSE) # by propagation
  #otherDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, relativeHeight, basalAreaTaller, bootstrapStandBasalAreaPerHectare, slope, bs = "ts", by = as.factor(isPlantation), k = 57), data = other2022, distinct = FALSE) # min k = 57 > 27-38
  #otherDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, elevation, bs = "ts", k = 3), data = other2022, distinct = FALSE) # collapses to linear
  #otherDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, slope, bs = "ts", k = 10), data = other2022) # 10x10 RMSE 8.61, AIC = 1036
  #otherDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, sin(3.141593/180*aspect), bs = "ts", k = 4), data = other2022, distinct = FALSE) # min k = 4 > ~1
  #otherDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, cos(3.141593/180*aspect), bs = "ts", k = 9), data = other2022) # 10x10 RMSE 8.54, AIC 1001
  #otherDiameterFromHeight$gamRoughness = fit_gam("REML GAM roughness", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, terrainRoughness, bs = "ts", k = 8), data = other2022) # 10x10 RMSE 8.58, AIC 1003
  #otherDiameterFromHeight$gamWetness = fit_gam("REML GAM westness", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, topographicWetnessFD8f, bs = "ts", k = 8), data = other2022) # 10x10 RMSE 8.60, AIC 985
  #otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, cos(3.141593/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", k = 16), data = other2022, distinct = FALSE) # min k = 16 > ~8
  #otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, slope, cos(3.141593/180*aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", k = 57), data = other2022, distinct = FALSE) # min k = 57 > ~27
  #otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, cos(3.141593/180*aspect), terrainRoughness, bs = "ts", k = 11), data = other2022, distinct = FALSE) # min k = 11 > ~7
  otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, cos(3.141593/180*aspect), bs = "ts", k = 9), data = other2022) # 10x10 RMSE 8.54, AIC 1001
  otherDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, relativeHeight, bs = "ts", k = 4), data = other2022, distinct = FALSE) # min k = 4, 10x10 RMSE 8.98, AIC 1007
  otherDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, cos(3.141593/180*aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 11), data = other2022, distinct = FALSE) # by propagation
  #otherDiameterFromHeight$gamRelHtPhysio$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  #lapply(otherDiameterFromHeight$gam$fit, k.check)
  #lapply(otherDiameterFromHeight$gam$fit, summary)

  otherDiameterFromHeight$randomForest = fit_ranger("random forest", c("dbh", "height"),
                                                    other2022, mtry = 1, minNodeSize = 75, sampleFraction = 0.688) # from setup.R tuneRanger @ 500 trees
  otherDiameterFromHeight$randomForestAbatRelHtPhysio = fit_ranger("random forest RelHt ABA+T physio", c("dbh", "height", "relativeHeight", "topHeight", "basalAreaTaller", "bootstrapStandBasalAreaPerHectare", "slope", "terrainRoughness", "elevation"), # from setup.R VSURF
                                                                   other2022, mtry = 4, minNodeSize = 3, sampleFraction = 0.733) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(otherDiameterFromHeight, paste0("trees/height-diameter/data/other DBH ", get_cross_validation_data_suffix(), ".Rds"))
}


if (otherOptions$fitDbhMixed)
{
  otherDiameterFromHeightMixed = list(linear = fit_lme("linear", fixedFormula = dbh ~ I(height - 1.37), other2022, randomFormula = ~1|stand/plot))
  otherDiameterFromHeightMixed$parabolic = fit_lme("parabolic", fixedFormula = dbh ~ I(height - 1.37) + I((height - 1.37)^2), other2022, randomFormula = ~1|stand/plot)

  #otherDiameterFromHeightMixed$chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1)*(exp((b1)*(height - 1.37)) - 1)^(b2 + b2ra), other2022, # a1ra, b1ra, b2ra singularity in backsolve
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot, # |stand, |uniquePlotID singularity in backsolve
  #                                                       start = list(fixed = c(a1 = 4, b1 = 0.3, b2 = 1.0)))
  otherDiameterFromHeightMixed$chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ (a1)*(1 + a1rm)*(exp((b1)*(height - 1.37)) - 1)^(b2), other2022, # a1rm, b1rm, b2rm singularity in backsolve
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2rm ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 4, b1 = 0.3, b2 = 1.0)), distinct = FALSE)
  otherDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(exp((b1)*(height - 1.37)) - 1)^(b2), other2022,
                                                             fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                             start = list(fixed = c(a1 = 4, a1aba = 0.0, b1 = 0.3, b2 = 1.0)), distinct = FALSE)
  otherDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1rh * relativeHeight)*(exp((b1)*(height - 1.37)) - 1)^(b2 + b2ra), other2022, # a1ra, b1ra, b2ra singularity in backsolve
                                                              fixedFormula = a1 + a1rh + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 1.0, a1rh = 0, b1 = 1.3, b2 = 0.35)), distinct = FALSE)
  otherDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2p * isPlantation + b2ra), 0.9999)), other2022, # a1ra max iterations, b1ra, b2ra step halving
                                                          fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = b2ra ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = -6.5, a1p = 2, b1 = 0.3, b2 = 0.4, b2p = 0.08)))
  otherDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*log(1 - pmin((b1)*(height - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2022, 
                                                              fixedFormula = a1 + a1aba + b1 + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = -7, a1aba = 0.0, b1 = 0.35, b2 = 0.37, b2p = -0.025)), distinct = FALSE)
  #otherDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2as * sin(pi/180 * aspect) + b2e * elevation + b2ra), 0.9999)), other2022, # a1ra, b1ra, b2ra step halving
  #                                                                     fixedFormula = a1 + b1 + b2 + b2as + b2e ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                                     start = list(fixed = c(a1 = -5, b1 = 0.3, b2 = 0.4, b2as = 0.0, b2e = 0.0001)))
  #otherDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1)*(1 + a1rm)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2as * sin(pi/180 * aspect) + b2e * elevation), 0.9999)), other2022, # a1rm, b1rm, b2rm step halving
  #                                                                     fixedFormula = a1 + b1 + b2 + b2as + b2e ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                                     start = list(fixed = c(a1 = -5, b1 = 0.3, b2 = 0.4, b2as = 0.0, b2e = 0.0001)))
  otherDiameterFromHeightMixed$chapmanRichardsPhysio = create_fit_statistics("Chapman-Richards inverse physio", fitting = "nlme")
  otherDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1ra)*log(1 - pmin(b1*(height - 1.37)^(b2 + b2rh * relativeHeight), 0.9999)), other2022, # a1ra job step halving, b1ra, b2ra step halving
                                                               fixedFormula = a1 + b1 + b2 + b2rh ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = -2.5, b1 = 0.5, b2 = 0.4, b2rh = 2)), distinct = FALSE)
  #otherDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation + a1ra) * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, # a1ra, a2ra, b1ra step halving
  #                                                               fixedFormula = a1 + a1p + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot, # |stand max iterations, |uniquePlotID singularity in backsolve
  #                                                               start = list(fixed = c(a1 = -100, a1p = 0, a2 = -200, b1 = 1.5)))
  #otherDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ (a1 + a1p * isPlantation)*(1 + a1rm)*(height - 1.37)^b1 / (a2 - (height - 1.37)^b1), other2022, # a1rm, b1rm step halving
  #                                                               fixedFormula = a1 + a1p + a2 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                               start = list(fixed = c(a1 = -100, a1p = 0, a2 = -200, b1 = 1.5)))
  otherDiameterFromHeightMixed$michaelisMentenReplace = create_fit_statistics("Michaelis-Menten replace", fitting = "nlme")
  otherDiameterFromHeightMixed$michaelisMentenReplaceAbat = fit_nlme("Michaelis-Menten replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), other2022, 
                                                                     fixedFormula = a1 + a1aba + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                     start = list(a1 = -100, a1aba = 0, a2 = -200, b1 = 1.5), distinct = FALSE)
  otherDiameterFromHeightMixed$michaelisMentenReplaceAbatRelHt = fit_nlme("Michaelis-Menten replace RelHt ABA+T", dbh ~ (a1 + a1rh * relativeHeight + a1aba * bootstrapStandBasalAreaPerHectare + a1ra) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), other2022, 
                                                                          fixedFormula = a1 + a1rh + a1aba + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                          start = list(a1 = -100, a1rh = 0, a1aba = 10, a2 = -200, b1 = 1.5), distinct = FALSE)
  otherDiameterFromHeightMixed$michaelisMentenReplaceRelHt = fit_nlme("Michaelis-Menten replace RelHt", dbh ~ (a1 + a1rh * relativeHeight + a1ra) * (height - 1.37)^(b1) / (a2 - (height - 1.37)^(b1)), other2022, 
                                                                      fixedFormula = a1 + a1rh + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                      start = list(a1 = -100, a1rh = 0, a2 = -200, b1 = 1.5), distinct = FALSE)
  otherDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + (a2 + a2ra)*sqrt(height - 1.37)), other2022, # a1ra max iterations
                                                  fixedFormula = a1 + a2 ~ 1, randomFormula = a2ra ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.8, a2 = -0.2))) # 2x250 step halving, 2x500, 5x200, 10x100 max iterations
  #otherDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ a1*sqrt(height - 1.37) / (1 + (a2)*(1 + a2rm)*sqrt(height - 1.37)), other2022, # a1rm step halving, a1 not significant with a2rm -> max iterations
  #                                                fixedFormula = a1 + a2 ~ 1, randomFormula = a2rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.8, a2 = -0.2)))
  otherDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ (a1)*(height - 1.37)^(b1 + b1ra), other2022, # a1ra step halving, 10x10 RMSE 9.18, AIC 1001
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.3, b1 = 1.1)))
  #otherDiameterFromHeightMixed$powerA1 = fit_nlme("power", dbh ~ (a1)*(1 + a1rm)*(height - 1.37)^(b1), other2022, # 10x10 RMSE 10.8, AIC 1061
  #                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 1.5, b1 = 1.2)))
  #otherDiameterFromHeightMixed$powerB1 = fit_nlme("power", dbh ~ (a1)*(height - 1.37)^(b1*(1 + b1rm)), other2022, # 10x10 RMSE 9.69, AIC 965
  #                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 1.0, b1 = 1.3)))
  #otherDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^b1, other2022, 
  #                                                  fixedFormula = a1 + a1aba + a2 + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 0.8, a1aba = 0, b1 = 1.3)), distinct = FALSE)
  otherDiameterFromHeightMixed$powerAbat = create_fit_statistics("power ABA+T", fitting = "nlme", distinct = FALSE)
  otherDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ (a1 + a1ra + a1e * elevation)*(height - 1.37)^b1, other2022,
                                                      fixedFormula = a1 + a1e + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 2.7, a1e = -0.001, b1 = 0.8)), distinct = FALSE)
  otherDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1ra + a1rh * relativeHeight)*(height - 1.37)^b1, other2022,
                                                     fixedFormula = a1 + a1rh + b1 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 3.0, a1rh = 2.3, b1 = 0.4)/0), distinct = FALSE)
  otherDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ a1*(height - 1.37)^((b1 + b1ra) * exp((b2) * (height - 1.37))), other2022, # a1ra step halving, b2ra job max iterations, 10x10 RMSE 9.45, AIC 907
                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.3, b1 = 1.1, b2 = 0.002)))
  #otherDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ a1*(height - 1.37)^((b1*(1 + b1rm)) * exp((b2) * (height - 1.37))), other2022, # a1rm step halving, b2rm max iterations, 10x10 RMSE 9.43, AIC 960
  #                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1rm ~ 1|stand/plot,
  #                                              start = list(fixed = c(a1 = 1.0, b1 = 1.4, b2 = -0.002)))
  otherDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^(b1) * exp((b2) * (height - 1.37)), other2022, 
                                                    fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
                                                          fixedFormula = a1 + a1aba + b1 + b2 + b2s ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0, b2s = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp((b2 + b2s * slope) * (height - 1.37)), other2022, 
                                                               fixedFormula = a1 + a1aba + a1e + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0, b2s = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, 
                                                         fixedFormula = a1 + a1aba + a1rh + b1 + b2 + b2s ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1e * elevation + a1ra)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022,
                                                      fixedFormula = a1 + a1e + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                      start = list(fixed = c(a1 = 2.67, a1e = 0, b1 = 0.813, b2 = 0.0067)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5), distinct = FALSE)
  otherDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1ra + a1rh * relativeHeight + a1ra)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed effects not significant
                                                     fixedFormula = a1 + a1rh + b1 + b2 + b2p ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 1, a1rh = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1ra + a1e * elevation + a1rh * relativeHeight + a1ra)*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), other2022, # fixed effects not significant 
                                                           fixedFormula = a1 + a1e + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                           start = list(fixed = c(a1 = 3.4, a1e = -0.002, a1rh = 1.2, b1 = 0.5, b2 = 0.0035)), distinct = FALSE)
  #otherDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/(a1 + a1ra) * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), other2022, # a1ra, a2ra, b2ra, Hara singularity in backsolve
  #                                                fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = a1ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.00003, a2 = 0.01, b1 = 1.4, Ha = 70)))
  #otherDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/(a1*(1 + a1rm)) * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), other2022, # a1rm, a2rm, b1rm, Harm singularity in backsolve
  #                                                fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = 0.00003, a2 = 0.01, b1 = 1.4, Ha = 70)))
  otherDiameterFromHeightMixed$schnute = create_fit_statistics("Schnute inverse", fitting = "nlme")
  otherDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^(b3 + b3ra)*(height - 1.37)) - 1), other2022, # a1ra, b1ra, b2ra step halving
                                                       fixedFormula = a1 + b1 + b2 + b3 ~ 1, randomFormula = b3ra ~ 1|stand/plot,
                                                       start = list(fixed = c(a1 = 35, b1 = -0.2, b2 = 0.1, b3 = -0.1)))
  #otherDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1*(1 + a1rm)*(height - 1.37)^(b1)*(exp(b2*(standTreesPerHectare/topHeight)^(b3)*(height - 1.37)) - 1), other2022, # a1rm, b1rm, b2rm, b3rm step halving
  #                                                     fixedFormula = a1 + b1 + b2 + b3 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                     start = list(fixed = c(a1 = 35, b1 = -0.2, b2 = 0.1, b3 = -0.1)))
  otherDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ a1*(height - 1.37)^((b1 + b1ra)*(height - 1.37)^(b2)), other2022, # a1ra step halving, b2ra singular precision, a1-b1 evaporation risk, 10x10 RMSE 11.4, AIC 930
                                                          fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                          start = list(fixed = c(a1 = 1, b1 = 1.5, b2 = 0))) # 2x500 max iterations
  #otherDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ a1*(1 + a1rm)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), other2022, # a1rm job step halving, b1rm step halving, b2 not significant, b2rm leading minor
  #                                                        fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 1.7, b1 = 1.0, b2 = 0)), distinct = FALSE) # collapses to power with a1rm
  otherDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1ra)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), other2022, 
                                                              fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                              start = list(fixed = c(a1 = 1, a1aba = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022,
                                                                    fixedFormula = a1 + a1aba + a1tr + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                    start = list(fixed = c(a1 = 1, a1aba = 0, a1tr = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1tr * terrainRoughness + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, #
                                                                         fixedFormula = a1 + a1aba + a1tr + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                         start = list(fixed = c(a1 = 1, a1aba = 0, a1tr = 0, a1rh = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, 
                                                                   fixedFormula = a1 + a1aba + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                   start = list(fixed = c(a1 = 1, a1aba = 0, a1rh = 0, b1 = 1.5, b2 = 0)), distinct = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1ra + a1tr * terrainRoughness)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022,
                                                                fixedFormula = a1 + a1tr + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 3.5, a1tr = -0.01, b1 = 0.3, b2 = 0.33)), distinct = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1ra + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022,
                                                               fixedFormula = a1 + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 3.0, a1rh = -0.7, b1 = 0.3, b2 = 0.348)), distinct = FALSE)
  otherDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1ra + a1tr * terrainRoughness + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), other2022, # fixed effects not significant 
                                                                     fixedFormula = a1 + a1tr + a1rh + b1 + b2 ~ 1, randomFormula = a1ra ~ 1|stand/plot,
                                                                     start = list(fixed = c(a1 = 3.5, a1tr = -0.01, a1rh = 0, b1 = 0.3, b2 = 0.33)), distinct = FALSE)
  #otherDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2 + b2ra), other2022, # a1ra, b1ra, b2ra singularity in backsolve
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2ra ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = -10000, b1 = 0.0001, b2 = 1.3)))
  #otherDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", dbh ~ (a1*(1 + a1rm)*log(1 - pmin(b1*(height - 1.37), 0.9999)))^(b2), other2022, # a1rm, b1rm, b2rm singularity in backsolve
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1rm ~ 1|stand/plot,
  #                                                start = list(fixed = c(a1 = -10000, b1 = 0.0001, b2 = 1.3)))
  otherDiameterFromHeightMixed$weibull = create_fit_statistics("Weibull inverse", fitting = "nlme")
  #to_fixed_coeffficients(otherDiameterFromHeightMixed$sibbesenReplace)

  otherDiameterFromHeightMixed$gam = fit_gam("REML GAM", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, bs = "ts", k = 5), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE) # collapses to linear
  otherDiameterFromHeightMixed$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherDiameterFromHeightMixed$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, basalAreaTaller, cos(3.141593/180*aspect), bs = "ts", k = 11), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherDiameterFromHeightMixed$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, relativeHeight, basalAreaTaller, cos(3.141593/180*aspect), bs = "ts", by = as.factor(isPlantation), k = 16), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherDiameterFromHeightMixed$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, cos(3.141593/180*aspect), bs = "ts", k = 9), random = list(stand = ~1, plot = ~1), data = other2022)
  otherDiameterFromHeightMixed$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, relativeHeight, bs = "ts", k = 4), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  otherDiameterFromHeightMixed$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ s(height, cos(3.141593/180*aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(stand = ~1, plot = ~1), data = other2022, distinct = FALSE)
  #lapply(otherDiameterFromHeightMixed$gamRelHt$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(otherDiameterFromHeightMixed$gamRelHt$fit, function(fit) { return(summary(fit$gam)) })
  
  saveRDS(otherDiameterFromHeightMixed, paste0("trees/height-diameter/data/other DBH mixed ", get_cross_validation_data_suffix(), ".Rds"))
}


## collect model parameters
if (otherOptions$fitHeight & otherOptions$fitHeightMixed & otherOptions$fitDbh & otherOptions$fitDbhMixed)
{
  if (exists("otherHeightFromDiameter") == FALSE)
  { 
    otherHeightFromDiameter = readRDS(paste0("trees/height-diameter/data/blocked by stand/other height ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("otherHeightFromDiameterMixed") == FALSE)
  { 
    otherHeightFromDiameterMixed = readRDS(paste0("trees/height-diameter/data/blocked by stand/other height mixed ", get_cross_validation_data_suffix(), ".Rds"))
  }
  if (exists("otherDiameterFromHeight") == FALSE)
  { 
    otherDiameterFromHeight = readRDS(paste0("trees/height-diameter/data/blocked by stand/other DBH ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  if (exists("otherDiameterFromHeightMixed") == FALSE)
  { 
    otherDiameterFromHeightMixed = readRDS(paste0("trees/height-diameter/data/blocked by stand/other DBH mixed ", get_cross_validation_data_suffix(), ".Rds")) 
  }
  
  otherCoefficients = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, unnest_coefficients)),
                                          bind_rows(lapply(otherHeightFromDiameterMixed, unnest_coefficients))) %>%
                                  mutate(responseVariable = "height"),
                                bind_rows(bind_rows(lapply(otherDiameterFromHeight, unnest_coefficients)),
                                          bind_rows(lapply(otherDiameterFromHeightMixed, unnest_coefficients))) %>%
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
  saveRDS(otherCoefficients, paste0("trees/height-diameter/data/other coefficients ", get_cross_validation_data_suffix(), ".Rds"))
  saveRDS(otherFitStats, paste0("trees/height-diameter/data/other fit stats ", get_cross_validation_data_suffix(), ".Rds"))
}


## preferred forms identified (results.R, Figure 9)
if (otherOptions$recalcPreferredModels)
{
  # 10x10 model selections
  # responseVariable species        folds  maeName                  rmseName                 aicName          nseName                  isGam isRandomForest maeFitting rmseFitting aicFitting nseFitting
  # height           other species  10x10  random forest            Sibbesen                 Sibbesen         Sibbesen                     0              1 ranger     nlme        gsl_nls    nlme      
  #                                 2x250  Sibbesen                 Sibbesen                 Korf             Sibbesen                     0              0 gsl_nls    nlme        nlme       nlme     
  # DBH              other species  10x10  REML GAM                 REML GAM                 REML GAM         REML GAM                     4              0 gam        gam         gam        gam 
  #                                 2x250  REML GAM                 REML GAM                 REML GAM         REML GAM                     4              0 gam        gam         gam        gam 
  #otherHeightFromDiameterPreferred = list(curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1*dbh / (1 + dbh)^b1, other2022, start = list(a1 = 2, b1 = 0.4), folds = 1, repetitions = 1))
  #otherHeightFromDiameterPreferred = list(korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1 * dbh^b2), other2022, start = list(a1 = 300, b1 = -5, b2 = -0.1), folds = 1, repetitions = 1))
  otherHeightFromDiameterPreferred = list(korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1 * dbh^b2), other2022, 
                                                             start = list(a1 = 300, b1 = -5, b2 = -0.1), folds = 1, repetitions = 1))
  otherHeightFromDiameterPreferred$randomForest = fit_ranger("random forest", c("height", "dbh", "isPlantation"), other2022, 
                                                             mtry = 2, minNodeSize = 76, sampleFraction = 0.869, folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$sharmaPartonB4 = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh^b4d)), other2022, 
                                                                start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4d = 1.0), folds = 1, repetitions = 1)
  #otherHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^(b4 + b4as * sin(pi/180 * aspect)), other2022, start = list(a1 = 5, b1 = 0.5, b2 = -0.01, b3 = 0.3, b4 = 1.0, b4as = -0.1), folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1 * dbh^b2), other2022, 
                                                          start = list(a1 = 1.2, b1 = 1.0, b2 = -0.08), folds = 1, repetitions = 1)
  otherHeightFromDiameterPreferred$sibbesenMixed = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^((b1 + b1ra) * dbh^b2), other2022,
                                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1ra ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 1.2, b1 = 0.8, b2 = -0.05)), folds = 1, repetitions = 1)
  saveRDS(otherHeightFromDiameterPreferred, "trees/height-diameter/data/other preferred height models.Rds")

  otherDiameterFromHeightPreferred = list(gam = fit_gam("REML GAM", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, bs = "ts", by = as.factor(isPlantation), k = 5), data = other2022, folds = 1, repetitions = 1))
  otherDiameterFromHeightPreferred$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.834 * (height - 1.37)^1.280) + s(height, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022, folds = 1, repetitions = 1)
  saveRDS(otherDiameterFromHeightPreferred, "trees/height-diameter/data/other preferred diameter models.Rds")
}


## error residual distributions
if (htDiaOptions$includeInvestigatory)
{
  otherDbhGam = gam(dbh ~ I(0.750 * height^1.263) + s(height, bs = "ts", by = as.factor(isPlantation), k = 7), data = other2022)
  ggplot() +
    geom_point(aes(x = relativeHeight, y = (dbh - residuals(otherDbhGam)) / dbh), other2022, alpha = 0.1, shape = 16) +
    geom_smooth(aes(x = relativeHeight, y = (dbh - residuals(otherDbhGam)) / dbh), other2022, formula = y ~ x, method = "loess") + # GAM linear, loess kinked linear
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "relative height", y = "scale error")
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

  otherDbhGam = gam(dbh ~ I(0.834 * (height - 1.37)^1.280) + 
                          s(height, bs = "ts", by = as.factor(isPlantation), k = 6) +
                          s(relativeHeight, bs = "ts", k = 3), # plantation not significant, linear
                          #s(bootstrapStandBasalAreaPerHectare, bs = "ts", k = 3) + # not significant
                          #s(basalAreaTaller, bs = "ts", k = 3) + # not significant
                          #s(elevation, bs = "ts", k = 3) + # not significant
                          #s(slope, bs = "ts", k = 5) + # not significant
                          #s(aspect, bs = "ts", k = 3) + # not significant
                          #s(terrainRoughness, bs = "ts", k = 3) + # not significant
                          #s(topographicWetnessFD8f, bs = "ts", k = 6), # not significant
                    data = other2022, folds = 1, repetitions = 1, method = "REML", select = TRUE, , weights = heightWeight)
  k.check(otherDbhGam)
  summary(otherDbhGam)
  plot_gam_smooths(otherDbhGam)

  otherHeightGam = gam(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11), data = other2022, family = Gamma, method = "REML", select = TRUE, weights = dbhWeight) # quasi(inverse, mu^2) (not as stable as gamma, major outlier)
  otherHeightGammRidge = gam(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11) + s(uniquePlotID, bs = "re"), data = other2022, control = gam.control(nthreads = htDiaOptions$rangerThreads), method = "REML", select = TRUE, weights = dbhWeight)
  sqrt(mean((other2022$height - otherHeightGammRidge$fitted.values)^2))
  otherHeightGammRidgePopulationPrediction = predict(otherHeightGammRidge, other2022, exclude = "s(uniquePlotID)", newdata.guaranteed = TRUE) # drop random effect to make population predictions
  #otherHeightGamm = gamm(height ~ I(1.523 * dbh^0.697) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), bs = "ts", k = 11), data = other2022, family = Gamma, method = "REML", random = list(stand = ~1, plot = ~1), weights = dbhWeight)
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
