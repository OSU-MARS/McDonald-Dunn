# load libraries, functions, and acma2022 from setup.R

acmaOptions = tibble(fitHeight = TRUE,
                     fitHeightMixed = TRUE,
                     fitDbh = TRUE,
                     fitDbhMixed = TRUE,
                     recalcPreferredModels = FALSE)

## bigleaf maple height regressions
if (acmaOptions$fitHeight)
{
  acmaHeightFromDiameter = list(linear = fit_lm("linear", height ~ dbh + I(isPlantation*dbh), acma2022))
  acmaHeightFromDiameter$parabolic = fit_lm("parabolic", height ~ dbh + I(dbh^2) + I(isPlantation*dbh) + I(isPlantation*dbh^2), acma2022)
  
  acmaHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, b1 = -0.06, b2 = 1.1)) # a1p, b1p, b2p not significant
  # reciprocal BA+L not significant
  acmaHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + (a1) * (1 - exp((b1)*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), acma2022, start = list(a1 = 27, b1 = -0.06, b2 = 1.1, b2bal = 0, b2balr = 10), distinct = FALSE) # a1ba+r, a1bal+r, b1ba+r, b1bal+r, b2ba+r, b2bal+r not significant
  acmaHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1s * slope) * (1 - exp(b1*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), acma2022, start = list(a1 = 27, a1s = 0, b1 = -0.06, b2 = 1.1, b2bal = 0, b2balr = 10), distinct = FALSE) # a1s not significant
  acmaHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1s * slope) * (1 - exp(b1*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), acma2022, start = list(a1 = 30, a1rd = 0, a1s = 0, b1 = -0.06, b2 = 1.1, b2bal = 0, b2balr = 10), distinct = FALSE) # a1s, a1rd, b2bal+r not significant
  acmaHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1rd * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2bal / (b2balr + basalAreaLarger)), acma2022, start = list(a1 = 27, a1rd = 0, b1 = -0.027, b2 = 1.08, b2bal = 0, b2balr = 10), distinct = FALSE) # a1rd, b2bal+r not significant
  # linear BA+L not significant
  #acmaHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", height ~ 1.37 + a1 * (1 - exp(b1*dbh))^(b2 + b2ba * standBasalAreaPerHectare), acma2022, start = list(a1 = 27, b1 = -0.06, b2 = 1.1, b2ba = 0), distinct = FALSE) # a1ba, a1bal, a1balp, b1ba, b1bal, b1balp, b2ba, b2bal, b2balp not significant
  #acmaHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1s * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 27, a1bal = 0, a1s = 0, b1 = -0.06, b2 = 1.1, b2p = -0.2), distinct = FALSE) # by propagation
  #acmaHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L physio", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter + a1s * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 30, a1bal = 0, a1rd = 0, a1s = 0, b1 = -0.06, b2 = 1.1, b2p = -0.2), distinct = FALSE) # by propagation
  #acmaHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards RelDbh BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger + a1rd * relativeDiameter) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, start = list(a1 = 27, a1bal = 0, a1rd = 0, b1 = -0.027, b2 = 1.08, b2p = -0.2), distinct = FALSE) # a1ba, a1rd not significant
  acmaHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", height ~ 1.37 + a1 * (1 - exp((b1)*dbh))^(b2 + b2p * isPlantation + a1tw * topographicWetnessFD8f), acma2022, start = list(a1 = 27, a1tw = 0, b1 = -0.06, b2 = 1.1, b2p = -0.2), distinct = FALSE) # { a1, b1, b2 } x { s, e, a, tr, tw } not significant, a1tw+e NaN-Inf
  acmaHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*(1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, a1rd = 0, b1 = -0.06, b2 = 1.1), distinct = FALSE) # a1rd, b2p, b2rd not significant, b1rd nan-inf
  acmaHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", height ~ 1.37 + (a1 + a1rd * relativeDiameter + a1as * sin(3.14159/180 * slope)) * (1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, a1rd = 0, a1as = 0, b1 = -0.06, b2 = 1.1), distinct = FALSE) # by propagation
  acmaHeightFromDiameter$curtis = fit_gsl_nls("Curtis", height ~ 1.37 + a1 * dbh / (1 + dbh)^b1, acma2022, start = list(a1 = 2.4, b1 = 0.4)) # a1p, b1p not significant
  acmaHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", height ~ 1.37 + a1 / (1 + a2 * dbh^b1), acma2022, start = list(a1 = 32, a2 = 33, b1 = -1.2)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$korf = fit_gsl_nls("Korf", height ~ 1.37 + a1*exp(b1*dbh^b2), acma2022, start = list(a1 = 70, b1 = -5.0, b2 = -0.35)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 33, a2 = 30, b1 = 1.1)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$michaelisMentenBal = fit_gsl_nls("Michaelis-Menten BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1ba / standBasalAreaPerHectare)), acma2022, start = list(a1 = 33, a2 = 30, b1 = 1.3, b1ba = -3.7)) # b1bap not significant, 10x10 MAE 2.38, AIC 1181
  #acmaHeightFromDiameter$michaelisMentenBaA1 = fit_gsl_nls("Michaelis-Menten BA", height ~ 1.37 + (a1 + a1ba * standBasalAreaPerHectare)*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 30, a1ba = 0.05, a2 = 30, b1 = 1.1)) # a2ba, b1ba not significant, 10x10 MAE 2.45, AIC 1209
  #acmaHeightFromDiameter$michaelisMentenBarA1 = fit_gsl_nls("Michaelis-Menten BA", height ~ 1.37 + (a1 + a1ba / standBasalAreaPerHectare)*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 35, a1ba = -140, a2 = 30, b1 = 1.2)) # 10x10 MAE 2.40, AIC 1181
  #acmaHeightFromDiameter$michaelisMentenBarA2 = fit_gsl_nls("Michaelis-Menten BA", height ~ 1.37 + (a1)*dbh^b1 / (a2 + a2ba / standBasalAreaPerHectare + dbh^b1), acma2022, start = list(a1 = 33, a2 = 25, a2ba = 175, b1 = 1.2)) # 10x10 MAE 2.45, AIC 1230
  #acmaHeightFromDiameter$michaelisMentenBalA1 = fit_gsl_nls("Michaelis-Menten BAL", height ~ 1.37 + (a1 + a1bal * basalAreaLarger)*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 30, a1bal = 0.08, a2 = 30, b1 = 1.1)) # a1bal+r, a2bal+r, b1bal+r not significant, 10x10 MAE 2.43, AIC 1224
  #acmaHeightFromDiameter$michaelisMentenBalB1 = fit_gsl_nls("Michaelis-Menten BAL", height ~ 1.37 + (a1)*dbh^(b1 + b1bal * basalAreaLarger) / (a2 + dbh^(b1 + b1bal * basalAreaLarger)), acma2022, start = list(a1 = 35, a2 = 32, b1 = 1.1, b1bal = 0.0002)) # 10x10 MAE 2.45, AIC 1220
  #acmaHeightFromDiameter$michaelisMentenBal = fit_gsl_nls("Michaelis-Menten BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1ba / standBasalAreaPerHectare)), acma2022, start = list(a1 = 33, a1bal = 0, a2 = 30, b1 = 1.3, b1ba = -3.7), distinct = FALSE) # a1bal not significant
  acmaHeightFromDiameter$michaelisMentenBalRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1rd * relativeDiameter + b1ba / standBasalAreaPerHectare)), acma2022, start = list(a1 = 33, a2 = 30, b1 = 1.3, b1rd = 0, b1ba = -3.7), distinct = FALSE) # b1ba-b1rd not mutually significant
  acmaHeightFromDiameter$michaelisMentenRelDbh = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + dbh^(b1 + b1rd * relativeDiameter)), acma2022, start = list(a1 = 37, a2 = 36, b1 = 1.2, b1rd = -0.05)) # b1rdp not significant, 10x10 MAE 2.43, AIC 1208
  #acmaHeightFromDiameter$michaelisMentenRelDbhA1 = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a1rd * relativeDiameter)*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 40, a1rd = -3, a2 = 35, b1 = 1.1)) # 10x10 MAE 2.44, AIC 1216
  #acmaHeightFromDiameter$michaelisMentenRelDbhA2 = fit_gsl_nls("Michaelis-Menten RelDbh", height ~ 1.37 + (a1)*dbh^b1 / (a2 + a2rd * relativeDiameter + dbh^b1), acma2022, start = list(a1 = 37, a2 = 35, a2rd = 12, b1 = 1.2)) # 10x10 MAE 2.49, AIC 1208
  acmaHeightFromDiameter$prodan = fit_gsl_nls("Prodan", height ~ 1.37 + dbh^2 / (a1 * dbh^2 + (a2 + a2p * isPlantation)*dbh + a3), acma2022, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19))
  acmaHeightFromDiameter$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, acma2022, start = list(a1 = 2.0, b1 = 0.66)) # a1p, b1p not significant, all data power = 0.666
  acmaHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", height ~ 1.37 + a1*exp(b1/(dbh + b2)), acma2022, start = list(a1 = 30, b1 = -13, b2 = 3.5)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", height ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * dbh)/d^(d/(1 - d))))^(1/(1 - d)), acma2022, start = list(Ha = 25, d = 0.7, kU = 0.035)) # Hap, kUp, dp not significant
  acmaHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh))^b4, acma2022, start = list(a1 = 12, b1 = 0.2, b2 = -0.04, b3 = 0.1, b4 = 1.1)) # a1p, b1p, b2p, b3p, b4p not significant, 10x10 RMSE 3.51, AIC 1204
  acmaHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4), acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.04, b3 = 0.1, b3bal = 0, b3balr = 10, b4 = 1.1), distinct = FALSE) # a1bal+r, b1bal+r, b3bal+r, b4bal+r not significant, b2bal+r NaN-inf
  acmaHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4 + b4tw * topographicWetnessFD8f), acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b3bal = 0, b3balr = 10, b4 = 1.1, b4tw = 0), distinct = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw } not significant
  acmaHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter + b3bal / (b3balr + basalAreaLarger) + b3tr * terrainRoughness)*dbh))^b4, acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b3rd = -0.04, b3bal = 0, b3balr = 10, b3tr = 0, b4 = 1.1), distinct = FALSE) # by propagation
  acmaHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter + b3bal / (b3balr + basalAreaLarger))*dbh))^b4, acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.05, b3 = 0.1, b3rd = -0.05, b3bal = 0, b3balr = 10, b4 = 1.1), distinct = FALSE) # by propagation
  acmaHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^(b4 + b4s * slope), acma2022, start = list(a1 = 12.4, b1 = 0.2, b2 = -0.05, b3 = 0.1, b4 = 1.1, b4s = 0), distinct = FALSE) # { a1, b1, b2, b3, b4 } x { e, s, a, tr, tw }, b4tw+e not significant
  acmaHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh))^b4, acma2022, start = list(a1 = 17, b1 = 0.14, b2 = -0.04, b3 = 0.12, b3rd = -0.04, b4 = 1.2)) # a1rd, b1rd, b4rd not significant, b2rd nan-inf, 10x10 RMSE 3.62, AIC 1207
  acmaHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", height ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^(b3 + b3rd * relativeDiameter + b3as * sin(3.14159/180 * aspect))*dbh))^b4, acma2022, start = list(a1 = 17, b1 = 0.14, b2 = -0.04, b3 = 0.12, b3as = 0, b3rd = -0.04, b4 = 1.2), distinct = FALSE) # by propagation
  #->acmaHeightFromDiameter$sharmaPartonB4 = fit_gsl_nls("Sharma-Parton b₄", height ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(standTreesPerHectare/standBasalAreaPerHectare)^b3*dbh^b4)), acma2022, start = list(a1 = 12, b1 = 0.2, b2 = -0.03, b3 = 0.1, b4 = 1)) # b3, b4 not significant, 10x10 RMSE 3.45, AIC 1222
  #->acmaHeightFromDiameter$sharmaPartonB4relDbh = fit_gsl_nls("Sharma-Parton b₄ RelDbh", height ~ 1.37 + (a1)*topHeight^(b1)*(1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3 + b3rd * relativeDiameter)*dbh^(b4))), acma2022, start = list(a1 = 15, b1 = 0.2, b2 = -0.02, b3 = 0.15, b3rd = -0.06, b4 = 1.1)) # a1rd, b1, b1rd, b2rd, b4rd not significant, 10x10 RMSE 3.49, AIC 1199
  #acmaHeightFromDiameter$sharmaPartonB4bal = fit_gsl_nls("Sharma-Parton b₄ BA+L", height ~ 1.37 + (a1)*topHeight^(b1) * (1 - exp((b2)*(standTreesPerHectare/standBasalAreaPerHectare)^(b3)*dbh^(b4 + b4bal / (b4balr + basalAreaLarger)))), acma2022, start = list(a1 = 13, b1 = 0.2, b2 = -0.04, b3 = 0.1, b4bal = 0, b4balr = 10, b4 = 1.1), distinct = FALSE) # a1bal, a1bal+r, b1bal, b1bal+r, b2bal, b2bal+r, b3bal, b3bal+r, b4bal, b4bal+r not significant
  acmaHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, acma2022, start = list(a1 = 17, b1 = 0.1, b2 = -0.05, b3 = -0.1, b4 = 1.0)) # a1p, b1p, b2p, b3p, b4p not significant, 10x10 RMSE 3.55, AIC 1183
  #acmaHeightFromDiameter$sharmaZhangB4 = fit_gsl_nls("Sharma-Zhang b₄", height ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*standTreesPerHectare^b3*dbh^b4)), acma2022, start = list(a1 = 16, b1 = 0.1, b2 = -0.03, b3 = 0, b4 = 1.1)) # b3 not significant, 10x10 RMSE 3.63, AIC 1222
  acmaHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^b3*dbh))^b4, acma2022, start = list(a1 = 10, b1 = 0.3, b2 = -0.06, b2ba = 0.0005, b3 = 0.04, b4 = 1.1)) # a1ba, b1bal, b2bal, b3ba, b3bal, b4ba, b4bal not significant, { a1bal, b1ba, b2ba } all significant
  #acmaHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1) * (1 - exp((b2)*standTreesPerHectare^(b3 + b3bal / (b3balr + basalAreaLarger))*dbh))^(b4), acma2022, start = list(a1 = 10, b1 = 0.1, b2 = -0.05, b3 = -0.1, b3bal = 0, b3balr = 10, b4 = 1.0)) # a1ba, a1ba+r, a1bal+r, b1ba, b1bal+r, b2ba, b2ba+r, b2bal+r, b3ba, b3ba+r, b3bal+r, b4ba, b4ba+r, b4bal+r NaN-inf, b3ba+r not significant
  acmaHeightFromDiameter$sharmaZhangBalRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^b3*dbh))^(b4 + b4rd * relativeDiameter), acma2022, start = list(a1 = 10, b1 = 0.3, b2 = -0.06, b2ba = 0.0005, b3 = 0.04, b4 = 1.1, b4rd = 0))
  acmaHeightFromDiameter$sharmaZhangRelDbh = fit_gsl_nls("Sharma-Zhang RelDbh", height ~ 1.37 + (a1)*standBasalAreaPerHectare^(b1)*(1 - exp((b2)*standTreesPerHectare^(b3)*dbh))^(b4 + b4rd * relativeDiameter), acma2022, start = list(a1 = 17, b1 = 0.1, b2 = -0.05, b3 = -0.1, b4 = 1.0, b4rd = 0), distinct = FALSE) # a1rd, b1rd, b4rd not significant, b2rd, b3rd NaN-Inf
  acmaHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", height ~ 1.37 + a1*dbh^(b1*dbh^b2), acma2022, start = list(a1 = 0.7, b1 = 1.5, b2 = -0.16)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$weibull = fit_gsl_nls("Weibull", height ~ 1.37 + (a1)*(1 - exp(b1*dbh^b2)), acma2022, start = list(a1 = 27, b1 = -0.04, b2 = 1.0)) # a1p, b1p, b2p not significant
  acmaHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger))), acma2022, start = list(a1 = 27, b1 = -0.04, b2 = 1.0, b2bal = 0.002)) # a1ba, b1ba, b1bal, b2ba not significant, a1bal also significant but typically less accurate
  #acmaHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1) * (1 - exp(b1*dbh^(b2 + b2bal / (b2balr + basalAreaLarger)))), acma2022, start = list(a1 = 27, b1 = -0.04, b2 = 1.0, b2bal = 0, b2balr = 10)) # a1bal+r, b1bal+r, b2bal+r not significant
  #acmaHeightFromDiameter$weibullBalA1 = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp(b1*dbh^(b2))), acma2022, start = list(a1 = 26, a1bal = 0.05, b1 = -0.04, b2 = 1.0))
  #acmaHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", height ~ 1.37 + (a1 + a1bal * basalAreaLarger) * (1 - exp(b1*dbh^(b2))), acma2022, start = list(a1 = 26, a1bal = 0.08, b1 = -0.04, b2 = 1.0)) # overlaps at 10x10 blocked cross validation but tends to be less accurate
  #print(to_parameter_confidence_intervals(acmaHeightFromDiameter$sharmaPartonBalRelDbh), n = 24)
  #to_fixed_coeffficients(acmaHeightFromDiameter$weibullBal)
  #acmaHeightFromDiameter$weibullBalA1$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  # all significant GAM forms NaN with Γ link
  # Slight increases in accuracy from power to Hossfeld IV and Michaelis-Menten bases, the two of which are mathematically identical for the regression coefficients reached in this case.
  # Fit with Hossfeld IV as it's slightly simpler.
  acmaHeightFromDiameter$gam = fit_gam("REML GAM", height ~ I(1.9315 * dbh^0.666) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022) # 10x10 MAE 2.52, AIC 1238, gamma NaN
  #acmaHeightFromDiameter$gamHossfeld = fit_gam("REML GAM", height ~ I(32.450/(1 + 30.977 * dbh^-1.169)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022) # 10x10 MAE 2.48, AIC 1197 with Hossfeld IV
  #acmaHeightFromDiameter$gamMichaelis = fit_gam("REML GAM", height ~ I(32.450 * dbh^1.169/(30.977 + dbh^1.169)) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022) # 10x10 MAE 2,46, AIC 1236 with Michaelis-Menten
  #acmaHeightFromDiameter$gamBa = fit_gam("REML GAM BA", height ~ I(1.9315 * dbh^0.666) + s(dbh, standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1264
  acmaHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", height ~ I(1.9315 * dbh^0.666) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 10), data = acma2022) # 10x10 AIC 1249
  #acmaHeightFromDiameter$gamBaBal = fit_gam("REML GAM BA+L", height ~ I(1.9315 * dbh^0.666) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 11), data = acma2022, distinct = FALSE) # k = 11 minimum, 10x10 AIC 1270
  acmaHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, basalAreaLarger, slope, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), data = acma2022, distinct = FALSE) # k = 16 minimum > ~10
  acmaHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, basalAreaLarger, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), data = acma2022, distinct = FALSE) # k = 16 minimum > ~10, also by propagation
  acmaHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I(1.9315 * dbh^0.666) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), data = acma2022, distinct = FALSE) # k = 16 minimum > ~10
  # slight MAE reductions with slope, sin(aspect), roughness, and FD8f, AIC neutral -> drop sin(aspect) on highest AIC + k = 57 > ~30 -> drop wetness on AIC + k = 16 > 10
  #acmaHeightFromDiameter$gamElevation = fit_gam("REML GAM elevation", height ~ I(1.9315 * dbh^0.666) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022) # 10x10 AIC 1242 (781-1912)
  #acmaHeightFromDiameter$gamSlope = fit_gam("REML GAM slope", height ~ I(1.9315 * dbh^0.666) + s(dbh, slope, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1268 (814-1938)
  #acmaHeightFromDiameter$gamSinAspect = fit_gam("REML GAM sin(aspect)", height ~ I(1.9315 * dbh^0.666) + s(dbh, sin(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 10), data = acma2022) # 10x10 AIC 1263 (806-1880)
  #acmaHeightFromDiameter$gamCosAspect = fit_gam("REML GAM cos(aspect)", height ~ I(1.9315 * dbh^0.666) + s(dbh, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1273 (753-2007)
  #acmaHeightFromDiameter$gamRoughness = fit_gam("REML GAM roughness", height ~ I(1.9315 * dbh^0.666) + s(dbh, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1264 (742-1984)
  #acmaHeightFromDiameter$gamWetness = fit_gam("REML GAM wetness", height ~ I(1.9315 * dbh^0.666) + s(dbh, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 9), data = acma2022) # 10x10 AIC 1264 (799-1985)
  acmaHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, slope, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 11), data = acma2022) # k = 11 minimum
  acmaHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(1.9315 * dbh^0.666) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 8), data = acma2022) # 10x10 AIC 1258 (785-2087)
  acmaHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", k = 16, by = as.factor(isPlantation)), data = acma2022, distinct = FALSE) # k = 16 minimum, 10x10 AIC 1485 (1028-2767)
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, summary)
  #acmaHeightFromDiameter$gamRelDbhPhysio$validation %>% summarize(maeMin = min(mae), maeMedian = median(mae), maeMean = mean(mae), maeMax = max(mae), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))

  acmaHeightFromDiameter$randomForest = fit_ranger("random forest", c("height", "dbh"),
                                                   acma2022, mtry = 1, minNodeSize = 43, sampleFraction = 0.255) # from setup.R tuneRanger @ 500 trees
  acmaHeightFromDiameter$randomForestBalRelDbhPhysio = fit_ranger("random forest RelDbh BA+L physio", c("height", "dbh", "relativeDiameter", "standQmd", "topHeight", "standBasalAreaPerHectare", "basalAreaLarger", "terrainRoughness", "slope", "elevation", "topographicWetnessFD8f"), # from setup.R VSURF
                                                                  acma2022, mtry = 3, minNodeSize = 2, sampleFraction = 0.751) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(acmaHeightFromDiameter, paste0("trees/height-diameter/data/ACMA3 height ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}

if (acmaOptions$fitHeightMixed)
{
  acmaHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", height ~ 1.37 + (a1 + a1r)*(1 - exp(b1*dbh))^b2, acma2022, # b1r, b2r singularity in backsolve
                                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                                start = list(fixed = c(a1 = 24, b1 = -0.06, b2 = 1.1))))
  #acmaHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", height ~ 1.37 + (a1 + a1r + a1bal * basalAreaLarger) * (1 - exp(b1*dbh))^b2, acma2022, # a1ba not significant
  #                                                          fixedFormula = a1 + a1bal + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 24, a1bal = 0, b1 = -0.06, b2 = 1.1)), distinct = FALSE)
  acmaHeightFromDiameterMixed$chapmanRichardsBal = create_fit_statistics(name = "Chapman-Richards BA+L", fitting = "nlme", distinct = FALSE)
  #acmaHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", height ~ 1.37 + (a1 + a1r + a1bal * basalAreaLarger + a1s * slope) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022,
  #                                                                fixedFormula = a1 + a1ba + a1s + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 30, a1bal = 0.06, a1s = -0.1, b1 = -0.03, b2 = 1.1, b2p = -0.19)))
  acmaHeightFromDiameterMixed$chapmanRichardsBalPhysio = create_fit_statistics(name = "Chapman-Richards BA+L physio", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$chapmanRichardsBalPhysioRelDbh = create_fit_statistics("Chapman-Richards RelDbh BA+L physio", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$chapmanRichardsBalRelDbh = create_fit_statistics("Chapman-Richards RelDbh BA+L", fitting = "nlme", distinct = FALSE)
  #acmaHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", height ~ 1.37 + (a1 + a1r + a1as * sin(3.14159/180 * slope)) * (1 - exp(b1*dbh))^(b2 + b2p * isPlantation), acma2022, 
  #                                                             fixedFormula = a1 + a1as + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 30, a1as = -7, b1 = -0.034, b2 = 1.1, b2p = -0.28)))
  acmaHeightFromDiameterMixed$chapmanRichardsPhysio = create_fit_statistics(name = "Chapman-Richards physio", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$chapmanRichardsRelDbh = create_fit_statistics("Chapman-Richards RelDbh", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$chapmanRichardsRelDbhPhysio = create_fit_statistics("Chapman-Richards RelDbh physio", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$curtis = fit_nlme("Curtis", height ~ 1.37 + (a1) * dbh / (1 + dbh)^(b1 + b1r), acma2022, # a1r and b1r both viable, b1r maybe more tractable
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                start = list(fixed = c(a1 = 1.5, b1 = 0.25)))
  acmaHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", height ~ 1.37 + (a1 + a1r) / (1 + a2 * dbh^b1), acma2022, # a2r, b1r max iterations
                                                  fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 27, b1 = 30, b2 = -1.3)))
  #acmaHeightFromDiameterMixed$korf = fit_nlme("Korf", height ~ 1.37 + (a1 + a1r)*exp(b1*dbh^b2), acma2022, # a1r max iterations, b1r, b2r step halving
  #                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                            start = list(fixed = c(a1 = 70, b1 = -5.0, b2 = -0.35)), control = nlmeControl(maxIter = 500))
  acmaHeightFromDiameterMixed$korf = create_fit_statistics("Korf", fitting = "nlme")
  acmaHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", height ~ 1.37 + (a1 + a1r)*dbh^b1 / (a2 + dbh^b1), acma2022, # a1r, a2r NaN, b1r max iterations
                                                         fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                         start = list(fixed = c(a1 = 33, a2 = 30,  b1 = 1.3)))
  acmaHeightFromDiameterMixed$michaelisMentenBal = fit_nlme("Michaelis-Menten BA+L", height ~ 1.37 + (a1 + a1r)*dbh^(b1 + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1ba / standBasalAreaPerHectare)), acma2022, # a2r, b1r max iterations
                                                            fixedFormula = a1 + a2 + b1 + b1ba ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 33, a2 = 30, b1 = 1.3, b1ba = -3.7)))
  #acmaHeightFromDiameterMixed$michaelisMentenBalRelDbh = fit_nlme("Michaelis-Menten RelDbh BA+L", height ~ 1.37 + (a1)*dbh^(b1 + b1rd * relativeDiameter + b1ba / standBasalAreaPerHectare) / (a2 + dbh^(b1 + b1rd * relativeDiameter + b1ba / standBasalAreaPerHectare)), acma2022, 
  #                                                                start = list(fixed = c(a1 = 33, a2 = 30, b1 = 1.3, b1rd = 0, b1ba = -3.7)), distinct = FALSE)
  acmaHeightFromDiameterMixed$michaelisMentenBalRelDbh = create_fit_statistics("Michaelis-Menten RelDbh BA+L", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$michaelisMentenRelDbh = fit_nlme("Michaelis-Menten RelDbh", height ~ 1.37 + (a1 + a1r)*dbh^(b1 + b1rd * relativeDiameter) / (a2 + dbh^(b1 + b1rd * relativeDiameter)), acma2022, # a2r max iterations, b1r computationally singular
                                                               fixedFormula = a1 + a2 + b1 + b1rd ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                               start = list(fixed = c(a1 = 37, a2 = 36, b1 = 1.2, b1rd = -0.05)))
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
  #acmaHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b3r max iterations
  #                                                       fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 11, b1 = 0.25, b2 = -0.04, b3 = 0.1, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaPartonBal = create_fit_statistics("Sharma-Parton BA+L", fitting = "nlme", distinct = FALSE)
  #acmaHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", height ~ 1.37 + (a1 + a1r + a1as * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^b3*dbh))^b4, acma2022, # fixed effects not significant
  #                                                             fixedFormula = a1 + a1as + b1 + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 12, a1as = 0.4, b1 = 0.2, b1p = 0.04, b2 = -0.037, b3 = -0.046, b4 = 1.06)))
  acmaHeightFromDiameterMixed$sharmaPartonBalPhysio = create_fit_statistics("Sharma-Parton BA+L physio", fitting = "nlme", distinct = FALSE)
  #acmaHeightFromDiameterMixed$sharmaPartonBalRelDbh = fit_nlme("Sharma-Parton RelDbh BA+L", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3rd * relativeDiameter)*dbh))^(b4), acma2022, # a1r, b4r singularity in backsolve, b2r step halving, b3r also tractable but slow to converge
  #                                                             fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 11, b1 = 0.25, b2 = -0.04, b3 = 0.1, b3rd = -0.05, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaPartonBalPhysioRelDbh = create_fit_statistics("Sharma-Parton RelDbh BA+L physio", fitting = "nlme", distinct = FALSE)
  #acmaHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", height ~ 1.37 + (a1 + a1r + a1as * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^b3*dbh))^b4, acma2022, # fixed effects not significant
  #                                                          fixedFormula = a1 + a1as + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 11, a1as = 0.3, b1 = 0.2, b1p = 0.06, b2 = -0.03, b3 = 0, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaPartonPhysio = create_fit_statistics("Sharma-Parton physio", fitting = "nlme", distinct = FALSE)
  acmaHeightFromDiameterMixed$sharmaPartonRelDbh = fit_nlme("Sharma-Parton RelDbh", height ~ 1.37 + a1*topHeight^(b1 + b1r) * (1 - exp(b2*(standTreesPerHectare/(standBasalAreaPerHectare))^(b3 + b3rd * relativeDiameter)*dbh))^b4, acma2022, # a1r, b4r singularity in backsolve, b2r step halving, b3r max iterations
                                                            fixedFormula = a1 + b1 + b2 + b3 + b3rd + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 14, b1 = 0.2, b2 = -0.04, b3 = 0.1, b3rd = -0.05, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaPartoRelDbhPhysio = create_fit_statistics("Sharma-Parton RelDbh physio", fitting = "nlme", distinct = FALSE)
  #acmaHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1r)*(1 - exp(b2*standTreesPerHectare^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b1r job singularity in backsolve, b3r max iterations
  #                                                   fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 17, b1 = 0.1, b2 = -0.04, b3 = 0.6, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sharmaZhang = create_fit_statistics("Sharma-Zhang", fitting = "nlme")
  acmaHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1r) * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^b3*dbh))^b4, acma2022, # a1r, b2r, b4r singularity in backsolve, b3r tractable but iterations may be an issue
                                                        fixedFormula = a1 + b1 + b2 + b2ba + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                        start = list(fixed = c(a1 = 7, b1 = 0.3, b2 = -0.05, b2ba = 0.0003, b3 = 0.1, b4 = 1.0)))
  acmaHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", height ~ 1.37 + a1*dbh^(b1*dbh^(b2 + b2r)), acma2022, # a1r singularity in backsolve, b1r also tractable
                                                  fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot,
                                                  start = list(fixed = c(a1 = 0.5, b1 =2.0, b2 = -0.18)))
  acmaHeightFromDiameterMixed$weibull = fit_nlme("Weibull", height ~ 1.37 + (a1 + a1r)*(1 - exp(b1*dbh^b2)), acma2022, # b1r step halving, b2r max iterations
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 24, b1 = -0.05, b2 = 1.1)))
  acmaHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", height ~ 1.37 + (a1 + a1r) * (1 - exp(b1*dbh^(b2 + b2bal * basalAreaLarger))), acma2022, #
                                                    fixedFormula = a1 + b1 + b2 + b2bal ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                    start = list(fixed = c(a1 = 24, b1 = -0.05, b2 = 1.0, b2bal = 0.001)))
  #to_fixed_coeffficients(acmaHeightFromDiameterMixed$weibullBal)
  #print(to_parameter_confidence_intervals(acmaHeightFromDiameterMixed$weibullBal), n = 16)

  acmaHeightFromDiameterMixed$gam = fit_gam("REML GAM", height ~ I(1.9315 * dbh^0.666) + s(dbh, bs = "ts", by = as.factor(isPlantation), k = 7), random = list(uniquePlotID = ~1), data = acma2022)
  acmaHeightFromDiameterMixed$gamBal = fit_gam("REML GAM BA+L", height ~ I(1.9315 * dbh^0.666) + s(dbh, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 10), random = list(uniquePlotID = ~1), data = acma2022)
  acmaHeightFromDiameterMixed$gamBalPhysio = fit_gam("REML GAM BA+L physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, basalAreaLarger, slope, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE)
  acmaHeightFromDiameterMixed$gamBalPhysioRelDbh = fit_gam("REML GAM RelDbh BA+L physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, basalAreaLarger, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE)
  acmaHeightFromDiameterMixed$gamBalRelDbh = fit_gam("REML GAM RelDbh BA+L", height ~ I(1.9315 * dbh^0.666) + s(dbh, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE)
  acmaHeightFromDiameterMixed$gamPhysio = fit_gam("REML GAM physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, slope, terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 11), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE) # k = 11 minimum > ~6
  acmaHeightFromDiameterMixed$gamRelDbh = fit_gam("REML GAM RelDbh", height ~ I(1.9315 * dbh^0.666) + s(dbh, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 8), random = list(uniquePlotID = ~1), data = acma2022)
  acmaHeightFromDiameterMixed$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", height ~ I(1.9315 * dbh^0.666) + s(dbh, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), relativeDiameter, bs = "ts", k = 16, by = as.factor(isPlantation)), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE)
  #lapply(acmaHeightFromDiameterMixed$gamRelDbhPhysio$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(acmaHeightFromDiameterMixed$gamRelDbhPhysio$fit, function(fit) { return(summary(fit$gam)) })
  
  saveRDS(acmaHeightFromDiameterMixed, paste0("trees/height-diameter/data/ACMA3 height mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## bigleaf maple diameter regressions
if (acmaOptions$fitDbh)
{
  acmaDiameterFromHeight = list(linear = fit_lm("linear", dbh ~ 0 + I(height - 1.37), acma2022)) # plantation not significant
  acmaDiameterFromHeight$parabolic = fit_lm("parabolic", dbh ~ 0 + I(height - 1.37) + I((height - 1.37)^2), acma2022) # plantation, plantation^2 not significant

  acmaDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", dbh ~ a1*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, start = list(a1 = 40, b1 = 0.02, b2 = 1.1)) # a1p, b1p, b2p not significant, a1-b1 evaporation
  acmaDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(exp((b1)*(height - 1.37)) - 1)^(b2), acma2022, start = list(a1 = 40, a1aba = 0, b1 = 0.02, b2 = 1.1), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  acmaDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", dbh ~ a1*(exp((b1 + b1rh * relativeHeight)*(height - 1.37)) - 1)^b2, acma2022, start = list(a1 = 50, b1 = 0.02, b1rh = 0, b2 = 1.1), distinct = FALSE) # a1rh, b1rh, b2rh not significant, a1-b1 evaporation
  acmaDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, start = list(a1 = 40, b1 = -0.02, b2 = 1.3)) # a1p, b1p, b2p not significant, a1-b1 evaporation
  acmaDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*log(1 - pmin((b1)*(height - 1.37)^(b2), 0.9999)), acma2022, start = list(a1 = 40, a1aba = 0, b1 = -0.02, b2 = 1.3), distinct = FALSE) # { a1, b1, b2 } x { aba, bat } not significant
  acmaDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2tw * topographicWetnessFD8f), 0.9999)), acma2022, start = list(a1 = 40, b1 = -0.023, b2 = 1.6, b2tw = 0), distinct = FALSE) # { a1, b1, b2 } x { e, s, a, tr, tw } not significant
  acmaDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", dbh ~ a1*log(1 - pmin(b1*(height - 1.37)^(b2 + b2rh * relativeHeight), 0.9999)), acma2022, start = list(a1 = 40, b1 = -0.02, b2 = 1.3, b2rh = 0), distinct = FALSE) # a1rh, b1rh, b2rh not significant
  acmaDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^b1 / (a2 - (height - 1.37)^b1), acma2022, start = list(a1 = -116, a2 = -128, b1 = 1.30)) # a1p, a2p, b1p not significant
  acmaDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", dbh ~ (a1 + a1p * isPlantation) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), acma2022, start = list(a1 = 3, a1p = -1, a2 = -0.12, a2p = -0.02))
  acmaDiameterFromHeight$power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, acma2022, start = list(a1 = 0.8, b1 = 1.2)) # a1p, b1p not significant
  acmaDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + (b1bat + b1batp * isPlantation) * basalAreaTaller), acma2022, start = list(a1 = 0.8, b1 = 1.2, b1bat = -0.003, b1batp = -0.002)) # a1aba not significant, 10x10 RMSE 6.45, AIC 1648
  #acmaDiameterFromHeight$powerAbatB1aba = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare), acma2022, start = list(a1 = 0.7, b1aba = -0.0002, b1 = 1.4)) # 10x10 RMSE 6.65, AIC 1677
  #acmaDiameterFromHeight$powerAbatA1bat = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1bat * basalAreaTaller)*(height - 1.37)^(b1), acma2022, start = list(a1 = 1.0, a1bat = -0.007, b1 = 1.2)) # 10x10 RMSE 6.48, AIC 1695
  #acmaDiameterFromHeight$powerAbatA1batB1aba = fit_gsl_nls("power ABA+T", dbh ~ (a1 + a1bat * basalAreaTaller)*(height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare), acma2022, start = list(a1 = 0.8, a1bat = -0.007, b1aba = -0.0002, b1 = 1.2)) # 10x10 RMSE 6.56, AIC 1629
  #acmaDiameterFromHeight$powerAbatB1abaBat = fit_gsl_nls("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + b1aba * bootstrapStandBasalAreaPerHectare + b1bat * basalAreaTaller), acma2022, start = list(a1 = 0.8, b1aba = -0.0002, b1bat = -0.003, b1 = 1.2)) # 10x10 RMSE 6.50, AIC 1640
  acmaDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", dbh ~ a1*(height - 1.37)^(b1 + b1tw * topographicWetnessFD8f), acma2022, start = list(a1 = 0.8, b1 = 1.1, b1tw = 0.02)) # { a1, b1 } x { e, s, a, tr }, a1tw not significant
  acmaDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", dbh ~ a1*(height - 1.37)^(b1 + b1rh * relativeHeight), acma2022, start = list(a1 = 0.8, b1 = 1.0, b1rh = 0), distinct = FALSE) # a1rh, b1rh not significant
  acmaDiameterFromHeight$ruark = fit_gsl_nls("Ruark", dbh ~ a1*(height - 1.37)^b1 * exp(b2 * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0), distinct = FALSE) # a1p, b1p, b2, b2p not significant -> collapses to power form
  acmaDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0, b2aba = 0), distinct = FALSE) # b2 not significant for all ABA+T placements
  acmaDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2e * elevation) * (height - 1.37)), acma2022, start = list(a1 = 0.9, a1aba = 0, b1 = 1.3, b2 = 0, b2e = 0), distinct = FALSE) # by propagation
  acmaDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight + b2e * elevation) * (height - 1.37)), acma2022, start = list(a1 = 0.9, a1aba = 0, b1 = 1.3, b2 = 0, b2rh = 0, b2e = 0), distinct = FALSE) # by propagation
  acmaDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), acma2022, start = list(a1 = 0.9, a1aba = 0, b1 = 1.3, b2 = 0, b2rh = 0), distinct = FALSE) # by propagation
  acmaDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2e * elevation) * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0, b2e = 0), distinct = FALSE) # { a1, b1 } x { e, s, a, tr, tw }, b2 not significant
  acmaDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0, b2rh = 0), distinct = FALSE) # a1rh, b1rh, b2rh not significant
  acmaDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", dbh ~ a1*(height - 1.37)^b1 * exp((b2 + b2e * elevation + b2rh * relativeHeight) * (height - 1.37)), acma2022, start = list(a1 = 0.9, b1 = 1.3, b2 = 0, b2e = 0, b2rh = 0), distinct = FALSE) # by propagation
  acmaDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", dbh ~ -1/a1 * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), acma2022, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.3, Ha = 161)) # a1p, b1p, Hap not significant, Ha-a1 evaporation increases with Levenberg
  acmaDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", dbh ~ a1*(height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1), acma2022, start = list(a1 = 20, b1 = 1, b2 = 0.1, b3 = -0.1, b4 = 0.3), control = gsl_nls_control(scale = "levenberg")) # a1p, b1p, b2p, b3p, b4p all NaN-inf, a1-b2 evaporation, job NaN-inf with More
  acmaDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b2 = 0), distinct = FALSE) # a1p, b1p, b2, b2p not significant, 10x10 RSMSE 10.7, AIC 1690
  acmaDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2bat * basalAreaTaller)), acma2022, start = list(a1 = 1.2, b1 = 0.9, b2 = 0, b2bat = -0.001)) # 10x10 RMSE 6.39, AIC 1668
  #acmaDiameterFromHeight$sibbesenReplaceAbatB1aba = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1 + b1aba * bootstrapStandBasalAreaPerHectare)*(height - 1.37)^(b2)), acma2022, start = list(a1 = 1.5, b1 = 0.7, b1aba = -0.001, b2 = 0), distinct = FALSE) # a1aba not significant, ABA+T significant in all other positions but b2 always not significant -> collapses to power unless modulated by b2aba or b2bat, 10x10 RMSE 6.69, AIC 1683
  #acmaDiameterFromHeight$sibbesenReplaceAbatB2aba = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2aba * bootstrapStandBasalAreaPerHectare)), acma2022, start = list(a1 = 1.5, b1 = 0.7, b2 = 0, b2aba = -0.0007)) # 10x10 RMSE 6.51, AIC 1645
  #acmaDiameterFromHeight$sibbesenReplaceAbatA1bat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1bat * basalAreaTaller)*(height - 1.37)^((b1)*(height - 1.37)^(b2)), acma2022, start = list(a1 = 1.4, a1bat = -0.01, b1 = 1.0, b2 = 0), distinct = FALSE) # 10x10 RMSE 6.35, AIC 1658
  #acmaDiameterFromHeight$sibbesenReplaceAbatB1bat = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1 + b1bat * basalAreaTaller)*(height - 1.37)^(b2)), acma2022, start = list(a1 = 1.1, b1 = 1.2, b1bat = -0.003, b2 = 0), distinct = FALSE) # 10x10 RMSE 6.40, AIC 1673
  #acmaDiameterFromHeight$sibbesenReplaceAbatA1batB2aba = fit_gsl_nls("Sibbesen replace ABA+T", dbh ~ (a1 + a1bat * basalAreaTaller)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2aba * bootstrapStandBasalAreaPerHectare)), acma2022, start = list(a1 = 1.4, a1bat = -0.01, b1 = 1.0, b2 = 0, b2aba = -0.0007), distinct = FALSE) # a1bat-b2aba not mutually significant, 10x10 RMSE 6.30, AIC 1698
  acmaDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", dbh ~ (a1)*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^(b2 + b2bat * basalAreaTaller)), acma2022, start = list(a1 = 1.3, b1 = 1, b1tw = 0.02, b2 = 0, b2bat = -0.001), distinct = FALSE) # a1bat-b2aba not mutually significant, b2 not significant
  acmaDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^(b2 + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare)), acma2022, start = list(a1 = 1.3, b1 = 0.7, b1tw = 0.01, b2 = 0.2, b2aba = -0.001, b2rh = -0.05)) # all parameters significant
  acmaDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace RelHt ABA+T", dbh ~ (a1)*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight + b2bat * basalAreaTaller)), acma2022, start = list(a1 = 1.3, b1 = 0.7, b2 = 0.2, b2bat = -0.001, b2rh = -0.05)) # b2aba-b2bat, b2aba-b2rh not mutually significant, 10x10 RMSE 10.9, AIC 1627
  acmaDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", dbh ~ a1*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b1tw = 0.02, b2 = 0), distinct = FALSE) # b2 not significant
  acmaDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2rh * relativeHeight)), acma2022, start = list(a1 = 1.3, b1 = 1, b2 = 0, b2rh = 0), distinct = FALSE) # a1rh, b1rh, b2rh not significant
  #acmaDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ a1*(height - 1.37)^(b1*(relativeHeight*(height - 1.37))^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b2 = 0), distinct = FALSE) # b2 not significant
  #acmaDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", dbh ~ a1*(relativeHeight*(height - 1.37))^(b1*(height - 1.37)^b2), acma2022, start = list(a1 = 5, b1 = 0.5, b2 = 0.1)) # 10x10 RMSE 12.0, AIC 1679
  acmaDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", dbh ~ a1*(height - 1.37)^((b1 + b1rh * relativeHeight + b1tw * topographicWetnessFD8f)*(height - 1.37)^b2), acma2022, start = list(a1 = 1.3, b1 = 1, b1rh = 0, b1tw = 0.02, b2 = 0), distinct = FALSE) # b1rh not significant, also by propagation
  acmaDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", dbh ~ (a1*log(1 - pmin(b1*(height - 1.37), 0.9999)))^b2, acma2022, start = list(a1 = -30, b1 = 0.07, b2 = 0.63)) # a1p, b1p, b2p not significant
  #print(to_parameter_confidence_intervals(acmaDiameterFromHeight$weibull), n = 12)
  #to_fixed_coeffficients(acmaDiameterFromHeight$weibull)
  #acmaDiameterFromHeight$sibbesenReplaceAbatRelHt$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))
  
  # all significant GAM forms NaN with Γ link
  acmaDiameterFromHeight$gam = fit_gam("REML GAM", dbh ~ I(0.750 * height^1.263) + s(height, bs = "ts", by = as.factor(isPlantation), k = 7), data = acma2022) # 10x10 median RMSE 10.6, AIC 1644
  #acmaDiameterFromHeight$gamSchnute = fit_gam("REML GAM", dbh ~ I(-1/3.265e-06 * log(1 - (1 - exp(-0.001811)) * (height^1.343 - 1.37^1.343)/(172.2^1.343 - 1.37^1.343))) + s(height, bs = "ts", k = 7), data = acma2022) # 10x10 median RMSE 10.6, AIC 1654
  #acmaDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.750 * height^1.263) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 13), data = acma2022)
  #acmaDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(0.750 * height^1.263) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, slope, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  #acmaDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.750 * height^1.263) + s(height, bootstrapStandBasalAreaPerHectare, basalAreaTaller, slope, terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 85), data = acma2022)
  # modest RMSE reductions with sin(aspect), cos(aspect), roughness, wetness -> drop elevation, slope -> drop cos(aspect) on k = 56 > ~23
  #acmaDiameterFromHeight$gamElevation = fit_gam("REML GAM elevation", dbh ~ I(0.750 * height^1.263) + s(height, elevation, bs = "ts", k = 4), data = acma2022, distinct = FALSE) # typically collapses to linear, plantation not significant
  #acmaDiameterFromHeight$gamSlope = fit_gam("REML GAM slope", dbh ~ I(0.750 * height^1.263) + s(height, slope, bs = "ts", k = 8), data = acma2022, distinct = FALSE) # plantation not significant, 10x10 AIC 1694
  #acmaDiameterFromHeight$gamSinAspect = fit_gam("REML GAM sin(aspect)", dbh ~ I(0.750 * height^1.263) + s(height, sin(3.14159/180 * aspect), bs = "ts", k = 6), data = acma2022) # plantation not significant, 10x10 AIC 1631
  #acmaDiameterFromHeight$gamCosAspect = fit_gam("REML GAM cos(aspect)", dbh ~ I(0.750 * height^1.263) + s(height, cos(3.14159/180 * aspect), bs = "ts", k = 7), data = acma2022) # plantation not significant, 10x10 AIC 1694
  #acmaDiameterFromHeight$gamRoughness = fit_gam("REML GAM physio", dbh ~ I(0.750 * height^1.263) + s(height, terrainRoughness, bs = "ts", k = 7), data = acma2022) # plantation not significant, 10x10 AIC 1626
  #acmaDiameterFromHeight$gamWetness = fit_gam("REML GAM physio", dbh ~ I(0.750 * height^1.263) + s(height, topographicWetnessFD8f, bs = "ts", k = 8), data = acma2022) # plantation not significant, 10x10 AIC 1676
  acmaDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.750 * height^1.263) + s(height, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 16), data = acma2022, threads = 8) # k = 16 minimum
  acmaDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.750 * height^1.263) + s(height, relativeHeight, bs = "ts", k = 4), data = acma2022) # plantation not significant
  acmaDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(0.750 * height^1.263) + s(height, relativeHeight, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 57), data = acma2022, distinct = FALSE) # k = 57 minimum > ~22
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, k.check)
  #lapply(acmaHeightFromDiameter$gamRelDbhPhysio$fit, summary)
  #acmaDiameterFromHeight$gamRelHt$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))

  acmaDiameterFromHeight$randomForest = fit_ranger("random forest", c("dbh", "height"),
                                                   acma2022, mtry = 1, minNodeSize = 112, sampleFraction = 0.359) # from setup.R tuneRanger @ 500 trees
  acmaDiameterFromHeight$randomForestAbatRelHtPhysio = fit_ranger("random forest RelHt ABA+T physio", c("dbh", "height", "relativeHeight", "topHeight", "basalAreaTaller", "bootstrapStandBasalAreaPerHectare", "topographicWetnessFD8f", "terrainRoughness", "aspect"), # from setup.R VSURF
                                                                  acma2022, mtry = 2, minNodeSize = 3, sampleFraction = 0.586) # from setup.R tuneRanger @ 500 trees
  
  saveRDS(acmaDiameterFromHeight, paste0("trees/height-diameter/data/ACMA3 DBH ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}

if (acmaOptions$fitDbhMixed)
{
  #acmaDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", dbh ~ a1 * (exp(b1*(height - 1.37)) - 1)^(b2 + b2r), acma2022, # a1r, b1r, b2r step halving
  #                                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b2r ~ 1|stand/plot, # |stand, |uniquePlotID step halving
  #                                                             start = list(fixed = c(a1 = 40, b1 = 0.001, b2 = 1.1))))
  acmaDiameterFromHeightMixed = list(chapmanReplace = create_fit_statistics(name = "Chapman-Richards replace", fitting = "nlme"))
  #acmaDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare)*(exp((b1)*(height - 1.37)) - 1)^(b2), acma2022, 
  #                                                          fixedFormula = a1 + a1bat + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                          start = list(fixed = c(a1 = 40, a1aba = 0, b1 = 0.02, b2 = 1.1)), distinct = FALSE)
  acmaDiameterFromHeightMixed$chapmanReplaceAbat = create_fit_statistics(name = "Chapman-Richards replace ABA+T", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", dbh ~ (a1 + a1rh * relativeHeight + a1r)*(exp(b1*(height - 1.37)) - 1)^b2, acma2022, # relative height not significant
  #                                                           fixedFormula = a1 + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 50, b1 = 0.02, b1rh = 0, b2 = 1.1))), distinct = FALSE)
  acmaDiameterFromHeightMixed$chapmanReplaceRelHt = create_fit_statistics(name = "Chapman-Richards replace RelHt", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", dbh ~ a1 * log(1 - pmin(b1*(height - 1.37)^(b2 + b1r), 0.9999)), acma2022, # a1r, b1r, b2r step halving
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot, # |stand, |uniquePlotID step halving
  #                                                       start = list(fixed = c(a1 = 40, b1 = -0.02, b2 = 1.3))))
  acmaDiameterFromHeightMixed$chapmanRichards = create_fit_statistics(name = "Chapman-Richards inverse", fitting = "nlme")
  #acmaDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, 
  #                                                           fixedFormula = a1 + a1bat + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 40, a1aba = 0, b1 = -0.02, b2 = 1.3)), distinct = FALSE)
  acmaDiameterFromHeightMixed$chapmanRichardsAbat = create_fit_statistics(name = "Chapman-Richards inverse ABA+T", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", dbh ~ (a1 + a1r + a1tr * terrainRoughness)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, # physio not significant
  #                                                             fixedFormula = a1 + a1tr + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 40, a1tr = 0, b1 = -0.023, b2 = 1.6)))
  acmaDiameterFromHeightMixed$chapmanRichardsPhysio = create_fit_statistics(name = "Chapman-Richards inverse physio", fitting = "nlme")
  #acmaDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight)*log(1 - pmin(b1*(height - 1.37)^b2, 0.9999)), acma2022, # relative height not significant
  #                                                            fixedFormula = a1 + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 80, a1rh = -30, b1 = -0.012, b2 = 1.6)))
  acmaDiameterFromHeightMixed$chapmanRichardsRelHt = create_fit_statistics(name = "Chapman-Richards inverse RelHt", fitting = "nlme")
  #acmaDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", dbh ~ a1 * (height - 1.37)^(b1 + b1r) / (a2 - (height - 1.37)^(b1 + b1r)), acma2022, # a1r, b1r, b2r step halving
  #                                                              fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                                              start = list(fixed = c(a1 = -116, a3 = -128, b1 = 1.3)))
  acmaDiameterFromHeightMixed$michaelisMentenReplace = create_fit_statistics(name = "Michaelis-Menten replace", fitting = "nlme")
  acmaDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", dbh ~ (a1 + a1p * isPlantation + a1r) * sqrt(height - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(height - 1.37)), acma2022, # a2r max iterations
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
                                                 start = list(fixed = c(a1 = 3, a1p = -1, a2 = -0.12, a2p = -0.02)))
  acmaDiameterFromHeightMixed$power = fit_nlme("power", dbh ~ a1*(height - 1.37)^(b1 + b1r), acma2022, # a1r max iterations
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                               start = list(fixed = c(a1 = 0.8, b1 = 1.2)))
  acmaDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", dbh ~ (a1)*(height - 1.37)^(b1 + (b1bat + b1batp * isPlantation) * basalAreaTaller + b1r), acma2022, # a1r max iterations
                                                   fixedFormula = a1 + b1 + b1bat + b1batp ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                   start = list(fixed = c(a1 = 0.8, b1 = 1.2, b1bat = -0.003, b1batp = -0.002)))
  acmaDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", dbh ~ a1*(height - 1.37)^(b1 + b1r + b1tw * topographicWetnessFD8f), acma2022, # a1r job max iterations
                                                     fixedFormula = a1 + b1 + b1tw ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                     start = list(fixed = c(a1 = 0.8, b1 = 1.1, b1tw = 0.02)))
  #acmaDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight)*(height - 1.37)^b1, acma2022, # relative height not significant
  #                                                  fixedFormula = a1 + a1rh + b1 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 2.5, a1rh = -0.6, b1 = 1.0)))
  acmaDiameterFromHeightMixed$powerRelHt = create_fit_statistics(name = "power RelHt", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruark = fit_nlme("Ruark", dbh ~ (a1) * (height - 1.37)^(b1 + b1r) * exp((b2) * (height - 1.37)), acma2022, # b2 not significant -> collapses to power, a1r step halving, b1r, b2r tractable
  #                                             fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
  #                                             start = list(fixed = c(a1 = 0.9, b1 = 1.3, b2 = 0)))
  acmaDiameterFromHeightMixed$ruark = create_fit_statistics(name = "Ruark", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", dbh ~ (a1)*(height - 1.37)^(b1) * exp((b2 + b2aba * bootstrapStandBasalAreaPerHectare) * (height - 1.37)), acma2022, 
  #                                                 fixedFormula = a1 + b1 + b2 + b2aba ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                 start = list(fixed = c(a1 = 0.9, b1 = 1.3, b2 = 0, b2aba = 0)), distinct = FALSE)
  acmaDiameterFromHeightMixed$ruarkAbat = create_fit_statistics(name = "Ruark ABA+T", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^(b1) * exp((b2 + b2e * elevation) * (height - 1.37)), acma2022, 
  #                                                       fixedFormula = a1 + a1aba + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 0.9, a1aba = 0, b1 = 1.3, b2 = 0, b2e = 0)), distinct = FALSE)
  acmaDiameterFromHeightMixed$ruarkAbatPhysio = create_fit_statistics(name = "Ruark ABA+T physio", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark RelHt ABA+T physio", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight + b2e * elevation) * (height - 1.37)), acma2022,
  #                                                            fixedFormula = a1 + a1aba + b1 + b2 + b2rh + b2e ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 0.9, a1aba = 0, b1 = 1.3, b2 = 0, b2rh = 0, b2e = 0)), distinct = FALSE)
  acmaDiameterFromHeightMixed$ruarkAbatPhysioRelHt = create_fit_statistics(name = "Ruark RelHt ABA+T physio", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark RelHt ABA+T", dbh ~ (a1 + a1aba * bootstrapStandBasalAreaPerHectare + a1r)*(height - 1.37)^(b1) * exp((b2 + b2rh * relativeHeight) * (height - 1.37)), acma2022, 
  #                                                      fixedFormula = a1 + a1bat + a1rh + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                      start = list(fixed = c(a1 = 0.9, a1aba = 0, b1 = 1.3, b2 = 0, b2rh = 0)), distinct = FALSE)
  acmaDiameterFromHeightMixed$ruarkAbatRelHt = create_fit_statistics(name = "Ruark RelHt ABA+T", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", dbh ~ (a1 + a1r + a1as * sin(3.14159/180 * aspect))*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, # fixed effects not significant
  #                                                   fixedFormula = a1 + a1as + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                   start = list(fixed = c(a1 = 1.2, a1as = -0.055, b1 = 1.45, b1p = -0.064, b2 = 0.03)))
  acmaDiameterFromHeightMixed$ruarkPhysio = create_fit_statistics(name = "Ruark physio", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, # fixed effects not significant
  #                                                  fixedFormula = a1 + a1rh + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                  start = list(fixed = c(a1 = 1.4, a1rh = 0, b1 = 1.45, b1p = -0.3, b2 = -0.03)))
  acmaDiameterFromHeightMixed$ruarkRelHt = create_fit_statistics(name = "Ruark RelHt", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", dbh ~ (a1 + a1r + a1as * sin(3.14159/180 * aspect) + a1rh * relativeHeight)*(height - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (height - 1.37)), acma2022, # fixed effects not significant
  #                                                        fixedFormula = a1 + a1as + a1rh + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                        start = list(fixed = c(a1 = 1.2, a1as = -0.06, a1rh = 0, b1 = 1.45, b1p = -0.064, b2 = 0.03)))
  acmaDiameterFromHeightMixed$ruarkRelHtPhysio = create_fit_statistics(name = "Ruark RelHt physio", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", dbh ~ -1/(a1 + a1r) * log(1 - (1 - exp(-a2))*(height^b1 - 1.37^b1)/(Ha^b1 - 1.37^b1)), acma2022, # a1r, b1r, Har step halving
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1|stand/plot,
  #                                               start = list(fixed = c(a1 = 0.000003, a2 = 0.002, b1 = 1.3, Ha = 161)))
  acmaDiameterFromHeightMixed$schnute = create_fit_statistics(name = "Schnute inverse", fitting = "nlme")
  #acmaDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", dbh ~ a1 * (height - 1.37)^b1*(exp(b2*(standTreesPerHectare/topHeight)^b3*(height - 1.37)) - 1)^(b4 + b4r), acma2022, # a1r, b1r, b2r, b3r, b4r singularity in backsolve
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = b4r ~ 1|stand/plot,
  #                                                    start = list(fixed = c(a1 = 0.23, b1 = 1.14, b2 = 0.00001, b3 = 0.7, b4 = -0.1)))
  acmaDiameterFromHeightMixed$sharmaParton = create_fit_statistics(name = "modified Sharma-Parton", fitting = "nlme")
  #acmaDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", dbh ~ (a1 + a1r) * (height - 1.37)^(b1*(height - 1.37)^b2), acma2022, # b2 not significant -> collapses to power
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                       start = list(fixed = c(a1 = 1.3, b1 = 1, b2 = 0)))
  acmaDiameterFromHeightMixed$sibbesenReplace = create_fit_statistics(name = "Sibbesen replace", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2bat * basalAreaTaller + b2r)), acma2022, # a1r step halving, b1r singularity in backsolve, b2r job singular
  #                                                           fixedFormula = a1 + b1 + b2 + b2bat ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                           start = list(fixed = c(a1 = 1.2, b1 = 0.9, b2 = 0, b2bat = -0.001)))
  acmaDiameterFromHeightMixed$sibbesenReplaceAbat = create_fit_statistics(name = "Sibbesen replace ABA+T", fitting = "nlme")
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", dbh ~ (a1 + a1r)*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^(b2 + b2bat * basalAreaTaller)), acma2022,
  #                                                                 fixedFormula = a1 + b1 + b1tw + b2 + b2bat ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                 start = list(fixed = c(a1 = 1.3, b1 = 1, b1tw = 0.02, b2 = 0, b2bat = -0.001)), distinct = FALSE)
  acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = create_fit_statistics(name = "Sibbesen replace ABA+T physio", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace RelHt ABA+T physio", dbh ~ (a1)*(height - 1.37)^((b1 + b1tw * topographicWetnessFD8f)*(height - 1.37)^(b2 + b2rh * relativeHeight + b2aba * bootstrapStandBasalAreaPerHectare + b2r)), acma2022, # a1r, b1r step halving, b2r singular
  #                                                                      fixedFormula = a1 + b1 + b1tw + b2 + b2aba + b2rh ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                                      start = list(fixed = c(a1 = 1.3, b1 = 0.7, b1tw = 0.01, b2 = 0.2, b2aba = -0.001, b2rh = -0.05)))
  acmaDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = create_fit_statistics(name = "Sibbesen replace RelHt ABA+T physio", fitting = "nlme")
  #acmaDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace RelHt ABA+T", dbh ~ (a1)*(height - 1.37)^((b1)*(height - 1.37)^(b2 + b2rh * relativeHeight + b2bat * basalAreaTaller + b2r)), acma2022, # a1r, b1r step halving, b2r job singular
  #                                                                fixedFormula = a1 + b1 + b2 + b2bat + b2rh ~ 1, randomFormula = b2r ~ 1|stand/plot,
  #                                                                start = list(fixed = c(a1 = 1.3, b1 = 0.7, b2 = 0.2, b2bat = -0.001, b2rh = -0.05)))
  acmaDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = create_fit_statistics(name = "Sibbesen replace RelHt ABA+T", fitting = "nlme")
  #acmaDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", dbh ~ (a1 + a1r + a1as * sin(3.14159/180 * slope))*(height - 1.37)^(b1*(height - 1.37)^(b2 + b2p * isPlantation)), acma2022, # fixed effects not significant
  #                                                             fixedFormula = a1 + a1as + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                             start = list(fixed = c(a1 = 0.5, a1as = -0.2, b1 = 2.6, b2 = -0.18, b2p = 0.02)))
  acmaDiameterFromHeightMixed$sibbesenReplacePhysio = create_fit_statistics(name = "Sibbesen replace physio", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", dbh ~ (a1 + a1r + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, # fixed effects not significant
  #                                                            fixedFormula = a1 + a1rh +b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                            start = list(fixed = c(a1 = 0.65, a1rh = -0.25, b1 = 2.0, b2 = -0.1)))
  acmaDiameterFromHeightMixed$sibbesenReplaceRelHt = create_fit_statistics(name = "Sibbesen replace RelHt", fitting = "nlme", distinct = FALSE)
  #acmaDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", dbh ~ (a1 + a1r + a1as * sin(3.14159/180 * slope) + a1rh * relativeHeight)*(height - 1.37)^(b1*(height - 1.37)^b2), acma2022, # fixed effects not significant
  #                                                                  fixedFormula = a1 + a1as + a1rh + b1 + b2 ~ 1, randomFormula = a1r ~ 1|stand/plot,
  #                                                                  start = list(fixed = c(a1 = 0.5, a1as = -0.2, a1rh = -0.1, b1 = 2.6, b2 = -0.14)))
  acmaDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = create_fit_statistics(name = "Sibbesen replace RelHt physio", fitting = "nlme", distinct = FALSE)
  acmaDiameterFromHeightMixed$weibull = tryCatch({ fit_nlme("Weibull inverse", dbh ~ (a1 * log(1 - pmin((b1 + b1r)*(height - 1.37), 0.9999)))^b2, acma2022, # a1r singularity in backsolve, b1r 2x50 job step halving, b2r step halving
                                                            fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = -30, b1 = 0.07, b2 = 0.63))) },
                                                 error = function(e) { return(create_fit_statistics("Weibull inverse", fittingMethod = "nlme")) })
  #acmaDiameterFromHeightMixed$powerPhysio$validation %>% summarize(rmseMin = min(rmse), rmseMedian = median(rmse), rmseMean = mean(rmse), rmseMax = max(rmse), aicMin = min(aic), aicMedian = median(aic), aicMean = mean(aic), aicMax = max(aic))

  acmaDiameterFromHeightMixed$gam = fit_gam("REML GAM", dbh ~ I(0.750 * height^1.263) + s(height, bs = "ts", by = as.factor(isPlantation), k = 7), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE) # collapses to linear
  #acmaDiameterFromHeightMixed$gamAbat = fit_gam("REML GAM ABA+T", dbh ~ I(0.750 * height^1.263) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 13), random = list(uniquePlotID = ~1), data = acma2022)
  #acmaDiameterFromHeightMixed$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", dbh ~ I(0.750 * height^1.263) + s(height, basalAreaTaller, bootstrapStandBasalAreaPerHectare, slope, sin(3.14159/180 * aspect), terrainRoughness, bs = "ts", by = as.factor(isPlantation), k = 85), random = list(uniquePlotID = ~1), data = acma2022)
  #acmaDiameterFromHeightMixed$gamAbatPhysioRelHt = fit_gam("REML GAM RelHt ABA+T physio", dbh ~ I(0.750 * height^1.263) + s(height, bootstrapStandBasalAreaPerHectare, basalAreaTaller, slope, terrainRoughness, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 85), random = list(uniquePlotID = ~1), data = acma2022)
  acmaDiameterFromHeightMixed$gamPhysio = fit_gam("REML GAM physio", dbh ~ I(0.750 * height^1.263) + s(height, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", k = 16), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE) # k = 16 minimum > ~9, plantation not significant
  acmaDiameterFromHeightMixed$gamRelHt = fit_gam("REML GAM RelHt", dbh ~ I(0.750 * height^1.263) + s(height, relativeHeight, bs = "ts", k = 4), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE) # collapses to linear
  acmaDiameterFromHeightMixed$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", dbh ~ I(0.750 * height^1.263) + s(height, relativeHeight, sin(3.14159/180 * aspect), terrainRoughness, topographicWetnessFD8f, bs = "ts", by = as.factor(isPlantation), k = 57), random = list(uniquePlotID = ~1), data = acma2022, distinct = FALSE)
  #lapply(acmaDiameterFromHeightMixed$gamRelHt$fit, function(fit) { return(k.check(fit$gam)) })
  #lapply(acmaDiameterFromHeightMixed$gamRelHt$fit, function(fit) { return(summary(fit$gam)) })
  
  saveRDS(acmaDiameterFromHeightMixed, paste0("trees/height-diameter/data/ACMA3 DBH mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## collect model parameters
if (acmaOptions$fitHeight & acmaOptions$fitHeightMixed & acmaOptions$fitDbh & acmaOptions$fitDbhMixed)
{
  if (exists("acmaHeightFromDiameter") == FALSE)
  { 
    acmaHeightFromDiameter = readRDS(paste0("trees/height-diameter/data/ACMA3 height ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds")) 
  }
  if (exists("acmaHeightFromDiameterMixed") == FALSE) 
  { 
    acmaHeightFromDiameterMixed = readRDS(paste0("trees/height-diameter/data/ACMA3 height mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  }
  if (exists("acmaDiameterFromHeight") == FALSE)
  { 
    acmaDiameterFromHeight = readRDS(paste0("trees/height-diameter/data/ACMA3 DBH ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  }
  if (exists("acmaDiameterFromHeightMixed") == FALSE)
  { 
    acmaDiameterFromHeightMixed = readRDS(paste0("trees/height-diameter/data/ACMA3 DBH mixed ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  }
  
  acmaFixedEffects = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, unnest_fixed_effects)),
                                         bind_rows(lapply(acmaHeightFromDiameterMixed, unnest_fixed_effects))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(acmaDiameterFromHeight, unnest_fixed_effects)),
                                         bind_rows(lapply(acmaDiameterFromHeightMixed, unnest_fixed_effects))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "ACMA3")
  acmaFitStats = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, unnest_validation_statistics)),
                                    bind_rows(lapply(acmaHeightFromDiameterMixed, unnest_validation_statistics))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(acmaDiameterFromHeight, unnest_validation_statistics)),
                                    bind_rows(lapply(acmaDiameterFromHeightMixed, unnest_validation_statistics))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "ACMA3")
  
  check_plot_fit_stats(acmaFitStats)
  saveRDS(acmaFixedEffects, paste0("trees/height-diameter/data/ACMA3 coefficients ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
  saveRDS(acmaFitStats, paste0("trees/height-diameter/data/ACMA3 fit stats ", htDiaOptions$folds, "x", htDiaOptions$repetitions, ".Rds"))
}


## preferred forms identified (results.R, Figure 7)
if (acmaOptions$recalcPreferredModels)
{
  acmaHeightFromDiameterPreferred = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", height ~ 1.37 + a1*(1 - exp(b1*dbh))^b2, acma2022, start = list(a1 = 27, b1 = -0.06, b2 = 1.1), folds = 1, repetitions = 1))
  acmaHeightFromDiameterPreferred$michaelisMenten = fit_gsl_nls("Michaelis-Menten", height ~ 1.37 + a1*dbh^b1 / (a2 + dbh^b1), acma2022, start = list(a1 = 33, a2 = 30, b1 = 1.1), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$power = fit_gsl_nls("power", height ~ 1.37 + a1*dbh^b1, acma2022, start = list(a1 = 2.0, b1 = 0.66), folds = 1, repetitions = 1)
  acmaHeightFromDiameterPreferred$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", height ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1r) * (1 - exp((b2 + b2ba * standBasalAreaPerHectare)*standTreesPerHectare^b3*dbh))^b4, acma2022,
                                                            fixedFormula = a1 + b1 + b2 + b2ba + b3 + b4 ~ 1, randomFormula = b1r ~ 1|stand/plot,
                                                            start = list(fixed = c(a1 = 7, b1 = 0.3, b2 = -0.05, b2ba = 0.0003, b3 = 0.1, b4 = 1.0)), folds = 1, repetitions = 1)
  saveRDS(acmaHeightFromDiameterPreferred, "trees/height-diameter/data/ACMA3 preferred height models.Rds")

  acmaDiameterFromHeightPreferred = list(power = fit_gsl_nls("power", dbh ~ a1*(height - 1.37)^b1, acma2022, start = list(a1 = 0.8, b1 = 1.2), folds = 1, repetitions = 1))
  acmaDiameterFromHeightPreferred$powerPhysio = fit_gsl_nls("power physio", dbh ~ a1*(height - 1.37)^(b1 + b1tw * topographicWetnessFD8f), acma2022, start = list(a1 = 0.8, b1 = 1.1, b1tw = 0.02), folds = 1, repetitions = 1)
  acmaDiameterFromHeightPreferred$randomForest = fit_ranger("random forest", c("dbh", "height"), acma2022, mtry = 1, minNodeSize = 112, sampleFraction = 0.359, folds = 1, repetitions = 1)
  saveRDS(acmaDiameterFromHeightPreferred, "trees/height-diameter/data/ACMA3 preferred diameter models.Rds")
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  acmaHeightGam = gam(height ~ I(1.9315 * dbh^0.666) +
                               s(dbh, bs = "ts", by = as.factor(isPlantation), k = 9) +
                               s(relativeDiameter, bs = "ts", k = 4) + # plantation not significant, linear
                               s(standBasalAreaPerHectare, bs = "ts", k = 4) + # plantation not significant, broadly linear
                               s(basalAreaLarger, bs = "ts", k = 6) + # plantation not significant, effect 0-20 m²/ha
                               #s(elevation, bs = "ts", k = 5) + # not significant
                               #s(slope, bs = "ts", k = 4) + # not significant
                               #s(aspect, bs = "ts", k = 3) + # not significant
                               #s(terrainRoughness, bs = "ts", k = 4) + # not significant
                               s(topographicWetnessFD8f, bs = "ts", k = 4), # hyperbolic
                      data = acma2022, method = "REML", select = TRUE, weights = dbhWeight)
  k.check(acmaHeightGam)
  summary(acmaHeightGam)
  plot_gam_smooths(acmaHeightGam)
  
  acmaDbhGam = gam(dbh ~ I(0.750 * height^1.263) +
                         s(height, bs = "ts", by = as.factor(isPlantation), k = 5) + # no effect for isPlantation = TRUE
                         s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 5) + # broadly linear, opposite signs for plantation and natural regen
                         #s(bootstrapStandBasalAreaPerHectare, bs = "ts", k = 4) + # not significant
                         s(basalAreaTaller, bs = "ts", k = 3), # broadly linear
                         #s(elevation, bs = "ts", k = 3) + # not significant
                         #s(slope, bs = "ts", k = 4) + # not significant
                         #s(aspect, bs = "ts", k = 3) + # not significant
                         #s(terrainRoughness, bs = "ts", k = 3) + # not significant
                         #s(topographicWetnessFD8f, bs = "ts", k = 4), # not significant
                    data = acma2022, method = "REML", select = TRUE, weights = heightWeight)
  k.check(acmaDbhGam)
  summary(acmaDbhGam)
  plot_gam_smooths(acmaDbhGam)
}
