
prodiff <- function(t, x, params, scenario) {
  
  f_lin1 <- input_lin1(scenario)
  f_mls2 <- input_mls2(scenario)
  
  with(as.list(c(x, params)), {
    
    MLS2 <- f_mls2(t)
    LIN1 <- f_lin1(t)
    
    dHLH1 <- tau_hlh1 * (B_hlh1(MLS2, HLH1, FOS1, k_mls2_hlh1, k_hlh1_hlh1, k_fos1_hlh1, HLH1LIN35, HLH1CYD1) - HLH1)
    dFOS1 <- tau_fos1 * (B_fos1(LIN1, HLH1, k_lin1_fos1, k_hlh1_fos1) - FOS1)
    dCYD1 <- tau_cyd1 * (B_cyd1(FOS1, k_fos1_cyd1, HLH1CYD1) - CYD1)
    dCYE1 <- tau_cye1 * (B_cye1(E2F, k_e2f_cye1) - CYE1)
    dCKI1 <- tau_cki1 * (B_cki1(HLH1, k_hlh1_cki1) - CKI1)
    dLIN35 <- tau_lin35 * (B_lin35(HLH1, k_hlh1_lin35, HLH1LIN35, E2FLIN35, CYD1, CYE1, CKI1) - LIN35)
    dE2F <- tau_E2F * (phosphorylate_E2FLIN35(E2FLIN35, CYD1, CYE1, CKI1) - E2F)
    dMEF2 <- tau_mef2 * (B_mef2(HLH1LIN35, MEF2, FOS1, k_hlh1_mef2, k_mef2_mef2, k_fos1_mef2) - MEF2)
    dHLH1LIN35 <- tau_HLH1LIN35 * (HLH1 * LIN35 - HLH1LIN35)
    dHLH1CYD1 <- tau_HLH1CYD1 * (HLH1 * CYD1 - HLH1CYD1)
    dE2FLIN35 <- -dE2F
    dRNR1 <- tau_rnr1 * (B_rnr1(E2F, k_e2f_rnr1) - RNR1)
    dUNC15 <- tau_unc15 * (B_unc15(HLH1LIN35, MEF2, FOS1, k_hlh1_unc15, k_mef2_unc15, k_fos1_unc15) - UNC15)
    
    return(list(c(dHLH1, dFOS1, dCYD1, dCYE1, dCKI1, dLIN35, dE2F, dMEF2, dHLH1LIN35, dHLH1CYD1, dE2FLIN35, dRNR1, dUNC15), 
                LIN1 = LIN1, MLS2 = MLS2))
  }) 
}


solve_prodiff <- function(params, t_start = 0, t_end = 20, scenario) {
  x0 <- c(HLH1 = 0, FOS1 = 0, CYD1 = 0, CYE1 = 0, CKI1 = 0, LIN35 = 0, 
          E2F = 0, MEF2 = 0, HLH1LIN35 = 0, HLH1CYD1 = 0, E2FLIN35 = 1, 
          RNR1 = 0, UNC15 = 0)
  t_evl <- seq(from = t_start, to = t_end, by = 10/60)
  return(ode(times = t_evl, y = x0, func = prodiff, parms = params, scenario = scenario))
}
