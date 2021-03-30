
prodiff <- function(t, x, params) {
  with(as.list(c(x, params)), {
    MLS2 <- f_mls2(t)
    LIN1 <- f_lin1(t)
    dHLH1 <- tau_hlh1 * (B_hlh1(MLS2, HLH1, FOS1, k_mls2_hlh1, k_hlh1_hlh1, k_fos1_hlh1) - HLH1)
    dFOS1 <- tau_fos1 * (B_fos1(LIN1, HLH1, k_lin1_fos1, k_hlh1_fos1) - FOS1)
    dCYD1 <- tau_cyd1 * (B_cyd1(FOS1, k_fos1_cyd1) - CYD1)
    dCYE1 <- tau_cye1 * (B_cye1(E2F, k_e2f_cye1) - CYE1)
    dCKI1 <- tau_cki1 * (B_cki1(HLH1, k_hlh1_cki1) - CKI1)
    dLIN35 <- tau_lin35 * (B_lin35(HLH1, E2F, CYD1, CYE1, CKI1, k_hlh1_lin35) - LIN35)
    dE2F <- tau_E2F * (phosphorylate_E2FLIN35(E2F, CYD1, CYE1, CKI1, HLH1) - E2F)
    dMEF2 <- tau_mef2 * (B_mef2(HLH1, LIN35, MEF2, FOS1, k_hlh1_mef2, k_mef2_mef2, k_fos1_mef2) - MEF2)
    proliferation_rate <- proliferation(E2F, k_e2f_proliferation)
    differentiation_rate <- differentiation(HLH1, LIN35, MEF2, FOS1, k_hlh1_differentiation, k_mef2_differentiation, k_fos1_differentiation)
    
    return(list(c(dHLH1, dFOS1, dCYD1, dCYE1, dCKI1, dLIN35, dE2F, dMEF2), 
                LIN1 = LIN1, MLS2 = MLS2, Proliferation = proliferation_rate, Differentiation = differentiation_rate))
  }) 
}


solve_prodiff <- function(params, t_start = 0, t_end = 20) {
  x0 <- c(HLH1 = 0, FOS1 = 0, CYD1 = 0, CYE1 = 0, CKI1 = 0, LIN35 = 0, 
          E2F = 0, MEF2 = 0)
  t_evl <- seq(from = t_start, to = t_end, by = 10/60)
  return(ode(times = t_evl, y = x0, func = prodiff, parms = params))
}
