
HillCube <- function(x, kd, n = 4) {
  return(x^n / (x^n + kd^n) / (1^n / (1^n + kd^n)))
}

B_hlh1 <- function(MLS2, HLH1, FOS1, k_mls2_hlh1, k_hlh1_hlh1, k_fos1_hlh1) {
  induce_by_mls2 <- HillCube(MLS2, k_mls2_hlh1)
  induce_by_hlh1 <- HillCube(HLH1, k_hlh1_hlh1)
  repress_by_fos1 <- HillCube(FOS1, k_fos1_hlh1)
  
  # (MLS2 OR HLH1) AND (NOT FOS1)
  return((1 - (1 - induce_by_mls2)*(1 - induce_by_hlh1))*(1 - repress_by_fos1))
}

B_fos1 <- function(LIN1, HLH1, k_lin1_fos1, k_hlh1_fos1) {
  induce_by_lin1 <- HillCube(LIN1, k_lin1_fos1)
  repress_by_hlh1 <- HillCube(HLH1, k_hlh1_fos1)
  # LIN1 AND NOT HLH1
  return(induce_by_lin1 * (1 - repress_by_hlh1))
}

B_cyd1 <- function(fos1, k_fos1_cyd1) {
  return(HillCube(fos1, k_fos1_cyd1))
}

B_cye1 <- function(E2F, k_e2f_cye1) {
  return(HillCube(E2F, k_e2f_cye1))
}

B_cki1 <- function(HLH1, k_hlh1_cki1) {
  return(HillCube(HLH1, k_hlh1_cki1))
}

phosphorylate_E2FLIN35 <- function(E2F, CYD1, CYE1, CKI1, HLH1) {
  E2F_tot = 1
  E2FLIN35 = E2F_tot - E2F
  #E2F:LIN35 AND ((CYD1 AND NOT HLH1) OR CYE1) AND NOT CKI1
  phosphorylation_by_cyd1 <- CYD1 * (1 - HLH1)
  phosphorylation_by_cyclins <- 1 - (1 - phosphorylation_by_cyd1)*(1 - CYE1)
  return(E2FLIN35 * phosphorylation_by_cyclins * (1 - CKI1))
}

B_lin35 <- function(HLH1, E2F, CYD1, CYE1, CKI1, k_hlh1_lin35) {
  # HLH1 OR E2FLIN35
  induce_by_hlh1 <- HillCube(HLH1, k_hlh1_lin35)
  release_from_E2FLIN35 <- phosphorylate_E2FLIN35(E2F, CYD1, CYE1, CKI1, HLH1)
  return(1 - (1 - induce_by_hlh1)*(1 - release_from_E2FLIN35))
}

proliferation <- function(E2F, k_e2f_proliferation) {
  return(HillCube(E2F, k_e2f_proliferation))
}

B_mef2 <- function(HLH1, LIN35, MEF2, FOS1, k_hlh1_mef2, k_mef2_mef2, k_fos1_mef2) {
  # ((HLH1 AND LIN35) OR MEF2) AND NOT FOS1
  HLH1LIN35 <- HLH1 * LIN35
  induce_by_hlh1 <- HillCube(HLH1LIN35, k_hlh1_mef2)
  induce_by_mef2 <- HillCube(MEF2, k_mef2_mef2)
  repress_by_fos1 <- HillCube(FOS1, k_fos1_mef2)
  
  return((1 - (1 - induce_by_hlh1) * (1 - induce_by_mef2)) * (1 - repress_by_fos1))
}

differentiation <- function(HLH1, LIN35, MEF2, FOS1, 
                            k_hlh1_differentiation, k_mef2_differentiation, k_fos1_differentiation) {
  # ((HLH1 AND LIN35) AND MEF2) AND NOT FOS1
  HLH1LIN35 <- HLH1 * LIN35
  induce_by_hlh1 <- HillCube(HLH1LIN35, k_hlh1_differentiation)
  induce_by_mef2 <- HillCube(MEF2, k_mef2_differentiation)
  repress_by_fos1 <- HillCube(FOS1, k_fos1_differentiation)
  
  return(induce_by_hlh1 * induce_by_mef2 * (1 - repress_by_fos1))
}
