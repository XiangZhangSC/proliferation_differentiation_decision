
HillCube <- function(x, kd, n = 4, normalized = TRUE) {
  if (normalized == TRUE) {
    return(x^n / (x^n + kd^n) / (1^n / (1^n + kd^n)))
  } else {
    return(x^n / (x^n + kd^n))
  }
  
}

B_hlh1 <- function(MLS2, HLH1, FOS1, k_mls2_hlh1, k_hlh1_hlh1, k_fos1_hlh1, 
                   HLH1LIN35, kd_HLH1LIN35, CYE1, CKI1, k_cye1_hlh1, k_cki1_hlh1) {
  induce_by_mls2 <- HillCube(MLS2, k_mls2_hlh1)
  induce_by_hlh1 <- HillCube(HLH1, k_hlh1_hlh1)
  repress_by_fos1 <- HillCube(FOS1, k_fos1_hlh1)
  
  # (MLS2 OR HLH1) AND (NOT FOS1)
  transcriptional_induction <- (1 - (1 - induce_by_mls2)*(1 - induce_by_hlh1))*(1 - repress_by_fos1)
  
  # HLH1LIN35
  deassociation_from_complexes <- HillCube(HLH1LIN35, kd_HLH1LIN35)
  
  # CKI1 AND (NOT CYE1)
  degrade_by_cye1 <- HillCube(CYE1, k_cye1_hlh1)
  stabilize_by_cki1 <- HillCube(CKI1, k_cki1_hlh1)
  degradation_by_phophosrylation <- degrade_by_cye1 * (1 - stabilize_by_cki1)
  
  # transcriptional induction OR release from complex AND (NOT phosphorylation)
  production <- 1 - (1 - transcriptional_induction)*(1 - deassociation_from_complexes)
  return(production * (1 - degradation_by_phophosrylation))
}

B_fos1 <- function(LIN1, HLH1, k_lin1_fos1, k_hlh1_fos1) {
  induce_by_lin1 <- HillCube(LIN1, k_lin1_fos1)
  repress_by_hlh1 <- HillCube(HLH1, k_hlh1_fos1)
  # LIN1 AND NOT HLH1
  return(induce_by_lin1 * (1 - repress_by_hlh1))
}

B_cyd1 <- function(FOS1, k_fos1_cyd1) {
  return(HillCube(FOS1, k_fos1_cyd1))
}

B_cye1 <- function(E2F, k_e2f_cye1) {
  return(HillCube(E2F, k_e2f_cye1))
}

B_cki1 <- function(HLH1, k_hlh1_cki1) {
  return(HillCube(HLH1, k_hlh1_cki1))
}

phosphorylate_E2FLIN35 <- function(E2FLIN35, CYD1, HLH1, CYE1, CKI1) {
  #E2F:LIN35 AND (CYD1 AND (NOT HLH1)) OR CYE1) AND (NOT CKI1)
  
  phosphorylation_by_cyd1 <- CYD1 * (1 - HLH1)
  phosphorylation_by_cye1 <- CYE1
  phosphorylation_by_cyclins <- 1 - (1 - phosphorylation_by_cyd1) * (1 - phosphorylation_by_cye1)
  
  return(E2FLIN35 * phosphorylation_by_cyclins * (1 - CKI1))
}

B_lin35 <- function(HLH1, k_hlh1_lin35, HLH1LIN35, kd_HLH1LIN35, E2FLIN35, CYD1, CYE1, CKI1) {
  # HLH1 OR HLH1LIN35 OR E2FLIN35
  induce_by_hlh1 <- HillCube(HLH1, k_hlh1_lin35)
  release_from_HLH1LIN35 <- HillCube(HLH1LIN35, kd_HLH1LIN35)
  release_from_E2FLIN35 <- phosphorylate_E2FLIN35(E2FLIN35, CYD1, HLH1, CYE1, CKI1)
  
  return(1 - (1 - induce_by_hlh1)*(1 - release_from_E2FLIN35)*(1 - release_from_HLH1LIN35))
}

B_rnr1 <- function(E2F, k_e2f_rnr1) {
  return(HillCube(E2F, k_e2f_rnr1))
}

B_unc120 <- function(HLH1LIN35, FOS1, CYD1, k_hlh1_unc120, k_fos1_unc120) {
  # ((HLH1 AND LIN35)) AND (NOT FOS1) AND (NOT CYD1)
  induce_by_hlh1 <- HillCube(HLH1LIN35, k_hlh1_unc120)
  repress_by_fos1 <- HillCube(FOS1, k_fos1_unc120)
  
  return(induce_by_hlh1 * (1 - repress_by_fos1) * (1 - CYD1))
}

B_unc15 <- function(HLH1LIN35, UNC120, FOS1, CYD1, k_hlh1_unc15, k_unc120_unc15, k_fos1_unc15) {
  # ((HLH1 AND LIN35) AND UNC120) AND (NOT FOS1)
  induce_by_hlh1 <- HillCube(HLH1LIN35, k_hlh1_unc15)
  induce_by_unc120 <- HillCube(UNC120, k_unc120_unc15)
  repress_by_fos1 <- HillCube(FOS1, k_fos1_unc15)
  
  return(induce_by_hlh1 * induce_by_unc120 * (1 - repress_by_fos1) * (1 - CYD1))
}
