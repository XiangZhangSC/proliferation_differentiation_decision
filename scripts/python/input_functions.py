from toolbox import HillCube
import numpy as np

def rate_hlh1_prod(mls2, myod, 
                   k_mls2_hlh1, k_myod_hlh1, 
                   n_myod_hlh1=4.0, n_mls2_hlh1=4.0):
  
  p_mls2_bound_hlh1 = HillCube(mls2, k_mls2_hlh1, n_mls2_hlh1, normalized=True)
  p_myod_bound_hlh1 = HillCube(myod, k_myod_hlh1, n_myod_hlh1, normalized=True)
  ## binding MyoD or MLS-2 will induce hlh-1 expression
  return 1.0 - (1.0 - p_mls2_bound_hlh1) * (1.0 - p_myod_bound_hlh1)

def rate_jun1_prod(myod, jun1, fos1, k_myod_jun1, n_myod_jun1=4.0):
  return HillCube(myod, k_myod_jun1, n_myod_jun1, normalized=True)

def rate_fos1_prod(myod, fos1, jun1, k_myod_fos1, n_myod_fos1=4.0):
  
  myod_bound_fos1 = HillCube(myod, k_myod_fos1, n_myod_fos1, normalized=True)
  
  return 1.0 - myod_bound_fos1

def rate_cyd1_prod(ap1, k_ap1_cyd1, n_ap1_cyd1=4.0):
  """
  AP1 = jun-1:fos-1
  """
  
  return HillCube(ap1, k_ap1_cyd1, n_ap1_cyd1, normalized=True)

def rate_cki1_prod(myod, k_myod_cki1, n_myod_cki1=4.0):
  
  p_myod_bound_cki1 = HillCube(myod, k_myod_cki1, n_myod_cki1, normalized=True)
  
  return p_myod_bound_cki1

def rate_lin35_phos(e2f, cyd1, myod, km_e2f, e2f_tot=1.0):
  """
  Cyclin D:cdk has active kinase activity and phosphrylate E2F:pRb complex.
  
  The E2F within E2F:pRb complex is inactive, whereas free E2F is active. 
  Total E2F = [E2F] + [E2F:pRb]
  
  MyoD-cdk4 interaction directly inhibits the phosphorylation of pRb
  
  Cyclin D1 is rate limiting in the formation of active cdk4
  """
  km_e2f_star = km_e2f * (1 + myod)
  return cyd1 * (e2f_tot - e2f) / ((e2f_tot - e2f) + km_e2f_star)

def rate_e2f_prod(e2f, cyd1, cki1, lin35, myod, km_e2f, k_cki1):
  """
  Active E2F is released from E2F:pRb by adding phosphorylation group
  by cyclin D:cdk4 complex
  """
  cki1_absent = 1.0 - HillCube(cki1, k_cki1, 4.0, True)
  phos_lin35 = rate_lin35_phos(e2f, cyd1, myod, km_e2f)
  return cki1_absent * phos_lin35

def rate_lin35_prod(myod, e2f, cyd1, cki1, lin35,  
                    k_myod_lin35, km_e2f, k_cki1, n_myod_lin35=4.0):
  
  myod_induction = HillCube(myod, k_myod_lin35, n_myod_lin35, normalized=True)
  release_from_e2flin35_complex = rate_e2f_prod(e2f, cyd1, cki1, lin35, myod, km_e2f, k_cki1)
  
  return 1.0 - (1.0 - myod_induction)*(1.0 - release_from_e2flin35_complex)

def rate_proliferation(e2f, kd_e2f, n_e2f=4.0):
  
  return HillCube(e2f, kd_e2f, n_e2f, normalized=True)

def rate_differentiation(myod, cyd1, kd_myod, n_myod=4.0):
  """
  Excess cyclin D1 activates more cdk4 translocated to the nucleus where 
  cdk4 interacts with MyoD and inhibits the activation of the myogenic program
  """
  
  myod_cdk4 = np.minimum(myod, cyd1)
  myod_free = myod - myod_cdk4
  
  return HillCube(myod_free, kd_myod, n_myod, normalized=True) 
