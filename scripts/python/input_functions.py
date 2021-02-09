from toolbox import HillCube

def rate_hlh1_prod(mls2, myod, 
                   k_mls2_hlh1, k_myod_hlh1, 
                   n_myod_hlh1=4.0, n_mls2_hlh1=4.0):
  
  p_mls2_bound_hlh1 = HillCube(mls2, k_mls2_hlh1, n_mls2_hlh1, normalized=True)
  p_myod_bound_hlh1 = HillCube(myod, k_myod_hlh1, n_myod_hlh1, normalized=True)
  ## binding MyoD or MLS-2 will induce hlh-1 expression
  return 1.0 - (1.0 - p_mls2_bound_hlh1) * (1.0 - p_myod_bound_hlh1)

def rate_fos1_prod(myod, k_myod_fos1, n_myod_fos1=4.0):
  
  p_myod_bound_fos1 = HillCube(myod, k_myod_fos1, n_myod_fos1, normalized=True)
  
  return 1.0 - p_myod_bound_fos1

def rate_cyd1_prod(fos1, k_fos1_cyd1, n_fos1_cyd1=4.0):
  
  p_fos1_bound_cyd1 = HillCube(fos1, k_fos1_cyd1, n_fos1_cyd1, normalized=True)
  
  return p_fos1_bound_cyd1

def rate_cki1_prod(myod, k_myod_cki1, n_myod_cki1=4.0):
  
  p_myod_bound_cki1 = HillCube(myod, k_myod_cki1, n_myod_cki1, normalized=True)
  
  return p_myod_bound_cki1

def rate_lin35_phos(e2f, cyd1, km_e2f, e2f_tot=1.0):
  """
  Cyclin D:cdk has active kinase activity and phosphrylate E2F:pRb complex.
  
  The E2F within E2F:pRb complex is inactive, whereas free E2F is active. 
  Total E2F = [E2F] + [E2F:pRb]
  """
  return cyd1 * (e2f_tot - e2f) / ((e2f_tot - e2f) + km_e2f)

def rate_e2f_prod(e2f, cyd1, cki1, lin35, km_e2f, k_on, k_cki1):
  """
  Active E2F is released from E2F:pRb by adding phosphorylation group
  by cyclin D:cdk4 complex
  """
  cki1_absent = 1.0 - HillCube(cki1, k_cki1, 4.0, True)
  phos_lin35 = rate_lin35_phos(e2f, cyd1, km_e2f)
  lin35_not_binds_e2f = 1.0 - k_on * e2f * lin35
  return cki1_absent * phos_lin35 * lin35_not_binds_e2f

def rate_lin35_prod(myod, e2f, cyd1, lin35, 
                    k_myod_lin35, km_e2f, k_on, n_myod_lin35=4.0):
  
  myod_induction = HillCube(myod, k_myod_lin35, n_myod_lin35, normalized=True)
  phos_lin35 = rate_lin35_phos(e2f, cyd1, km_e2f)
  lin35_not_binds_e2f = 1.0 - k_on * e2f * lin35
  
  return (1 - (1 - myod_induction)*(1 - phos_lin35))*lin35_not_binds_e2f

def rate_proliferation(e2f, kd_e2f, n_e2f=4.0):
  
  return HillCube(e2f, kd_e2f, n_e2f, normalized=True)

def rate_differentiation(myod, lin35, kd_myod, n_myod=4.0):
  """
  
  """
  return HillCube(myod*lin35, kd_myod, n_myod, normalized=True)
  
