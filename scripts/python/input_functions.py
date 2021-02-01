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

def rate_lin35_phos(e2f, km_e2f, e2f_tot=1.0):
  return (e2f_tot - e2f) / ((e2f_tot - e2f) + km_e2f)

def rate_proliferation(e2f, kd_e2f, n_e2f=4.0):
  
  return HillCube(e2f, kd_e2f, n_e2f, normalized=True)

def rate_differentiation(myod, kd_myod, n_myod=4):
  """
  Hypophosphorylated Rb is required for myogenesis and muscle-specific gene 
  expression through its interaction with MyoD
  """
  
  return HillCube(myod, kd_myod, n_myod, normalized=True)
  
