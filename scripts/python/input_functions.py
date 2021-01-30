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

def rate_lin35_prod(myod, k_myod_lin35, n_myod_lin35=4.0):
  
  return HillCube(myod, k_myod_lin35, n_myod_lin35, normalized=True)

def rate_lin35_phos(lin35, cyd1, cki1, km_lin35, kd_cyd1, kd_cki1, n_cyd1=4.0, n_cki1=4.0):
  # probability of positive cell cycle regulator binding
  pos_cycle_reg_present = HillCube(cyd1, kd_cyd1, n_cyd1, normalized=True)
  
  # probability of negative cell cycle regulaotr unbinding
  neg_cycle_reg_absent = 1.0 - HillCube(cki1, kd_cki1, n_cki1, normalized=True)
  
  lin35_phos_active = pos_cycle_reg_present * neg_cycle_reg_absent
  
  lin35_phos = cyd1 * lin35 / (lin35 + km_lin35)
  
  # postive regulaotr ON AND NOT negative regulator
  return lin35_phos_active * lin35_phos


def rate_e2f_prod(lin35, kd_lin35_e2f, n_lin35_e2f=4.0):
  
  return 1 - HillCube(lin35, kd_lin35_e2f, n_lin35_e2f, normalized=True)

def rate_proliferation(e2f, kd_e2f, n_e2f=4.0):
  
  return HillCube(e2f, kd_e2f, n_e2f, normalized=True)

def rate_differentiation(myod, kd, n=4):
  return HillCube(myod, kd, n, normalized=True)
