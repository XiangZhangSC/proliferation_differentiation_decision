from input_functions import rate_hlh1_prod, rate_fos1_prod, rate_cyd1_prod, \
  rate_cki1_prod, rate_lin35_prod, rate_lin35_phos, rate_e2f_prod

def pdd(x,t, 
        tau_mls2, 
        tau_hlh1, k_mls2_hlh1, k_myod_hlh1, 
        tau_fos1, k_myod_fos1, 
        tau_cyd1, k_fos1_cyd1, 
        tau_cki1, k_myod_cki1, 
        tau_lin35, k_myod_lin35, 
        km_lin35, k_cyd1, k_cki1, 
        tau_e2f, k_lin35_e2f):
  # mls-2
  tau_mls2 = tau_mls2
  if 6 <= t <= 14:
    mls2_in = 1
  else:
    mls2_in = 0
  
  # hlh-1
  tau_hlh1 = tau_hlh1
  k_myod_hlh1 = k_myod_hlh1 + x[3]
  hlh1_in = rate_hlh1_prod(x[0], x[1], k_mls2_hlh1, k_myod_hlh1)
  
  # fos-1
  tau_fos1 = tau_fos1
  k_myod_fos1 = k_myod_fos1 + x[3]
  fos1_in = rate_fos1_prod(x[1], k_myod_fos1)
  
  # cyd-1
  tau_cyd1 = tau_cyd1
  cyd1_in = rate_cyd1_prod(x[2], k_fos1_cyd1)
  
  # cki-1
  tau_cki1 = tau_cki1
  k_myod_cki1 = k_myod_cki1 + x[3]
  cki1_in = rate_cki1_prod(x[1], k_myod_cki1)
  
  # lin-35 (hypophosphorylated, inhibiting proliferation)
  tau_lin35 = tau_lin35
  k_myod_lin35 = k_myod_lin35 + x[3]
  lin35_in = rate_lin35_prod(x[1], k_myod_lin35)
  
  lin35_phos = rate_lin35_phos(x[5], x[3], x[4], km_lin35, k_cyd1, k_cki1)
  
  # E2F
  tau_e2f = tau_e2f
  e2f_in = rate_e2f_prod(x[5], k_lin35_e2f)
  
  # ODEs
  dxdt = [0, 0, 0, 0, 0, 0, 0]
  
  dxdt[0] = 1.0/tau_mls2 * (mls2_in - x[0]) # mls-2
  dxdt[1] = 1.0/tau_hlh1 * (hlh1_in - x[1]) # hlh-1/MyoD
  dxdt[2] = 1.0/tau_fos1 * (fos1_in - x[2]) # fos-1
  dxdt[3] = 1.0/tau_cyd1 * (cyd1_in - x[3]) # cyd-1
  dxdt[4] = 1.0/tau_cki1 * (cki1_in - x[4]) # cki-1
  dxdt[5] = 1.0/tau_lin35 * (lin35_in - lin35_phos - x[5]) # lin-35
  dxdt[6] = 1.0/tau_e2f * (e2f_in - x[6]) # E2F
  
  return dxdt