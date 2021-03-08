from input_functions import rate_hlh1_prod, rate_fos1_prod, rate_cyd1_prod, \
  rate_cki1_prod, rate_e2f_prod, prod_lin35_by_myod, rate_cye1_prod

def pdd(x,t, 
        tau_mls2, f_mls2_in, 
        tau_hlh1, k_mls2_hlh1, k_myod_hlh1, ko_hlh1,
        tau_fos1, k_myod_fos1, 
        tau_cyd1, k_fos1_cyd1,
        tau_cki1, k_myod_cki1, ko_cki1,
        km_e2f, 
        tau_lin35, k_myod_lin35, ko_lin35, 
        tau_cye1, k_e2f_cye1):
  # mls-2
  tau_mls2 = tau_mls2
  mls2_in = f_mls2_in(t)
  
  # hlh-1
  tau_hlh1 = tau_hlh1
  k_mls2_hlh1 = k_mls2_hlh1
  k_myod_hlh1 = k_myod_hlh1
  hlh1_in = rate_hlh1_prod(x[0], x[1], k_mls2_hlh1, k_myod_hlh1) * (1.0-ko_hlh1)
  
  # fos-1
  tau_fos1 = tau_fos1
  k_myod_fos1 = k_myod_fos1
  fos1_in = rate_fos1_prod(x[1], k_myod_fos1)
  
  # cyd-1
  tau_cyd1 = tau_cyd1
  k_fos1_cyd1 = k_fos1_cyd1
  cyd1_in = rate_cyd1_prod(x[2], k_fos1_cyd1) 
  
  # cki-1
  tau_cki1 = tau_cki1
  k_myod_cki1 = k_myod_cki1
  cki1_in = rate_cki1_prod(x[1], k_myod_cki1) * (1.0-ko_cki1)
  
  # E2F
  km_e2f = km_e2f
  e2f_in = rate_e2f_prod(e2f=x[5], cyd1=x[3], cki1=x[4], myod=x[1], cye1=x[7],
                         km_e2f=km_e2f)
  
  # lin-35
  tau_lin35 = tau_lin35
  k_myod_lin35 = k_myod_lin35
  lin35_in = 1.0 - (1.0 - prod_lin35_by_myod(x[1], k_myod_lin35)*(1.0 - ko_lin35)) * (1.0 - e2f_in)
  
  # cye-1
  tau_cye1 = tau_cye1
  k_e2f_cye1 = k_e2f_cye1
  cye1_in = rate_cye1_prod(x[5], k_e2f_cye1)
  
  # ODEs
  dxdt = [0, 0, 0, 0, 0, 0, 0, 0]
  
  dxdt[0] = tau_mls2 * (mls2_in - x[0]) # mls-2
  dxdt[1] = tau_hlh1 * (hlh1_in - x[1]) # hlh-1/MyoD
  dxdt[2] = tau_fos1 * (fos1_in - x[2]) # fos-1
  dxdt[3] = tau_cyd1 * (cyd1_in - x[3]) # cyd-1
  dxdt[4] = tau_cki1 * (cki1_in - x[4]) # cki-1
  dxdt[5] = e2f_in - x[5]*x[6] # E2F (E2F either is free or is bound with pRb, no degredation)
  dxdt[6] = tau_lin35 * (lin35_in - x[6]) # lin-35
  dxdt[7] = tau_cye1 * (cye1_in - x[7])
  
  return dxdt