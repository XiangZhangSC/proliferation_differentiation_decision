from input_functions import rate_hlh1_prod, rate_fos1_prod, rate_cyd1_prod, \
  rate_cki1_prod, rate_e2f_prod, rate_lin35_prod, rate_jun1_prod

def pdd(x,t, 
        tau_mls2, 
        tau_hlh1, k_mls2_hlh1, k_myod_hlh1, 
        tau_fos1, k_myod_fos1, 
        tau_cyd1, k_ap1_cyd1, 
        tau_cki1, k_myod_cki1, 
        km_e2f, k_cki1, 
        tau_lin35, k_myod_lin35, 
        tau_jun1, k_myod_jun1):
  # mls-2
  tau_mls2 = tau_mls2
  if 4.0 <= t <= 10.0:
    mls2_in = 1.0
  else:
    mls2_in = 0.0
  
  # hlh-1
  tau_hlh1 = tau_hlh1
  k_mls2_hlh1 = k_mls2_hlh1
  k_myod_hlh1 = k_myod_hlh1
  hlh1_in = rate_hlh1_prod(x[0], x[1], k_mls2_hlh1, k_myod_hlh1)
  
  # fos-1
  tau_fos1 = tau_fos1
  k_myod_fos1 = k_myod_fos1
  fos1_in = rate_fos1_prod(x[1], x[2], x[7], k_myod_fos1)
  
  # cyd-1
  tau_cyd1 = tau_cyd1
  ap1 = min(x[2], x[7])
  k_ap1_cyd1 = k_ap1_cyd1
  cyd1_in = rate_cyd1_prod(ap1, k_ap1_cyd1)
  
  # cki-1
  tau_cki1 = tau_cki1
  k_myod_cki1 = k_myod_cki1
  cki1_in = rate_cki1_prod(x[1], k_myod_cki1)
  
  # E2F
  k_cki1 = k_cki1
    
  km_e2f = km_e2f
  e2f_in = rate_e2f_prod(e2f=x[5], cyd1=x[3], cki1=x[4], lin35=x[6], myod=x[1],
                         km_e2f=km_e2f, k_cki1=k_cki1)
  
  # lin-35
  tau_lin35 = tau_lin35
  k_myod_lin35 = k_myod_lin35
  lin35_in = rate_lin35_prod(myod=x[1], e2f=x[5], cyd1=x[3], cki1=x[4], lin35=x[6],  
                    k_myod_lin35=k_myod_lin35, km_e2f=km_e2f, k_cki1=k_cki1)
  
  # jun-1
  tau_jun1 = tau_jun1
  k_myod_jun1 = k_myod_jun1
  jun1_in = rate_jun1_prod(x[1], x[7], x[2], k_myod_jun1)
  
  # ODEs
  dxdt = [0, 0, 0, 0, 0, 0, 0, 0]
  
  dxdt[0] = tau_mls2 * (mls2_in - x[0]) # mls-2
  dxdt[1] = tau_hlh1 * (hlh1_in - x[1]) # hlh-1/MyoD
  dxdt[2] = tau_fos1 * (fos1_in - x[2]) # fos-1
  dxdt[3] = tau_cyd1 * (cyd1_in - x[3]) # cyd-1
  dxdt[4] = tau_cki1 * (cki1_in - x[4]) # cki-1
  dxdt[5] = e2f_in - x[5]*x[6] # E2F (E2F either is free or is bound with pRb, no degredation)
  dxdt[6] = tau_lin35 * (lin35_in - x[6]) # lin-35
  dxdt[7] = tau_jun1 * (jun1_in - x[7]) # jun-1
  
  return dxdt