from input_functions import rate_hlh1_prod, rate_fos1_prod, rate_cyd1_prod, rate_cki1_prod

def pdd(x,t, 
        a_mls2, 
        a_hlh1, max_hlh1_prod, k_mls2_hlh1, k_myod_hlh1, 
        a_fos1, max_fos1_prod, k_myod_fos1, 
        a_cyd1, max_cyd1_prod, k_fos1_cyd1, 
        a_cki1, max_cki1_prod, k_myod_cki1):
  # mls-2
  a_mls2 = a_mls2
  
  # hlh-1
  a_hlh1 = a_hlh1
  hlh1_in = rate_hlh1_prod(x[0], x[1], max_hlh1_prod, k_mls2_hlh1, k_myod_hlh1)
  
  # fos-1
  a_fos1 = a_fos1
  fos1_in = rate_fos1_prod(x[1], max_fos1_prod, k_myod_fos1)
  
  # cyd-1
  a_cyd1 = a_cyd1
  cyd1_in = rate_cyd1_prod(x[2], max_cyd1_prod, k_fos1_cyd1)
  
  # cki-1
  a_cki1 = a_cki1
  cki1_in = rate_cki1_prod(x[1], max_cki1_prod, k_myod_cki1)
  
  # ODEs
  dxdt = [0, 0, 0, 0, 0]
  
  dxdt[0] = -a_mls2*x[0] # mls-2
  dxdt[1] = hlh1_in - a_hlh1*x[1] # hlh-1/MyoD
  dxdt[2] = fos1_in - a_fos1*x[2] # fos-1
  dxdt[3] = cyd1_in - a_cyd1*x[3] # cyd-1
  dxdt[4] = cki1_in - a_cki1*x[4] # cki-1
  
  return dxdt