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

def rate_proliferation(pos_reg, neg_reg, kd_pos, kd_neg, n_pos=4.0, n_neg=4.0):
  """
  rate is between 0 and 1 (activity)

  Parameters
  ----------
  pos_reg : TYPE
    Positive cell cycle regulator quantity
  kd_pos : TYPE
    Kd for positive cell cycle regulator
  n_pos : TYPE
    Hill coefficient for positive cell cycle regualtor
  neg_reg : TYPE
    Negative cell cycle regulator quantity
  kd_neg : TYPE
    Kd for negative cell cycle regulator
  n_neg : TYPE
    Hill coefficient for negative cell cycle regulator

  Returns
  -------
  None.

  """
  
  # probability of positive cell cycle regulator binding
  pos_cycle_reg_present = HillCube(pos_reg, kd_pos, n_pos, normalized=True)
  
  # probability of negative cell cycle regulaotr unbinding
  neg_cycle_reg_absent = 1.0 - HillCube(neg_reg, kd_neg, n_neg, normalized=True)
  
  # AND gate for positive and negative cell cycle regulators
  return pos_cycle_reg_present * neg_cycle_reg_absent

def rate_differentiation(myod, kd, n=4):
  return HillCube(myod, kd, n, normalized=True)
