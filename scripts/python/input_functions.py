from toolbox import binding_prob

def rate_hlh1_prod(mls2, myod, max_hlh1_prod, 
                   k_mls2_hlh1, k_myod_hlh1, 
                   n_myod_hlh1=2, n_mls2_hlh1=2):
  
  p_mls2_bound_hlh1 = binding_prob(mls2, k_mls2_hlh1, n_mls2_hlh1)
  p_myod_bound_hlh1 = binding_prob(myod, k_myod_hlh1, n_myod_hlh1)
  ## binding MyoD or MLS-2 will induce hlh-1 expression
  return max_hlh1_prod * (1 - (1 - p_mls2_bound_hlh1) * (1 - p_myod_bound_hlh1))

def rate_fos1_prod(myod, max_fos1_prod, k_myod_fos1, n_myod_fos1=2):
  
  p_myod_bound_fos1 = binding_prob(myod, k_myod_fos1, n_myod_fos1)
  
  return max_fos1_prod * (1 - p_myod_bound_fos1)

def rate_cyd1_prod(fos1, max_cyd1_prod, k_fos1_cyd1, n_fos1_cyd1=2):
  
  p_fos1_bound_cyd1 = binding_prob(fos1, k_fos1_cyd1, n_fos1_cyd1)
  
  return max_cyd1_prod * p_fos1_bound_cyd1

def rate_cki1_prod(myod, max_cki1_prod, k_myod_cki1, n_myod_cki1=2):
  
  p_myod_bound_cki1 = binding_prob(myod, k_myod_cki1, n_myod_cki1)
  
  return max_cki1_prod * p_myod_bound_cki1

def rate_proliferation(pos_reg, neg_reg, kd_pos, kd_neg, n_pos=4, n_neg=4):
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
  pos_cycle_reg_present = binding_prob(pos_reg, kd_pos, n_pos)
  
  # probability of negative cell cycle regulaotr unbinding
  neg_cycle_reg_absent = 1 - binding_prob(neg_reg, kd_neg, n_neg)
  
  # AND gate for positive and negative cell cycle regulators
  return pos_cycle_reg_present * neg_cycle_reg_absent

def rate_differentiation(myod, kd, n):
  return binding_prob(myod, kd, n)
