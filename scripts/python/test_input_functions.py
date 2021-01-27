#%%
from input_functions import rate_fos1_prod

def test_rate_fos1_prod():
  """
  Test that rate_fos1_prod is a decreasing function

  Returns
  -------
  None.

  """
  max_fos1_prod = 1; k_myod_fos1 = 5
  low_myod = rate_fos1_prod(0.1, max_fos1_prod, k_myod_fos1)
  mid_myod = rate_fos1_prod(5, max_fos1_prod, k_myod_fos1)
  high_myod = rate_fos1_prod(10, max_fos1_prod, k_myod_fos1)
  
  tol = 1E-14
  
  assert low_myod > mid_myod > high_myod
  assert abs(mid_myod - 0.5) < tol

#%%  
from input_functions import rate_cyd1_prod  
  
def test_rate_cyd1_prod():
  """
  Test that rate_cyd1_prod is a increasing function

  Returns
  -------
  None.

  """
  max_cyd1_prod = 1; k_fos1_cyd1 = 5
  low_fos1 = rate_cyd1_prod(0.1, max_cyd1_prod, k_fos1_cyd1)
  mid_fos1 = rate_cyd1_prod(5, max_cyd1_prod, k_fos1_cyd1)
  high_fos1 = rate_cyd1_prod(10, max_cyd1_prod, k_fos1_cyd1)
  
  tol=1E-14
  
  assert low_fos1 < mid_fos1 < high_fos1
  assert abs(mid_fos1 - 0.5) < tol