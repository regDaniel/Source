''' Plot hydraulic conductivity as a function of porosity and soil water
content. '''

import numpy as np
import matplotlib.pyplot as plt

plt.style("ggplot") # ggplot style 

# Parameter (for clay)
air_dp =
k_sat =
k1 = 
porosity_0
# root depth (forest) [m]
root_dpt = 1.0
#tuning shape parameter
f = 1.0 


# Conductivity
def conductivity(w_content, porosity):
    
    frac = (porosity - w_content)/(porosity - air_dp)
    k = np.exp(frac)
    
    return k
    

# porosity (exponential...)
z = np.arange(0, 10, 100)
porosity = np.zeros(100, dtype=float)
porosity =  porosity_0 * np.exp(-f_tun*(z - root_dp))

# water content [air_dp, porosity]
w_content_max = porosity # array
w_content_min = air_dp # scalar



hydraulic_conductivity = conductivity(w_content, porosity)