import os
from cdo import *

os.chdir("/home/redaniel/Data")

cdo = Cdo()
# model
cdo.seldate("1989-06-01,1998-08-31", input="cosmo5_validation/fivefold_gamma_20y/monmean/monmean_out04.nc", output="tmp1.nc")
cdo.selseas("JJA",input="tmp1.nc", output="tmp2.nc")
cdo.selname("ALHFL_S", input="tmp2.nc", output="tmp3.nc")
cdo.remapbil("grid_1deg.txt", input="tmp3.nc", output="cosmo5_validation/fivefold_gamma_20y/monmean04_1989-1998_JJA_regridded.nc")

# observations
cdo.selseas("JJA", input="landflux_dataset/LandFluxEVAL.merged.89-05.monthly.diagnostic.nc", output="tmp4.nc")
cdo.seldate("1989-06-01,1998-08-31", input="tmp4.nc", output="tmp5.nc")
cdo.sellonlatbox("-20.5,44.5,25.5,70.5", input="tmp5.nc", output="landflux_dataset/LandFluxEVAL.merged.89-05.monthly.diagnostic.regridded.nc")

os.system("rm tmp*.nc")

