load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; order of loading is important
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------------
begin

   ; time control
   season = ("JJA")
   year1 = "1989"
   year2 = "1998"

   ; file paths
   simulation  = "fivefold_gamma_20y"
   datafolder  = "$HOME/Data/landflux_dataset/"
   datafile    = datafolder+"LandFluxEVAL.merged.89-05.monthly.diagnostic.regridded.nc"
   modelfolder = "$HOME/Data/cosmo5_validation/"+simulation+"/"
   modelfile   = modelfolder +"monmean04_"+year1+"-"+year2+"_"+season+"_regridded.nc"

   ; open ps output file
   outfolder = "$HOME/Results/landflux/"
   figname = outfolder+"bias_"+simulation+"_"+year1+"-"+year2
   wks     = gsn_open_wks("ps",figname)

   
   in_model     = addfile(modelfile, "r")
   lat= (/in_model->lat/)
   lat@units = "degrees_north"
   lon     = (/in_model->lon/)
   lon@units = "degrees_east"
   latent_heatflux = in_model->ALHFL_S
;   latent_heatflux!0 = "time"
   

   lhevap = 2.501e06 ; latent heat of evaporation (wikipedia)
   conversion_fac = -60.*60.*24./lhevap ; conversion factor (-time/lhevap)
   evap_model = latent_heatflux*conversion_fac

   in_obs = addfile(datafile, "r")
   evap_obs = in_obs->ET_mean

   evap_diff = evap_model - evap_obs
   printVarSummary(evap_diff)
   evap_bias = dim_avg_n(evap_diff,0)
   printVarSummary(evap_bias)
   evap_bias!1 = "lon"
   evap_bias!0 = "lat"
   evap_bias&lon = lon
   evap_bias&lat = lat

       
    
   ; Set resources for plot
   res2d                     = True
   
   res2d@pmLegendDisplayMode = "Always"              ; turn on legend
   
   res2d                      = True          ; plot mods desired
   res2d@gsnFrame             = False          ; False for overlay
   res2d@gsnDraw              = False          ; False for overlay
   res2d@gsnMaximize          = True          ; Maximize plot in frame.
   res2d@gsnSpreadColors      = True
   
   res2d@cnFillOn             = True          ; Turn on contour fill
   ;res2d@cnFillMode           = "RasterFill"  ; Turn on raster fill
   res2d@cnFillMode           = "CellFill"  ; Turn on raster fill
   res2d@cnLinesOn            = False         ; Turn off contour lines
   res2d@cnLineLabelsOn       = False    ; turns off contour line labels
   res2d@cnLevelSelectionMode = "ManualLevels"
     ;  res2d@cnMinLevelValF       =  0
     ;  res2d@cnMaxLevelValF       =  1
     ;  res2d@cnLevelSpacingF      =  0.1
    
     ;  res2d@lbLabelBarOn         = True
   res2d@lbOrientation        = "vertical"     ; vertical label bar
   
   
   res2d@mpDataBaseVersion = "MediumRes"
   res2d@mpOutlineBoundarySets = "Geophysical"
   ;res2d@mpDataBaseVersion     = "Ncarg4_1"
   res2d@mpDataSetName         = "Earth..4"
   
   ;res2d@mpProjection          = "CylindricalEquidistant"
   res2d@mpProjection          = "Mercator"
   res2d@mpFillOn              = False
   res2d@mpGridAndLimbOn       = True             ; turn on lat/lon lines
   res2d@mpPerimOn             = False            ; turn off box around plot
   res2d@mpGeophysicalLineThicknessF = 2
   res2d@mpNationalLineThicknessF    = 2
   
   res2d@mpCenterLonF          = 10
   res2d@mpCenterLatF          = 51
   res2d@mpLimitMode          = "LatLon"
   
   res2d@mpMinLatF            = 24.
   res2d@mpMaxLatF            = 75.
   res2d@mpMinLonF            = -12.
   res2d@mpMaxLonF            = 36.
   
   res2d@gsnAddCyclic = False 
   
   res2d@tiMainString          = "Evapotranspiration Bias"
   res2d@mpMinLatF             = 28.
   res2d@mpMaxLatF             = 71.
   res2d@mpMinLonF             = -10.
   res2d@mpMaxLonF             = 30.
   
   gsn_define_colormap(wks,"BlWhRe")  ; choose a colormap
   gsn_reverse_colormap(wks)           ; Reverse the color map.
   
   res2d@cnMinLevelValF       =  -3.
   res2d@cnMaxLevelValF       =  3.
   res2d@cnLevelSpacingF      = 0.1
   res2d@gsnLeftString        = simulation +" "+season+" "+year1+"-"+year2
   res2d@lbTitleOn        = True                  ; turn on title
   res2d@lbTitleString    = "mm/d"
   res2d@lbTitleFontHeightF= .015                 ; make title smaller
;   labels                     = fspan(-15., 15., 101)
;   res2d@lbLabelStrings       = (/"-15","0","15"/)
   res2d@lbLabelAutoStride          = True
;   res2d@gsnRightString       = "mm/d"
   ;delete(res2d@cnFillColors)
   
   
   plot = gsn_csm_contour_map(wks,evap_bias,res2d)
   draw(plot)
   frame(wks)
    

   print(figname+".ps")

end
