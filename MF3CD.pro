PRO MF3CD

  ; written by Subhadip Dey
  ; Research scholar
  ; Centre for Studies in Resources Engineering
  ; Indian Institute of Technology Bombay
  ; email: sdey2307@gmail.com
  
  path = "E:\urban_mapping\RS-2_SF_FP\T3\" ; change this path for your dataset
  
  T11_path = path + "T11.bin"
  T12_imag_path = path + "T12_imag.bin"
  T12_real_path = path + "T12_real.bin"
  T22_path = path + "T22.bin"
  
  e = ENVI()
  T11 = e.OpenRaster(T11_path)
  T12_imag = e.OpenRaster(T12_imag_path)
  T12_real = e.OpenRaster(T12_real_path)
  T22 = e.OpenRaster(T22_path)
  
  T11_dat = T11.GetData(BANDS=[0])
  T12_imag_dat = T12_imag.GetData(BANDS=[0])
  T12_real_dat = T12_real.GetData(BANDS=[0])
  T22_dat = T22.GetData(BANDS=[0])
  
  T12_dat = COMPLEX(T12_real_dat, T12_imag_dat)
  T21_dat = CONJ(T12_dat)
  
  dims = SIZE(T11_dat)
  nrow = dims[1] ; reads row of the raster
  ncol = dims[2] ; reads column of the raster
  
  ; to make an average of image using 7 x 7 window size
  ; change the window size if required
  
  T11_dat_mean = MEAN_FILTER(T11_dat, 7, $
     /ARITHMETIC)
  T12_dat_mean = MEAN_FILTER(T12_dat, 7, $
     /ARITHMETIC)
  T21_dat_mean = MEAN_FILTER(T21_dat, 7, $
     /ARITHMETIC)
  T22_dat_mean = MEAN_FILTER(T22_dat, 7, $
     /ARITHMETIC)
  
  det_T = REAL_PART(T11_dat_mean*T22_dat_mean-T12_dat_mean*T21_dat_mean)
  tr_T = T11_dat_mean + T22_dat_mean
  
  dop_b = SQRT(1-((4*det_T)/(tr_T*tr_T))) ; computation of Barakat degree of polarization
  dop_b_raster = e.CreateRaster("E:\urban_mapping\Dop_Barakat.bin", dop_b, NBANDS=1)
  dop_b_raster.Save
  
  val = (dop_b*tr_T*(T11_dat_mean - T22_dat_mean))/(T11_dat_mean*T22_dat_mean+dop_b*dop_b*tr_T*tr_T)
  theta_d = 180/!PI*ATAN(val)
  theta_raster = e.CreateRaster("E:\urban_mapping\theta_dcp.bin", theta_d, NBANDS=1)
  theta_raster.Save
  
  ps = ((dop_b*tr_T)/2)*(1 + SIN(2*(theta_d*!PI/180)))
  ps_raster = e.CreateRaster("E:\urban_mapping\ps_dcp.bin", ps, NBANDS=1)
  ps_raster.Save
  
  pd = ((dop_b*tr_T)/2)*(1 - SIN(2*(theta_d*!PI/180)))
  pd_raster = e.CreateRaster("E:\urban_mapping\pd_dcp.bin", pd, NBANDS=1)
  pd_raster.Save
  
  pv = tr_T*(1-dop_b)
  pv_raster = e.CreateRaster("E:\urban_mapping\pv_dcp.bin", pv, NBANDS=1)
  pv_raster.Save

  e.Close
  print, "Done"
END
