# Model free 3 component scattering power decomposition for dual co-polarimetric (HH--VV) SAR data

Model free 3-component scattering power decomposition using dual co-polarimetric SAR data. The program file types available are: ENVI IDL, MATLAB and C

## C usage

To compile and run with T2 matrix data in T2 folder and analysis window size 3:
```
   gcc MF3CD.c -o MF3CD.exe -lm
   ./MF3CD.exe T2 3
```

# Reference

 Subhadip Dey, Narayanarao Bhogapurapu, Avik Bhattacharya, Dipankar Mandal, Juan M. Lopez-Sanchez, Heather McNairn & Alejandro C. Frery (2021) Rice phenology mapping using novel target characterization parameters from polarimetric SAR data, International Journal of Remote Sensing, 42:14, 5515-5539, DOI: [10.1080/01431161.2021.1921876](https://doi.org/10.1080/01431161.2021.1921876)
