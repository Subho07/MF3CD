# Model free 3 component scattering power decomposition for dual co-polarimetric (HH--VV) SAR data

Model free 3-component scattering power decomposition using dual co-polarimetric SAR data. The program file types available are: ENVI IDL, MATLAB and C

## C usage

To compile and run with T2 matrix data in T2 folder and analysis window size 3:
```
   gcc MF3CD.c -o MF3CD.exe -lm
   ./MF3CD.exe T2 3
```
