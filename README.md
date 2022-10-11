# COARE-algorithm
Repository of the TOGA-COARE bulk air-sea flux algorithm in:
- Python
- MATLAB
- FORTRAN

## Origins
The international TOGA-COARE field program which took place in the western Pacific warm pool over 4 months from November 1992 to February 1993 ([Fairall et al. 1996a](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/Fairall%20et%20al.%201996a%20-%20cool%20skin%20warm%20layer.pdf), [1996b](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/Fairall%20et%20al.%201996b%20-%20bulk%20fluxes%20of%20variables.pdf) and [1997](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/Fairall%20et%20al.%201997%20-%20ship%20measurements%20MABL.pdf)) spurred the development of the COARE model. The algorithm is intended to provide estimates of `momentum, sensible heat`, and `latent heat fluxes` using inputs of bulk atmospheric variables (`wind speed, SST, air temperature, air humidity`). The algorithm contains subroutines/functions to handle near-surface gradients of temperature in the ocean.

## Versions
- **Version 2.5** was published in 1996. 
- **Version 3.0** was published in 2003; it was a major update from Version 2.5. This update was based on new observations at higher wind speeds (10 to 20 m/s). Available in MATLAB and FORTRAN only.
- **Version 3.5** was released in 2013 following the publication of [Edson et al. 2013](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/Edson%20et%20al.%202013%20-%20momentum%20flux.pdf), which made adjustments to the wind speed dependence of the Charnock parameter based on a large database of direct covariance stress observations (principally from a buoy). This led to an increase in stress for wind speeds greater than about 18 m/s. The roughness Reynolds number formulation of the scalar roughness length was tuned slightly to give the same values of `Ch` and `Ce` as Version 3.0. The diurnal warm layer model was structured as a separate routine instead of embedded in a driver program. COARE 3.5 was based on Edson’s buoy data ([Edson et al. 2013](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/Edson%20et%20al.%202013%20-%20momentum%20flux.pdf)) and was compared to a large database (a total of 16,000 hours of observations) combining observations from NOAA, WHOI, and U. Miami ([Fairall et al. 2011](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/Fairall%20et%20al.%202011%20-%20COAREG.pdf)). It is available in Python and MATLAB.
- **Version 3.6** is slightly restructured and built around improvements in the representation of the effects of waves on fluxes. This includes improved relationships of `surface roughness`, $z_o$, and `whitecap fraction`, $W_f$, on wave parameters.  More details can be found in [coare3.6\_readme\_1.pdf](https://github.com/noaa-psd/COARE-algorithm/blob/master/References/coare36_readme_1.pdf). This version is available in Python, MATLAB and FORTRAN.


# References:

*    Edson, J.B., J. V. S. Raju, R.A. Weller, S. Bigorre, A. Plueddemann, C.W. Fairall, S. Miller, L. Mahrt, Dean Vickers, and Hans Hersbach, 2013: On the Exchange of momentum over the open ocean. J. Phys. Oceanogr., 43, 1589–1610. doi: http://dx.doi.org/10.1175/JPO-D-12-0173.1

*    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson, and G.S. Young, 1996a: The cool skin and the warm layer in bulk flux calculations. J. Geophys. Res. 101, 1295-1308. https://doi.org/10.1029/95JC03190

*    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson, G.S. Young, 1996b: Bulk parameterization of air-sea fluxes for TOGA COARE. J. Geophys. Res. 101, 3747-3764. https://doi.org/10.1029/95JC03205

*    Fairall, C. W., White, A. B., Edson, J. B., and Hare, J. E.: Integrated Shipboard Measurements of the Marine Boundary Layer, J. Atmos. Ocean. Tech., 14, 338–359, 1997. https://doi.org/10.1175/1520-0426(1997)014<0338:ISMOTM>2.0.CO;2

*    Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson, 2003: Bulk parameterization of air-sea fluxes: Updates and verification for the COARE algorithm. J. Climate 16, 571-591. https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2

*    Fairall, C.W., Mingxi Yang, Ludovic Bariteau, J.B. Edson, D. Helmig, W. McGillis, S. Pezoa, J.E. Hare, B. Huebert, and B. Blomquist, 2011: Implementation of the COARE flux algorithm with CO2, DMS, and O3. J. Geophys. Res., 116, C00F09, https://doi.org/10.1029/2010JC006884
