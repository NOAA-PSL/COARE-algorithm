# Python implementation of COARE air/sea flux algorithm version 3.6.

This includes functions for bulk flux calculations:
- Without warm layer computations [coare36vn\_zrf\_et.py](https://github.com/noaa-psd/COARE-algorithm/blob/feature/sanAkel/decorated_doc/Python/COARE3.6/coare36vn_zrf_et.py).
- With warm layer computations [coare36vnWarm\_et.py](https://github.com/noaa-psd/COARE-algorithm/blob/feature/sanAkel/decorated_doc/Python/COARE3.6/coare36vnWarm_et.py).

The python codes were translated from the MATLAB scripts. They can be run over the same input data set [test\_36\_data.txt](https://github.com/noaa-psd/COARE-algorithm/blob/feature/sanAkel/decorated_doc/Python/COARE3.6/test_36_data.txt) that is used to exercise the MATLAB code. Output with and without wave effects is included in [test\_36\_output\_withwavesinput\_withwarmlayer.txt](https://github.com/noaa-psd/COARE-algorithm/blob/feature/sanAkel/decorated_doc/Python/COARE3.6/test_36_output_withwavesinput_withwarmlayer.txt) and [test\_36\_output\_withnowavesinput\_withwarmlayer.txt](https://github.com/noaa-psd/COARE-algorithm/blob/feature/sanAkel/decorated_doc/Python/COARE3.6/test_36_output_withnowavesinput_withwarmlayer.txt) respectively.

## Instructions
- For the bulk flux calculations without warm layer computations, run: `coare36vn_zrf_et.py` from the iPython command line. Edit line `959` to set path to test data file: `test_36_data.txt`. 
- With warm layer computations, run: `coare36vnWarm_et.py`. Edit line 403 to set path to test data file: `test_36_data.txt`. 

Depending if the waves parameters: `cp` and `sigH` are used as input to COARE3.6, this will output a file of results that you can compare to the ones provided (`test_36_output_withnowavesinput_withwarmlayer.txt`, `test_36_output_withwavesinput_withwarmlayer.txt`).  

## Sample output
This file contains a time series of following variables.


| Variable Name | Description | Notes |
| :------------ | :---------: | ----: |
| usr           | friction velocity that includes gustiness (m/s) | Denoted by u* |
| tau           | wind stress that includes gustiness (N/m^2)| |
| hsb           | sensible heat flux (W/m^2) | positive for $T_{air}$ < $T_{skin}$ |
| hlb           | latent heat flux (W/m^2) | positive for $q_{air} < q_s$
| hbb           | atmospheric buoyany flux (W/m^2) | positive when `hlb` and `hsb` heat the atmosphere |
| hsbb          | atmospheric buoyancy flux from sonic | as above, computed with sonic anemometer `T` |
| hlwebb        | Webb factor to be added to `hl` | covariance and ID latent heat fluxes |
| tsr           | temperature scaling parameter (K) | Denoted by t* |             
| qsr           | specific humidity scaling parameter (kg/kg) | Denoted by q* |
| zo | momentum roughness length (m)| |
| zot | thermal roughness length (m)| |
| zoq | moisture roughness length (m)| |
| Cd | wind stress transfer (drag) coefficient at height zu (unitless)| |
| Ch |sensible heat transfer coefficient (Stanton number) at height zu (unitless)| |
| Ce | latent heat transfer coefficient (Dalton number) at height zu (unitless)| |
| L | Monin-Obukhov length scale (m)
| zeta | Monin-Obukhov stability parameter zu/L (dimensionless)| |
| dT_skin | cool-skin temperature depression (degC)| positive value means skin is cooler than subskin|
| dq_skin | cool-skin humidity depression (g/kg)| |
| dz_skin | cool-skin thickness (m)| |
| Urf | wind speed at reference height |user can select height at input |
| Trf | air temperature at reference height| |
| Qrf |air specific humidity at reference height| |
| RHrf | air relative humidity at reference height| |
| UrfN | neutral value of wind speed at reference height| |
| TrfN | neutral value of air temp at reference height | |
| qarfN | neutral value of air specific humidity at reference height | |
| lw_net | Net IR radiation computed by COARE (W/m2) | positive heating ocean |
| sw_net | Net solar radiation computed by COARE (W/m2) | positive heating ocean | 
| Le | latent heat of vaporization (J/K)| |
| rhoa | density of air at input parameter height zt (kg/m3)|typically same as zq |
| UN | neutral value of wind speed at zu (m/s)| |
| U10 | wind speed adjusted to 10 m (m/s)| |
| UN10 | neutral value of wind speed at 10m (m/s)| |
| Cdn_10 | neutral value of drag coefficient at 10m (unitless)| |
| Chn_10 | neutral value of Stanton number at 10m (unitless)| |
| Cen_10 | neutral value of Dalton number at 10m (unitless)| |
| hrain | rain heat flux (W/m^2) | positive cooling ocean |
| Qs | sea surface specific humidity, assuming saturation (g/kg)| |
| Evap | evaporation rate (mm/h)| |
| T10 | air temperature at 10m (deg C)| |
| T10N | neutral air temperature at 10m (deg C) | |
| Q10 | air specific humidity at 10m (g/kg) | |
| Q10N | neutral air specific humidity at 10m (g/kg) | |
| RH10 | air relative humidity at 10m (unitless) | |
| P10 | air pressure at 10m (mb) | |
| rhoa10 | air density at 10m (kg/m3) | |
| gust | gustiness velocity (m/s) | |
| wc_frac | whitecap fraction (ratio) | |
| Edis | energy dissipated by wave breaking (W/m^2) | |
| dT_warm | dT from base of warm layer to skin, i.e. warming across entire warm layer depth (deg C) | |
| dz_warm | warm layer thickness (m) | |
| dT\_warm\_to\_skin | dT from measurement depth to skin due to warm layer, such that $T_{skin} = T_{sea} + dT_{warm_to_skin} - dT_{skin}$ | |
| du_warm | total current accumulation in warm layer (m/s) | |
