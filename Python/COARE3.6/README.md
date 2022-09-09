Python implementation of COARE air/sea flux algorithm version 3.6.

This includes the functions for COARE model version 3.6 bulk flux calculations, coare36vn_zrf_et.py,
and the functions for the bulk flux calculations with warm layer computations coare36vnWarm_et.py

The python codes are translated from the MATLAB scripts and run the same input data set 'test_36_data.txt' used to exercise the matlab code. 

For the bulk flux calculations, execute 'run coare36vn_zrf_et.py' from the iPython command line. Edit line 959 to indicate path to test data file 'test_36_data.txt'. 
For the bulk flux calculations with warm layer computations, execute 'run coare36vnWarm_et.py' from the iPython command line. Edit line 403 to indicate path to test data file 'test_36_data.txt'. 
Depending if the waves parameters cp and sigH are used as input to COARE3.6, this will output a file of results that you can compare to the 'test_36_output_withnowavesinput_withwarmlayer.txt' and 'test_36_output_withwavesinput_withwarmlayer.txt' provided.  
The file contains a time series of flux variables: 
*    usr = friction velocity that includes gustiness (m/s), u*
*    tau = wind stress that includes gustiness (N/m^2)
*    hsb = sensible heat flux (W/m^2) ... positive for Tair < Tskin
*    hlb = latent heat flux (W/m^2) ... positive for qair < qs
*    hbb = atmospheric buoyany flux (W/m^2)... positive when hlb and hsb heat the atmosphere
*    hsbb = atmospheric buoyancy flux from sonic ... as above, computed with sonic anemometer T
*    hlwebb = webb factor to be added to hl covariance and ID latent heat fluxes
*    tsr = temperature scaling parameter (K), t*
*    qsr = specific humidity scaling parameter (kg/kg), q*
*    zo = momentum roughness length (m)
*    zot = thermal roughness length (m)
*    zoq = moisture roughness length (m)
*    Cd = wind stress transfer (drag) coefficient at height zu (unitless)
*    Ch = sensible heat transfer coefficient (Stanton number) at height zu (unitless)
*    Ce = latent heat transfer coefficient (Dalton number) at height zu (unitless)
*    L = Monin-Obukhov length scale (m)
*    zeta = Monin-Obukhov stability parameter zu/L (dimensionless)
*    dT_skin = cool-skin temperature depression (degC), pos value means skin is cooler than subskin
*    dq_skin = cool-skin humidity depression (g/kg)
*    dz_skin = cool-skin thickness (m)
*    Urf = wind speed at reference height (user can select height at input)
*    Trf = air temperature at reference height
*    Qrf = air specific humidity at reference height
*    RHrf = air relative humidity at reference height
*    UrfN = neutral value of wind speed at reference height
*    TrfN = neutral value of air temp at reference height
*    qarfN = neutral value of air specific humidity at reference height
*    lw_net = Net IR radiation computed by COARE (W/m2)... positive heating ocean
*    sw_net = Net solar radiation computed by COARE (W/m2)... positive heating ocean
*    Le = latent heat of vaporization (J/K)
*    rhoa = density of air at input parameter height zt, typically same as zq (kg/m3)
*    UN = neutral value of wind speed at zu (m/s)
*    U10 = wind speed adjusted to 10 m (m/s)
*    UN10 = neutral value of wind speed at 10m (m/s)
*    Cdn_10 = neutral value of drag coefficient at 10m (unitless)
*    Chn_10 = neutral value of Stanton number at 10m (unitless)
*    Cen_10 = neutral value of Dalton number at 10m (unitless)
*    hrain = rain heat flux (W/m^2)... positive cooling ocean
*    Qs = sea surface specific humidity, i.e. assuming saturation (g/kg)
*    Evap = evaporation rate (mm/h)
*    T10 = air temperature at 10m (deg C)
*    T10N = neutral air temperature at 10m (deg C)
*    Q10 = air specific humidity at 10m (g/kg)
*    Q10N = neutral air specific humidity at 10m (g/kg)
*    RH10 = air relative humidity at 10m (#)
*    P10 = air pressure at 10m (mb)
*    rhoa10 = air density at 10m (kg/m3)
*    gust = gustiness velocity (m/s)
*    wc_frac = whitecap fraction (ratio)
*    Edis = energy dissipated by wave breaking (W/m^2)
*    dT_warm = dT from base of warm layer to skin, i.e. warming across entire warm layer depth (deg C)
*    dz_warm = warm layer thickness (m)
*    dT_warm_to_skin = dT from measurement depth to skin due to warm layer, such that Tskin = tsea + dT_warm_to_skin - dT_skin
*    du_warm = total current accumulation in warm layer (m/s?)
