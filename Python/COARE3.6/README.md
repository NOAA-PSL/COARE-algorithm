Python implementation of COARE air/sea flux algorithm version 3.6.

This includes the functions for COARE model version 3.6 bulk flux calculations, coare36vn_zrf_et.py,
and the functions for the bulk flux calculations with warm layer computations coare36vnWarm_et.py

The python codes are translated from the MATLAB scripts and run the same input data set 'test_36_data.txt' used to exercise the matlab code. 

For the bulk flux calculations, execute 'run coare36vn_zrf_et.py' from the iPython command line. Edit line 959 to indicate path to test data file 'test_36_data.txt'. 
For the bulk flux calculations with warm layer computations, execute 'run coare36vnWarm_et.py' from the iPython command line. Edit line 403 to indicate path to test data file 'test_36_data.txt'. 

Depending if the waves parameters cp and sigH are used, this will output a file of results that you can compare to the 'test_36_output_withnowavesinput_withwarmlayer.txt' and 'test_36_output_withwavesinput_withwarmlayer.txt' provided.  
The file contains a time series of flux variables: 
usr	tau	hsb	hlb	hbb	hlwebb	tsr	qsr	zo	zot	zoq	Cd	Ch	Ce	L	zeta	dT_skinx	dq_skinx	dz_skin	Urf	Trf	Qrf	RHrf	UrfN	TrfN	QrfN	lw_net	sw_net	Le	rhoa	UN	U10	U10N	Cdn_10	Chn_10	Cen_10	hrain	Qs	Evap	T10	T10N	Q10	Q10N	RH10	P10	rhoa10	gust	wc_frac	Edis	dT_warm	dz_warm	dT_warm_to_skin	du_warm


