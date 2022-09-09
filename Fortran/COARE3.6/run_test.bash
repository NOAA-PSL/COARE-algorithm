#!/bin/bash

echo "Running test_coare_36, driven by test cases in test_35_data.txt"

rm -f *mod *.o test_coare_36
ifort -no-wrap-margin -o test_coare_36 machine.f90 physcons.F90 module_coare36_parameters.f90 module_coare36_model.f90 test_coare_36.f90
nice ./test_coare_36 > log_coare_36_output


