"""
Utility functions gathered from various sources or custom.

V.1 Oct 2014, BWB
"""

def load_txt_file(path,fil,cols,sep='\t',hdr=True):
    """
    loads delimited text file of column variables
    header, if present, is assumed to be one line of var names only

    usage: data, varNames = load_txt_file(path,fil,cols)

    inputs: path: string variable for path to file
            fil: string variable of file name
            cols: number of data columns in file
            sep: delimiter character, tab is default
            hdr: boolean, True if file has header (default)

    output: data: 2-D array of data as column variables
            varNames: list of variable names
    """
    import os
    import numpy as np

    fpath = os.path.join(path,fil)
    fin = open(fpath)
    first = True
    data = np.zeros((0,cols))
    temp = np.zeros((1,cols))
    for line in fin:
        if first and hdr:
            cleanLine = line.strip()
            varNames = cleanLine.split(sep)
            first = False
            continue
        elif first and not hdr:
            varNames = []
            for jj in range(cols): # default var names
                varNames.append('var'+str(jj))
            first = False
        cleanLine = line.strip()
        fields = cleanLine.split(sep)
        for jj in range(cols):
            temp[0,jj] = np.float(fields[jj])
        data = np.row_stack((data,temp))
    fin.close()
    return (data, varNames)


def noaa_tcode(ddd, hh, code_str):
    """
    Converts NOAA serial data acquisition system time code to decimal
        julian date.

    Usage: jd = noaa_tcode(ddd, hh, code_str)

    Input strings: ddd = julian date, hh = hour string, and code_str = time code
        (e.g.1230423), where first two chr are minutes, next two are sec, and
        last three are msec

    Output is decimal julian date/time for given year

    """

    import numpy as np

    day = np.float(ddd)
    hr = np.float(hh)
    min = np.float(code_str[0:2])
    sec = np.float(code_str[2:4])
    msec = np.float(code_str[4:])/1000
    jd = day + (hr + (min + (sec + msec)/60)/60)/24;
    return jd


def find(b):
    """
    Usage: idx = find(b) - Returns sorted array of indices where boolean
    input array b is true.  Similar to MATLAB find function.

    Input may be a 1-D boolean array or any expression that evaluates to
    a 1-D boolean: e.g. ii = find(x < 3), where x is a 1-D ndarray.
    This syntax is similar to MATLAB usage.

    2-D or higher arrays could be flattened prior to calling find() and
    then reconstituted with reshape.  This modification could be added to
    this function as well to support N-D arrays.

    Returns 1-D ndarray of int64.
    """
    from numpy import sum, argsort, sort, ndarray

    if type(b) != ndarray:
        raise ValueError ('find: Input should be ndarray')
    if b.dtype != 'bool':
        raise ValueError ('find: Input should be boolean array')
    if b.ndim > 1:
        raise ValueError ('find: Input should be 1-D')

    F = b.size - sum(b)    # number of False in b
    idx = argsort(b)[F:]   # argsort puts True at the end, so select [F:]
    idx = sort(idx)        # be sure values in idx are ordered low to high
    return idx


def md2jd(YYYY, MM, DD):
    """
    Converts Month and Day into day-of-year (julian date)

    Usage: yday = md2jd(YYYY, MM, DD)

    Inputs YYYY, MM & DD may be strings or numbers.

    Returns integer day of year

    """
    import numpy as np

    day_tab = [[31,28,31,30,31,30,31,31,30,31,30,31],
               [31,29,31,30,31,30,31,31,30,31,30,31]]
    days = np.array(day_tab)
    yr = np.int(YYYY)
    mo = np.int(MM)
    day = np.int(DD)

    leap = np.logical_and((np.remainder(yr,4) == 0),
                          (np.remainder(yr,100) != 0))
    leap = np.logical_or(leap, (np.remainder(yr,400) == 0))

    yday = day
    for ii in np.arange(1,mo):
        yday += days[leap,ii-1]

    return np.int(yday)


def jd2md(YYYY, DOY):
    """
    Converts day-of-year (julian date) to Month and Day of given Year

    Usage: mo, da = jd2md(YYYY, JD)

    Inputs YYYY, & DOY may be strings or numbers.

    Returns tuple of integers, month and day

    """
    import numpy as np

    day_tab = [[31,28,31,30,31,30,31,31,30,31,30,31],
               [31,29,31,30,31,30,31,31,30,31,30,31]]
    days = np.array(day_tab)
    yr = np.int(YYYY)
    jd = np.int(DOY)

    leap = np.logical_and((np.remainder(yr,4) == 0),
                          (np.remainder(yr,100) != 0))
    leap = np.logical_or(leap, (np.remainder(yr,400) == 0))

    ii = 1
    while jd > days[leap,ii-1]:
        jd -= days[leap,ii-1]
        ii += 1

    mo = np.int(ii)
    day = np.int(jd)
    return (mo, day)


def replace_nan(x):
    """
    Replaces NaNs in 1D array with nearest finite value.

    Usage: y = replace_nan(x)

    Returns filled array y without altering input array x.
    Assumes input is numpy array.

    3/2015 BWB
    """
    import numpy as np
#
    x2 = np.zeros(len(x))
    np.copyto(x2,x)
#
    bads = find(np.isnan(x)) # indices of NaNs
    if bads.size == 0:
        return x2
    else:
        fins = find(np.isfinite(x))        # indices for all finites
        for ii in np.arange(0,bads.size):  # for all NaNs
            # locate index of nearest finite
            diffs = np.abs(fins-bads[ii])
            idx = diffs.argmin()
            # replace NaN with nearest finite
            x2[bads[ii]] = x[fins[idx]]
    return x2


def invert_dict(d):
    """Inverts a dictionary, returning a map from val to a list of keys.

    If the mapping key->val appears in d, then in the new dictionary
    val maps to a list that includes key.

    d: dict

    Returns: dict

    Copyright 2012 Allen B. Downey
    License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
    """
    inverse = {}
    for key, val in d.iteritems():
        inverse.setdefault(val, []).append(key)
    return inverse


def find_module(module):
    """
    Searches the PYTHONPATH for a module
    """
    import sys
    import os
    import glob

    result = []
    # Loop over the list of paths in sys.path
    for subdir in sys.path:
        # Join the subdir path with the module we're searching for
        pth = os.path.join(subdir, module)
        # Use glob to test if the pth is exists
        res = glob.glob(pth)
        # glob returns a list, if it is not empty, the pth exists
        if len(res) > 0:
            result.append(res)
    return result


def print_attributes(obj):
    """
    Prints attributes for given object
    """
    for attr in obj.__dict__:
        print(attr), getattr(obj,attr)


def walk(dirname):
    """
    Recursively traverse all files in given directory, printing
    file names.
    """
    import os

    for name in os.listdir(dirname):
        path = os.path.join(dirname,name)

        if os.path.isfile(path):
            print(path)
        else:
            walk(path)

# Test code
if __name__ == '__main__':
    print('test')
#     print 'testing invert_dict'
#     d = dict(a=1, b=2, c=3, z=1)
#     print d
#     inverse = invert_dict(d)
#     for val, keys in inverse.iteritems():
#         print val, keys
#
#     print "testing find_module('site.py')"
#     result = find_module('site.py')
#     print result
#
#     print 'testing walk(os.getcwd())'
#     walk(os.getcwd())

