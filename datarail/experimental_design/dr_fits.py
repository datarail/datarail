import numpy as np
from scipy.optimize import curve_fit


def biphasic_fit_function(x, a, b, c, d, e, f):
    """ Function for biphasic fit
    Parameters
    ----------
    x : 1d array
      array/list of doses
    a : float
      Gr inf for term 1
    b : float
      GE50 for term 1
    c : float
      HS for term 1
    d : float
      Gr inf for term 2
    e : float
      GE50 for term 2
    f : float
      HS for term 2 
   
    Returns
    -------
    biphasic_function

    """
    term1 = 1 + a + ((1 - a)/(1 + (x * 10 ** b) ** c))
    term2 = 1 + d + ((1 - d)/(1 + (x * 10 ** e) ** f))

    biphasic_function = (2 ** 0.5 * (np.log2(term1) + np.log2(term2))) - 1
    return biphasic_function


def sigmoidal_fit_function(x, a, b, c):
    sigmoidal_function = a + ((1 - a)/(1 + (x * 10 ** b) ** c ))
    return sigmoidal_function


def fit(xdata, ydata):
    """ Scipy's curve fit uses non-linear least square to fit
    function "biphasic_fit_function" to data

    Parameters
    ----------
    xdata : 1d array
       list/array of doses
    ydata : 1daray
       list/array of GR values

    Returns
    -------
    yfit: 1d array
       array of GR values estimated based on fit
    """
    if type(xdata) == list:
        xdata = np.array(xdata)
    if type(ydata) == list:
        ydata = np.array(ydata)
    ge50_low = np.max((np.min(xdata) * 1e-4, 1e-7))
    ge50_high = np.min((np.max(xdata) * 1e2, 1e2))
    lower_bounds = [-.05, ge50_low, .025, -1, 0.3, 0.025]
    upper_bounds = [1, 1, 5, .5, ge50_high, 10]

    priors = [.1, np.median(xdata), 2, -0.1, 1, 2]
   
    popt_bp, pcov_bp = curve_fit(biphasic_fit_function, xdata, ydata,
                                 bounds=(lower_bounds, upper_bounds),
                                 p0=priors)
    yfit_bp = biphasic_fit_function(xdata, *popt_bp)

    popt_sig1, pcov_sig1 = curve_fit(sigmoidal_fit_function, xdata, ydata,
                                     bounds=(lower_bounds[:3], upper_bounds[:3]),
                                     p0=priors[:3])
    yfit_sig1 = sigmoidal_fit_function(xdata, *popt_sig1)

    popt_sig2, pcov_sig2 = curve_fit(sigmoidal_fit_function, xdata, ydata,
                                     bounds=(lower_bounds[3:], upper_bounds[3:]),
                                     p0=priors[3:])
    yfit_sig2 = sigmoidal_fit_function(xdata, *popt_sig2)

    
    return yfit_bp, yfit_sig1, yfit_sig2
    


