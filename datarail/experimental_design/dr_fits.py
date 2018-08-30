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
    trm1 = 1 + a + (1 - a)/(1 + (x * 10 ** b) ** c)
    trm2 = 1 + d + (1 - d)/(1 + (x * 10 ** e) ** f)

    biphasic_function = (2 ** 0.5 * (np.log2(trm1) + np.log2(trm2))) - 1
    return biphasic_function


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
       list/array of GR values estimated based on fit
    """
    ge50_low = np.max(np.min(xdata)*1e-4, 1e-7)
    ge50_high = np.min(np.max(xdata)*1e2, 1e2)
    lower_bounds = [-.05, ge50_low, .025, -1, 0.3, 0.025]
    upper_bounds = [1, 1, 5, -.5, ge50_high, 10]

    priors = [.1, np.median(xdata), 2, -0.1, 1, 2]
   
    popt, pcov = curve_fit(biphasic_fit_function, xdata, ydata,
                           bounds=(lower_bounds, upper_bounds),
                           p0=priors)
    yfit = biphasic_fit_function(xdata, *popt)
    return yfit
    


