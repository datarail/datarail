import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42



def biphasic_fit_function(x, a, b, c, d, e, f):
    """Function for biphasic fit

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
    term1 = 1 + (a + (1 - a)/(1 + (x * (10 ** b)) ** c))
    term2 = 1 + (d + (1 - d)/(1 + (x * (10 ** e)) ** f))

    biphasic_function = 2 ** (0.5 * (np.log2(term1) + np.log2(term2))) - 1
    return biphasic_function


def sigmoidal_fit_function(x, a, b, c):
    sigmoidal_function = a + ((1 - a)/(1 + (x * 10 ** b) ** c ))
    return sigmoidal_function


def fit(xdata, ydata, cap=1, extrapolrange=10, ax=None, fig_title=None):
    """Scipy's curve fit uses non-linear least square to fit
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

    # Cap on GR values
    # ----------------
    if cap > 0:
        ydata = np.array([np.min((yd, cap)) for yd in ydata])

    
    ge50_low = np.max((np.min(xdata) * 1e-4, 1e-7))
    ge50_high = np.min((np.max(xdata) * 1e2, 1e2))
    lower_bounds = [-.05, -np.log10(1), .025,
                    -1, -np.log10(ge50_high), 0.025]
    upper_bounds = [1, -np.log10(ge50_low), 5,
                    .5, -np.log10(0.3), 10]

    priors = [.1, -np.log10(np.median(xdata)), 2,
              -0.1, -np.log10(1), 2]

    cmin = np.log10(np.min(xdata)/extrapolrange)
    cmax = np.log10(np.max(xdata) * extrapolrange)
    xc = 10 ** (np.arange(cmin, cmax, 0.05))

    # Compute Biphasic fit
    # --------------------
    popt_bp, pcov_bp = curve_fit(biphasic_fit_function, xdata, ydata,
                                 bounds=(lower_bounds, upper_bounds),
                                 p0=priors)
    yfit_bp = biphasic_fit_function(xc, *popt_bp)
    #popt_bp[1] = 10 ** -popt_bp[1]
    #popt_bp[4] = 10 ** -popt_bp[4]
 
    # Compute Sigmoidal fit 1
    # ------------------------
    popt_sig1, pcov_sig1 = curve_fit(sigmoidal_fit_function, xdata, ydata,
                                     bounds=(lower_bounds[:3], upper_bounds[:3]),
                                     p0=priors[:3])
    sig1_rsquared = get_rsquare(sigmoidal_fit_function(xdata, *popt_sig1), ydata)
    yfit_sig1 = sigmoidal_fit_function(xc, *popt_sig1)
    popt_sig1[1] = 10 ** -popt_sig1[1]

    # Compute Sigmoidal fit 2
    # ------------------------
    popt_sig2, pcov_sig2 = curve_fit(sigmoidal_fit_function, xdata, ydata,
                                     bounds=(lower_bounds[3:], upper_bounds[3:]),
                                     p0=priors[3:])
    sig2_rsquared = get_rsquare(sigmoidal_fit_function(xdata, *popt_sig2), ydata)
    yfit_sig2 = sigmoidal_fit_function(xc, *popt_sig2)
    popt_sig2[1] = 10 ** -popt_sig2[1]
    
    if sig1_rsquared > sig2_rsquared:
        print('1st phase sigmoidal fit is the better of the 2 sigmoidal fits ')
        best_sig_fit = yfit_sig1
        sigmoidal_params = np.array(list(popt_sig1)+[1, -np.inf, .01])
    else:
        best_sig_fit = yfit_sig2
        print('2nd phase sigmoidal fit is the better of the 2 sigmoidal fits')
        sigmoidal_params = np.array([1, -np.inf, .01] + list(popt_sig2))

    # Plot data, biphasic and best sigmoidal fits
    # -------------------------------------------
    if ax is not None:
        ax.semilogx(xdata, ydata, 'ob', label='Measured GR value')    
        ax.semilogx(xc, yfit_bp, 'lightblue', label='Biphasic fit')
        ax.semilogx(xc, best_sig_fit, '-k', label='Best sigmoidal fit')
        ax.set_ylim((-0.5, 1))
        xlim = (10 ** cmin, 10 ** cmax)
        ax.set_xlim(xlim)
        ax.plot(xlim, [0, 0], '--k')
        ax.set_title(fig_title)
        
    return yfit_bp, popt_bp, best_sig_fit, sigmoidal_params
    


def get_sse(ypred, ytrue):
    sse = np.sum([(yt - yp) ** 2 for yp, yt in zip(ypred, ytrue)])
    return sse


def get_rsquare(ypred, ytrue):
    sst = np.sum([(yt - np.mean(ytrue))**2  for yt in ytrue])
    sse = get_sse(ypred, ytrue)
    rsquare = 1 - (sse/sst)
    return rsquare




# area = np.sum((1 - (ydata[1:] + ydata[:-1])/2) * (np.log10(xdata[1]) - np.log10(xdata[0])))/\
#     (np.log10(xdata[-1]) - np.log10(xdata[0]))
    

# grmax = np.min(ydata[-2:])
