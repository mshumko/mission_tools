import numpy as np
import scipy.stats

def fit_exponent(xData, yData):
    """
    NAME: fit_exponent(self, xData, yData)
    USE:  Fits an exponent to xData and yData using linear regression after
          the exponent has been transformed into a linear function.
    RETURNS: A dictionary of E0 and J0 for that fit. 
    MOD: 2017-10-30
    AUTHOR: Mykhaylo Shumko
    """
    wData = np.log(yData)
    slope, intercept, r_value, p_value, std_err = \
        scipy.stats.linregress(xData, y=wData)
    fitParams = {'J0': np.exp(intercept), 'E0':(-1/slope)}
    return fitParams
