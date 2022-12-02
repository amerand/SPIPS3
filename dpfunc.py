"""
collection of functions parametrized by dictionnaries, to allow easy fitting
with dpfit.py

author: amerand@eso.org
"""
import numpy as np

def polyN(x, params):
    """
    Polynomial function. e.g. params={'A0':1.0, 'A2':2.0} returns
    x->1+2*x**2. The coefficients do not have to start at 0 and do not
    need to be continuous. Actually, the way the function is written,
    it will accept any '_x' where '_' is any character and x is a float.
    """
    res = 0
    for k in list(params.keys()):
        res += params[k]*np.array(x)**float(k[1:])
    return res

def power(x,params):
    """
    return params['AMP']*(X - params['X0'])**params['POW'] + params['OFFSET']

    only params['POW'] is mandatory
    """
    res = x
    if 'X0' in params:
        res -= params['X0']
    res = res**params['POW']
    if 'AMP' in params:
        res *= params['AMP']
    if 'OFFSET' in params:
        res += params['OFFSET']
    return res

def gaussian(x,params):
    """
    1D gaussian function of moment 'MU' and variance 'SIGMA'**2

    params: {'MU':, 'SIGMA':, ('OFFSET':,) ('AMP':)} 'AMP' and
    'OFFSET' are optional: if AMP is not given, the amplitude is set
    to 1/(sigma*sqrt(2*pi)).
    """
    res = np.exp(-(x-params['MU'])**2/(2*params['SIGMA']**2))
    if 'AMP' in params:
        res *= params['AMP']
    else:
        res *= 1./(params['SIGMA']*np.sqrt(2*np.pi))

    if 'OFFSET' in params:
        res += params['OFFSET']
    return res

def lorentzian(x,params):
    """
    1D lorentzian function of moment 'MU' and parameter 'GAMMA'

    params: {'MU':, 'GAMMA':, ('OFFSET':,) ('AMP':)}. 'AMP' and
    'OFFSET' are optional: if AMP is not given, the amplitude is set to
    1/(sigma*sqrt(2*pi)). """
    res = 1/(1+(x-params['MU'])**2/(params['GAMMA']**2))
    if 'AMP' in params:
        res *= params['AMP']
    else:
        res *= 1./(params['GAMMA']*np.pi)
    if 'OFFSET' in params:
        res += params['OFFSET']
    return res

def sin(x, params):
    """
    sinusoidal wave
    params: {'AMP':amplitude, 'WAV':wavelength, ('OFFSET':offset), ('PHI':phase
    in radian)}

    returns AMP*sin(2*pi*x/WAV + PHI) + OFFSET
    """
    xx = 2*np.pi*x/params['WAV']
    if 'PHI' in params:
        xx += params['PHI']
    res = params['AMP']*np.sin(xx)
    if 'OFFSET' in params:
        res += params['OFFSET']
    return res

def cos(x, params):
    """
    cosinusoidal wave
    params: {'AMP':amplitude, 'WAV':wavelength, ('OFFSET':offset), ('PHI':phase
    in radian)}

    returns AMP*cos(2*pi*x/WAV + PHI) + OFFSET
    """
    xx = 2*np.pi*x/params['WAV']
    if 'PHI' in params:
        xx += params['PHI']
    res = params['AMP']*np.cos(xx)
    if 'OFFSET' in params:
        res += params['OFFSET']
    return res

def fourier(x, params):
    """
    fourier serie:
    (A0) + sum(Ak*cos(k*2*pi*x/WAV +PHIk))_k=1...

    A0 is optional. Parameters should contain at least WAV and one (Ak,PHIk),
    for example:

    {'WAV':1.0, 'A1':1.0, 'PHI1':0} -> cos(2*pi*x)
    {'A0':1.0, 'WAV':2*np.pi, 'A2':1.0, 'PHI2':np.pi/4} -> 1+cos(2*x+np.pi/4)

    warning, will crash if dictionnary is ill formed, for example:
    {'WAV':1.0, 'A1':1.0, 'PHI2':0}
    """
    phi = [k for k in list(params.keys()) if k[:3]=='PHI']
    res = np.copy(x)
    res *=0
    for f in phi:
        res += params['A'+f[3:]]*np.cos(float(f[3:])*2*np.pi*x/params['WAV']
                                                + params[f])
    if 'A0' in params:
        res += params['A0']
    return res

def gauss2d(xy, params):
    """
    xy = [x, y]
    """
    if 'SIGMAX' in list(params.keys()):
        sigmax = params['SIGMAX']
    if 'SIGMAY' in list(params.keys()):
        sigmay = params['SIGMAY']
    if 'SIGMA' in list(params.keys()):
        sigmax = params['SIGMA']
        sigmay = params['SIGMA']

    res = np.exp(-(xy[0]-params['X0'])**2/(2*sigmax)**2)
    res *= np.exp(-(xy[1]-params['Y0'])**2/(2*sigmay)**2)
    if 'AMP' in list(params.keys()):
        res *= params['AMP']
    if 'OFFSET' in list(params.keys()):
        res *= params['OFFSET']
    return res
