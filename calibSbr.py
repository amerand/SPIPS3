import numpy as np
import matplotlib.pyplot as plt
import dpfit

Sv_0 = 5.5

def model(data, param):
    if isinstance(data, list):
        res = []
        for x in data:
            res.append(model(x, param))
        return res
    # assumes data is a tuple!
    if 'Z_'+data[0] in param and 'S_'+data[0] in param:
        return (param['T_'+str(data[1])]-Sv_0)*param['S_'+data[0]]+ param['Z_'+data[0]]
    elif 'S_'+data[0] in param and 'Sv_0' in param:
        return param['Sv_0'] + (param['T_'+str(data[1])]-param['Sv_0'])*param['S_'+data[0]]
    else:
        return param['T_'+str(data[1])]

def model_2(data, p):
    """
    data = [('D', hip),
            ('V', hip), ('Ic', hip), ('H', hip), ('K', hip), etc.]
    """
    if isinstance(data, list):
        res = []
        for x in data:
            res.append(model_2(x, p))
        return res
    # -- assumes data is a tuple!
    D = p['D_'+str(data[1])]
    if data[0] == 'D':
        return D
    V = p['V_'+str(data[1])]
    if data[0] == 'V':
        return V

    # -- compute logTheta Mag0
    Sv = V + 5*np.log10(D)
    if 'Z_'+data[0] in p:
        S = p['Z_'+data[0]] + p['S_'+data[0]]*(Sv - Sv_0)
        return S - 5*np.log10(D)
    elif 'Sv_0' in p:
        return V - p['S_'+data[0]]*(Sv - p['Sv_0'])

    # -- compute mag from surface brightness in V
    Sv = V + 5*np.log10(D)
    if 'Z_'+data[0] in p:
        S = p['Z_'+data[0]] + p['S_'+data[0]]*(Sv - Sv_0)
        return S - 5*np.log10(D)
    elif 'Sv_0' in p:
        return V - p['S_'+data[0]]*(Sv - p['Sv_0'])


# -- dwarfs and giants:
hip = [3765, 3821, 4151, 4436, 5336, 7513, 7981, 8102, 12114, 12777, 16537,
       16852, 19849, 22449, 24813, 27435, 27913, 32349, 32362, 35350, 36366,
       37279, 40843, 43587, 45343, 46733, 46853, 47080, 51459, 53910, 56997,
       57757, 57939, 64394, 64924, 65721, 66249, 67927, 71284, 72567, 72659,
       78459, 81300, 91262, 92043, 93747, 96100, 96441, 96895, 98505, 102422,
       108870, 112447, 113368, 114570, 114622, 116771, 120005]

tmp = [3092,7607,7884,9884,13328,20205,20455,20885,20889,21421,22453,37826,42527,
       45860,46390,49637,53229,54539,55219,56343,57399,57477,59746,60202,63608,
       67459,68594,69673,72607,74666,74793,75260,75458,77070,79882,80331,80816,
       81833,82611,86182,87833,90344,93194,94376,96837,97938,98337,99663,102488,
       104732,110538,111944]
hip.extend(tmp)

# spType = [
# K3V,F9V,F9V,A6V,K1V(Fe-2),F9V,K1V,G8.5V,K3V,F7V,K2V(k),F9IV-V,K0.5V,F6IV-V,G1V,
# G2V,G0V(CH-0.3),A0V(mA1Va),F5IV-V,A3V,F1V,F5IV-V,F6V,K0IV-V,M0.0V,F0V,F7V,G8IV,
# F8V,A1IV(spSr),G8V,F8.5IV-V,G8V(P),G0V,G7V,G5V,A2V(an),G0IV,F4V(kF2mF1),F9IV-V,
# G7V,G0V,K0V(k),A1V,F5IV-V,A1V,G9V,F3V,G2V,K2V,K0IV,K5V,F6V,A4V,F1V,K3V,F7V,K7V
# ]

diam= [(0.868, 0.004), (1.623, 0.004), (0.865, 0.010), (0.708, 0.013),
        (0.972, 0.009), (1.143, 0.010), (1.000, 0.004), (2.080, 0.030),
        (1.030, 0.007), (1.103, 0.009), (2.126, 0.014), (1.081, 0.014),
        (1.446, 0.022), (1.419, 0.027), (0.981, 0.015), (0.572, 0.009),
        (1.051, 0.009), (5.959, 0.059), (1.401, 0.009), (0.835, 0.013),
        (0.853, 0.014), (5.434, 0.050), (0.706, 0.013), (0.711, 0.004),
        (0.871, 0.015), (1.133, 0.009), (1.632, 0.005), (0.821, 0.013),
        (0.794, 0.014), (1.149, 0.014), (0.910, 0.009), (1.431, 0.006),
        (0.686, 0.006), (1.127, 0.011), (1.073, 0.005), (1.010, 0.020),
        (0.852, 0.009), (2.252, 0.036), (0.841, 0.013), (0.569, 0.011),
        (1.196, 0.014), (0.735, 0.014), (0.724, 0.011), (3.280, 0.010),
        (1.000, 0.006), (0.895, 0.017), (1.254, 0.012), (0.844, 0.009),
        (0.554, 0.011), (0.385, 0.006), (2.650, 0.040), (1.881, 0.017),
        (1.091, 0.008), (2.230, 0.020), (0.648, 0.008), (1.106, 0.007),
        (1.082, 0.009), (0.856, 0.016)]

tmp = [(4.168, 0.047), (3.760, 0.070), (2.810, 0.030), (6.847, 0.071), (4.060, 0.040),
        (2.520, 0.030), (2.302, 0.040), (2.310, 0.040), (2.572, 0.046), (20.297, 0.384),
        (2.727, 0.013), (8.177, 0.130), (2.225, 0.020), (8.025, 0.142), (9.700, 0.100),
        (3.330, 0.040), (2.540, 0.030), (4.107, 0.053), (4.745, 0.060), (2.386, 0.021),
        (3.230, 0.020), (1.606, 0.006), (1.498, 0.028), (1.651, 0.016), (3.254, 0.037),
        (4.720, 0.050), (0.948, 0.012), (20.877, 0.277), (10.300, 0.100), (2.744, 0.036),
        (2.336, 0.020), (1.690, 0.031), (3.596, 0.015), (4.828, 0.062), (2.961, 0.007),
        (3.633, 0.066), (3.492, 0.050), (2.529, 0.050), (1.440, 0.004), (1.515, 0.010),
        (9.978, 0.180), (2.120, 0.020), (0.753, 0.009), (3.268, 0.054), (1.765, 0.012),
        (1.726, 0.008), (6.821, 0.098), (1.859, 0.003), (4.610, 0.050), (2.820, 0.030),
        (1.920, 0.020), (2.731, 0.024)]
diam.extend(tmp)
diam, ediam = np.array([d[0] for d in diam]), np.array([d[1] for d in diam])

V = [5.74, 3.460, 4.800, 3.860, 5.170, 4.100, 5.240, 3.490, 5.790, 4.100, 3.720,
     4.290, 4.430, 3.190, 4.690, 5.970, 4.390, -1.440, 3.350, 3.580, 4.160, 0.400,
     5.130, 5.960, 7.640, 3.650, 3.170, 5.400, 4.820, 2.340, 5.310, 3.590, 6.420,
     4.230, 4.740, 4.970, 3.380, 2.680, 4.470, 5.860, 4.540, 5.390, 5.770, 0.030,
     4.190, 2.990, 4.670, 4.490, 5.990, 7.670, 3.410, 4.690, 4.200, 1.170, 4.530,
     5.570, 4.130, 7.700]
tmp = [3.270, 3.590, 4.450, 2.010, 4.560, 3.650, 3.770, 3.840, 3.530, 0.870, 4.890,
       1.160, 4.590, 3.140, 1.990, 4.390, 3.790, 3.000, 3.490, 3.540, 3.690, 5.270,
       5.720, 4.720, 2.850, 4.050, 6.180, -0.050, 2.070, 3.460, 5.020, 5.720, 3.290,
       2.630, 3.230, 2.730, 2.780, 3.480, 5.990, 5.350, 2.240, 4.820, 3.250, 3.070,
       4.390, 4.710, 3.510, 5.810, 2.480, 3.210, 4.420, 4.500]
V.extend(tmp)
V = np.array(V)
eV = np.ones(len(V))*0.02

Ic = [4.780, np.nan, 4.210, 3.700, 4.360, 3.500, 4.360, 2.630, 4.740, 3.530, np.nan,
    3.640, 3.530, 2.650, 4.040, np.nan, np.nan, -1.430, 2.870, 3.450, 3.780, -0.140,
    np.nan, np.nan, np.nan, 3.270, 2.610, np.nan, 4.240, 2.380, 4.580, 3.000, 5.570,
    3.620, 3.990, 4.190, 3.280, 2.080, 4.020, np.nan, np.nan, np.nan, np.nan, 0.080,
    3.660, 2.990, 3.850, 4.020, 5.440, 6.680, 2.510, 3.530, 3.590, 1.090, 4.160,
    4.470, 3.520,  np.nan]
tmp = [2.040, 2.310, 3.050, 0.860, 2.820, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
      0.160, 3.420, 1.460, 0.550, 2.890, 2.770, 1.920, 2.110, 2.620, 2.580, np.nan, np.nan,
      3.720, 1.970, 2.440, np.nan, -1.330, 0.590, 2.480, np.nan, np.nan, 2.200,
      1.560, 2.290, 1.890, 1.880, 2.580, np.nan, np.nan, 0.630, 3.630, 3.260, 2.120, 3.410,
      3.630, 1.790, np.nan, 1.450, 2.270, 3.400, 3.190]
Ic.extend(tmp)
Ic = np.array(Ic)
eIc = np.ones(len(Ic))*0.027
eIc[Ic==6.680] = 0.008


H = [np.nan, 2.020, 3.560, 3.370, np.nan, 2.990, 3.345, 1.727, 3.542, 3.070, 1.749,
    np.nan, np.nan, 2.148, 3.330, 4.499, 3.050, -1.387, np.nan,  np.nan, np.nan,
    -0.569, 3.940, 4.140, 4.253, np.nan, 2.025, 3.770, np.nan, np.nan, np.nan,
    2.345, np.nan, 2.923, np.nan, 3.320, 3.050, 1.390, 3.516, 4.530, 3.000, 3.945,
    3.910, 0.004, np.nan, np.nan, np.nan, np.nan, 4.731, np.nan, np.nan, np.nan,
    np.nan, 1.054, np.nan, 3.400, np.nan, 4.253]
tmp = [0.551, np.nan, 1.409, -0.558, np.nan, 1.500, np.nan, np.nan, np.nan,
      -2.653, np.nan, -1.003, 1.941, -0.475, -1.074, 1.190, np.nan, 0.539,
      0.415, 1.577, 1.020, np.nan, np.nan, np.nan, 0.770, np.nan, 3.775,
      -2.951, np.nan, 1.260, np.nan, np.nan, np.nan, 0.197, np.nan, np.nan, 0.690,
      np.nan, np.nan, np.nan, -1.160, np.nan, 3.195, np.nan, 2.210, np.nan, -0.042, np.nan,
      0.206, 1.155, 2.217, 1.600 ]
H.extend(tmp)
H = np.array(H)
eH = np.ones(len(H))*0.050
eH[H==3.370] = 0.060
eH[H==1.600] = 0.007

K = [np.nan, 1.821, np.nan, 3.365, np.nan, 2.841, np.nan, 1.631, np.nan, 2.761,
    1.601, 2.871, np.nan, 2.031, 3.255, np.nan, 2.971, np.nan, 2.111, np.nan,
    np.nan, -0.669, np.nan, np.nan, np.nan, 2.711, 1.951, np.nan, np.nan, 2.361,
    np.nan, 2.301, np.nan, 2.851, np.nan, np.nan, np.nan, 1.291, np.nan, np.nan,
    2.651, 3.901, np.nan, -0.079, 2.941, 2.921, 2.811, np.nan, 4.569, np.nan,
    1.201, np.nan, 2.851, 0.981, np.nan, np.nan, 2.731, np.nan]
tmp = [0.421, 0.771, 1.241, -0.649, 0.721, 1.481, 1.581, 1.621, 1.291, np.nan,
    1.441, -1.139, 1.901, -0.699, -1.379, 1.011, 1.351, 0.371, 0.251, 1.491,
    0.901, 2.531, 2.921, 2.321, 0.731, 0.391, 3.666, np.nan, -1.259, 1.121,
    1.901, 2.721, 0.701, 0.041, 0.961, 0.601, 0.621, 1.281, 2.811, 2.651, -1.319,
    1.931, np.nan, 0.741, 2.071, 2.351, -0.309, 2.311, 0.101, 1.051, 1.961, 1.391]
K.extend(tmp)
K = np.array(K)
eK =  [np.nan, 0.060 , np.nan, 0.071 , np.nan, 0.080 , np.nan, 0.060 , np.nan,
       0.090 , 0.060 , 0.100 , np.nan, 0.060 , 0.045 , np.nan, 0.070 , np.nan,
       0.060, np.nan, np.nan, 0.051 , np.nan, np.nan, np.nan, 0.090 , 0.070,
       np.nan, np.nan,0.060 , np.nan, 0.060 , np.nan, 0.100 , np.nan, np.nan, np.nan,
       0.051 , np.nan, np.nan, 0.080 , 0.045, np.nan, 0.060 , 0.090 , 0.080 ,
       0.080 , np.nan, 0.045 , np.nan, 0.051 , np.nan, 0.080 , 0.051 , np.nan,
       np.nan, 0.080, np.nan]
eK.extend([0.051 ,0.041 ,0.031 ,0.051 ,0.051 ,0.041 ,0.051 ,0.060 ,0.051 ,
           np.nan,0.041 ,0.051 ,0.070 ,0.031 ,0.060 ,0.070 ,0.041 ,0.041 ,
           0.041 ,0.060 ,0.031 ,0.060 ,0.080 ,0.051 ,0.051 ,0.041 ,0.015 ,
           np.nan,0.070 ,0.031 ,0.051 ,0.060 ,0.041 ,0.051 ,0.051 ,0.031 ,
           0.041 ,0.031 ,0.090 ,0.070 ,0.041 ,0.051 ,np.nan,0.051 ,0.060 ,
           0.090 ,0.041 ,0.070 ,0.070 ,0.031 ,0.051 ,0.070])
eK = np.array(eK)
# ---------------------------------------------------


def fitSb():
    data = []
    logD = np.log10(diam)
    elogD = logD - np.log10(diam-ediam)
    for i in range(len(hip)):

        Sv = V[i]  + 5*logD[i]
        eSv = np.sqrt(eV[i]**2 + 25*elogD[i]**2)
        Si = Ic[i] + 5*logD[i]
        eSi = np.sqrt(eIc[i]**2 + 25*elogD[i]**2)
        Sh = H[i]  + 5*logD[i]
        eSh = np.sqrt(eH[i]**2 + 25*elogD[i]**2)
        Sk = K[i]  + 5*logD[i]
        eSk = np.sqrt(eK[i]**2 + 25*elogD[i]**2)

        data.append(('V', hip[i], Sv, eSv))
        if not np.isnan(Si):
            data.append(('Ic', hip[i], Si, eSi))
        if not np.isnan(Sh):
            data.append(('H', hip[i], Sh, eSh))
        if not np.isnan(Sk):
            data.append(('K', hip[i], Sk, eSk))

    param = {'T_'+str(hip[i]):V[i]+5*logD[i] for i in range(len(hip))}

    # -- linear
    if True:
        param.update({'Z_'+b:0.0 for b in ['Ic', 'H', 'K']})
        param.update({'S_'+b:1.0 for b in ['Ic', 'H', 'K']})
        fit = dpfit.leastsqFit(model, data, param, [d[2] for d in data],
                                [d[3] for d in data],
                                doNotFit=['S_V','Z_V'], # T==Sv basically
                                verbose=0, maxfev=1000)
        show = [x for x in list(param.keys()) if not x.startswith('T_') and not '_V' in x]
        show.sort()
    else:
        param.update({'S_'+b:1.0 for b in ['Ic', 'H', 'K']})
        param['Sv_0'] = 2.6
        fit = dpfit.leastsqFit(model, data, param, [d[2] for d in data],
                                [d[3] for d in data],
                                doNotFit=['S_V', 'Sv_0'], # T==Sv basically
                                verbose=0, maxfev=1000)
        show = ['S_Ic', 'S_H', 'S_K']
        #fit = dpfit.randomParam(fit)

    print('chi2=%.2f'%fit['chi2'])
    dpfit.dispCor(fit, params=show)

    plt.figure(0,figsize=(9,7))
    plt.clf()
    plt.subplots_adjust(hspace=0.3, wspace=0.25, top=0.93, right=0.98, left=0.08,
                        bottom=0.08)
    plt.suptitle('Adams et al. 2017')

    ax1 = plt.subplot(221)
    ax1.set_ylabel('Sb')
    ax2 = plt.subplot(223, sharex=ax1)

    colors = {'V':'g', 'Ic':'y', 'H':'orange', 'K':'r'}
    for b in ['V','Ic','H','K']:
        _data = [x for x in data if x[0]==b]
        X = [fit['best']['T_'+str(d[1])] for d in _data]
        eX = [fit['uncer']['T_'+str(d[1])] for d in _data]
        Y = np.array([d[2] for d in _data])
        eY = np.array([d[3] for d in _data])
        M = np.array(model(_data, fit['best']))
        R = Y-M
        eR = eY
        ax1.errorbar(X,Y, marker=',', color=colors[b],
                    xerr = eX, yerr=eY, linestyle='',
                    alpha=0.5, label='data: '+b)
        ax1.plot(X, M, 'o', color=colors[b], markersize=3, alpha=0.5,
                label='model: '+b)
        ax2.errorbar(X, R, marker='.', color=colors[b], alpha=0.5,
                    linestyle='', xerr=eX, yerr=eR,
                    label=b+' $\sigma$=%.3f'%np.std(R))
    ax2.set_xlabel('fitted Sb(V)')
    ax1.set_xlabel('fitted Sb(V)')

    ax2.set_ylabel('residuals')
    ax1.legend(fontsize=6)
    ax1.set_ylim(2, 8)
    ax2.legend(fontsize=6)
    ax2.set_ylim(-0.35, 0.35)
    ax2.set_xlim(2, 8)

    # -- data:
    Sv_d = V + 5*logD
    eSv_d = np.sqrt(eV**2 + 25*elogD**2)
    # -- fitted:
    Sv_f = np.array([fit['best']['T_'+str(h)] for h in hip])

    S = {'Ic':Ic + 5*logD, 'H': H + 5*logD,
         'K': K + 5*logD }
    eS = {'Ic':np.sqrt(eIc**2 + 25*elogD**2),
          'H': np.sqrt(eH**2 + 25*elogD**2),
          'K': np.sqrt(eK**2 + 25*elogD**2)}
    eM = {'Ic':eIc, 'H': eH, 'K': eK}

    _Sv = np.linspace(2.3,7.9,10)
    print('Sb(%2s)_data = %.3f + %.3f[Sb(V)_fitted-%.1f]'%('V', Sv_0, 1.0, Sv_0), end=' ')
    print('[FIXED]')

    for i,b in enumerate(['Ic', 'H', 'K']):
        print('Sb(%2s) = %.3f + %.3f[Sb(V)-%.1f]'%(b, fit['best']['Z_'+b],
                                                fit['best']['S_'+b], Sv_0))
        print('   +-    %.3f   %.3f'%(fit['uncer']['Z_'+b], fit['uncer']['S_'+b]))
        plt.subplot(3,2,2+2*i)
        plt.ylim(2,8)
        w = ~np.isnan(S[b])
        # -- data points
        plt.errorbar(Sv_d[w]-S[b][w], Sv_d[w], marker=',', label='data',
                    alpha=0.5, linestyle='', yerr=eSv_d[w],
                    xerr = np.sqrt(eV[w]**2+eM[b][w]**2))
        # -- extrapolated
        S2 = fit['best']['Z_'+b] + (Sv_d-Sv_0)*fit['best']['S_'+b]
        plt.plot(Sv_f[w]-S2[w], Sv_f[w],  '.k', label='fitted', alpha=0.2)
        plt.xlabel('V - '+b)
        plt.ylabel('Sb(V)')
        plt.legend()

        tmp = []
        # for p in fit['r_param']:
        #     tmp.append(p['Z_'+b] + (_Sv-Sv_0)*p['S_'+b])
        # _S2 = np.mean(tmp, axis=0)
        # plt.plot(_Sv-_S2, _Sv, '-b')
        # -- uncertainty is too small
        #plt.plot(_Sv-_S2-np.std(tmp, axis=0), _Sv, '--b')
        #plt.plot(_Sv-_S2+np.std(tmp, axis=0), _Sv, '--b')

def fitSb_2():
    data = []
    param = {}
    for i in range(len(hip)):
        data.append(('D', hip[i], diam[i], ediam[i]))
        data.append(('V', hip[i], V[i], eV[i]))
        param['D_'+str(hip[i])] = round(diam[i],2)
        param['V_'+str(hip[i])] = round(V[i],2)
        if not np.isnan(Ic[i]):
            data.append(('Ic', hip[i], Ic[i], eIc[i]))
        if not np.isnan(H[i]):
            data.append(('H', hip[i], H[i], eH[i]))
        if not np.isnan(K[i]):
            data.append(('K', hip[i], K[i], eK[i]))
    #multiZero = False
    multiZero = True

    if multiZero:
        param.update({'Z_Ic':4.56, 'S_Ic':0.67,
                      'Z_H':3.45, 'S_H':0.28,
                      'Z_K':3.32, 'S_K':0.25})
        show = ['S_Ic', 'Z_Ic', 'S_H', 'Z_H', 'S_K', 'Z_K']
    else:
        param.update({'S_Ic':0.67,'S_H':0.28, 'S_K':0.25, 'Sv_0':2.6})
        show = ['S_Ic', 'S_H', 'S_K', 'Sv_0']

    print(len(data), len(param))
    doNotFit = [x for x in list(param.keys()) if x.startswith('D_')]
    fit = dpfit.leastsqFit(model_2, data, param,
                            [d[2] for d in data],
                            [d[3] for d in data],
                            #doNotFit=doNotFit,
                            verbose=0, maxfev=1000)
    print('chi2=%.2f'%fit['chi2'])
    dpfit.dispCor(fit, params=show)
    plt.figure(0)
    plt.clf()
    plt.subplots_adjust(left=0.1)
    # -- _d for data, _f for fitted/modeled
    D_d = np.array([d[2] for d in data if d[0]=='D'])
    eD_d = np.array([d[3] for d in data if d[0]=='D'])
    D_m = np.array([fit['best']['D_'+str(d[1])] for d in data if d[0]=='D'])
    eD_m = np.array([fit['uncer']['D_'+str(d[1])] for d in data if d[0]=='D'])

    V_d = np.array([d[2] for d in data if d[0]=='V'])
    eV_d = np.array([d[3] for d in data if d[0]=='V'])
    V_m = np.array([fit['best']['V_'+str(d[1])] for d in data if d[0]=='V'])
    eV_m = np.array([fit['uncer']['V_'+str(d[1])] for d in data if d[0]=='V'])

    plt.subplot(221)
    plt.errorbar(0.2*V_d+np.log10(D_d), V_d - V_m,
                #xerr=eV_d,
                yerr=np.sqrt(eV_d**2+0*eV_m**2),
                marker='o', linestyle='')
    plt.xlabel(r'log$\theta_{V=0}$ [data] ')
    plt.ylabel(r'Vmag data - model')
    plt.grid()
    plt.ylim(-0.06,0.06)

    rD = D_d/D_m
    erD = np.sqrt(eD_d**2/D_m**2+ eD_m**2*D_d**2/D_m**4)
    print(np.array(hip)[(rD-1)>1.5*erD])
    plt.subplot(223)
    plt.errorbar(0.2*V_d+np.log10(D_d), 100*(rD-1),
                #xerr=eV_d,
                yerr=100*erD, marker='o', linestyle='')
    plt.xlabel(r'log$\theta_{V=0}$ [data] ')
    plt.ylabel(r'$\theta$ data/model (%)')
    plt.grid()
    plt.ylim(-11,11)

    plt.subplot(1,2,2)
    for i, b in enumerate(['Ic', 'H', 'K']):
        if multiZero:
            print('Sb(%2s) = %.3f + %.3f*[Sb(V) - %.1f]'%(b, fit['best']['Z_'+b],
                                                    fit['best']['S_'+b], Sv_0))
            print('   +-    %.3f   %.3f'%(fit['uncer']['Z_'+b], fit['uncer']['S_'+b]))
        else:
            print('%2s = V - %.3f*[Sb(V) - %.3f]'%(b, fit['best']['S_'+b], fit['best']['Sv_0']))
        D_d, eD_d, V_d, eV_d = [], [], [], []
        hb = [d[1] for d in [x for x in data if x[0]==b]]
        for h in hb:
            for d in data:
                if d[1]==h and d[0]=='V':
                    V_d.append(d[2])
                    eV_d.append(d[3])
                if d[1]==h and d[0]=='D':
                    D_d.append(d[2])
                    eD_d.append(d[3])

        V_d = np.array(V_d)
        eV_d = np.array(eV_d)
        D_d = np.array(D_d)
        eD_d = np.array(eD_d)

        V_m = np.array([fit['best']['V_'+str(d[1])] for d in data if d[0]==b])
        eV_m = np.array([fit['uncer']['V_'+str(d[1])] for d in data if d[0]==b])
        D_m = np.array([fit['best']['D_'+str(d[1])] for d in data if d[0]==b])
        eD_m = np.array([fit['uncer']['D_'+str(d[1])] for d in data if d[0]==b])

        M_d = np.array([d[2] for d in data if d[0]==b])
        eM_d = np.array([d[3] for d in data if d[0]==b])
        M_m = np.array([fit['model'][j] for j in range(len(data)) if data[j][0]==b])

        Sv_d = V_d + 5*np.log10(D_d)
        Sv_m = V_m + 5*np.log10(D_m)
        elD_d = np.log10(D_d) - np.log10(D_d-eD_d)
        eSv_d = np.sqrt(eV_d**2 + 25*elD_d**2)
        elD_m = np.log10(D_m) - np.log10(D_m-eD_m)
        eSv_m = np.sqrt(eV_m**2 + 25*elD_m**2)

        S_d = M_d + 5*np.log10(D_d)
        S_m = M_m + 5*np.log10(D_m)
        eS_d = np.sqrt(eM_d**2 + 25*elD_d**2)
        eS_m = np.sqrt(25*elD_m**2)

        lD0_d = np.log10(D_d) + 0.2*M_d
        elD0_d = np.sqrt(elD_d**2 + 0.04*eM_d**2)
        lD0_m = np.log10(D_m) + 0.2*M_m
        elD0_m = np.sqrt(elD_m**2 )

        # plt.subplot(3,3,2+3*i)
        # plt.errorbar(Sv_d - S_d, Sv_d, xerr=eSv_d, yerr=eS_d, linestyle='',
        #             label='data', color='k')
        # plt.errorbar(Sv_m - S_m, Sv_m, xerr=eSv_m, yerr=eS_m, linestyle='',
        #             label='fitted', color='y')
        # for j in range(len(Sv_d)):
        #     plt.plot([(Sv_d - S_d)[j], (Sv_m - S_m)[j]],
        #             [Sv_d[j], Sv_m[j]],
        #             color='k', alpha=0.2)
        # plt.legend()
        # plt.xlabel('V - '+b)
        # plt.ylabel('SbV')
        # plt.ylim(2,8)

        #plt.subplot(3,3,3+3*i)
        #plt.subplot(3,2,2+2*i)

        plt.errorbar(Sv_d - S_d, lD0_d, marker=',', linestyle='',
            yerr=elD0_d, xerr=np.sqrt(eV_d**2+eM_d**2),
            label='Mag = '+b)
        plt.errorbar(Sv_m - S_m, lD0_m, marker=',', color='k', linestyle='',
            yerr=elD0_m, xerr=np.sqrt(eV_m**2))
        for j in range(len(Sv_d)):
            plt.plot([Sv_d[j]-S_d[j], Sv_m[j]-S_m[j]],
                    [lD0_d[j], lD0_m[j]],
                      color='k', alpha=0.2)
    plt.legend()
    plt.grid()
    plt.xlabel('V - Mag')
    plt.ylabel(r'log$\theta_{Mag=0}$')

    return
