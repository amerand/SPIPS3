#!/home/amerand/anaconda2/bin/python
import os
from astropy.io import fits

# -- list directories:
dirs = os.listdir('./')
dirs = [x for x in dirs if os.path.isdir(x)]

filename = 'index.html'
f = open(filename, 'w')
f.write('<!DOCTYPE html>\n<html><body>\n')
f.write('SPIPS models:\n')
f.write('<table style="width=100%" border="1">\n')
f.write('<tr>\n')
f.write('<th>Star Name</th>')
f.write('<th>Period<BR>(d)</th>')
f.write('<th>distance<BR>(pc)</th>')
f.write('<th>P-factor</th>')
f.write('<th>[ang. diam.]<BR>(mas)</th>')
f.write('<th>[Teff]<BR>(K)</th>')
f.write('<th>E(B-V)</th>')
f.write('<th>Ex 2um<BR>(mag)</th>')
f.write('<th>Ex 6um<BR>(mag)</th>')
f.write('<th>Ex 12um<BR>(mag)</th>')
f.write('</tr>\n')
for d in dirs:
    f.write('<tr>\n')
    # -- find FITS files
    fitsfile = [x for x in os.listdir(d) if x.endswith('.fits')][0]
    h = fits.open(os.path.join(d, fitsfile))
    f.write("<td> <a href='./%s'>%s</td>\n"%(d, h[0].header['STARNAME']))
    f.write("<td> %10.7f <BR>+- %9.7f</td>"%(h[0].header['PARAM PERIOD'],
                                          h[0].header['UNCER PERIOD']))
    if h[0].header['UNCER d_kpc']>0:
        f.write("<td> %7.1f <BR>+- %6.1f</td>"%(1000*h[0].header['PARAM d_kpc'],
                                          1000*h[0].header['UNCER d_kpc']))
    else:
        f.write("<td> %7.1f </td>"%(1000*h[0].header['PARAM d_kpc']))
    if h[0].header['UNCER P-FACTOR']>0:
        f.write("<td> %6.4f <BR>+- %6.4f</td>"%(h[0].header['PARAM P-FACTOR'],
                                          h[0].header['UNCER P-FACTOR']))
    else:
        f.write("<td> %6.4f </td>"%(h[0].header['PARAM P-FACTOR']))
    f.write("<td> %6.4f <BR>+- %6.4f</td>"%(h[0].header['PARAM DIAMAVG'],
                                          h[0].header['UNCER DIAMAVG']))
    f.write("<td>%4.0f</td>"%(h[0].header['MODEL AVG_TEFF']))
    f.write("<td>%.3f</td>"%(h[0].header['PARAM E(B-V)']))
    f.write("<td>%.3f</td>"%(h[0].header['MODEL IR EXCESS  2.0UM']))
    f.write("<td>%.3f</td>"%(h[0].header['MODEL IR EXCESS  6.0UM']))
    f.write("<td>%.3f</td>"%(h[0].header['MODEL IR EXCESS 12.0UM']))
    f.write('</tr>\n')
    h.close()
f.write('<hr>\n')
f.write('</html></body>\n')
f.close()

