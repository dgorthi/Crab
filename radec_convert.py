import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time

time_init = '2017-09-27 22:39:00'    # TIME OF OBSERVATION
rasc = 282                           # RA IN DEGREES
decl = -1                            # DEC IN DEGREES
pdt_summer_istrue = 1                # IF 1, DAYLIGHT SAVINGS TIME IN EFFECT

obj = SkyCoord(ra=rasc*u.degree, dec=decl*u.degree, frame='icrs')

leuschner = EarthLocation(lat = 37.91934*u.degree, lon = -122.15385*u.degree, height = 304*u.m)

utcoffset = -(8-pdt_summer_istrue)*u.hour
time = Time(time_init) - utcoffset

objaltaz = obj.transform_to(AltAz(obstime=time,location=leuschner))
print(objaltaz)
