import astropy
from astropy.coordinates import AltAz
from astropy.coordinates import get_sun
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

observing_time = Time('2017-10-25T22:50:00', format = 'isot', scale = 'utc')#use UTC time
observing_location = EarthLocation(lat='37d55.1m', lon='-122d9.4m', height=304*u.m) #leuschner location#
AA = AltAz(location=observing_location, obstime=observing_time)
sun =  SkyCoord(get_sun(observing_time))
#Change to local Az and Alt#
LA_sun_at_obs = sun.transform_to(AA)
print LA_sun_at_obs
print '#######Postion for source in Observatory########'
print 'Azimuth:',LA_sun_at_obs.az,'Altitude:',LA_sun_at_obs.alt


time_init = '2017-10-04 10:45:00'    # TIME OF OBSERVATION
rasc = LA_sun_at_obs                           # RA IN DEGREES
decl = -1                            # DEC IN DEGREES
pdt_summer_istrue = 1                # IF 1, DAYLIGHT SAVINGS TIME IN EFFECT

obj = SkyCoord(ra=(rasc*u.degree), dec=(decl*u.degree), frame='icrs')

leuschner = EarthLocation(lat = 37.91934*u.degree, lon = -122.15385*u.degree, height = 304*u.m)

utcoffset = -(8-pdt_summer_istrue)*u.hour
time = Time(time_init) - utcoffset

objaltaz = obj.transform_to(AltAz(obstime=time,location=leuschner))
print(objaltaz)
