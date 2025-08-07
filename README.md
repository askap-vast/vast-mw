# vast-mw
Tools for multi-wavelength searches of VAST objects

Requirements:
* [astroquery](https://astroquery.readthedocs.io/en/latest/)
* [psrqpy](https://psrqpy.readthedocs.io/en/latest/#)
* astropy
* numpy
* loguru

Install via:
```
git clone git@github.com:askap-vast/vast-mw.git
pip install .
```
or:
`pip install git+ssh://git@github.com:askap-vast/vast-mw.git`
if you don't want to check out the repository.

---
## Tools
* [check_gaia](#check_gaia-look-for-matches-in-gaia-currently-dr3): check for matches in Gaia
* [check_simbad](#check_simbad-look-for-matches-in-simbad-currently-dr3): check for matches in Simbad
* [check_pulsarscraper](#check_pulsarscraper-search-for-pulsars-in-atnf-or-unpublished-catalogs): check for matches in pulsar survey scraper
* [check_atnf](#check_atnf-search-for-pulsars-in-atnf-catalog): check for matches in pulsar catalog
* [check_casda](#check_casda-check-for-askap-observations): check for ASKAP observations
* [check_vla](#check_vla-check-for-vla-or-evla-observations): check for VLA/EVLA observations
* [check_planets](#check_planets-look-for-solar-system-planets): check for solar system planets
* [check_tgss](#check_tgss): check for matches in TGSSADR1
* [check_nvss](#check_nvss): check for matches in NVSS
* [check_first](#check_first): check for matches in FIRST
* [check_all](#check_all-check-against-all-available-services): query all available services

---
## `check_gaia`: look for matches in Gaia (currently DR3)
### Search for a single source, with position corrected to a given time:
```
check_gaia -c "11h05m21.52536s,+43d31m34.9932s" -t "2023-12-21 21:07:40" --radius=45 -vv
DEBUG   : Input time is '2023-12-21 21:07:40.000'
INFO    : For source at '11h05m21.53s, +43d31m35.0s' = '166.340d, +43.526d', found 2 Gaia matches within 45.0 arcsec
VAST J1105.3+4332	Gaia DR3 778947608243864320:  6.3 arcsec
VAST J1105.3+4332	Gaia DR3 778947814402602752: 37.5 arcsec
```

### Search through XML files produced by Emil with CP sources.
By default, only looks for sources with name `unknown`:
```
check_gaia -x ~/Downloads/RACS-low3_cp_summary.xml                                      
VAST J1706.5-0824[unknown]	Gaia DR3 4336579344345914240:  0.6 arcsec
VAST J1706.5-0824[unknown]	Gaia DR3 4336579348644101120:  5.3 arcsec
VAST J1706.5-0824[unknown]	Gaia DR3 4336579348644340224: 13.7 arcsec
```

### Search through XML files, look for a specific source
```
check_gaia -x ~/Downloads/RACS-low3_cp_summary.xml -k "* alf For B"
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349158319689344:  1.2 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059348952161258624:  5.9 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349192679586304: 10.4 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349158319689344:  0.9 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059348952161258624:  6.0 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349192679586304: 10.1 arcsec
```

### Search through XML files, look for all sources
May be slow, as searches are not vectorized
```
check_gaia -x ~/Downloads/RACS-low3_cp_summary.xml -k all          
VAST J1554.6-2515[*   3 Sco]	Gaia DR3 6235747125966268928:  0.7 arcsec
VAST J1554.6-2515[*   3 Sco]	Gaia DR3 6235747125967029120: 10.2 arcsec
VAST J1554.6-2515[*   3 Sco]	Gaia DR3 6235747881879907712: 13.6 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349158319689344:  1.2 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059348952161258624:  5.9 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349192679586304: 10.4 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349158319689344:  0.9 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059348952161258624:  6.0 arcsec
VAST J0312.0-2859[* alf For B]	Gaia DR3 5059349192679586304: 10.1 arcsec
VAST J0734.5+3153[* alf Gem A]	Gaia DR3 892348694913501952:  8.0 arcsec
VAST J1955.3+0624[* bet Aql]	Gaia DR3 4296708789289490816:  0.6 arcsec
VAST J1955.3+0624[* bet Aql]	Gaia DR3 4296708789290712064: 12.8 arcsec
VAST J0743.3+2853[* sig Gem]	Gaia DR3 878467085735262720:  4.4 arcsec
...
```

### API
```python
from vast_mw import vast_mw
results = vast_mw.check_gaia(source, t=..., radius=...)
```
`source` is a `astropy.coordinates.SkyCoord`.  If an `obstime` is supplied as part of that object, then the Gaia coordinates are corrected (proper motion only) to that time.  Otherwise the time can be specified with `t` (`astropy.time.Time`). 

The returned object is a dictionary containing pairs of Gaia ID, angular separation.

---
## `check_simbad`: look for matches in Simbad (currently DR3)
### Search for a single source, with position corrected to a given time:
```
check_simbad -r 24.771674208211856 -d -17.947682860008488 -t 2016 -vv --radius=60                
DEBUG   : Input time is '2016-01-01 00:00:00.000'
INFO    : For source at '1h39m05.20s, -17d56m51.7s' = '24.772d, -17.948d', found 3 Gaia matches within 60.0 arcsec
VAST J0139.0-1757	G 272-61B:  0.0 arcsec	https://simbad.u-strasbg.fr/simbad/sim-id?Ident=G+272-61B&submit=SIMBAD+search&NbIdent=1
VAST J0139.0-1757	G 272-61:  1.6 arcsec	https://simbad.u-strasbg.fr/simbad/sim-id?Ident=G+272-61&submit=SIMBAD+search&NbIdent=1
VAST J0139.0-1757	G 272-61A:  2.3 arcsec	https://simbad.u-strasbg.fr/simbad/sim-id?Ident=G+272-61A&submit=SIMBAD+search&NbIdent=1
```
Note that the positions in Simbad are epoch 2000.  This means you need a large search radius in cases like these (which is UV Ceti).  The query first looks for matches with the epoch 2000 position and the requested position, and only for those matches are the separations at the desired epoch (2016 in this case) computed.  (I think this is a bug in the Simbad query interface but don't have confirmation yet.)

Note that the Simbad query also prints out the single-source URLs.

### API
```python
from vast_mw import vast_mw
results = vast_mw.check_simbad(source, t=..., radius=...)
```
`source` is a `astropy.coordinates.SkyCoord`.  If an `obstime` is supplied as part of that object, then the Gaia coordinates are corrected (proper motion only) to that time.  Otherwise the time can be specified with `t` (`astropy.time.Time`). 

The returned object is a dictionary containing pairs of Simbad ID, angular separation. The URL can be further constructed by `vast_mw.simbad_url(name)`.

---
## `check_pulsarscraper`: Search for pulsars in ATNF or unpublished catalogs

### Search through XML files, look for a specific source
```
check_pulsarscraper -x ~/Downloads/RACS-low3_cp_summary.xml -k "PSR B0301+19"
VAST J0304.5+1933[PSR B0301+19]	B0301+19[ATNF]:  0.0 deg
```

### API
```python
from vast_mw import vast_mw
results = vast_mw.check_pulsarscraper(source, radius=...)
```
`source` is a `astropy.coordinates.SkyCoord`.  

The returned object is a dictionary containing pairs of Pulsar name (including survey), angular separation.

---
## `check_atnf`: Search for pulsars in ATNF catalog

### Search through XML files, look for a specific source
```
check_atnf -x ~/Downloads/RACS-low3_cp_summary.xml -k "PSR B0301+19"
VAST J0304.5+1933[PSR B0301+19]	J0304+1932:  4.5 arcsec
```

### API
```python
from vast_mw import vast_mw
results = vast_mw.check_atnf(source, radius=...)
```
`source` is a `astropy.coordinates.SkyCoord`.  

The returned object is a dictionary containing pairs of Pulsar name, angular separation.

---

## `check_planets`: look for solar system planets
### Search for all solar system planets, with positions corrected to a given time:
```
check_planets -c "347.498197, -7.741003" -t "2024-09-10 15:24:05"  -r 60 -vv
DEBUG   : Input time is '2024-09-10 15:24:05.000'
INFO    : For source at '23h09m59.57s, -07d44m27.6s' = '347.498d, -7.741d', found 1 planets within 60 arcsec
VAST J2309.9-0744	saturn: 0.01 deg
```
Note that you might use a larger radius here, since the planet ephemerides can be a bit uncertain.

### API
```python
from vast_mw import vast_mw
results = vast_mw.check_planets(source, t=..., radius=..., obs=...)
```
`source` is a `astropy.coordinates.SkyCoord`.  An `obstime` must be supplied as part of that object or as a separate argument, or with `t` (`astropy.time.Time`). `obs` is the name of an observatory.

The returned object is a dictionary containing pairs of planet name, angular separation.

---
## `check_casda`: Check for ASKAP observations

Note that this interface is different from the others, and is not included in `check_all`

### Find ASKAP observations covering a source and after a particular time
```
check_casda -c "12:34:56,+39:00:00" --tstart=60000                            
   obs_id         t_min                 start          t_exptime     Frequency              obs_collection         
                    d                                      s            MHz                                        
----------- ------------------ ----------------------- --------- ----------------- --------------------------------
ASKAP-55588      60298.9521875 2023-12-20 22:51:09.000   905.748 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55589  60298.96336226852 2023-12-20 23:07:14.500   905.748 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55590  60298.97453587963 2023-12-20 23:23:19.900   895.795 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55591 60298.985594907404 2023-12-20 23:39:15.400   905.748 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55665  60299.95082175926 2023-12-21 22:49:11.000   895.795 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55666 60299.961880787036 2023-12-21 23:05:06.500   905.748 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55667  60299.97305555556 2023-12-21 23:21:12.000   895.795 943.4909994999999 The Rapid ASKAP Continuum Survey
ASKAP-55668  60299.98434490741 2023-12-21 23:37:27.400   905.748 943.4909994999999 The Rapid ASKAP Continuum Survey
```

### API
```
from vast_mw import vast_mw
result = vast_mw.check_casda(source,allcolumns=True)
```
`source` is a `astropy.coordinates.SkyCoord`.  

The returned object is a `astropy.table.Table`.

---
## `check_vla`: Check for VLA or EVLA observations

Note that this interface is different from the others, and is not included in `check_all`

### Find VLA/EVLA observations covering a source and after a particular time
```
check_vla -r 18:32:48.41 --dec="-09:11:15.8" -v --group --tstart=59000 
INFO    : For source at '18h32m48.41s, -09d11m15.8s' = '278.202d, -9.188d', found 3 VLA/EVLA matches within None arcsec
               obs_publisher_did                  target_name       Separation           t_min              t_max        t_exptime   freq_max   configuration
                                                                      arcmin                                                                                 
------------------------------------------------ -------------- ------------------ ------------------ ------------------ --------- ------------ -------------
 20A-087.sb38658939.eb38659647.59112.17612017361    G22.90+0.00  16.36948764293314  59112.17612847222 59112.266087951386      44.9  464000000.0         {'B'}
 20A-087.sb37821358.eb38688477.59125.91483408565    G22.90+0.00  16.36948764293314 59125.914842592596  59126.00370717593     44.85  464000000.0         {'B'}
21A-285.sb39768958.eb39990570.59405.256784178244 HESS_J1832-093 10.675953504240203 59405.256791087966  59405.31905960648     839.7 1968000000.0         {'C'}
```
These are additionally grouped by scheduling block.

### API
```
from vast_mw import vast_mw
result = vast_mw.check_vla(source,allcolumns=True)
```
`source` is a `astropy.coordinates.SkyCoord`.  

The returned object is a `astropy.table.Table`.

---

## `check_all`: Check against all available services

### Search for a specific source against all services
```
check_all -r 24.771674208211856 -d -17.947682860008488 -t 2016 -vv --radius=60
DEBUG   : Input time is '2016-01-01 00:00:00.000'
INFO    : For source at '1h39m05.20s, -17d56m51.7s' = '24.772d, -17.948d', found 4 Gaia matches within 60.0 arcsec
VAST J0139.0-1757	Gaia DR3 5140693571158946048:  0.0 arcsec
VAST J0139.0-1757	Gaia DR3 5140693571158739712:  0.5 arcsec
VAST J0139.0-1757	Gaia DR3 5140693571158739840:  2.3 arcsec
VAST J0139.0-1757	Gaia DR3 5140693773021476224: 33.1 arcsec
INFO    : For source at '1h39m05.20s, -17d56m51.7s' = '24.772d, -17.948d', found 3 Simbad matches within 60.0 arcsec
VAST J0139.0-1757	G 272-61B:  0.0 arcsec
VAST J0139.0-1757	G 272-61:  1.6 arcsec
VAST J0139.0-1757	G 272-61A:  2.3 arcsec
WARNING : For source at '1h39m05.20s, -17d56m51.7s' = '24.772d, -17.948d', found 0 Pulsar Survey Scraper matches within 60.0 arcsec
WARNING : For source at '1h39m05.20s, -17d56m51.7s' = '24.772d, -17.948d', found 0 ATNF Pulsar Catalog matches within 60.0 arcsec
```
