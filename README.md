# vast-mw
Tools for multi-wavelength searches of VAST objects

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

### Search through XML files, look for a all sources
May be slow, as searches are not vecetorized
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
```

