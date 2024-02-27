import numpy as np
from astropy import units as u, constants as c
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from loguru import logger as log
import sys
from typing import Dict
import requests
import warnings

logformat = "<level>{level: <8}</level>: <level>{message}</level>"
log.add(sys.stderr, format=logformat, level="WARNING")

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
scraper_url = "https://pulsar.cgca-hub.org/api"
simbad_url = "https://simbad.u-strasbg.fr/simbad/sim-id"
cSimbad = Simbad()
cSimbad.add_votable_fields("pmra", "pmdec")


def format_radec_decimal(coord: SkyCoord) -> str:
    """Return coordinates as 'ddd.ddd, ddd.ddd'

    Parameters
    ----------
    coord: SkyCoord

    Returns
    -------
    str
    """
    return f"{coord.icrs.ra.to_string(decimal=True,precision=3)}d, {coord.icrs.dec.to_string(decimal=True, alwayssign=True, precision=3)}d"


def format_radec(coord: SkyCoord) -> str:
    """Return coordinates as 'HHhMMmSS.SSs DDdMMmSS.Ss'

    Parameters
    ----------
    coord: SkyCoord

    Returns
    -------
    str
    """
    sra = coord.icrs.ra.to_string(u.hour, decimal=False, sep="hms", precision=2)
    sdec = coord.icrs.dec.to_string(
        u.degree, decimal=False, sep="dms", precision=1, pad=True, alwayssign=True
    )
    return f"{sra}, {sdec}"


def format_name(coord: SkyCoord, prefix: str = "VAST") -> str:
    """Return coordinates as 'VAST JHHMM.M+DDMM'

    Parameters
    ----------
    coord: SkyCoord
    prefix: str, optional
        Name for output

    Returns
    -------
    str
    """
    hms = coord.ra.hms
    return f"{prefix} J{int(hms[0]):02d}{int(hms[1]):02d}.{int(10*hms[2]/60)}{coord.dec.to_string(u.deg,alwayssign=True,sep='',fields=2,pad=True)}"


def check_gaia(
    source: SkyCoord, t: Time = None, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against Gaia, correcting for proper motion

    Parameters
    ----------
    source : SkyCoord
    t : Time, optional
        Will override ``source.obstime`` is supplied, or if ``source.obstime`` is not supplied
    radius : u.Quantity, optional
        Search radius for cone search

    Returns
    -------
    dict
        Pairs of Gaia identifier and angular separation
    """
    if t is None:
        if source.obstime is None:
            log.error(
                "Must supply either SkyCoord with obstime or separate time for coordinate check"
            )
            return {}
        t = source.obstime
    q = Gaia.cone_search(coordinate=source, radius=radius)
    r = q.get_results()
    separations = {}
    for i in range(len(r)):
        gaia_source = SkyCoord(
            r[i]["ra"] * u.deg,
            r[i]["dec"] * u.deg,
            pm_ra_cosdec=r[i]["pmra"] * u.mas / u.yr,
            pm_dec=r[i]["pmdec"] * u.mas / u.yr,
            distance=(
                (r[i]["parallax"] * u.mas).to(u.kpc, equivalencies=u.parallax())
                if r[i]["parallax"] > 0
                else 1 * u.kpc
            ),
            obstime=Time(r[0]["ref_epoch"], format="decimalyear"),
        )
        separations[r[i]["DESIGNATION"]] = (
            gaia_source.apply_space_motion(t).separation(source).arcsec * u.arcsec
        )
    return separations


def check_pulsarscraper(
    source: SkyCoord, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against the Pulsar survey scraper

    Parameters
    ----------
    source : SkyCoord
    radius : u.Quantity, optional
        Search radius for cone search

    Returns
    -------
    dict
        Pairs of pulsar identifier and angular separation
    """
    response = requests.get(
        scraper_url,
        params={
            "type": "search",
            "ra": source.ra.deg,
            "dec": source.dec.deg,
            "radius": radius.to_value(u.deg),
        },
    )
    if not response.ok:
        log.error(
            f"Unable to query pulsarsurveyscraper: received code={response.status_code} ({response.reason})"
        )
        return {}
    out = {}
    for k in response.json():
        if k.startswith("search") or k.startswith("nmatches"):
            continue
        out[f"{k}[{response.json()[k]['survey']['value']}]"] = (
            response.json()[k]["distance"]["value"] * u.deg
        )
    return out


def check_simbad(
    source: SkyCoord, t: Time = None, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against Simbad, correcting for proper motion

    Parameters
    ----------
    source : SkyCoord
    t : Time, optional
        Will override ``source.obstime`` is supplied, or if ``source.obstime`` is not supplied
    radius : u.Quantity, optional
        Search radius for cone search

    Returns
    -------
    dict
        Pairs of Simbad identifier and angular separation
    """
    if t is None:
        if source.obstime is None:
            log.error(
                "Must supply either SkyCoord with obstime or separate time for coordinate check"
            )
            return {}
        t = source.obstime
    # do the query at the requested time although that doesn't seem to work
    # but still need to update positions to get real separations
    # r = cSimbad.query_criteria(
    #     f"region(CIRCLE, ICRS, J{t.jyear}, {t.jyear}, {source.ra.deg} {source.dec.deg},{radius.to_value(u.arcmin)})"
    # )
    r = cSimbad.query_region(source, radius=radius)
    if r is None:
        return {}
    separations = {}
    for i in range(len(r)):
        # simbad gives positions in epoch 2000
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            simbad_source = SkyCoord(
                r[i]["RA"],
                r[i]["DEC"],
                unit=("hour", "deg"),
                pm_ra_cosdec=r[i]["PMRA"] * u.mas / u.yr,
                pm_dec=r[i]["PMDEC"] * u.mas / u.yr,
                obstime=Time(2000, format="decimalyear"),
            )
            separations[r[i]["MAIN_ID"]] = (
                simbad_source.apply_space_motion(t).separation(source).arcsec * u.arcsec
            )
    return separations
