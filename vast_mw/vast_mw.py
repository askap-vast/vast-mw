import numpy as np
from astropy import units as u, constants as c
from astropy.table import Table, Column
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.casda import Casda
import psrqpy
from loguru import logger as log
import sys
from typing import Dict
import requests
import urllib
import warnings

logformat = "<level>{level: <8}</level>: <level>{message}</level>"
log.add(sys.stderr, format=logformat, level="WARNING")

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
_scraper_url = "https://pulsar.cgca-hub.org/api"
_simbad_url = "https://simbad.u-strasbg.fr/simbad/sim-id"
_gaia_url = "https://gaia.ari.uni-heidelberg.de/singlesource.html"
cSimbad = Simbad()
cSimbad.add_votable_fields("pmra", "pmdec")


services = {
    "Gaia": ["check_gaia", "gaia_url"],
    "Simbad": ["check_simbad", "simbad_url"],
    "Pulsar Survey Scraper": ["check_pulsarscraper", None],
    "ATNF Pulsar Catalog": ["check_atnf", None],
}


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
        _scraper_url,
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


def gaia_url(name: str) -> str:
    """Return the single-source Gaia URL for an object

    Parameters
    ----------
    name : str

    Returns
    -------
    str

    """
    if " " in name:
        name = name.split()[-1]
    return f"{_gaia_url}#{urllib.parse.urlencode({'gaiadr3_id': name})}"


def simbad_url(name: str) -> str:
    """Return the single-source Simbad URL for an object

    Parameters
    ----------
    name : str

    Returns
    -------
    str

    """
    return f"{_simbad_url}?{urllib.parse.urlencode({'Ident': name, 'submit': 'SIMBAD search','NbIdent': 1})}"


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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
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


def check_atnf(
    source: SkyCoord, t: Time = None, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against ATNF pulsar catalog, correcting for proper motion

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
        Pairs of Pulsar identifier and angular separation
    """
    if t is None:
        if source.obstime is None:
            log.error(
                "Must supply either SkyCoord with obstime or separate time for coordinate check"
            )
            return {}
        t = source.obstime
    r = psrqpy.QueryATNF(
        params=["JNAME", "RAJD", "DECJD", "POSEPOCH", "PMRA", "PMDEC"],
        circular_boundary=[
            source.ra.to_string(u.hour, decimal=False, sep="hms", precision=2),
            source.dec.to_string(
                u.degree,
                decimal=False,
                sep="dms",
                precision=1,
                pad=True,
                alwayssign=True,
            ),
            radius.to_value(u.deg),
        ],
    )
    if r is None or len(r) == 0:
        return {}
    separations = {}
    for i in range(len(r)):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            atnf_source = SkyCoord(
                r.table[i]["RAJD"] * u.deg,
                r.table[i]["DECJD"] * u.deg,
                pm_ra_cosdec=r.table[i]["PMRA"] * u.mas / u.yr,
                pm_dec=r.table[i]["PMDEC"] * u.mas / u.yr,
                obstime=Time(r.table[i]["POSEPOCH"], format="mjd"),
            )
            separations[r.table[i]["JNAME"]] = (
                atnf_source.apply_space_motion(t).separation(source).arcsec * u.arcsec
            )
    return separations


def check_casda(
    source: SkyCoord,
    radius: u.Quantity = 15 * u.arcsec,
    tstart: Time = None,
    tstop: Time = None,
    vastonly: bool = False,
    allcolumns: bool = False,
) -> Table:
    """Check a source against ATNF pulsar catalog, correcting for proper motion

    Parameters
    ----------
    source : SkyCoord
    radius : u.Quantity, optional
        Search radius for cone search
    tstart : Time, optional
        Will only return observations >= this time if supplied
    tstop : Time, optional
        Will only return observations <= this time if supplied
    vastonly : bool, optional
        Will only return VAST observations if True
    allcolumns : bool, optional
        Will only return a subset of columns unless this is True

    Returns
    -------
    Table
        Table of matching observations
    """
    result = Casda.query_region(source, radius=radius)
    # try to filter to get only a single Stokes I entry per observation
    filter = (
        (result["dataproduct_type"] == "cube")
        & (result["pol_states"] == "/I/")
        & (result["facility_name"] == "ASKAP")
        & (result["dataproduct_subtype"] == "cont.restored.t0")
        & (
            (np.array(["conv" in x["filename"] for x in result]))
            | (np.array(["restored.fcor" in x["filename"] for x in result]))
        )
    )
    if vastonly:
        filter = filter & (result["obs_collection"] == "VAST")
    filtered_result = result[filter]
    freq = (result[filter]["em_max"].data.data * u.m).to(
        u.MHz, equivalencies=u.spectral()
    )
    filtered_result.add_column(Column(freq, name="Frequency"))
    filtered_result.add_column(
        Column(Time(filtered_result["t_min"], format="mjd").iso, name="start")
    )
    filtered_result.sort("t_min")
    if (tstart is not None) and (len(filtered_result) > 0):
        filtered_result = filtered_result[filtered_result["start"] >= tstart]
    if (tstop is not None) and (len(filtered_result) > 0):
        filtered_result = filtered_result[filtered_result["start"] <= tstop]
    if allcolumns:
        return filtered_result
    else:
        return filtered_result[
            "obs_id", "t_min", "start", "t_exptime", "Frequency", "obs_collection"
        ]
