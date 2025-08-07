import numpy as np
from astropy import units as u, constants as c
from astropy.table import Table, Column
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.casda import Casda
from astroquery.vizier import Vizier
import pyvo as vo
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
_vizier_url = "https://vizier.cds.unistra.fr/viz-bin/VizieR-5"
cSimbad = Simbad()
cSimbad.add_votable_fields("pmra", "pmdec")


services = {
    "Gaia": ["check_gaia", "gaia_url"],
    "Simbad": ["check_simbad", "simbad_url"],
    "Pulsar Survey Scraper": ["check_pulsarscraper", None],
    "ATNF Pulsar Catalog": ["check_atnf", None],
    "Planets": ["check_planets", None],
    "TGSS": ["check_tgss", None],
    "FIRST": ["check_first", "vizier_url"],
    "NVSS": ["check_nvss", None],
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


def vizier_url(catalog: str, idname: str, sourcename: str) -> str:
    return f"{_vizier_url}?{urllib.parse.urlencode({'-source': catalog, idname: sourcename})}"


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
    """Check a source against public ASKAP observations on CASDA

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


def check_vla(
    source: SkyCoord,
    radius: u.Quantity = None,
    tstart: Time = None,
    tstop: Time = None,
    allcolumns: bool = False,
    group: bool = False,
) -> Table:
    """Check a source against ATNF pulsar catalog, correcting for proper motion

    Parameters
    ----------
    source : SkyCoord
    radius : u.Quantity, optional
        Search radius for cone search.  If None, will attempt to check against actual FOV.
    tstart : Time, optional
        Will only return observations >= this time if supplied
    tstop : Time, optional
        Will only return observations <= this time if supplied
    allcolumns : bool, optional
        Will only return a subset of columns unless this is True
    group : bool, optional
        Group together individual scans within a block

    Returns
    -------
    Table
        Table of matching observations
    """
    checkfov = False
    if radius is None:
        checkfov = True
        # use a pretty large radius to start
        radius = 1.5 * (0.2 * u.m / (25 * u.m)).to(
            u.degree, equivalencies=u.dimensionless_angles()
        )
    service = vo.dal.TAPService("https://data-query.nrao.edu/tap")
    query = (
        "SELECT * FROM tap_schema.obscore WHERE CONTAINS(POINT('ICRS',s_ra,s_dec),CIRCLE('ICRS',"
        + str(source.ra.degree)
        + ","
        + str(source.dec.degree)
        + ","
        + str(radius.value)
        + "))=1 and (instrument_name='EVLA' OR instrument_name='VLA')"
    )
    result = service.search(query)
    output = result.to_table()
    wl = ((output["freq_max"] + output["freq_min"]) / 2 * u.Hz).to(
        u.m, equivalencies=u.spectral()
    )
    # compare to https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/fov
    hpbw = 1.2 * (wl / (25 * u.m)).to(u.arcmin, equivalencies=u.dimensionless_angles())
    obs_coords = SkyCoord(ra=output["s_ra"] * u.degree, dec=output["s_dec"] * u.degree)
    output.add_column(
        Column(source.separation(obs_coords).to(u.arcmin), name="Separation")
    )
    good = np.ones(len(output), dtype=bool)
    if checkfov:
        good = good & (source.separation(obs_coords) < hpbw / 2)
    if tstart is not None:
        good = good & (Time(output["t_min"], format="mjd") >= tstart)
    if tstop is not None:
        good = good & (Time(output["t_max"], format="mjd") <= tstop)
    if np.any(~good):
        output = output[good]
    output.sort("t_min")
    if group:
        unique_obs = set(
            [(x[0], x[1]) for x in output["obs_publisher_did", "target_name"]]
        )
        freq_max = np.zeros(len(unique_obs))
        t_min = np.zeros(len(unique_obs))
        t_max = np.zeros(len(unique_obs))
        t_exptime = np.zeros(len(unique_obs))
        separation = np.zeros(len(unique_obs)) * u.arcmin
        configuration = []
        obs_publisher_did = []
        target_name = []
        for i, obs in enumerate(unique_obs):
            obs_publisher_did.append(obs[0])
            target_name.append(obs[1])
            match = np.where(
                (output["obs_publisher_did"] == obs[0])
                & (output["target_name"] == obs[1])
            )[0]
            freq_max[i] = output[match]["freq_max"].max()
            t_min[i] = output[match]["t_min"].min()
            t_max[i] = output[match]["t_max"].max()
            t_exptime[i] = output[match]["t_exptime"].sum()
            separation[i] = output[match]["Separation"].quantity.min()
            configuration.append(set(output[match]["configuration"]))
        grouped_output = Table(
            [
                Column(obs_publisher_did, name="obs_publisher_did"),
                Column(target_name, name="target_name"),
                Column(separation, name="Separation"),
                Column(t_min, name="t_min"),
                Column(t_max, name="t_max"),
                Column(t_exptime, name="t_exptime"),
                Column(freq_max, name="freq_max"),
                Column(configuration, name="configuration"),
            ]
        )
        grouped_output.sort("t_min")
        return grouped_output
    if allcolumns:
        return output
    else:
        return output[
            "obs_publisher_did",
            "target_name",
            "Separation",
            "t_min",
            "t_max",
            "t_exptime",
            "freq_max",
            "configuration",
        ]


def check_planets(
    source: SkyCoord,
    t: Time = None,
    radius: u.Quantity = 1 * u.arcmin,
    obs: str = "mwa",
) -> Dict[str, u.Quantity]:
    """Check a source against solar system planets, correcting for proper motion

    Parameters
    ----------
    source : SkyCoord
    t : Time, optional
        Will override ``source.obstime`` is supplied, or if ``source.obstime`` is not supplied
    radius : u.Quantity, optional
        Search radius for cone search
    obs : str, optional
        Observatory name

    Returns
    -------
    dict
        Pairs of planet name and angular separation
    """
    if t is None:
        if source.obstime is None:
            log.error(
                "Must supply either SkyCoord with obstime or separate time for coordinate check"
            )
            return {}
        t = source.obstime
    loc = EarthLocation.of_site(obs)
    separations = {}
    with solar_system_ephemeris.set("builtin"):
        for planet_name in solar_system_ephemeris.bodies:
            planet = get_body(planet_name, t, loc)
            if planet.separation(source) < radius:
                separations[planet_name] = planet.separation(source)

    return separations


def check_tgss(
    source: SkyCoord, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against the TGSS ADR1

    Parameters
    ----------
    source : SkyCoord
    radius : u.Quantity, optional
        Search radius for cone search

    Returns
    -------
    dict
        Pairs of TGSSADR1 identifier and angular separation
    """
    result = Vizier().query_region(source, radius=radius, catalog="J/A+A/598/A78")
    out = {}
    for r in result:
        names = r["TGSSADR"]
        matchpos = SkyCoord(r["RAJ2000"], r["DEJ2000"])
        for i in range(len(r)):
            out[names[i]] = matchpos[i].separation(source)
    return out


def check_first(
    source: SkyCoord, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against FIRST

    Parameters
    ----------
    source : SkyCoord
    radius : u.Quantity, optional
        Search radius for cone search

    Returns
    -------
    dict
        Pairs of FIRST identifier and angular separation
    """
    result = Vizier().query_region(source, radius=radius, catalog="VIII/92/first14")
    out = {}
    for r in result:
        names = r["FIRST"]
        matchpos = SkyCoord(r["RAJ2000"], r["DEJ2000"], unit=("hour", "deg"))
        for i in range(len(r)):
            out[f"{names[i]}"] = matchpos[i].separation(source)
    return out


def check_nvss(
    source: SkyCoord, radius: u.Quantity = 15 * u.arcsec
) -> Dict[str, u.Quantity]:
    """Check a source against NVSS

    Parameters
    ----------
    source : SkyCoord
    radius : u.Quantity, optional
        Search radius for cone search

    Returns
    -------
    dict
        Pairs of NVSS identifier and angular separation
    """
    result = Vizier().query_region(source, radius=radius, catalog="VIII/65/nvss")
    out = {}
    for r in result:
        names = r["NVSS"]
        matchpos = SkyCoord(r["RAJ2000"], r["DEJ2000"], unit=("hour", "deg"))
        for i in range(len(r)):
            out[names[i]] = matchpos[i].separation(source)
    return out
