import numpy as np
from astropy import units as u, constants as c
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.table import Table
import argparse
from loguru import logger as log
import sys
import warnings

from vast_mw import vast_mw


def parse_time(t):
    time_format = "iso" if "-" in t else "decimalyear"
    if time_format == "decimalyear":
        t = float(t)
        if t > 50000:
            time_format = "mjd"
    try:
        return Time(t, format=time_format)
    except ValueError:
        log.error(f"Cannot parse input time '{t}' with input format '{time_format}'")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Query CASDA for a position",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--coord",
        type=str,
        default=None,
        help="Input coordinates (separated by ,)",
    )
    parser.add_argument(
        "-r",
        "--ra",
        type=str,
        default=None,
        help="RA to search (hms or deg) - only if not specified with -c",
    )
    parser.add_argument(
        "-d",
        "--dec",
        type=str,
        default=None,
        help="Dec to search (deg) - only if not specified with -c",
    )
    parser.add_argument(
        "--tstart",
        type=str,
        default=None,
        help="Start time for query",
    )
    parser.add_argument("--tstop", type=str, default=None, help="Stop time for query")
    parser.add_argument(
        "--radius",
        default=None,
        type=float,
        help="Search radius (arcsec): if None, use FOV",
    )
    parser.add_argument(
        "--allcolumns", default=False, action="store_true", help="Return all columns"
    )
    parser.add_argument(
        "--group",
        default=False,
        action="store_true",
        help="Group observations",
    )
    parser.add_argument(
        "-v", "--verbosity", default=0, action="count", help="Increase output verbosity"
    )

    args = parser.parse_args()
    log.remove()
    log.add(sys.stderr, format=vast_mw.logformat, level="WARNING")
    if args.verbosity == 1:
        log.remove()
        log.add(sys.stderr, format=vast_mw.logformat, level="INFO")
    elif args.verbosity >= 2:
        log.remove()
        log.add(sys.stderr, format=vast_mw.logformat, level="DEBUG")

    if args.coord is not None:
        ra, dec = args.coord.split(",")
    elif args.ra is not None and args.dec is not None:
        ra, dec = args.ra, args.dec

    ra_units = "hour" if any(x in ra for x in [" ", ":", "h"]) else "deg"
    dec_units = "deg"
    tstart = None
    if args.tstart is not None:
        tstart = parse_time(args.tstart)
        if tstart is None:
            sys.exit(1)
        log.debug(f"Input start time is '{tstart.iso}'")
    tstop = None
    if args.tstop is not None:
        tstop = parse_time(args.tstop)
        if tstop is None:
            sys.exit(1)
        log.debug(f"Input stop time is '{tstop.iso}'")
    try:
        source = SkyCoord(ra, dec, unit=(ra_units, dec_units))
    except:
        log.error(
            f"Cannot parse input coordinates '{ra}, {dec}' with input units '{ra_units}, {dec_units}'"
        )
        sys.exit(1)
    results = vast_mw.check_vla(
        source,
        radius=args.radius * u.arcsec if args.radius is not None else args.radius,
        tstart=tstart,
        tstop=tstop,
        allcolumns=args.allcolumns,
        group=args.group,
    )
    level = log.info if len(results) > 0 else log.warning
    level(
        f"For source at '{vast_mw.format_radec(source)}' = '{vast_mw.format_radec_decimal(source)}', found {len(results)} VLA/EVLA matches within {args.radius} arcsec"
    )
    print(results)
