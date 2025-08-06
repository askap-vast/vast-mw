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


def main():
    parser = argparse.ArgumentParser(
        description="Query planet positions",
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
        "-t", "--time", type=str, default=None, help="Time for input coordinates"
    )
    parser.add_argument(
        "-x", "--xml", default=None, type=str, help="XML table for input"
    )
    parser.add_argument(
        "-k", "--key", default="unknown", help="Source name in XML to search (or 'all')"
    )
    parser.add_argument("-o","--obs",default="mwa",help="Observatory name")
    parser.add_argument(
        "--radius", default=60, type=float, help="Search radius (arcsec)"
    )
    parser.add_argument(
        "-v", "--verbosity", default=0, action="count", help="Increase output verbosity"
    )

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()    
    log.remove()
    log.add(sys.stderr, format=vast_mw.logformat, level="WARNING")
    if args.verbosity == 1:
        log.remove()
        log.add(sys.stderr, format=vast_mw.logformat, level="INFO")
    elif args.verbosity >= 2:
        log.remove()
        log.add(sys.stderr, format=vast_mw.logformat, level="DEBUG")
    if args.xml is not None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=u.UnitsWarning)
            data = Table.read(args.xml)
        if args.key != "all":
            data = data[data["src"] == args.key]
        sources = SkyCoord(
            data["ra_deg"] * u.deg, data["dec_deg"] * u.deg, obstime=data["scan_start"]
        )
        names = data["src"]
        log.debug(f"Found {len(sources)} sources in '{args.xml}'")
    else:
        ra,dec=None,None
        if args.coord is not None:
            ra, dec = args.coord.split(",")
        elif args.ra is not None and args.dec is not None:
            ra, dec = args.ra, args.dec
        if ra is None or dec is None:
            raise ValueError("Must supply a XML file or source coordinate")

        ra_units = "hour" if any(x in ra for x in [" ", ":", "h"]) else "deg"
        dec_units = "deg"
        time_format = "iso" if "-" in args.time else "decimalyear"
        if time_format == "decimalyear":
            args.time = float(args.time)
            if args.time > 50000:
                time_format = "mjd"
        try:
            t = Time(args.time, format=time_format)
        except ValueError:
            log.error(
                f"Cannot parse input time '{args.time}' with input format '{time_format}'"
            )
            sys.exit(1)
        log.debug(f"Input time is '{t.iso}'")
        try:
            source = SkyCoord(ra, dec, unit=(ra_units, dec_units), obstime=t)
        except:
            log.error(
                f"Cannot parse input coordinates '{ra}, {dec}' with input units '{ra_units}, {dec_units}'"
            )
            sys.exit(1)
        sources = [source]
        names = [None]
    for name, source in zip(names, sources):
        results = vast_mw.check_planets(source, radius=args.radius * u.arcsec, obs=args.obs)
        level = log.info if len(results) > 0 else log.warning
        level(
            f"For source at '{vast_mw.format_radec(source)}' = '{vast_mw.format_radec_decimal(source)}', found {len(results)} planets within {args.radius} arcsec"
        )
        for k, v in sorted(results.items(), key=lambda x: x[1]):
            s = vast_mw.format_name(source)
            if name is not None:
                s += f"[{name}]"
            out = f"{s}\t{k}: {v:4.2f}"
            print(out)
