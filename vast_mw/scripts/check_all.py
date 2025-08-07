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
        description="Query all services for a number of positions",
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
    parser.add_argument(
        "--radius", default=15, type=float, help="Search radius (arcsec)"
    )
    parser.add_argument("-u", "--url", action="store_true", help="Return URL")
    parser.add_argument(
        "-v", "--verbosity", default=0, action="count", help="Increase output verbosity"
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    log.remove()
    log.add(sys.stderr, format=vast_mw.logformat, level="WARNING")

    sources, names = vast_mw._parse_input(args, require_time=True)
    if sources is None or len(sources) == 0:
        sys.exit(1)
    for name, source in zip(names, sources):
        for service in vast_mw.services:
            function_to_call = getattr(vast_mw, vast_mw.services[service][0])
            results = function_to_call(source, radius=args.radius * u.arcsec)
            level = log.info if len(results) > 0 else log.warning
            level(
                f"For source at '{vast_mw.format_radec(source)}' = '{vast_mw.format_radec_decimal(source)}', found {len(results)} {service} matches within {args.radius} arcsec"
            )
            for k, v in sorted(results.items(), key=lambda x: x[1]):
                s = vast_mw.format_name(source)
                if name is not None:
                    s += f"[{name}]"
                out = f"{s}\t{k}: {v:4.1f}"
                if vast_mw.services[service][1] is not None and args.url:
                    out += f"\t{getattr(vast_mw, vast_mw.services[service][1])(k)}"
                print(out)
