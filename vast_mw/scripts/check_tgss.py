import numpy as np
from astropy import units as u, constants as c
from astropy.coordinates import SkyCoord
from astropy.table import Table
import argparse
from loguru import logger as log
import sys
import warnings

from vast_mw import vast_mw


def main():
    parser = argparse.ArgumentParser(
        description="Query TGSSADR1 catalog for a number of positions",
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
    if args.verbosity == 1:
        log.remove()
        log.add(sys.stderr, format=vast_mw.logformat, level="INFO")
    elif args.verbosity >= 2:
        log.remove()
        log.add(sys.stderr, format=vast_mw.logformat, level="DEBUG")

    sources, names = vast_mw._parse_input(args, require_time=False)
    if sources is None or len(sources) == 0:
        sys.exit(1)
    for name, source in zip(names, sources):
        results = vast_mw.check_tgss(source, radius=args.radius * u.arcsec)
        level = log.info if len(results) > 0 else log.warning
        level(
            f"For source at '{vast_mw.format_radec(source)}' = '{vast_mw.format_radec_decimal(source)}', found {len(results)} TGSSADR1 matches within {args.radius} arcsec"
        )
        for k, v in sorted(results.items(), key=lambda x: x[1]):
            s = vast_mw.format_name(source)
            if name is not None:
                s += f"[{name}]"
            out = f"{s}\t{k}: {v:4.1f}"
            if args.url:
                out += f"\t{vast_mw.vizier_url(catalog='J/A+A/598/A78/table3',idname='TGSSADR',sourcename=k)}"
            print(out)
