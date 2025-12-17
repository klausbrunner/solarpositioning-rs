#!/usr/bin/env python3
"""
Helper to fetch rise/transit/set times from JPL Horizons for a given site.

Uses the Horizons JSON API with a configurable step size and standard RTS markers.
Time zone handling here is simple: a fixed offset in hours added to UTC for display.
"""

import argparse
import datetime as dt
import json
import sys
import urllib.parse
import urllib.request


def build_url(lat: float, lon: float, start: str, stop: str, step_minutes: int) -> str:
    params = {
        "format": "json",
        "COMMAND": "'10'",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": "coord@399",
        "COORD_TYPE": "GEODETIC",
        "SITE_COORD": f"'{lon},{lat},0'",
        "START_TIME": f"'{start}'",
        "STOP_TIME": f"'{stop}'",
        "STEP_SIZE": f"'{step_minutes} m'",
        "QUANTITIES": "'4,9,10'",
    }
    return "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(params)


def fetch(url: str) -> str:
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    payload = json.loads(data)
    if "result" not in payload:
        raise RuntimeError("Unexpected Horizons response")
    return payload["result"]


def parse_events(result: str) -> dict[str, dt.datetime]:
    events: dict[str, dt.datetime] = {}
    in_table = False
    for line in result.splitlines():
        if line.startswith("$$SOE"):
            in_table = True
            continue
        if line.startswith("$$EOE"):
            break
        if not in_table:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        date_str, time_str, marker = parts[0], parts[1], parts[2]
        try:
            timestamp = dt.datetime.strptime(f"{date_str} {time_str}", "%Y-%b-%d %H:%M")
        except ValueError:
            continue
        if "r" in marker:
            events.setdefault("rise", timestamp)
        if "t" in marker:
            events.setdefault("transit", timestamp)
        if "s" in marker:
            events.setdefault("set", timestamp)
    return events


def format_time(ts: dt.datetime, offset_hours: float) -> str:
    local = ts + dt.timedelta(hours=offset_hours)
    return f"{ts.isoformat()}Z  | local {local.isoformat()} (offset {offset_hours:+g}h)"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Fetch sunrise/transit/set from JPL Horizons."
    )
    parser.add_argument("--lat", type=float, required=True, help="Latitude in degrees")
    parser.add_argument(
        "--lon", type=float, required=True, help="Longitude in degrees (east+)"
    )
    parser.add_argument(
        "--date", required=True, help="UTC date YYYY-MM-DD for the request window start"
    )
    parser.add_argument(
        "--tz-offset",
        type=float,
        default=0.0,
        help="Fixed offset hours to add for display (e.g. +5 or -7). Default: 0 (UTC only).",
    )
    parser.add_argument(
        "--step-min",
        type=int,
        default=10,
        help="Step size in minutes for Horizons table (affects RTS resolution). Default: 10.",
    )
    args = parser.parse_args(argv)

    start = f"{args.date} 00:00"
    # Stop at next UTC midnight so we bracket the local day.
    stop_dt = dt.datetime.strptime(args.date, "%Y-%m-%d") + dt.timedelta(days=1)
    stop = stop_dt.strftime("%Y-%m-%d 00:00")

    url = build_url(args.lat, args.lon, start, stop, args.step_min)
    result = fetch(url)
    events = parse_events(result)

    if not events:
        print("No events found in Horizons output", file=sys.stderr)
        return 1

    print(f"Site: lat={args.lat}, lon={args.lon}, start={start} UTC")
    for name in ("rise", "transit", "set"):
        ts = events.get(name)
        if ts:
            print(f"{name.capitalize():7}: {format_time(ts, args.tz_offset)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
