# Reference data

`spa_reference_testdata.csv` checks the unchanged position algorithm against NREL's
reference positions. `test/deltat/` contains the existing delta-T reference series.

`usno/` contains sunrise/sunset and civil-twilight tables from the
[US Naval Observatory](https://aa.usno.navy.mil/data/RS_OneYear), including their
original text under `raw/`. The comparison uses UTC dates and a 90-second tolerance
for these minute-resolution tables. Event counts and polar states must match.

`jpl_events.csv` is the same independent fixture as the Java library: 14 UTC
dates/locations at four horizons, including polar transitions, exact poles,
one-event dates, repeated events and dates without transit. Regenerate it with:

```sh
python3 scripts/generate-jpl-events.py /path/to/de440s.bsp > tests/data/jpl_events.csv
```

The recorded versions are Skyfield 1.55, jplephem 2.24, NumPy 2.5.3 and sgp4 2.27.
Download [DE440s](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp).
(SHA-256: `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`.)
Both sides use fixed TT−UT1 = 69.184 s, label UT1 instants as UTC, and use apparent,
unrefracted topocentric solar-centre positions on WGS84 at sea level.

Counts must match exactly. The rise/set tolerance is SPA's stated 0.0003° angular
uncertainty divided by local elevation speed, plus 2 ms. Transit tolerance is one
second. Very shallow polar crossings can differ by seconds despite small angular
errors. The generator's five-minute sampling is sufficient for these selected
cases; analytically known curves separately test close crossings and tangencies.
