# Present-day Earth data

`noncompiled/generators/EarthData.js` holds the real Earth used by the **Earth** preset. These scripts rebuild it from
the original datasets. You only need them to change the data; the site itself just loads `EarthData.js`.

## Sources

| File | What | Source |
|---|---|---|
| `etopo1_ice_0p1.nc` | Elevation and bathymetry (ice sheet surface), 1 arc-minute grid sampled every 0.1° | ETOPO1, Amante & Eakins (2009), NOAA NGDC, via the [ERDDAP etopo180](https://coastwatch.pfeg.noaa.gov/erddap/info/etopo180/index.html) dataset |
| `age.2020.1.GTS2012.6m.grd` | Age of oceanic crust, 0.1° | Seton et al. (2020), G-cubed, doi:10.1029/2020GC009214, [EarthByte](https://www.earthbyte.org/webdav/ftp/earthbyte/agegrid/2020/Grids/), CC BY 4.0 |
| `PB2002_plates.dig.txt` | Outlines of the 52 plates | Bird (2003), G-cubed 4(3), 1027, doi:10.1029/2001GC000252, [peterbird.name](http://peterbird.name/oldFTP/PB2002/) |
| `PB2002_poles.dat.txt` | Euler poles of the 52 plates relative to the Pacific plate | Bird (2003), same as above |

## Rebuilding

Needs Python 3 with numpy (Pillow is optional, for preview images). No netCDF library is needed; `nc3.py` reads the
netCDF-3 files directly.

```sh
mkdir -p raw && cd raw
curl -o etopo1_ice_0p1.nc "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude%5B(-90):6:(90)%5D%5B(-180):6:(180)%5D"
curl -o age.2020.1.GTS2012.6m.grd "https://www.earthbyte.org/webdav/ftp/earthbyte/agegrid/2020/Grids/age.2020.1.GTS2012.6m.grd"
curl -o PB2002_plates.dig.txt "http://peterbird.name/oldFTP/PB2002/PB2002_plates.dig.txt"
curl -o PB2002_poles.dat.txt  "http://peterbird.name/oldFTP/PB2002/PB2002_poles.dat.txt"
cd ..
python preprocess.py   # writes out/earth_0p5.npz and preview PNGs
python export_js.py    # writes ../../noncompiled/generators/EarthData.js
```

## What the scripts do

* **Elevation** is averaged over each 0.5° pixel (area-weighted), and the game averages those pixels again over each
  grid cell, so cells carry their true mean elevation rather than a point sample.
* **Continental fraction** is the share of each pixel *not* covered by dated oceanic crust in the age grid.
* **Plates** are rasterized with a spherical winding-number test, which handles plates that cross the antimeridian or
  contain a pole. The 52 outlines tile the sphere exactly (their areas sum to 4π with no gaps or overlaps).
* **Euler poles** are converted from Pacific-fixed to the no-net-rotation (NNR) frame by requiring the area-weighted
  net rotation of all plates to vanish. The resulting Pacific pole (63.2°S, 107.7°E, 0.649°/Myr) agrees with the
  published NNR-NUVEL-1A Pacific pole (63.0°S, 107.4°E, 0.641°/Myr).
* **Depth vs. age**: the median seafloor depth in 5 Myr age bins over 0–105 Myr gives
  `depth = -3520 m - 18.9 m/Myr × age`. The game fits its effective oceanic plate thickness and sea level datum to this
  line, which also puts continental crust at a realistic ~41 km on average.
