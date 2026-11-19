# Tectonics.js
a 3d planet simulator in your web browser

## End user quickstart guide:
Just go [here](http://davidson16807.github.io/tectonics.js/) if you want to mess around with the model.

If you want to run the model locally, download the zipfile [here](https://github.com/davidson16807/tectonics.js/archive/gh-pages.zip), unzip the contents, and open "index.html" in the browser of your choice. If you want a bleeding edge version, just download the zipfile [here](https://github.com/davidson16807/tectonics.js/archive/master.zip)

## Present-day Earth
This fork opens on the real, present-day Earth instead of a random planet:

* **Elevation** comes from ETOPO1 (NOAA), averaged over each grid cell. The crust is built so the simulated surface
  reproduces those elevations exactly, with continental crust (from the Seton et al. 2020 seafloor age grid) and
  oceanic crust whose depth implies a realistic age.
* **Plates** are the 52 plates of Bird's PB2002 model, each rotating about its real present-day Euler pole
  (converted to the no-net-rotation frame). The smallest microplates appear once the grid is fine enough
  (49 at the default resolution, 51 at the highest).
* The real plates keep their present-day motions until the regular supercontinent cycle resegments the lithosphere
  (150 million years).

Use **Tune → Earth** to return to present-day Earth, **Tune → Randomize** for a random planet, or open
`index.html?world=random` to start on one. Time can now run at up to **10 million years per second**.

Data sources, credits, and the scripts that build `noncompiled/generators/EarthData.js` are in
[tools/earth-data](tools/earth-data/README.md).
