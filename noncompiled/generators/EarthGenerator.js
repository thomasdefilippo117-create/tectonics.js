'use strict';

// EarthGenerator turns a world into present-day Earth using the datasets in EarthData.js:
//  * real elevations (ETOPO1) become rock columns whose isostatic height reproduces them exactly
//  * the 52 real plates (PB2002) become plates that rotate about their real present-day Euler poles
// It works on any grid, so the preset can be rebuilt at every resolution.
const EarthGenerator = (function() {

const EarthGenerator = {};

// the model's mafic density rises linearly over this many seconds, see Crust.get_thickness()
const MAFIC_LIFETIME = 250 * Units.MEGAYEAR;
// fully continental columns carry no mafic rock; CrustGenerator uses the same nominal age for land
const CONTINENTAL_AGE = 1000 * Units.MEGAYEAR;

let decoded = void 0;

function decode_base64(text) {
    const binary = atob(text);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0, li = binary.length; i < li; i++) {
        bytes[i] = binary.charCodeAt(i);
    }
    return bytes;
}
// inverse of rle_uint8() in tools/earth-data/export_js.py: (value, LEB128 run length) pairs
function decode_rle(bytes, length) {
    const result = new Uint8Array(length);
    let i = 0;
    let pos = 0;
    while (pos < bytes.length) {
        const value = bytes[pos++];
        let run = 0;
        let scale = 1;
        let byte = 0;
        do {
            byte = bytes[pos++];
            run += (byte & 0x7F) * scale;
            scale *= 128;
        } while (byte & 0x80);
        if (i + run > length) { throw 'EarthData: run-length data overflows the raster'; }
        result.fill(value, i, i + run);
        i += run;
    }
    if (i !== length) { throw 'EarthData: run-length data does not fill the raster'; }
    return result;
}
function get_decoded() {
    if (decoded !== void 0) { return decoded; }
    const width = EarthData.width;
    const height = EarthData.height;
    const length = width * height;

    const elevation_bytes = decode_base64(EarthData.elevation);
    if (elevation_bytes.length !== 2 * length) { throw 'EarthData: elevation raster has the wrong size'; }
    const elevation_view = new DataView(elevation_bytes.buffer);
    const elevation = new Float32Array(length);
    for (let i = 0; i < length; i++) {
        elevation[i] = elevation_view.getInt16(2 * i, true);
    }

    const continental_bytes = decode_rle(decode_base64(EarthData.continental), length);
    const continental = new Float32Array(length);
    for (let i = 0; i < length; i++) {
        continental[i] = continental_bytes[i] / 255;
    }

    const plate = decode_rle(decode_base64(EarthData.plate), length);
    for (let i = 0; i < length; i++) {
        if (plate[i] >= EarthData.plates.length) { throw 'EarthData: plate raster references an unknown plate'; }
    }

    decoded = { width: width, height: height, elevation: elevation, continental: continental, plate: plate };
    return decoded;
}

function smoothstep(edge0, edge1, x) {
    const t = Math.min(Math.max((x - edge0) / (edge1 - edge0), 0), 1);
    return t * t * (3 - 2 * t);
}
function lerp_clamped(x, x0, x1, y0, y1) {
    const t = Math.min(Math.max((x - x0) / (x1 - x0), 0), 1);
    return y0 + (y1 - y0) * t;
}

// Resamples the equirectangular rasters onto a grid.
// Each cell takes the area-weighted mean of the pixels whose centres lie closest to it (its Voronoi cell),
// so elevations are true cell averages instead of point samples, and each cell joins the plate covering most of it.
EarthGenerator.sample = function(grid) {
    const data = get_decoded();
    const width = data.width;
    const height = data.height;
    const length = width * height;
    const plate_count = EarthData.plates.length;
    const cell_count = grid.vertices.length;
    const Δlat = Math.PI / height;
    const Δlon = 2 * Math.PI / width;

    // pixel centres in the model's coordinate system (see SphericalGeometry.spherical_to_cartesian)
    const pos = VectorRaster.OfLength(length, undefined);
    const area = new Float32Array(length);
    for (let row = 0; row < height; row++) {
        const lat = Math.PI / 2 - (row + 0.5) * Δlat;
        const band_area = Δlon * (Math.sin(Math.PI / 2 - row * Δlat) - Math.sin(Math.PI / 2 - (row + 1) * Δlat));
        const cos_lat = Math.cos(lat);
        const sin_lat = Math.sin(lat);
        for (let col = 0; col < width; col++) {
            const lon = -Math.PI + (col + 0.5) * Δlon;
            const i = row * width + col;
            pos.x[i] =  cos_lat * Math.cos(lon);
            pos.y[i] =  sin_lat;
            pos.z[i] = -cos_lat * Math.sin(lon);
            area[i] = band_area;
        }
    }
    const nearest = new Uint16Array(length);
    grid.getNearestIds(pos, nearest);

    const weight = new Float64Array(cell_count);
    const elevation_sum = new Float64Array(cell_count);
    const continental_sum = new Float64Array(cell_count);
    const votes = new Float64Array(cell_count * plate_count);
    for (let i = 0; i < length; i++) {
        const cell = nearest[i];
        const w = area[i];
        weight[cell] += w;
        elevation_sum[cell] += w * data.elevation[i];
        continental_sum[cell] += w * data.continental[i];
        votes[cell * plate_count + data.plate[i]] += w;
    }

    const elevation = Float32Raster(grid);
    const continental = Float32Raster(grid);
    const plate = Uint8Raster(grid);
    const lat = SphericalGeometry.get_latitudes(grid.pos.y, Float32Raster(grid));
    const lon = SphericalGeometry.get_longitudes(grid.pos.x, grid.pos.z, Float32Raster(grid));
    for (let cell = 0; cell < cell_count; cell++) {
        if (weight[cell] > 0) {
            elevation[cell] = elevation_sum[cell] / weight[cell];
            continental[cell] = continental_sum[cell] / weight[cell];
            const offset = cell * plate_count;
            let best = 0;
            for (let p = 1; p < plate_count; p++) {
                if (votes[offset + p] > votes[offset + best]) { best = p; }
            }
            plate[cell] = best;
        } else {
            // no pixel centre fell in this cell (only possible on grids finer than the data): take the nearest pixel
            const row = Math.min(height - 1, Math.max(0, Math.floor((Math.PI / 2 - lat[cell]) / Δlat)));
            const col = Math.min(width - 1, Math.max(0, Math.floor((lon[cell] + Math.PI) / Δlon)));
            const i = row * width + col;
            elevation[cell] = data.elevation[i];
            continental[cell] = data.continental[i];
            plate[cell] = data.plate[i];
        }
    }
    return { elevation: elevation, continental: continental, plate: plate };
};

// Writes rock columns into `crust` whose isostatic displacement reproduces `elevation` (metres above sea level).
// The model derives height from thickness and density (FluidMechanics.get_isostatic_displacements), so per cell:
//  * continental crust: felsic rock (85% plutonic, 15% volcanic, as in CrustGenerator) thick enough for the elevation
//  * oceanic crust: mafic rock whose age sets its density. The effective plate thickness and the sea level datum are
//    fitted to the real depth-age curve of the seafloor (EarthData.depth_age), so depth implies a realistic age.
//  * margins mix the two in proportion to the continental fraction from the seafloor age grid.
// Returns the constants that were used so the caller can configure sea level and new crust to match.
EarthGenerator.get_crust_from_elevations = function(elevation, continental, material_density, crust) {
    const ρm = material_density.mantle;
    const ρ_min = material_density.mafic_volcanic_min;
    const ρ_max = material_density.mafic_volcanic_max;
    // displacement gained per kg/m² of rock: thickness * (1 - ρ/ρm) = mass * (1/ρ - 1/ρm)
    const felsic_factor = 0.85 * (1 / material_density.felsic_plutonic - 1 / ρm)
                        + 0.15 * (1 / material_density.felsic_volcanic - 1 / ρm);
    const sediment_factor = 1 / material_density.sediment - 1 / ρm;
    const mafic_density = age => ρ_min + (ρ_max - ρ_min) * Math.min(Math.max(age / MAFIC_LIFETIME, 0), 1);
    const mafic_factor = age => 1 - mafic_density(age) / ρm;   // displacement gained per metre of mafic rock

    // depth = intercept + slope * age  <=>  displacement = T0 * mafic_factor(age) - S, with S the sea level datum
    const slope = EarthData.depth_age.slope / Units.MEGAYEAR;                  // metres per second
    const dfactor_dage = -(ρ_max - ρ_min) / (MAFIC_LIFETIME * ρm);            // per second
    const oceanic_thickness = slope / dfactor_dage;                           // T0, metres
    const sealevel = oceanic_thickness * mafic_factor(0) - EarthData.depth_age.intercept; // S, metres
    if (!(oceanic_thickness > 0) || !(sealevel > 0)) {
        throw 'EarthGenerator: depth-age fit is incompatible with the material densities';
    }

    for (let i = 0, li = elevation.length; i < li; i++) {
        const e = elevation[i];
        const target = sealevel + e;
        // seafloor deeper than 3-4 km is oceanic even where the age grid has no data (e.g. Caribbean, trenches):
        // the weight fades from 1 above -3000m to 0 below -4000m
        const w = continental[i] * smoothstep(-4000, -3000, e);
        const sediment = lerp_clamped(e, -3200, -1500, 0, 2500 * 5);  // same sediment template as CrustGenerator
        const residual = target - sediment * sediment_factor;
        const felsic_needed = Math.max(residual, 0) / felsic_factor;

        let felsic = 0;
        let mafic = 0;
        let age = CONTINENTAL_AGE;
        if (w >= 0.999 && residual > 0) {
            felsic = felsic_needed;
        } else {
            felsic = w * felsic_needed;
            const r = residual - felsic * felsic_factor;          // displacement left for the mafic rock
            const nominal = Math.max(1 - w, 1e-3) * oceanic_thickness;
            const required_density = ρm * (1 - r / nominal);
            age = MAFIC_LIFETIME * (required_density - ρ_min) / (ρ_max - ρ_min);
            let thickness = nominal;
            if (age < 0) {                     // shallower than a ridge: young crust, thicker than nominal
                age = 0;
                thickness = r / mafic_factor(0);
            } else if (age > MAFIC_LIFETIME) { // deeper than the oldest seafloor: old crust, thicker than nominal
                age = MAFIC_LIFETIME;
                thickness = r / mafic_factor(MAFIC_LIFETIME);
            }
            mafic = thickness * mafic_density(age);
        }

        crust.sediment[i]        = sediment;
        crust.sedimentary[i]     = 0;
        crust.metamorphic[i]     = 0;
        crust.felsic_plutonic[i] = 0.85 * felsic;
        crust.felsic_volcanic[i] = 0.15 * felsic;
        crust.mafic_volcanic[i]  = mafic;
        crust.mafic_plutonic[i]  = 0;
        crust.age[i]             = age;
    }
    return { sealevel: sealevel, oceanic_thickness: oceanic_thickness };
};

// Turns `world` into present-day Earth: crust, ocean volume, plates and their motions.
// Call before the world is initialized; afterwards call world.lithosphere.mergePlates().
EarthGenerator.generate = function(world) {
    const grid = world.grid;
    const lithosphere = world.lithosphere;
    const crust = lithosphere.total_crust;
    const sample = EarthGenerator.sample(grid);
    const fit = EarthGenerator.get_crust_from_elevations(sample.elevation, sample.continental, world.material_density, crust);

    // ocean volume that makes Hydrology.solve_sealevel() land exactly on the datum used above
    // (the solver averages ocean depth over cells, so we do the same)
    let ocean_depth_sum = 0;
    for (let i = 0, li = sample.elevation.length; i < li; i++) {
        ocean_depth_sum += Math.max(-sample.elevation[i], 0);
    }
    world.hydrosphere.average_ocean_depth = ocean_depth_sum / sample.elevation.length;

    // seafloor created at rifts gets the fitted effective thickness at age 0, so new ridges stand at ridge depth
    lithosphere.rifting_crust = new RockColumn({
        mafic_volcanic: fit.oceanic_thickness * world.material_density.mafic_volcanic_min,
    });

    // one model plate per PB2002 plate that covers at least one cell at this resolution
    const plates = EarthData.plates;
    const degrees = Math.PI / 180;
    const mask = Uint8Raster(grid);
    lithosphere.plates = [];
    for (let p = 0; p < plates.length; p++) {
        Uint8Field.eq_scalar(sample.plate, p, mask);
        if (Uint8Dataset.sum(mask) === 0) { continue; }
        const pole = SphericalGeometry.spherical_to_cartesian(plates[p].pole_lat * degrees, plates[p].pole_lon * degrees);
        const rate = plates[p].rate * degrees / Units.MEGAYEAR;   // radians per second
        lithosphere.plates.push(new Plate(grid, {
            mask:             mask.buffer.slice(0),
            crust:            crust.buffer.slice(0),
            angular_velocity: [pole.x * rate, pole.y * rate, pole.z * rate],
            code:             plates[p].code,
            name:             plates[p].name,
        }));
    }

    // keep the real plates until the regular supercontinent cycle resegments the lithosphere
    lithosphere.supercontinentCycle.age = 0;

    return {
        plate_count: lithosphere.plates.length,
        total_plate_count: plates.length,
        sealevel: fit.sealevel,
        oceanic_thickness: fit.oceanic_thickness,
    };
};

return EarthGenerator;
})();
