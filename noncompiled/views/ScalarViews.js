'use strict';

var scalarViews = {};
scalarViews.satellite = new RealisticWorldView();

scalarViews.npp     = new ScalarWorldView( 
        new ColorscaleRasterView( { min_color: 0xffffff, max_color: 0x00ff00, min: 0., max: 1. }),
        world => world.biosphere.npp.value()
    );
scalarViews.alt     = new ScalarWorldView( 
        new ColorscaleRasterView( { min_color: 0x000000, max_color: 0xffffff, scaling: true }),
        (world, result, scratch, options) => (options.ocean_visibility > 0.5? world.lithosphere.surface_height.value() : world.hydrosphere.elevation.value())
    );
scalarViews.surface_normal_map     = new ScalarWorldView( 
        new SurfaceNormalMapRasterView({}),
        (world, result, scratch, options) => (options.ocean_visibility > 0.5? world.lithosphere.surface_height.value() : world.lithosphere.displacement.value())
    );
scalarViews.topographic    = new ScalarWorldView( 
        new TopographicRasterView( { scaling: true }),
        world => world.hydrosphere.elevation.value()
    );
scalarViews.plates     = new ScalarWorldView(
        new HeatmapRasterView( { min: 0., max: 1. }),
        (world, result) => {
            // A color per plate index cannot keep dozens of plates (Earth has 52) apart, so plates are colored
            // like a political map: each takes the first palette color none of its neighbors use (greedy coloring,
            // in plate order so colors stay put from frame to frame).
            const palette = [0.15, 0.55, 0.85, 0.35, 1.0, 0.7, 0.25, 0.45, 0.95, 0.6];
            const ids = world.lithosphere.top_plate_map;
            const neighbor_lookup = world.grid.neighbor_lookup;
            const plate_count = world.lithosphere.plates.length;
            const neighbors = [];
            for (let p = 0; p < plate_count; p++) { neighbors.push(new Set()); }
            for (let i = 0, li = ids.length; i < li; i++) {
                const a = ids[i];
                if (a >= plate_count) { continue; }
                const cell_neighbors = neighbor_lookup[i];
                for (let j = 0, lj = cell_neighbors.length; j < lj; j++) {
                    const b = ids[cell_neighbors[j]];
                    if (b !== a && b < plate_count) { neighbors[a].add(b); }
                }
            }
            const color = new Int8Array(plate_count).fill(-1);
            for (let p = 0; p < plate_count; p++) {
                const taken = new Set();
                neighbors[p].forEach(q => { if (color[q] >= 0) { taken.add(color[q]); } });
                let c = 0;
                while (taken.has(c) && c < palette.length - 1) { c++; }
                color[p] = c;
            }
            for (let i = 0, li = ids.length; i < li; i++) {
                result[i] = ids[i] < plate_count? palette[color[ids[i]]] : 0;
            }
            return result;
        }
    );
scalarViews.plate_count     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 3. }),
        world => world.lithosphere.plate_count
    );
scalarViews.temp     = new ScalarWorldView( 
        new HeatmapRasterView( { min: -20., max: 30. }),
        (world, result) => ScalarField.sub_scalar(world.atmosphere.surface_temperature, 273.15, result)
    );
scalarViews.precipitation     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 2000., max: 1. }),
        world => world.atmosphere.precipitation.value()
    );
scalarViews.age     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 250., max: 0. }),
        (world, result) => ScalarField.div_scalar(world.lithosphere.top_crust.age, Units.MEGAYEAR, result)
    );
scalarViews.mafic_volcanic     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 7000. }),
        world => world.lithosphere.top_crust.mafic_volcanic
    );
scalarViews.felsic_plutonic     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 70000. }),
        world => world.lithosphere.top_crust.felsic_plutonic
    );
scalarViews.felsic_plutonic_erosion     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 100. }),
        world => world.lithosphere.top_crust.felsic_plutonic_erosion
    );
scalarViews.sediment     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 5. }),
        world => world.lithosphere.top_crust.sediment
    );
scalarViews.sediment_erosion     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 5. }),
        world => world.lithosphere.top_crust.sediment_erosion
    );
scalarViews.sediment_weathering     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 5. }),
        world => world.lithosphere.top_crust.sediment_weathering
    );
scalarViews.sedimentary     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 10000. }),
        world => world.lithosphere.top_crust.sedimentary
    );
scalarViews.metamorphic     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 10000. }),
        world => world.lithosphere.top_crust.metamorphic
    );
scalarViews.thickness     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 70000. }),
        world => world.lithosphere.thickness.value()
    );
scalarViews.density     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 2500, max: 3300 }),
        world => world.lithosphere.density.value()
    );
scalarViews.elevation     = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 10000., scaling: true }),
        function(world, result, scratch, options){
            return options
                world.lithosphere.elevation.value();
        }
    );
scalarViews.snow_coverage = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 1. }),
        world => world.hydrosphere.snow_coverage.value()
    );
scalarViews.land_coverage = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 1. }),
        world => world.lithosphere.land_coverage.value()
    );
scalarViews.ocean_coverage = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 1. }),
        world => world.hydrosphere.ocean_coverage.value()
    );
scalarViews.plant_coverage = new ScalarWorldView( 
        new HeatmapRasterView( { min: 0., max: 1. }),
        world => world.biosphere.plant_coverage.value()
    );
scalarViews.asthenosphere_pressure = new ScalarWorldView( 
        new HeatmapRasterView(  { min: 1., max: 0.  }),
        function (world, output, scratch, iterations) {return TectonicsModeling.get_asthenosphere_pressure(world.lithosphere.density, output, scratch);}
    );

scalarViews.surface_air_pressure = new ScalarWorldView( 
        new HeatmapRasterView( { min: 980000., max: 1030000. }),
        world => world.atmosphere.surface_pressure.value()
    );
