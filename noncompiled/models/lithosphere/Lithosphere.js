'use strict';


function Lithosphere(grid, parameters) {
    parameters = parameters || {}
    grid = grid || stop('missing parameter: "grid"');
    // NOTE: the cycle's state is saved under "supercontinent_cycle"; it was read from the wrong place before,
    //  so every loaded world reset its plates on the first update
    this.supercontinentCycle = new SupercontinentCycle(this, parameters['supercontinent_cycle'] || {});
    this.plates = (parameters['plates'] || []).map(plate_parameters => new Plate(grid, plate_parameters));

    // "rifting_crust" is the rock column given to new crust that forms where plates pull apart
    this.rifting_crust = new RockColumn(parameters['rifting_crust'] || {
        mafic_volcanic:         7100,     // +/- 800, White McKenzie and O'nions 1992
        // felsic_plutonic:     100, // This can be set above zero to "cheat" on felsic mass conservation
     });

    this.getParameters = function() {
        const rifting_crust = this.rifting_crust;
        return {
            //grid:                 grid. // TODO: add grid
            plates:                 this.plates.map(plate => plate.getParameters()),
            total_crust:             this.total_crust.buffer.slice(0),
            supercontinent_cycle:     this.supercontinentCycle.getParameters(),
            rifting_crust: {
                sediment:         rifting_crust.sediment,
                sedimentary:      rifting_crust.sedimentary,
                metamorphic:      rifting_crust.metamorphic,
                felsic_plutonic:  rifting_crust.felsic_plutonic,
                felsic_volcanic:  rifting_crust.felsic_volcanic,
                mafic_volcanic:   rifting_crust.mafic_volcanic,
                mafic_plutonic:   rifting_crust.mafic_plutonic,
                age:              rifting_crust.age,
            },
        };
    }

    let material_viscosity = undefined;
    let material_density = undefined;
    let sealevel = undefined;
    let surface_gravity = undefined;


    // The following are fields that are derived from other fields:
    // "displacement is the height of the crust relative to an arbitrary datum level
    // It is not called "elevation" because we want to emphasize that it is not relative to sea level
    const self = this; 
    // "surface_height" is the height of the surface relative to sealevel - if elevation < 0, then surface_height = 0
    this.surface_height = new Memo(
        Float32Raster(grid),  
        result => Hydrology.get_surface_heights(this.displacement.value(), sealevel.value(), result)
    ); 
    this.displacement = new Memo(  
        Float32Raster(grid),  
        result => FluidMechanics.get_isostatic_displacements(self.thickness.value(), self.density.value(), material_density, result) 
    ); 
    // the thickness of the crust in km
    this.thickness = new Memo(  
        Float32Raster(grid),  
        result => Crust.get_thickness(self.total_crust, material_density, result)
    ); 
    // total mass of the crust in kg
    this.total_mass = new Memo(  
        Float32Raster(grid),  
        result => Crust.get_total_mass(self.total_crust, result)
    ); 
    // the average density of the crust, in kg/m^3
    this.density = new Memo(
        Float32Raster(grid),  
        result => Crust.get_density(self.total_mass.value(), self.thickness.value(),    material_density.mafic_volcanic_min, result)
    ); 
    this.buoyancy = new Memo(
        Float32Raster(grid),  
        result => Crust.get_buoyancy(self.density.value(), material_density, surface_gravity, result)
    ); 

    this.top_plate_map             = Uint8Raster(grid);
    this.plate_count         = Uint8Raster(grid);
    this.asthenosphere_velocity = VectorRaster(grid);

    this.total_crust     = new Crust({grid: grid, buffer: parameters['total_crust']});
    this.top_crust         = new Crust({grid: grid});
    this.erosion         = new Crust({grid: grid});
    this.weathering     = new Crust({grid: grid});
    this.lithification     = new Crust({grid: grid});
    this.metamorphosis     = new Crust({grid: grid});
    this.accretion         = new Crust({grid: grid});
    this.crust_delta     = new Crust({grid: grid});
    this.crust_scratch     = new Crust({grid: grid});

    this.average_conserved_per_cell = 0;





    function move_plates(plates, seconds) {
        for (let i=0, li=plates.length; i<li; ++i) {
             plates[i].move(seconds, material_density, material_viscosity, surface_gravity);
         }
    }

    // Scratch space reused every frame by the functions below.
    // The plate loops only do real work on the cells a plate covers and on its borders, instead of running
    // several whole-world raster passes per plate. With dozens of plates (e.g. Earth's 52) that was the main cost.
    // Every function reproduces the original raster-based implementation exactly, including float32 rounding,
    // except for the gap filling at the end of update_rifting(), which is new.
    const cell_count = grid.vertices.length;
    const neighbor_lookup = grid.neighbor_lookup;
    const master_density = Float32Raster(grid);
    const visit_stamps = new Uint32Array(cell_count);   // cells already examined in the current pass
    let visit_stamp = 0;
    const claim_stamps = new Uint32Array(cell_count);   // world cells that received new crust this frame
    let claim_stamp = 0;
    const is_subducted_scratch = new Uint8Array(cell_count);
    const accretion_scratch = new Float32Array(cell_count);
    const selected_cells = new Uint32Array(cell_count);
    function next_visit_stamp() {
        visit_stamp++;
        if (visit_stamp === 0xFFFFFFFF) { visit_stamps.fill(0); visit_stamp = 1; }
        return visit_stamp;
    }

    // stitches the plates into the world maps: total_crust, top_crust, top_plate_map, plate_count
    // A plate is on top of a cell if it is the least dense plate there (ties go to the earlier plate).
    function merge_plates_to_master(plates, master) {
        const UINT8_NULL = 255;

        //WIPE MASTER RASTERS CLEAN
        Crust.reset(master.total_crust);
        Crust.reset(master.top_crust);
        Uint8Raster.fill(master.top_plate_map, UINT8_NULL);
        Uint8Raster.fill(master.plate_count, 0);
        Float32Raster.fill(master_density, 9999);

        const top_plate_map = master.top_plate_map;
        const plate_count = master.plate_count;
        const total_conserved = master.total_crust.conserved_pools;
        const total_mafic_volcanic = master.total_crust.mafic_volcanic;
        const total_mafic_plutonic = master.total_crust.mafic_plutonic;
        const total_age = master.total_crust.age;
        const top_pools = master.top_crust.all_pools;
        const conserved_pool_count = total_conserved.length;
        const pool_count = top_pools.length;

        for (let i=0, li=plates.length; i<li; ++i) {
            const plate = plates[i];
            const local_ids_of_global_cells = plate.local_ids_of_global_cells;
            const mask = plate.mask;
            const density = plate.density.value();
            const conserved = plate.crust.conserved_pools;
            const pools = plate.crust.all_pools;
            const mafic_volcanic = plate.crust.mafic_volcanic;
            const mafic_plutonic = plate.crust.mafic_plutonic;
            const age = plate.crust.age;

            for (let g = 0; g < cell_count; g++) {
                const l = local_ids_of_global_cells[g];
                if (mask[l] !== 1) { continue; }

                // add 1 to plate_count where current plate exists
                plate_count[g] += 1;
                // conserved rock from every plate that is present adds up
                for (let k = 0; k < conserved_pool_count; k++) {
                    total_conserved[k][g] += conserved[k][l];
                }
                // the plate seen from space supplies the rest
                if (density[l] < master_density[g]) {
                    master_density[g] = density[l];
                    top_plate_map[g] = i;
                    total_mafic_volcanic[g] = mafic_volcanic[l];
                    total_mafic_plutonic[g] = mafic_plutonic[l];
                    total_age[g] = age[l];
                    for (let k = 0; k < pool_count; k++) {
                        top_pools[k][g] = pools[k][l];
                    }
                }
            }
        }
    }

    // adds new crust just outside a plate's border wherever that region and all its neighbors are either
    // empty or covered only by this plate, i.e. where plates are pulling apart
    function update_rifting(world, plates) {
        const top_plate_map = world.top_plate_map;
        const plate_count = world.plate_count;
        const rifting_pools = world.rifting_crust.all_pools;
        const pool_count = rifting_pools.length;

        claim_stamp++;
        if (claim_stamp === 0xFFFFFFFF) { claim_stamps.fill(0); claim_stamp = 1; }

        for (let i=0, li=plates.length; i<li; ++i) {
            const plate = plates[i];
            const mask = plate.mask;
            const global_ids_of_local_cells = plate.global_ids_of_local_cells;
            const pools = plate.crust.all_pools;
            // is_riftable: count == 0 or (count = 1 and top_plate = i)
            const is_riftable = function(l) {
                const g = global_ids_of_local_cells[l];
                const count = plate_count[g];
                return count === 0 || (count === 1 && top_plate_map[g] === i);
            };

            // candidates are the cells just outside the border: not in the mask, but next to it
            const stamp = next_visit_stamp();
            let rifting_count = 0;
            for (let l = 0; l < cell_count; l++) {
                if (mask[l] !== 1) { continue; }
                const neighbors = neighbor_lookup[l];
                for (let j = 0, lj = neighbors.length; j < lj; j++) {
                    const n = neighbors[j];
                    if (mask[n] !== 0 || visit_stamps[n] === stamp) { continue; }
                    visit_stamps[n] = stamp;
                    // the cell will stay riftable: it and all of its neighbors are riftable
                    if (!is_riftable(n)) { continue; }
                    const neighbors_n = neighbor_lookup[n];
                    let will_stay_riftable = true;
                    for (let m = 0, lm = neighbors_n.length; m < lm; m++) {
                        if (!is_riftable(neighbors_n[m])) { will_stay_riftable = false; break; }
                    }
                    if (will_stay_riftable) { selected_cells[rifting_count++] = n; }
                }
            }
            // apply only after every candidate was judged against the unmodified mask
            for (let r = 0; r < rifting_count; r++) {
                const l = selected_cells[r];
                mask[l] = 1;
                for (let k = 0; k < pool_count; k++) {
                    pools[k][l] = rifting_pools[k];
                }
                claim_stamps[global_ids_of_local_cells[l]] = claim_stamp;
            }
        }

        // Gap filling: the rule above only rifts where a cell and all of its neighbors are free, so the one or two
        // cell wide gaps that open between diverging plates every frame stayed empty for good. Empty cells have no
        // crust and sat deep below the seafloor, so spreading ridges showed up as dashed trenches.
        // Each gap cell next to a plate now joins the plate that covers most of its neighbors (ties go to the earlier
        // plate), the way new seafloor accretes to the plates on either side of a ridge. One plate per cell, so no
        // cell is claimed twice.
        const plate_total = plates.length;
        for (let g = 0; g < cell_count; g++) {
            if (plate_count[g] !== 0 || claim_stamps[g] === claim_stamp) { continue; }
            const neighbors = neighbor_lookup[g];
            let best = -1;
            let best_votes = 0;
            for (let j = 0, lj = neighbors.length; j < lj; j++) {
                const p = top_plate_map[neighbors[j]];
                if (p >= plate_total) { continue; }
                let votes = 0;
                for (let m = 0; m < lj; m++) {
                    if (top_plate_map[neighbors[m]] === p) { votes++; }
                }
                if (votes > best_votes || (votes === best_votes && p < best)) {
                    best = p;
                    best_votes = votes;
                }
            }
            if (best < 0) { continue; }
            const plate = plates[best];
            const l = plate.local_ids_of_global_cells[g];
            if (plate.mask[l] === 1) { continue; }
            plate.mask[l] = 1;
            const pools = plate.crust.all_pools;
            for (let k = 0; k < pool_count; k++) {
                pools[k][l] = rifting_pools[k];
            }
            claim_stamps[g] = claim_stamp;
        }
    }

    // finds crust that lies under another plate: metamorphoses its felsic rock, and detaches dense crust
    // just inside the plate's border, handing its felsic mass to the overriding plate as accretion
    function update_subducted(world, plates) {
        const top_plate_map = world.top_plate_map;
        const plate_count = world.plate_count;
        const mantle_density = material_density.mantle;

        // WARNING: unfortunate side effect!
        // we calculate accretion delta during detachment for performance reasons
        Crust.reset(world.accretion);
        const accretion_felsic_plutonic = world.accretion.felsic_plutonic;
        const accretion_felsic_volcanic = world.accretion.felsic_volcanic;

        const fround = Math.fround;
        const is_subducted = is_subducted_scratch;

        for (let i=0, li=plates.length; i<li; ++i) {
            const plate = plates[i];
            const mask = plate.mask;
            const global_ids_of_local_cells = plate.global_ids_of_local_cells;
            const local_ids_of_global_cells = plate.local_ids_of_global_cells;
            const crust = plate.crust;
            const sediment = crust.sediment;
            const sedimentary = crust.sedimentary;
            const metamorphic = crust.metamorphic;
            const felsic_plutonic = crust.felsic_plutonic;
            const felsic_volcanic = crust.felsic_volcanic;

            // detect subduction: is_subducted: count > 1 and top_plate != i
            // WARNING: unfortunate side effect!
            // we metamorphose stuff here because trying to do it using deltas was causing problems with conservation
            // (every local cell is processed, as before, since cells outside the mask still feed plate velocity estimates)
            for (let l = 0; l < cell_count; l++) {
                const g = global_ids_of_local_cells[l];
                const subducted = plate_count[g] !== 1 && top_plate_map[g] !== i;
                is_subducted[l] = subducted? 1 : 0;
                if (subducted) {
                    metamorphic[l] = metamorphic[l] + sediment[l];
                    metamorphic[l] = metamorphic[l] + sedimentary[l];
                    metamorphic[l] = metamorphic[l] + felsic_plutonic[l];
                    metamorphic[l] = metamorphic[l] + felsic_volcanic[l];
                    sediment[l] = 0;
                    sedimentary[l] = 0;
                    felsic_plutonic[l] = 0;
                    felsic_volcanic[l] = 0;
                }
            }

            // detaching cells: just inside the border, subducted along with all their neighbors, and denser than mantle
            // NOTE: density is the value memoized before metamorphosis, exactly as before
            const density = plate.density.value();
            let detaching_count = 0;
            for (let l = 0; l < cell_count; l++) {
                if (mask[l] !== 1 || is_subducted[l] !== 1 || !(density[l] > mantle_density)) { continue; }
                const neighbors = neighbor_lookup[l];
                let is_just_inside_border = false;
                let will_stay_detachable = true;
                for (let j = 0, lj = neighbors.length; j < lj; j++) {
                    const n = neighbors[j];
                    if (mask[n] !== 1) { is_just_inside_border = true; }
                    if (is_subducted[n] !== 1) { will_stay_detachable = false; break; }
                }
                if (is_just_inside_border && will_stay_detachable) { selected_cells[detaching_count++] = l; }
            }
            if (detaching_count === 0) { continue; }

            // calculate accretion delta: the conserved mass of each detaching cell
            const stamp = next_visit_stamp();
            for (let r = 0; r < detaching_count; r++) {
                const l = selected_cells[r];
                mask[l] = 0;
                let mass = fround(0 + sediment[l]);
                mass = fround(mass + sedimentary[l]);
                mass = fround(mass + metamorphic[l]);
                mass = fround(mass + felsic_plutonic[l]);
                mass = fround(mass + felsic_volcanic[l]);
                accretion_scratch[l] = mass;
                visit_stamps[l] = stamp;
            }
            // every world cell that maps onto a detaching cell receives its mass
            for (let g = 0; g < cell_count; g++) {
                const l = local_ids_of_global_cells[g];
                if (visit_stamps[l] !== stamp) { continue; }
                accretion_felsic_plutonic[g] = accretion_felsic_plutonic[g] + accretion_scratch[l] * 0.85;
                accretion_felsic_volcanic[g] = accretion_felsic_volcanic[g] + accretion_scratch[l] * 0.15;
            }
        }
    }
    function calculate_deltas(lithosphere, seconds) {

           // CALCULATE DELTAS
        Crust.model_erosion(
            self.surface_height.value(), seconds,
            material_density, surface_gravity,
            lithosphere.top_crust, lithosphere.erosion, lithosphere.crust_scratch
        );

           // CALCULATE DELTAS
        Crust.model_weathering(
            self.surface_height.value(), seconds,
            material_density, surface_gravity,
            lithosphere.top_crust, lithosphere.weathering, lithosphere.crust_scratch
        );

           // CALCULATE DELTAS
        Crust.model_lithification(
            self.surface_height.value(), seconds,
            material_density, surface_gravity,
            lithosphere.top_crust, lithosphere.lithification, lithosphere.crust_scratch
        );

           // CALCULATE DELTAS
        Crust.model_metamorphosis(
            self.surface_height.value(), seconds,
            material_density, surface_gravity,
            lithosphere.top_crust, lithosphere.metamorphosis, lithosphere.crust_scratch
        );

        // COMPILE DELTAS
        const globalized_deltas = lithosphere.crust_delta;
        Crust.reset             (globalized_deltas);
        Crust.add_delta         (globalized_deltas, lithosphere.erosion,                 globalized_deltas);
        Crust.add_delta         (globalized_deltas, lithosphere.weathering,             globalized_deltas);
        Crust.add_delta         (globalized_deltas, lithosphere.lithification,            globalized_deltas);
        Crust.add_delta         (globalized_deltas, lithosphere.metamorphosis,            globalized_deltas);
        Crust.add_delta         (globalized_deltas, lithosphere.accretion,                globalized_deltas);
    }

    function integrate_deltas(world, plates, seconds) {
        // INTEGRATE DELTAS
        Crust.add_delta(world.total_crust, world.crust_delta, world.total_crust);

        // Each world cell's delta goes to the local cell of the plate on top there (adding together effects when a
        // local cell is mentioned twice). This method retains positive mass, and appears to give the best results
        // of any method attempted so far. One pass over the world does what used to take one pass per plate.
        const top_plate_map = world.top_plate_map;
        const delta_pools = world.crust_delta.mass_pools;
        const pool_count = delta_pools.length;
        const plate_count = plates.length;
        for (let g = 0; g < cell_count; g++) {
            const i = top_plate_map[g];
            if (i >= plate_count) { continue; } // no plate on top (255), or a stale index
            const plate = plates[i];
            const l = plate.local_ids_of_global_cells[g];
            const pools = plate.crust.mass_pools;
            for (let k = 0; k < pool_count; k++) {
                pools[k][l] += delta_pools[k][g];
            }
        }
        // The only problem with the aforementioned method is it applies the delta in a "patchwork" manner.
        // This isn't a big problem for most fields, but it is very noticeable for the "age" and "density" field,
        // which ought to have very smooth gradients.
        // These fields are very important since they are used to derive the motions and borders of the plates.
        // To get around this, we give the age field special treatment, and simply increment its values everywhere by the same timestep.
        for (let i = 0; i < plate_count; i++) {
            ScalarField.add_scalar(plates[i].crust.age, seconds, plates[i].crust.age);
        }
    }

    this.resetPlates = function() {
        // get plate masks from image segmentation of asthenosphere velocity
        const pressure = FluidMechanics.get_fluid_pressures(this.buoyancy.value());
        FluidMechanics.get_fluid_velocities(pressure, this.asthenosphere_velocity);
        const angular_velocity = VectorField.cross_vector_field(this.asthenosphere_velocity, grid.pos);
        const top_plate_map = Tectonophysics.guess_plate_map(angular_velocity, 7, 200);
        const plate_ids = Uint8Dataset.unique(top_plate_map);
        this.plates = [];

        let plate, mask;
        // TODO: overwrite plates instead of creating new ones, create separate function for plate initialization
        for (let i = 0, li = plate_ids.length; i < li; ++i) {
            mask = Uint8Field.eq_scalar(top_plate_map, plate_ids[i]);

            plate = new Plate(
                grid,
                {
                    mask:     mask.buffer,
                    crust:  this.total_crust.buffer.slice(0)
                },
            )
            plate.setDependencies({
                material_density:     material_density,
                material_viscosity: material_viscosity,
                surface_gravity: surface_gravity,
            });
            this.plates.push(plate);
        }
    };

    // stitches the plates into the world maps (total_crust, top_plate_map, etc.) without advancing time
    // use this after plates are created outside of the simulation loop so the maps are valid immediately
    this.mergePlates = function() {
        assert_dependencies();
        this.invalidate();
        merge_plates_to_master(this.plates, this);
        this.invalidate();
    };

    function assert_dependencies() {
        if (sealevel === void 0)             { throw '"sealevel" not provided'; }
        if (surface_gravity === void 0)    { throw '"surface_gravity" not provided'; }
        if (material_density === void 0)    { throw '"material_density" not provided'; }
        if (material_viscosity === void 0)    { throw '"material_viscosity" not provided'; }
    }

    this.setDependencies = function(dependencies) {
        for (let i = 0; i < this.plates.length; i++) {
            this.plates[i].setDependencies(dependencies);
        }

        sealevel             = dependencies['sealevel']             !== void 0?     dependencies['sealevel']                 : sealevel;
        surface_gravity     = dependencies['surface_gravity']     !== void 0?     dependencies['surface_gravity']         : surface_gravity;
        material_density     = dependencies['material_density']     !== void 0?     dependencies['material_density']         : material_density;
        material_viscosity     = dependencies['material_viscosity']!== void 0?     dependencies['material_viscosity']         : material_viscosity;
    };

    this.initialize = function() { 
        assert_dependencies();

        this.average_conserved_per_cell = Crust.get_average_conserved_per_cell(this.total_crust);
    }

    this.invalidate = function() {
        this.displacement.invalidate();
        this.surface_height.invalidate();
        this.thickness    .invalidate();
        this.total_mass    .invalidate();
        this.density    .invalidate();
        this.buoyancy    .invalidate();
        const plates = this.plates;
        for (let i=0, li=plates.length; i<li; ++i) {
            plates[i]    .invalidate();
        }
    }

    const mean_supercontinent_cycle_duration = 150 * Units.MEGAYEAR;

    this.calcChanges = function(seconds) {
        const max_perceivable_duration = 60*60*24*30 * seconds; // 1 day worth of real time at 30fps
        if (mean_supercontinent_cycle_duration > max_perceivable_duration) {
            return;
        }
        
        assert_dependencies();

        calculate_deltas        (this, seconds);             // this creates a world map of all additions and subtractions to crust (e.g. from erosion, accretion, etc.)
    };

    this.applyChanges = function(seconds){
        const max_perceivable_duration = 60*60*24*30 * seconds; // 1 day worth of real time at 30fps
        if (mean_supercontinent_cycle_duration > max_perceivable_duration) {
            return;
        }

        assert_dependencies();

        integrate_deltas         (this, this.plates, seconds);// this uses the map above in order to add and subtract crust

        move_plates             (this.plates, seconds);     // this performs the actual plate movement
        this.supercontinentCycle.update(seconds);             // this periodically splits the world into plates
        merge_plates_to_master    (this.plates, this);         // this stitches plates together to create a world map
        update_rifting            (this, this.plates);         // this identifies rifting regions on the world map and adds crust to plates where needed
        update_subducted         (this, this.plates);         // this identifies detaching regions on the world map and then removes crust from plates where needed
    };
}
