'use strict';

function Plate(grid, parameters)
{
    parameters = parameters || stop('missing parameter object')
    grid = grid || stop('missing parameter: "grid"');

    this.crust = new Crust({grid: grid, buffer: parameters['crust']});
    this.mask = Uint8Raster.FromBuffer(parameters['mask'], grid);
    this.local_to_global_matrix = Matrix3x3.Identity();
    if (parameters['local_to_global_matrix'] !== void 0) {
        this.local_to_global_matrix.set(parameters['local_to_global_matrix'])
    }
    this.global_to_local_matrix = Matrix3x3.invert(this.local_to_global_matrix);

    // "angular_velocity" is an optional Euler vector (radians per second, world frame, counter-clockwise positive).
    // When it is set, the plate rotates rigidly about that fixed axis instead of using the slab pull estimate.
    // This is how plates with known present-day motions (e.g. the PB2002 plates of the Earth preset) are driven.
    this.angular_velocity = void 0;
    if (parameters['angular_velocity'] !== void 0) {
        const ω = parameters['angular_velocity'];
        this.angular_velocity = Vector(ω[0], ω[1], ω[2]);
    }
    // optional identifiers for plates that correspond to real ones, e.g. code 'PA', name 'Pacific'
    this.code = parameters['code'];
    this.name = parameters['name'];

    this.getParameters = function() {
        const result = {
            //grid:         grid. // TODO: add grid
            crust:             this.crust.buffer.slice(0),
            mask:             this.mask.buffer.slice(0),
            local_to_global_matrix: Array.from(this.local_to_global_matrix)
        };
        if (this.angular_velocity !== void 0) {
            result.angular_velocity = [this.angular_velocity.x, this.angular_velocity.y, this.angular_velocity.z];
        }
        if (this.code !== void 0) { result.code = this.code; }
        if (this.name !== void 0) { result.name = this.name; }
        return result;
    }

    const self = this;
    const cell_count = grid.vertices.length;

    // "Sparse" plates: a plate with a prescribed motion never uses the slab pull estimate, so the lithosphere only
    // reads its crust on the plate itself and on a margin around it. For such plates, lookups and derived fields
    // are computed for those cells alone, so the cost follows the plate's size instead of the whole planet's:
    //  * active_cells:  local cells of the plate's mask plus two rings around it
    //                   (rifting judges the first ring by its neighbors, the second ring)
    //  * covered_cells: world cells close enough that their nearest local cell could be on the plate
    // Other plates keep a count of -1, meaning every cell is used. Their slab pull estimate reads crust
    // outside the mask, so they must keep processing the whole planet.
    this.active_cells = void 0;
    this.active_count = -1;
    this.covered_cells = void 0;
    this.covered_count = -1;
    this.isSparse = function() {
        return this.angular_velocity !== void 0 && Plate.sparse_prescribed_plates;
    }

    // Rasters are allocated on first use: sparse plates never touch most of them, which saves a lot of memory
    // on fine grids with many plates.
    // The following are fields that are derived from other fields:
    // "displacement is the height of the crust relative to an arbitrary datum level
    // It is not called "elevation" because we want to emphasize that it is not relative to sea level
    this.displacement = new Memo(
        void 0,
        result => FluidMechanics.get_isostatic_displacements(self.thickness.value(), self.density.value(), material_density, result || Float32Raster(grid))
    );
    // the thickness of the crust in km
    this.thickness = new Memo(
        void 0,
        result => Crust.get_thickness(self.crust, material_density, result || Float32Raster(grid))
    );
    // total mass of the crust in tons
    this.total_mass = new Memo(
        void 0,
        result => Crust.get_total_mass(self.crust, result || Float32Raster(grid))
    );
    // the average density of the crust, in kg/m^3
    // (for sparse plates only the active cells are computed; nothing reads the others)
    this.density = new Memo(
        void 0,
        result => self.active_count >= 0?
            Crust.get_density_of_cells(self.crust, material_density, self.active_cells, self.active_count, result || Float32Raster(grid)) :
            Crust.get_density(self.total_mass.value(), self.thickness.value(), material_density.mafic_volcanic_min, result || Float32Raster(grid))
    );
    this.buoyancy = new Memo(
        void 0,
        result => Crust.get_buoyancy(self.density.value(), material_density, surface_gravity, result || Float32Raster(grid))
    );
    this.boundary_normal = new Memo(
        void 0,
        result => {
            result = result || VectorRaster(grid);
            Uint8Field .gradient  (self.mask, result);
            VectorField.normalize (result,          result);
            return result;
        }
    );
    this.velocity = new Memo(
        void 0,
        result => Tectonophysics.guess_plate_velocity(self.boundary_normal.value(), self.buoyancy.value(), material_viscosity, result || VectorRaster(grid))
    );
    this.center_of_mass = new Memo(
        { x:0, y:0, z:0 },
        result => Tectonophysics.get_plate_center_of_mass    (self.total_mass.value(), self.mask)
    );



    this.global_ids_of_local_cells = Uint16Raster(grid);
    this.local_ids_of_global_cells = Uint16Raster(grid);

    // positions used by the full (non-sparse) lookups, allocated when first needed
    this.global_pos_of_local_cells = void 0;
    this.local_pos_of_global_cells = void 0;

    let material_density = undefined;
    let material_viscosity = undefined;
    let surface_gravity = undefined;

    function assert_dependencies() {
        if (material_density === void 0)    { throw '"material_density" not provided'; }
        if (material_viscosity === void 0)    { throw '"material_viscosity" not provided'; }
        if (surface_gravity === void 0)    { throw '"surface_gravity" not provided'; }
    }

    this.setDependencies = function(dependencies) {
        material_density     = dependencies['material_density']     !== void 0?     dependencies['material_density']         : material_density;
        material_viscosity     = dependencies['material_viscosity']!== void 0?     dependencies['material_viscosity']         : material_viscosity;
        surface_gravity     = dependencies['surface_gravity']     !== void 0?     dependencies['surface_gravity']         : surface_gravity;
    };

    this.invalidate = function(){
        this.displacement.invalidate();
        this.thickness.invalidate();
        this.total_mass.invalidate();
        this.density.invalidate();
        this.boundary_normal.invalidate();
        this.buoyancy.invalidate();
        this.velocity.invalidate();
        this.center_of_mass.invalidate();
    }

    // scratch space shared by every plate on the grid (plates update one at a time)
    function get_scratch() {
        let scratch = grid.plate_scratch;
        if (scratch === void 0) {
            scratch = grid.plate_scratch = {
                x:      new Float32Array(cell_count),
                y:      new Float32Array(cell_count),
                z:      new Float32Array(cell_count),
                ids:    new Uint16Array(cell_count),
                stamps: new Uint32Array(cell_count),
                stamp:  0,
            };
        }
        scratch.stamp++;
        if (scratch.stamp === 0xFFFFFFFF) { scratch.stamps.fill(0); scratch.stamp = 1; }
        return scratch;
    }
    // rotates the positions of `count` cells listed in `cells` by `matrix` (same arithmetic as VectorField.mult_matrix,
    // rounded to float32 the same way), finds the nearest grid cell of each, and writes those ids to `ids_out[cell]`
    function look_up_rotated(cells, count, matrix, scratch, ids_out) {
        const pos = grid.pos;
        const px = pos.x, py = pos.y, pz = pos.z;
        const x = scratch.x, y = scratch.y, z = scratch.z;
        const xx = matrix[0], xy = matrix[3], xz = matrix[6];
        const yx = matrix[1], yy = matrix[4], yz = matrix[7];
        const zx = matrix[2], zy = matrix[5], zz = matrix[8];
        for (let c = 0; c < count; c++) {
            const i = cells[c];
            const pxi = px[i], pyi = py[i], pzi = pz[i];
            x[c] = pxi * xx + pyi * xy + pzi * xz;
            y[c] = pxi * yx + pyi * yy + pzi * yz;
            z[c] = pxi * zx + pyi * zy + pzi * zz;
        }
        const ids = scratch.ids.subarray(0, count);
        grid.getNearestIds({ x: x.subarray(0, count), y: y.subarray(0, count), z: z.subarray(0, count) }, ids);
        for (let c = 0; c < count; c++) {
            ids_out[cells[c]] = ids[c];
        }
    }

    this.update_ids = function() {
        if (!this.isSparse()) {
            this.active_count = -1;
            this.covered_count = -1;
            this.local_pos_of_global_cells = this.local_pos_of_global_cells || VectorRaster(grid);
            this.global_pos_of_local_cells = this.global_pos_of_local_cells || VectorRaster(grid);
            // for each cell in the master's grid, this raster indicates the id of the corresponding cell in the plate's grid
            // this is used to convert between global and local coordinate systems
            VectorField.mult_matrix(grid.pos, this.global_to_local_matrix, this.local_pos_of_global_cells);
            grid.getNearestIds(this.local_pos_of_global_cells, this.local_ids_of_global_cells);

            // for each cell in the plate's grid, this raster indicates the id of the corresponding cell in the world's grid
            // this is used to convert between global and local coordinate systems
            VectorField.mult_matrix(grid.pos, this.local_to_global_matrix, this.global_pos_of_local_cells);
            grid.getNearestIds(this.global_pos_of_local_cells, this.global_ids_of_local_cells);
            return;
        }

        const scratch = get_scratch();
        const stamps = scratch.stamps;
        const stamp = scratch.stamp;
        const mask = this.mask;
        const pos = grid.pos;
        const neighbor_lookup = grid.neighbor_lookup;
        if (this.active_cells === void 0) {
            this.active_cells = new Uint32Array(cell_count);
            this.covered_cells = new Uint32Array(cell_count);
        }

        // active cells: the mask, then two rings around it
        const active = this.active_cells;
        let count = 0;
        let center_x = 0, center_y = 0, center_z = 0;
        for (let l = 0; l < cell_count; l++) {
            if (mask[l] !== 1) { continue; }
            active[count++] = l;
            stamps[l] = stamp;
            center_x += pos.x[l];
            center_y += pos.y[l];
            center_z += pos.z[l];
        }
        const mask_count = count;
        let ring_start = 0;
        for (let ring = 0; ring < 2; ring++) {
            const ring_end = count;
            for (let a = ring_start; a < ring_end; a++) {
                const neighbors = neighbor_lookup[active[a]];
                for (let j = 0, lj = neighbors.length; j < lj; j++) {
                    const n = neighbors[j];
                    if (stamps[n] === stamp) { continue; }
                    stamps[n] = stamp;
                    active[count++] = n;
                }
            }
            ring_start = ring_end;
        }
        this.active_count = count;

        // world ids of the active local cells
        look_up_rotated(active, count, this.local_to_global_matrix, scratch, this.global_ids_of_local_cells);

        // covered world cells: within a cap around the plate, padded by three cell spacings. Every world cell whose
        // nearest local cell lies on the plate is within one spacing of the unpadded cap, and the lithosphere also
        // looks up empty world cells next to those (gap filling), which lie within one more.
        if (mask_count === 0) {
            this.covered_count = 0;
            return;
        }
        const center_length = Math.sqrt(center_x*center_x + center_y*center_y + center_z*center_z);
        let cos_radius = 1;
        if (center_length > 1e-9) {
            center_x /= center_length; center_y /= center_length; center_z /= center_length;
            for (let a = 0; a < mask_count; a++) {
                const l = active[a];
                const d = pos.x[l]*center_x + pos.y[l]*center_y + pos.z[l]*center_z;
                if (d < cos_radius) { cos_radius = d; }
            }
        } else {
            cos_radius = -1; // the plate surrounds its own center: treat it as covering everything
        }
        const radius = Math.acos(Math.max(-1, Math.min(1, cos_radius))) + 3 * grid.average_distance;
        const cos_reach = radius >= Math.PI? -2 : Math.cos(radius);
        const M = this.local_to_global_matrix;
        const world_center_x = center_x*M[0] + center_y*M[3] + center_z*M[6];
        const world_center_y = center_x*M[1] + center_y*M[4] + center_z*M[7];
        const world_center_z = center_x*M[2] + center_y*M[5] + center_z*M[8];
        const covered = this.covered_cells;
        let covered_count = 0;
        for (let g = 0; g < cell_count; g++) {
            if (pos.x[g]*world_center_x + pos.y[g]*world_center_y + pos.z[g]*world_center_z >= cos_reach) {
                covered[covered_count++] = g;
            }
        }
        this.covered_count = covered_count;

        // local ids of the covered world cells
        look_up_rotated(covered, covered_count, this.global_to_local_matrix, scratch, this.local_ids_of_global_cells);
    }

    this.move = function(seconds){
        assert_dependencies();

        const ω = this.angular_velocity;
        if (ω !== void 0) {
            // rigid rotation about an axis that is fixed in the world frame, so it is applied on the world side
            const θ = Math.sqrt(ω.x*ω.x + ω.y*ω.y + ω.z*ω.z) * seconds;
            if (θ > 0) {
                const rotation_matrix = Matrix3x3.FromRotationVector(ω.x*seconds, ω.y*seconds, ω.z*seconds);
                Matrix3x3.mult_matrix(rotation_matrix, this.local_to_global_matrix, this.local_to_global_matrix);
            }
        } else {
            // the slab pull estimate is computed in the plate's own frame, so it is applied on the plate side
            const rotation_matrix = Tectonophysics.get_plate_rotation_matrix3x3(this.velocity.value(), this.center_of_mass.value(), seconds);
            Matrix3x3.mult_matrix(this.local_to_global_matrix, rotation_matrix, this.local_to_global_matrix);
        }
        Matrix3x3.invert(this.local_to_global_matrix, this.global_to_local_matrix);

        this.update_ids();
    }

    // Lookups must be valid as soon as the plate exists: the lithosphere merges plates (and integrates deltas
    // on loaded plates) before the plate's first move(). Left at zero, every cell would map to cell 0 for a frame.
    this.update_ids();
}
// switch for the sparse treatment of plates with prescribed motion (see Plate above); kept switchable for testing
Plate.sparse_prescribed_plates = true;
