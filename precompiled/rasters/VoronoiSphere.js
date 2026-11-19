var VoronoiSphere = (function() {
    
    const OCTAHEDRON_SIDE_COUNT = 8;    // number of sides on the data cube
    var OCTAHEDRON_SIDE_Z = VectorRaster.FromVectors([
            Vector(-1,-1,-1),
            Vector( 1,-1,-1),
            Vector(-1, 1,-1),
            Vector( 1, 1,-1),
            Vector(-1,-1, 1),
            Vector( 1,-1, 1),
            Vector(-1, 1, 1),
            Vector( 1, 1, 1)
        ]);
    VectorField.normalize(OCTAHEDRON_SIDE_Z, OCTAHEDRON_SIDE_Z);
    var OCTAHEDRON_SIDE_X = VectorField.cross_vector(OCTAHEDRON_SIDE_Z, Vector(0,0,1));
    VectorField.normalize(OCTAHEDRON_SIDE_X, OCTAHEDRON_SIDE_X);
    var OCTAHEDRON_SIDE_Y = VectorField.cross_vector_field(OCTAHEDRON_SIDE_Z, OCTAHEDRON_SIDE_X);
    VectorField.normalize(OCTAHEDRON_SIDE_Y, OCTAHEDRON_SIDE_Y);
    var OCTAHEDRON_SIDE_X = VectorRaster.ToArray(OCTAHEDRON_SIDE_X);
    var OCTAHEDRON_SIDE_Y = VectorRaster.ToArray(OCTAHEDRON_SIDE_Y);
    var OCTAHEDRON_SIDE_Z = VectorRaster.ToArray(OCTAHEDRON_SIDE_Z);
    // flat copies of the side axes so lookups avoid object property access and function calls
    var OCTAHEDRON_SIDE_X_FLAT = new Float64Array(OCTAHEDRON_SIDE_COUNT * 3);
    var OCTAHEDRON_SIDE_Y_FLAT = new Float64Array(OCTAHEDRON_SIDE_COUNT * 3);
    for (var side = 0; side < OCTAHEDRON_SIDE_COUNT; side++) {
        OCTAHEDRON_SIDE_X_FLAT[3*side+0] = OCTAHEDRON_SIDE_X[side].x;
        OCTAHEDRON_SIDE_X_FLAT[3*side+1] = OCTAHEDRON_SIDE_X[side].y;
        OCTAHEDRON_SIDE_X_FLAT[3*side+2] = OCTAHEDRON_SIDE_X[side].z;
        OCTAHEDRON_SIDE_Y_FLAT[3*side+0] = OCTAHEDRON_SIDE_Y[side].x;
        OCTAHEDRON_SIDE_Y_FLAT[3*side+1] = OCTAHEDRON_SIDE_Y[side].y;
        OCTAHEDRON_SIDE_Y_FLAT[3*side+2] = OCTAHEDRON_SIDE_Y[side].z;
    }


    var cell_count = function (dimensions_x, dimensions_y){
        return OCTAHEDRON_SIDE_COUNT * dimensions_x * dimensions_y;
    }
    var cell_id = function (side_id, xi, yi, dimensions_x, dimensions_y){
        return  side_id * dimensions_x * dimensions_y
              + xi      * dimensions_y 
              + yi;
    }

    //Data structure mapping coordinates on a sphere to the nearest neighbor
    //Retrievals from the map are of O(1) complexity. The result resembles a voronoi diagram, hence the name.
    function VoronoiSphere(pos, cell_width, farthest_distance){
        var dimension_x = Math.ceil(2./cell_width)+1;
        var dimension_y = Math.ceil(2./cell_width)+1;
        var cells = new Uint16Array(cell_count(dimension_x, dimension_y));

        this.cell_width = cell_width;
        this.dimension_x = dimension_x;
        this.dimension_y = dimension_y;
        this.cells = cells;

        //Feed locations into an integer lattice for fast lookups
        var points = [];
        var x = pos.x;
        var y = pos.y;
        var z = pos.z;
        for(var i=0, il = x.length; i<il; i++){
            points.push({x:x[i], y:y[i], z:z[i], i:i});
        }
        var getDistance = function(a,b) { 
                return (a.x - b.x)*(a.x - b.x) +  (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z); 
            };
        var lattice = new IntegerLattice(points, getDistance, farthest_distance);

        var side_x = OCTAHEDRON_SIDE_X[0];
        var side_y = OCTAHEDRON_SIDE_Y[0];
        var side_z = OCTAHEDRON_SIDE_Z[0];
        var cell_x = Vector();
        var cell_y = Vector();
        var cell_z = Vector();
        var cell_pos = Vector();
        var sqrt = Math.sqrt;
        var max = Math.max;

        var x2d = 0.;
        var y2d = 0.;
        var z2d = 0.;
        // populate cells using the slower IntegerLattice implementation
        for (var side_id = 0; side_id < OCTAHEDRON_SIDE_COUNT; side_id++)
        {
            side_x = OCTAHEDRON_SIDE_X[side_id];
            side_y = OCTAHEDRON_SIDE_Y[side_id];
            side_z = OCTAHEDRON_SIDE_Z[side_id];
            for (var xi2d = 0; xi2d < dimension_x; xi2d++)
            {
                for (var yi2d = 0; yi2d < dimension_y; yi2d++)
                {
                    // get position of the cell that's projected onto the 2d grid
                    // lookups floor() into a cell, so the cell is represented by its center, not its corner
                    x2d = (xi2d + 0.5) * cell_width - 1.;
                    y2d = (yi2d + 0.5) * cell_width - 1.;

                    // reconstruct the dimension omitted from the grid using pythagorean theorem
                    z2d = sqrt(max(1. - (x2d*x2d) - (y2d*y2d), 0.));

                    Vector.mult_scalar(side_x.x, side_x.y, side_x.z, x2d, cell_x);
                    Vector.mult_scalar(side_y.x, side_y.y, side_y.z, y2d, cell_y);
                    Vector.mult_scalar(side_z.x, side_z.y, side_z.z, z2d, cell_z);

                    // reset vector
                    cell_pos.x = 0;
                    cell_pos.y = 0;
                    cell_pos.z = 0;

                    Vector.add_vector(cell_pos.x, cell_pos.y, cell_pos.z, cell_x.x, cell_x.y, cell_x.z, cell_pos);
                    Vector.add_vector(cell_pos.x, cell_pos.y, cell_pos.z, cell_y.x, cell_y.y, cell_y.z, cell_pos);
                    Vector.add_vector(cell_pos.x, cell_pos.y, cell_pos.z, cell_z.x, cell_z.y, cell_z.z, cell_pos);
                    cells[cell_id(side_id, xi2d, yi2d, dimension_x, dimension_y)] =  lattice.nearest(cell_pos).i;
                }
            }
        }
    }
    VoronoiSphere.prototype.getNearestIds = function(pos_field, result) {
        result = result || new Uint16Array(pos_field.x.length);

        var cell_width = this.cell_width;
        var dimension_y = this.dimension_y;
        var side_size = this.dimension_x * this.dimension_y;
        var cells = this.cells;

        var side_x = OCTAHEDRON_SIDE_X_FLAT;
        var side_y = OCTAHEDRON_SIDE_Y_FLAT;
        var pos_x = pos_field.x;
        var pos_y = pos_field.y;
        var pos_z = pos_field.z;
        var floor = Math.floor;

        // same arithmetic as Vector.dot_vector() and cell_id(), inlined: this runs twice per plate every frame
        for (var i = 0, li = pos_x.length; i < li; i++)
        {
            var xi = pos_x[i];
            var yi = pos_y[i];
            var zi = pos_z[i];

            var side_id =
              (( xi > 0)     ) +
              (( yi > 0) << 1) +
              (( zi > 0) << 2) ;
            var s = 3 * side_id;

            var projection_x = (side_x[s] * xi + side_x[s+1] * yi + side_x[s+2] * zi);
            var projection_y = (side_y[s] * xi + side_y[s+1] * yi + side_y[s+2] * zi);

            var grid_pos_x = floor((projection_x + 1.) / cell_width);
            var grid_pos_y = floor((projection_y + 1.) / cell_width);

            result[i] = cells[side_id * side_size + grid_pos_x * dimension_y + grid_pos_y];
        }
        return result;
    }
    return VoronoiSphere;
})();