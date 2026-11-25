function Simulation (parameters) {
    parameters = parameters || {};
    let _model                 = void 0;
    let _focus_id              = void 0;
    this.paused             = parameters.paused || false;
    this.speed                 = parameters.speed || 1;
    this.elapsed_time        = parameters.elapsed_time || 0;
    this.seed                 = parameters.seed || 0;
    this.random             = new Random(parseSeed(this.seed));
    if (parameters.random !== void 0) {
        this.random.mt      = parameters.random.mt;
        this.random.mti     = parameters.random.mti;
    }
    this._last_update_timestamp = 0;

    function parseSeed(text) {
            let parsed = parseInt(text);
            if(isNaN(parsed)) {
                parsed = 0;
                for (let i = 0; i < Math.min(8, text.length); i++) {
                    parsed = (parsed * 256) + text.charCodeAt(i);
                }
            }
            return parsed;
        }
        
    // the "model" is the singular entity we are simulating
    this.model = function(model) {
        if (model === void 0) {
            return _model;
        };
        _model = model;
        model.initialize();
    };
    if (parameters.model !== void 0) {
        this.model(new Universe(parameters.model));
    }

    // the "focus" is component of the model on which we focus the camera and ui
    // it is currently set to always be the world we are simulating
    this.focus = function(value) {
        if (value === void 0) {
            return _model.bodies[_focus_id];
        };
        _focus_id = value.id;
    };
    if (parameters.focus !== void 0 && _model !== void 0) {
        _focus_id = parameters.focus;
    } 

    this.getParameters = function() {
        return {
            model:             _model !== void 0? _model.getParameters() : undefined,
            paused:         this.paused,
            speed:             this.speed,
            elapsed_time:     this.elapsed_time,
            seed:             this.seed,
            focus:             _focus_id, 
            random: {
                mt:  this.random.mt,
                mti: this.random.mti
            },
        };
    }

    // The longest step the model takes at once. Longer steps let plates jump several cells, opening gaps and
    // overlaps the lithosphere cannot resolve, so faster speeds are reached with several steps per frame.
    this.max_timestep = 2 * Units.MEGAYEAR;
    // Milliseconds a frame may spend on extra steps. Past that, the rest of the frame's simulated time is dropped:
    // the simulation runs slower than the chosen speed instead of freezing the page.
    this.frame_budget = 100;

    this.update = function() {
        const now = performance.now();
        // limit each frame to what a 5fps frame would advance. Frames are clamped rather than skipped: skipping
        // made slow updates (big grids, many plates) alternate between a skipped frame and a nearly empty one,
        // so the simulation barely advanced. Clamping also covers long gaps such as a backgrounded tab.
        const seconds = Math.min((now - this._last_update_timestamp)/1000, 1/5);
        this._last_update_timestamp = now;

        if (this.paused){
            return;
        }

        if (_model === void 0) {
            return;
        }

        const frame_time = this.speed * seconds;
        const step_count = Math.max(1, Math.ceil(frame_time / this.max_timestep));
        const timestep = frame_time / step_count;
        for (let i = 0; i < step_count; i++) {
            this.elapsed_time += timestep;
            _model.invalidate(timestep);
            _model.calcChanges(timestep);
            _model.applyChanges(timestep);
            if (performance.now() - now > this.frame_budget) {
                break;
            }
        }
    };

    this.toggle_pause = function () {
        this.paused = !this.paused;
    }
}
