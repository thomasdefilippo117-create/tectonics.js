'use strict';

function SupercontinentCycle(lithosphere, parameters){
    parameters = parameters || {};

    this.getRandomDuration = function() { return 150 * Units.MEGAYEAR; };
        // function() { return random.uniform(300, 500); };
        // from wikipedia
    
    this.lithosphere = lithosphere;
    // NOTE: saves store Infinity (breakup disabled) as null, since JSON has no Infinity
    // NOTE: an age of 0 is meaningful (a cycle that just started), so test for presence instead of truthiness
    this.duration = parameters['duration'] === null? Infinity
        : parameters['duration'] !== void 0? parameters['duration']
        : this.getRandomDuration();
    // without a saved age the cycle starts out ended, so plates are segmented on the first update
    this.age = typeof parameters['age'] === 'number'? parameters['age'] : this.duration;

    this.getParameters = function() {
        return { 
            duration:     this.duration,
            age:         this.age,
        };
    }
};
SupercontinentCycle.prototype.update = function(seconds) {
    this.age += seconds;
    if(this.isEnding() === true){
        this.restart();
    }
};
SupercontinentCycle.prototype.isEnding = function() {
    //return this.lithosphere.plates.length <= 2;
    return this.age >= this.duration;
};
SupercontinentCycle.prototype.restart = function() {
    this.age = 0;
    this.duration = this.getRandomDuration();
    this.lithosphere.resetPlates();
}; 
