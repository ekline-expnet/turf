'use strict';

//add any turf function needed by including them on the turf object
//to rebuild turf.js, run:
//  browserify main.js > turf.js
var turf = {
  "distance": require("@turf/distance")
}

window.turf = turf;
