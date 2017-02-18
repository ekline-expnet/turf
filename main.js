'use strict';

//add any turf function needed by including them on the turf object
//to rebuild turf.js, run:
//  browserify main.js > turf.js
var turf = {
  "along": require("@turf/along"),
  "area": require("@turf/area"),
  "bbox": require("@turf/bbox"),
  "bearing": require("@turf/bearing"),
  "centroid": require("@turf/centroid"),
  "circle": require("@turf/circle"),
  "distance": require("@turf/distance"),
  "envelope": require("@turf/envelope"),
  "midpoint": require("@turf/midpoint"),
  "square": require("@turf/square")
}

var helpers = require("@turf/helpers");
var invariant = require("@turf/invariant");


var toMergeWithTurf = [helpers, invariant];

toMergeWithTurf.forEach(function(obj) {
  var keys = Object.keys(obj);
  keys.forEach(function(k,v){
    turf[k] = obj[k];
  });
})

window.turf = turf;
