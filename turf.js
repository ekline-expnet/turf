(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
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

},{"@turf/along":3,"@turf/area":4,"@turf/bbox":6,"@turf/bearing":7,"@turf/centroid":8,"@turf/circle":9,"@turf/distance":11,"@turf/envelope":12,"@turf/helpers":13,"@turf/invariant":14,"@turf/midpoint":16,"@turf/square":17}],2:[function(require,module,exports){
var wgs84 = require('wgs84');

module.exports.geometry = geometry;
module.exports.ring = ringArea;

function geometry(_) {
    var area = 0, i;
    switch (_.type) {
        case 'Polygon':
            return polygonArea(_.coordinates);
        case 'MultiPolygon':
            for (i = 0; i < _.coordinates.length; i++) {
                area += polygonArea(_.coordinates[i]);
            }
            return area;
        case 'Point':
        case 'MultiPoint':
        case 'LineString':
        case 'MultiLineString':
            return 0;
        case 'GeometryCollection':
            for (i = 0; i < _.geometries.length; i++) {
                area += geometry(_.geometries[i]);
            }
            return area;
    }
}

function polygonArea(coords) {
    var area = 0;
    if (coords && coords.length > 0) {
        area += Math.abs(ringArea(coords[0]));
        for (var i = 1; i < coords.length; i++) {
            area -= Math.abs(ringArea(coords[i]));
        }
    }
    return area;
}

/**
 * Calculate the approximate area of the polygon were it projected onto
 *     the earth.  Note that this area will be positive if ring is oriented
 *     clockwise, otherwise it will be negative.
 *
 * Reference:
 * Robert. G. Chamberlain and William H. Duquette, "Some Algorithms for
 *     Polygons on a Sphere", JPL Publication 07-03, Jet Propulsion
 *     Laboratory, Pasadena, CA, June 2007 http://trs-new.jpl.nasa.gov/dspace/handle/2014/40409
 *
 * Returns:
 * {float} The approximate signed geodesic area of the polygon in square
 *     meters.
 */

function ringArea(coords) {
    var p1, p2, p3, lowerIndex, middleIndex, upperIndex, i,
    area = 0,
    coordsLength = coords.length;

    if (coordsLength > 2) {
        for (i = 0; i < coordsLength; i++) {
            if (i === coordsLength - 2) {// i = N-2
                lowerIndex = coordsLength - 2;
                middleIndex = coordsLength -1;
                upperIndex = 0;
            } else if (i === coordsLength - 1) {// i = N-1
                lowerIndex = coordsLength - 1;
                middleIndex = 0;
                upperIndex = 1;
            } else { // i = 0 to N-3
                lowerIndex = i;
                middleIndex = i+1;
                upperIndex = i+2;
            }
            p1 = coords[lowerIndex];
            p2 = coords[middleIndex];
            p3 = coords[upperIndex];
            area += ( rad(p3[0]) - rad(p1[0]) ) * Math.sin( rad(p2[1]));
        }

        area = area * wgs84.RADIUS * wgs84.RADIUS / 2;
    }

    return area;
}

function rad(_) {
    return _ * Math.PI / 180;
}
},{"wgs84":18}],3:[function(require,module,exports){
var measureDistance = require('@turf/distance');
var point = require('@turf/helpers').point;
var bearing = require('@turf/bearing');
var destination = require('@turf/destination');

/**
 * Takes a {@link LineString|line} and returns a {@link Point|point} at a specified distance along the line.
 *
 * @name along
 * @param {Feature<LineString>} line input line
 * @param {number} distance distance along the line
 * @param {string} [units=kilometers] can be degrees, radians, miles, or kilometers
 * @return {Feature<Point>} Point `distance` `units` along the line
 * @example
 * var line = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "LineString",
 *     "coordinates": [
 *       [-77.031669, 38.878605],
 *       [-77.029609, 38.881946],
 *       [-77.020339, 38.884084],
 *       [-77.025661, 38.885821],
 *       [-77.021884, 38.889563],
 *       [-77.019824, 38.892368]
 *     ]
 *   }
 * };
 *
 * var along = turf.along(line, 1, 'miles');
 *
 * var result = {
 *   "type": "FeatureCollection",
 *   "features": [line, along]
 * };
 *
 * //=result
 */
module.exports = function (line, distance, units) {
    var coords;
    if (line.type === 'Feature') coords = line.geometry.coordinates;
    else if (line.type === 'LineString') coords = line.coordinates;
    else throw new Error('input must be a LineString Feature or Geometry');

    var travelled = 0;
    for (var i = 0; i < coords.length; i++) {
        if (distance >= travelled && i === coords.length - 1) break;
        else if (travelled >= distance) {
            var overshot = distance - travelled;
            if (!overshot) return point(coords[i]);
            else {
                var direction = bearing(coords[i], coords[i - 1]) - 180;
                var interpolated = destination(coords[i], overshot, direction, units);
                return interpolated;
            }
        } else {
            travelled += measureDistance(coords[i], coords[i + 1], units);
        }
    }
    return point(coords[coords.length - 1]);
};

},{"@turf/bearing":7,"@turf/destination":10,"@turf/distance":11,"@turf/helpers":13}],4:[function(require,module,exports){
var geometryArea = require('@mapbox/geojson-area').geometry;

/**
 * Takes one or more features and returns their area
 * in square meters.
 *
 * @param {(Feature|FeatureCollection)} input input features
 * @return {number} area in square meters
 * @example
 * var polygons = {
 *   "type": "FeatureCollection",
 *   "features": [
 *     {
 *       "type": "Feature",
 *       "properties": {},
 *       "geometry": {
 *         "type": "Polygon",
 *         "coordinates": [[
 *           [-67.031021, 10.458102],
 *           [-67.031021, 10.53372],
 *           [-66.929397, 10.53372],
 *           [-66.929397, 10.458102],
 *           [-67.031021, 10.458102]
 *         ]]
 *       }
 *     }, {
 *       "type": "Feature",
 *       "properties": {},
 *       "geometry": {
 *         "type": "Polygon",
 *         "coordinates": [[
 *           [-66.919784, 10.397325],
 *           [-66.919784, 10.513467],
 *           [-66.805114, 10.513467],
 *           [-66.805114, 10.397325],
 *           [-66.919784, 10.397325]
 *         ]]
 *       }
 *     }
 *   ]
 * };
 *
 * var area = turf.area(polygons);
 *
 * //=area
 */
function area(input) {
    if (input.type === 'FeatureCollection') {
        for (var i = 0, sum = 0; i < input.features.length; i++) {
            if (input.features[i].geometry) {
                sum += geometryArea(input.features[i].geometry);
            }
        }
        return sum;
    } else if (input.type === 'Feature') {
        return geometryArea(input.geometry);
    } else {
        return geometryArea(input);
    }
}
module.exports = area;

},{"@mapbox/geojson-area":2}],5:[function(require,module,exports){
var polygon = require('@turf/helpers').polygon;

/**
 * Takes a bbox and returns an equivalent {@link Polygon|polygon}.
 *
 * @name bboxPolygon
 * @param {Array<number>} bbox extent in [minX, minY, maxX, maxY] order
 * @return {Feature<Polygon>} a Polygon representation of the bounding box
 * @example
 * var bbox = [0, 0, 10, 10];
 *
 * var poly = turf.bboxPolygon(bbox);
 *
 * //=poly
 */

module.exports = function (bbox) {
    var lowLeft = [bbox[0], bbox[1]];
    var topLeft = [bbox[0], bbox[3]];
    var topRight = [bbox[2], bbox[3]];
    var lowRight = [bbox[2], bbox[1]];

    return polygon([[
        lowLeft,
        lowRight,
        topRight,
        topLeft,
        lowLeft
    ]]);
};

},{"@turf/helpers":13}],6:[function(require,module,exports){
var each = require('@turf/meta').coordEach;

/**
 * Takes a set of features, calculates the bbox of all input features, and returns a bounding box.
 *
 * @name bbox
 * @param {(Feature|FeatureCollection)} geojson input features
 * @return {Array<number>} bbox extent in [minX, minY, maxX, maxY] order
 * @example
 * var pt1 = turf.point([114.175329, 22.2524])
 * var pt2 = turf.point([114.170007, 22.267969])
 * var pt3 = turf.point([114.200649, 22.274641])
 * var pt4 = turf.point([114.200649, 22.274641])
 * var pt5 = turf.point([114.186744, 22.265745])
 * var features = turf.featureCollection([pt1, pt2, pt3, pt4, pt5])
 *
 * var bbox = turf.bbox(features);
 *
 * var bboxPolygon = turf.bboxPolygon(bbox);
 *
 * //=bbox
 *
 * //=bboxPolygon
 */
module.exports = function (geojson) {
    var bbox = [Infinity, Infinity, -Infinity, -Infinity];
    each(geojson, function (coord) {
        if (bbox[0] > coord[0]) bbox[0] = coord[0];
        if (bbox[1] > coord[1]) bbox[1] = coord[1];
        if (bbox[2] < coord[0]) bbox[2] = coord[0];
        if (bbox[3] < coord[1]) bbox[3] = coord[1];
    });
    return bbox;
};

},{"@turf/meta":15}],7:[function(require,module,exports){
var getCoord = require('@turf/invariant').getCoord;
//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html

/**
 * Takes two {@link Point|points} and finds the geographic bearing between them.
 *
 * @name bearing
 * @param {Feature<Point>} start starting Point
 * @param {Feature<Point>} end ending Point
 * @returns {number} bearing in decimal degrees
 * @example
 * var point1 = {
 *   "type": "Feature",
 *   "properties": {
 *     "marker-color": '#f00'
 *   },
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [-75.343, 39.984]
 *   }
 * };
 * var point2 = {
 *   "type": "Feature",
 *   "properties": {
 *     "marker-color": '#0f0'
 *   },
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [-75.534, 39.123]
 *   }
 * };
 *
 * var points = {
 *   "type": "FeatureCollection",
 *   "features": [point1, point2]
 * };
 *
 * //=points
 *
 * var bearing = turf.bearing(point1, point2);
 *
 * //=bearing
 */
module.exports = function (start, end) {
    var degrees2radians = Math.PI / 180;
    var radians2degrees = 180 / Math.PI;
    var coordinates1 = getCoord(start);
    var coordinates2 = getCoord(end);

    var lon1 = degrees2radians * coordinates1[0];
    var lon2 = degrees2radians * coordinates2[0];
    var lat1 = degrees2radians * coordinates1[1];
    var lat2 = degrees2radians * coordinates2[1];
    var a = Math.sin(lon2 - lon1) * Math.cos(lat2);
    var b = Math.cos(lat1) * Math.sin(lat2) -
        Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1);

    var bearing = radians2degrees * Math.atan2(a, b);

    return bearing;
};

},{"@turf/invariant":14}],8:[function(require,module,exports){
var each = require('@turf/meta').coordEach;
var point = require('@turf/helpers').point;

/**
 * Takes one or more features and calculates the centroid using
 * the mean of all vertices.
 * This lessens the effect of small islands and artifacts when calculating
 * the centroid of a set of polygons.
 *
 * @name centroid
 * @param {(Feature|FeatureCollection)} features input features
 * @return {Feature<Point>} the centroid of the input features
 * @example
 * var poly = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Polygon",
 *     "coordinates": [[
 *       [105.818939,21.004714],
 *       [105.818939,21.061754],
 *       [105.890007,21.061754],
 *       [105.890007,21.004714],
 *       [105.818939,21.004714]
 *     ]]
 *   }
 * };
 *
 * var centroidPt = turf.centroid(poly);
 *
 * var result = {
 *   "type": "FeatureCollection",
 *   "features": [poly, centroidPt]
 * };
 *
 * //=result
 */
module.exports = function (features) {
    var xSum = 0, ySum = 0, len = 0;
    each(features, function (coord) {
        xSum += coord[0];
        ySum += coord[1];
        len++;
    }, true);
    return point([xSum / len, ySum / len]);
};

},{"@turf/helpers":13,"@turf/meta":15}],9:[function(require,module,exports){
var destination = require('@turf/destination');
var helpers = require('@turf/helpers');
var polygon = helpers.polygon;

/**
 * Takes a {@link Point} and calculates the circle polygon given a radius in degrees, radians, miles, or kilometers; and steps for precision.
 *
 * @name circle
 * @param {Feature<Point>} center center point
 * @param {number} radius radius of the circle
 * @param {number} [steps=64] number of steps
 * @param {string} [units=kilometers] miles, kilometers, degrees, or radians
 * @returns {Feature<Polygon>} circle polygon
 * @example
 * var center = turf.point([-75.343, 39.984]);
 * var radius = 5;
 * var steps = 10;
 * var units = 'kilometers';
 *
 * var circle = turf.circle(center, radius, steps, units);
 *
 * //=circle
 */
module.exports = function (center, radius, steps, units) {
    steps = steps || 64;
    var coordinates = [];

    for (var i = 0; i < steps; i++) {
        coordinates.push(destination(center, radius, i * 360 / steps, units).geometry.coordinates);
    }

    coordinates.push(coordinates[0]);

    return polygon([coordinates]);
};

},{"@turf/destination":10,"@turf/helpers":13}],10:[function(require,module,exports){
//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html
var getCoord = require('@turf/invariant').getCoord;
var helpers = require('@turf/helpers');
var point = helpers.point;
var distanceToRadians = helpers.distanceToRadians;

/**
 * Takes a {@link Point} and calculates the location of a destination point given a distance in degrees, radians, miles, or kilometers; and bearing in degrees. This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
 *
 * @name destination
 * @param {Feature<Point>} from starting point
 * @param {number} distance distance from the starting point
 * @param {number} bearing ranging from -180 to 180
 * @param {string} [units=kilometers] miles, kilometers, degrees, or radians
 * @returns {Feature<Point>} destination point
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {
 *     "marker-color": "#0f0"
 *   },
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [-75.343, 39.984]
 *   }
 * };
 * var distance = 50;
 * var bearing = 90;
 * var units = 'miles';
 *
 * var destination = turf.destination(point, distance, bearing, units);
 * destination.properties['marker-color'] = '#f00';
 *
 * var result = {
 *   "type": "FeatureCollection",
 *   "features": [point, destination]
 * };
 *
 * //=result
 */
module.exports = function (from, distance, bearing, units) {
    var degrees2radians = Math.PI / 180;
    var radians2degrees = 180 / Math.PI;
    var coordinates1 = getCoord(from);
    var longitude1 = degrees2radians * coordinates1[0];
    var latitude1 = degrees2radians * coordinates1[1];
    var bearing_rad = degrees2radians * bearing;

    var radians = distanceToRadians(distance, units);

    var latitude2 = Math.asin(Math.sin(latitude1) * Math.cos(radians) +
        Math.cos(latitude1) * Math.sin(radians) * Math.cos(bearing_rad));
    var longitude2 = longitude1 + Math.atan2(Math.sin(bearing_rad) *
        Math.sin(radians) * Math.cos(latitude1),
        Math.cos(radians) - Math.sin(latitude1) * Math.sin(latitude2));

    return point([radians2degrees * longitude2, radians2degrees * latitude2]);
};

},{"@turf/helpers":13,"@turf/invariant":14}],11:[function(require,module,exports){
var getCoord = require('@turf/invariant').getCoord;
var radiansToDistance = require('@turf/helpers').radiansToDistance;
//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html

/**
 * Calculates the distance between two {@link Point|points} in degrees, radians,
 * miles, or kilometers. This uses the
 * [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula)
 * to account for global curvature.
 *
 * @name distance
 * @param {Feature<Point>} from origin point
 * @param {Feature<Point>} to destination point
 * @param {string} [units=kilometers] can be degrees, radians, miles, or kilometers
 * @return {number} distance between the two points
 * @example
 * var from = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [-75.343, 39.984]
 *   }
 * };
 * var to = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [-75.534, 39.123]
 *   }
 * };
 * var units = "miles";
 *
 * var points = {
 *   "type": "FeatureCollection",
 *   "features": [from, to]
 * };
 *
 * //=points
 *
 * var distance = turf.distance(from, to, units);
 *
 * //=distance
 */
module.exports = function (from, to, units) {
    var degrees2radians = Math.PI / 180;
    var coordinates1 = getCoord(from);
    var coordinates2 = getCoord(to);
    var dLat = degrees2radians * (coordinates2[1] - coordinates1[1]);
    var dLon = degrees2radians * (coordinates2[0] - coordinates1[0]);
    var lat1 = degrees2radians * coordinates1[1];
    var lat2 = degrees2radians * coordinates2[1];

    var a = Math.pow(Math.sin(dLat / 2), 2) +
          Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);

    return radiansToDistance(2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a)), units);
};

},{"@turf/helpers":13,"@turf/invariant":14}],12:[function(require,module,exports){
var bbox = require('@turf/bbox');
var bboxPolygon = require('@turf/bbox-polygon');

/**
 * Takes any number of features and returns a rectangular {@link Polygon} that encompasses all vertices.
 *
 * @name envelope
 * @param {(Feature|FeatureCollection)} features input features
 * @return {Feature<Polygon>} a rectangular Polygon feature that encompasses all vertices
 * @example
 * var fc = {
 *   "type": "FeatureCollection",
 *   "features": [
 *     {
 *       "type": "Feature",
 *       "properties": {
 *         "name": "Location A"
 *       },
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [-75.343, 39.984]
 *       }
 *     }, {
 *       "type": "Feature",
 *       "properties": {
 *         "name": "Location B"
 *       },
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [-75.833, 39.284]
 *       }
 *     }, {
 *       "type": "Feature",
 *       "properties": {
 *         "name": "Location C"
 *       },
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [-75.534, 39.123]
 *       }
 *     }
 *   ]
 * };
 *
 * var enveloped = turf.envelope(fc);
 *
 * var resultFeatures = fc.features.concat(enveloped);
 * var result = {
 *   "type": "FeatureCollection",
 *   "features": resultFeatures
 * };
 *
 * //=result
 */

module.exports = function (features) {
    return bboxPolygon(bbox(features));
};

},{"@turf/bbox":6,"@turf/bbox-polygon":5}],13:[function(require,module,exports){
/**
 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
 *
 * @name feature
 * @param {Geometry} geometry input geometry
 * @param {Object} properties properties
 * @returns {FeatureCollection} a FeatureCollection of input features
 * @example
 * var geometry = {
 *      "type": "Point",
 *      "coordinates": [
 *        67.5,
 *        32.84267363195431
 *      ]
 *    }
 *
 * var feature = turf.feature(geometry);
 *
 * //=feature
 */
function feature(geometry, properties) {
    return {
        type: 'Feature',
        properties: properties || {},
        geometry: geometry
    };
}

module.exports.feature = feature;

/**
 * Takes coordinates and properties (optional) and returns a new {@link Point} feature.
 *
 * @name point
 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
 * @param {Object=} properties an Object that is used as the {@link Feature}'s
 * properties
 * @returns {Feature<Point>} a Point feature
 * @example
 * var pt1 = turf.point([-75.343, 39.984]);
 *
 * //=pt1
 */
module.exports.point = function (coordinates, properties) {
    if (!Array.isArray(coordinates)) throw new Error('Coordinates must be an array');
    if (coordinates.length < 2) throw new Error('Coordinates must be at least 2 numbers long');
    return feature({
        type: 'Point',
        coordinates: coordinates.slice()
    }, properties);
};

/**
 * Takes an array of LinearRings and optionally an {@link Object} with properties and returns a {@link Polygon} feature.
 *
 * @name polygon
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object=} properties a properties object
 * @returns {Feature<Polygon>} a Polygon feature
 * @throws {Error} throw an error if a LinearRing of the polygon has too few positions
 * or if a LinearRing of the Polygon does not have matching Positions at the
 * beginning & end.
 * @example
 * var polygon = turf.polygon([[
 *  [-2.275543, 53.464547],
 *  [-2.275543, 53.489271],
 *  [-2.215118, 53.489271],
 *  [-2.215118, 53.464547],
 *  [-2.275543, 53.464547]
 * ]], { name: 'poly1', population: 400});
 *
 * //=polygon
 */
module.exports.polygon = function (coordinates, properties) {

    if (!coordinates) throw new Error('No coordinates passed');

    for (var i = 0; i < coordinates.length; i++) {
        var ring = coordinates[i];
        if (ring.length < 4) {
            throw new Error('Each LinearRing of a Polygon must have 4 or more Positions.');
        }
        for (var j = 0; j < ring[ring.length - 1].length; j++) {
            if (ring[ring.length - 1][j] !== ring[0][j]) {
                throw new Error('First and last Position are not equivalent.');
            }
        }
    }

    return feature({
        type: 'Polygon',
        coordinates: coordinates
    }, properties);
};

/**
 * Creates a {@link LineString} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name lineString
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object=} properties an Object of key-value pairs to add as properties
 * @returns {Feature<LineString>} a LineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var linestring1 = turf.lineString([
 *	[-21.964416, 64.148203],
 *	[-21.956176, 64.141316],
 *	[-21.93901, 64.135924],
 *	[-21.927337, 64.136673]
 * ]);
 * var linestring2 = turf.lineString([
 *	[-21.929054, 64.127985],
 *	[-21.912918, 64.134726],
 *	[-21.916007, 64.141016],
 * 	[-21.930084, 64.14446]
 * ], {name: 'line 1', distance: 145});
 *
 * //=linestring1
 *
 * //=linestring2
 */
module.exports.lineString = function (coordinates, properties) {
    if (!coordinates) {
        throw new Error('No coordinates passed');
    }
    return feature({
        type: 'LineString',
        coordinates: coordinates
    }, properties);
};

/**
 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
 *
 * @name featureCollection
 * @param {Feature[]} features input features
 * @returns {FeatureCollection} a FeatureCollection of input features
 * @example
 * var features = [
 *  turf.point([-75.343, 39.984], {name: 'Location A'}),
 *  turf.point([-75.833, 39.284], {name: 'Location B'}),
 *  turf.point([-75.534, 39.123], {name: 'Location C'})
 * ];
 *
 * var fc = turf.featureCollection(features);
 *
 * //=fc
 */
module.exports.featureCollection = function (features) {
    return {
        type: 'FeatureCollection',
        features: features
    };
};

/**
 * Creates a {@link Feature<MultiLineString>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiLineString
 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
 * @param {Object=} properties an Object of key-value pairs to add as properties
 * @returns {Feature<MultiLineString>} a MultiLineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
 *
 * //=multiLine
 *
 */
module.exports.multiLineString = function (coordinates, properties) {
    if (!coordinates) {
        throw new Error('No coordinates passed');
    }
    return feature({
        type: 'MultiLineString',
        coordinates: coordinates
    }, properties);
};

/**
 * Creates a {@link Feature<MultiPoint>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPoint
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object=} properties an Object of key-value pairs to add as properties
 * @returns {Feature<MultiPoint>} a MultiPoint feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
 *
 * //=multiPt
 *
 */
module.exports.multiPoint = function (coordinates, properties) {
    if (!coordinates) {
        throw new Error('No coordinates passed');
    }
    return feature({
        type: 'MultiPoint',
        coordinates: coordinates
    }, properties);
};


/**
 * Creates a {@link Feature<MultiPolygon>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPolygon
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
 * @param {Object=} properties an Object of key-value pairs to add as properties
 * @returns {Feature<MultiPolygon>} a multipolygon feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
 *
 * //=multiPoly
 *
 */
module.exports.multiPolygon = function (coordinates, properties) {
    if (!coordinates) {
        throw new Error('No coordinates passed');
    }
    return feature({
        type: 'MultiPolygon',
        coordinates: coordinates
    }, properties);
};

/**
 * Creates a {@link Feature<GeometryCollection>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name geometryCollection
 * @param {Array<{Geometry}>} geometries an array of GeoJSON Geometries
 * @param {Object=} properties an Object of key-value pairs to add as properties
 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
 * @example
 * var pt = {
 *     "type": "Point",
 *       "coordinates": [100, 0]
 *     };
 * var line = {
 *     "type": "LineString",
 *     "coordinates": [ [101, 0], [102, 1] ]
 *   };
 * var collection = turf.geometryCollection([pt, line]);
 *
 * //=collection
 */
module.exports.geometryCollection = function (geometries, properties) {
    return feature({
        type: 'GeometryCollection',
        geometries: geometries
    }, properties);
};

var factors = {
    miles: 3960,
    nauticalmiles: 3441.145,
    degrees: 57.2957795,
    radians: 1,
    inches: 250905600,
    yards: 6969600,
    meters: 6373000,
    metres: 6373000,
    kilometers: 6373,
    kilometres: 6373,
    feet: 20908792.65
};

/*
 * Convert a distance measurement from radians to a more friendly unit.
 *
 * @name radiansToDistance
 * @param {number} distance in radians across the sphere
 * @param {string} [units=kilometers] can be degrees, radians, miles, or kilometers
 * inches, yards, metres, meters, kilometres, kilometers.
 * @returns {number} distance
 */
module.exports.radiansToDistance = function (radians, units) {
    var factor = factors[units || 'kilometers'];
    if (factor === undefined) {
        throw new Error('Invalid unit');
    }
    return radians * factor;
};

/*
 * Convert a distance measurement from a real-world unit into radians
 *
 * @name distanceToRadians
 * @param {number} distance in real units
 * @param {string} [units=kilometers] can be degrees, radians, miles, or kilometers
 * inches, yards, metres, meters, kilometres, kilometers.
 * @returns {number} radians
 */
module.exports.distanceToRadians = function (distance, units) {
    var factor = factors[units || 'kilometers'];
    if (factor === undefined) {
        throw new Error('Invalid unit');
    }
    return distance / factor;
};

/*
 * Convert a distance measurement from a real-world unit into degrees
 *
 * @name distanceToRadians
 * @param {number} distance in real units
 * @param {string} [units=kilometers] can be degrees, radians, miles, or kilometers
 * inches, yards, metres, meters, kilometres, kilometers.
 * @returns {number} degrees
 */
module.exports.distanceToDegrees = function (distance, units) {
    var factor = factors[units || 'kilometers'];
    if (factor === undefined) {
        throw new Error('Invalid unit');
    }
    return (distance / factor) * 57.2958;
};

},{}],14:[function(require,module,exports){
/**
 * Unwrap a coordinate from a Feature with a Point geometry, a Point
 * geometry, or a single coordinate.
 *
 * @param {*} obj any value
 * @returns {Array<number>} a coordinate
 */
function getCoord(obj) {
    if (Array.isArray(obj) &&
        typeof obj[0] === 'number' &&
        typeof obj[1] === 'number') {
        return obj;
    } else if (obj) {
        if (obj.type === 'Feature' &&
            obj.geometry &&
            obj.geometry.type === 'Point' &&
            Array.isArray(obj.geometry.coordinates)) {
            return obj.geometry.coordinates;
        } else if (obj.type === 'Point' &&
            Array.isArray(obj.coordinates)) {
            return obj.coordinates;
        }
    }
    throw new Error('A coordinate, feature, or point geometry is required');
}

/**
 * Enforce expectations about types of GeoJSON objects for Turf.
 *
 * @alias geojsonType
 * @param {GeoJSON} value any GeoJSON object
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function geojsonType(value, type, name) {
    if (!type || !name) throw new Error('type and name required');

    if (!value || value.type !== type) {
        throw new Error('Invalid input to ' + name + ': must be a ' + type + ', given ' + value.type);
    }
}

/**
 * Enforce expectations about types of {@link Feature} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @alias featureOf
 * @param {Feature} feature a feature with an expected geometry type
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} error if value is not the expected type.
 */
function featureOf(feature, type, name) {
    if (!name) throw new Error('.featureOf() requires a name');
    if (!feature || feature.type !== 'Feature' || !feature.geometry) {
        throw new Error('Invalid input to ' + name + ', Feature with geometry required');
    }
    if (!feature.geometry || feature.geometry.type !== type) {
        throw new Error('Invalid input to ' + name + ': must be a ' + type + ', given ' + feature.geometry.type);
    }
}

/**
 * Enforce expectations about types of {@link FeatureCollection} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @alias collectionOf
 * @param {FeatureCollection} featurecollection a featurecollection for which features will be judged
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function collectionOf(featurecollection, type, name) {
    if (!name) throw new Error('.collectionOf() requires a name');
    if (!featurecollection || featurecollection.type !== 'FeatureCollection') {
        throw new Error('Invalid input to ' + name + ', FeatureCollection required');
    }
    for (var i = 0; i < featurecollection.features.length; i++) {
        var feature = featurecollection.features[i];
        if (!feature || feature.type !== 'Feature' || !feature.geometry) {
            throw new Error('Invalid input to ' + name + ', Feature with geometry required');
        }
        if (!feature.geometry || feature.geometry.type !== type) {
            throw new Error('Invalid input to ' + name + ': must be a ' + type + ', given ' + feature.geometry.type);
        }
    }
}

module.exports.geojsonType = geojsonType;
module.exports.collectionOf = collectionOf;
module.exports.featureOf = featureOf;
module.exports.getCoord = getCoord;

},{}],15:[function(require,module,exports){
/**
 * Iterate over coordinates in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name coordEach
 * @param {Object} layer any GeoJSON object
 * @param {Function} callback a method that takes (value)
 * @param {boolean=} excludeWrapCoord whether or not to include
 * the final coordinate of LinearRings that wraps the ring in its iteration.
 * @example
 * var point = { type: 'Point', coordinates: [0, 0] };
 * turfMeta.coordEach(point, function(coords) {
 *   // coords is equal to [0, 0]
 * });
 */
function coordEach(layer, callback, excludeWrapCoord) {
    var i, j, k, g, l, geometry, stopG, coords,
        geometryMaybeCollection,
        wrapShrink = 0,
        isGeometryCollection,
        isFeatureCollection = layer.type === 'FeatureCollection',
        isFeature = layer.type === 'Feature',
        stop = isFeatureCollection ? layer.features.length : 1;

  // This logic may look a little weird. The reason why it is that way
  // is because it's trying to be fast. GeoJSON supports multiple kinds
  // of objects at its root: FeatureCollection, Features, Geometries.
  // This function has the responsibility of handling all of them, and that
  // means that some of the `for` loops you see below actually just don't apply
  // to certain inputs. For instance, if you give this just a
  // Point geometry, then both loops are short-circuited and all we do
  // is gradually rename the input until it's called 'geometry'.
  //
  // This also aims to allocate as few resources as possible: just a
  // few numbers and booleans, rather than any temporary arrays as would
  // be required with the normalization approach.
    for (i = 0; i < stop; i++) {

        geometryMaybeCollection = (isFeatureCollection ? layer.features[i].geometry :
        (isFeature ? layer.geometry : layer));
        isGeometryCollection = geometryMaybeCollection.type === 'GeometryCollection';
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (g = 0; g < stopG; g++) {
            geometry = isGeometryCollection ?
            geometryMaybeCollection.geometries[g] : geometryMaybeCollection;
            coords = geometry.coordinates;

            wrapShrink = (excludeWrapCoord &&
                (geometry.type === 'Polygon' || geometry.type === 'MultiPolygon')) ?
                1 : 0;

            if (geometry.type === 'Point') {
                callback(coords);
            } else if (geometry.type === 'LineString' || geometry.type === 'MultiPoint') {
                for (j = 0; j < coords.length; j++) callback(coords[j]);
            } else if (geometry.type === 'Polygon' || geometry.type === 'MultiLineString') {
                for (j = 0; j < coords.length; j++)
                    for (k = 0; k < coords[j].length - wrapShrink; k++)
                        callback(coords[j][k]);
            } else if (geometry.type === 'MultiPolygon') {
                for (j = 0; j < coords.length; j++)
                    for (k = 0; k < coords[j].length; k++)
                        for (l = 0; l < coords[j][k].length - wrapShrink; l++)
                            callback(coords[j][k][l]);
            } else if (geometry.type === 'GeometryCollection') {
                for (j = 0; j < geometry.geometries.length; j++)
                    coordEach(geometry.geometries[j], callback, excludeWrapCoord);
            } else {
                throw new Error('Unknown Geometry Type');
            }
        }
    }
}
module.exports.coordEach = coordEach;

/**
 * Reduce coordinates in any GeoJSON object into a single value,
 * similar to how Array.reduce works. However, in this case we lazily run
 * the reduction, so an array of all coordinates is unnecessary.
 *
 * @name coordReduce
 * @param {Object} layer any GeoJSON object
 * @param {Function} callback a method that takes (memo, value) and returns
 * a new memo
 * @param {*} memo the starting value of memo: can be any type.
 * @param {boolean=} excludeWrapCoord whether or not to include
 * the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {*} combined value
 */
function coordReduce(layer, callback, memo, excludeWrapCoord) {
    coordEach(layer, function (coord) {
        memo = callback(memo, coord);
    }, excludeWrapCoord);
    return memo;
}
module.exports.coordReduce = coordReduce;

/**
 * Iterate over property objects in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name propEach
 * @param {Object} layer any GeoJSON object
 * @param {Function} callback a method that takes (value)
 * @example
 * var point = { type: 'Feature', geometry: null, properties: { foo: 1 } };
 * turfMeta.propEach(point, function(props) {
 *   // props is equal to { foo: 1}
 * });
 */
function propEach(layer, callback) {
    var i;
    switch (layer.type) {
    case 'FeatureCollection':
        for (i = 0; i < layer.features.length; i++) {
            callback(layer.features[i].properties, i);
        }
        break;
    case 'Feature':
        callback(layer.properties, 0);
        break;
    }
}
module.exports.propEach = propEach;

/**
 * Reduce properties in any GeoJSON object into a single value,
 * similar to how Array.reduce works. However, in this case we lazily run
 * the reduction, so an array of all properties is unnecessary.
 *
 * @name propReduce
 * @param {Object} layer any GeoJSON object
 * @param {Function} callback a method that takes (memo, coord) and returns
 * a new memo
 * @param {*} memo the starting value of memo: can be any type.
 * @returns {*} combined value
 * @example
 * // an example of an even more advanced function that gives you the
 * // javascript type of each property of every feature
 * function propTypes (layer) {
 *   opts = opts || {}
 *   return turfMeta.propReduce(layer, function (prev, props) {
 *     for (var prop in props) {
 *       if (prev[prop]) continue
 *       prev[prop] = typeof props[prop]
 *     }
 *   }, {})
 * }
 */
function propReduce(layer, callback, memo) {
    propEach(layer, function (prop, i) {
        memo = callback(memo, prop, i);
    });
    return memo;
}
module.exports.propReduce = propReduce;

/**
 * Iterate over features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name featureEach
 * @param {Object} layer any GeoJSON object
 * @param {Function} callback a method that takes (value)
 * @example
 * var feature = { type: 'Feature', geometry: null, properties: {} };
 * turfMeta.featureEach(feature, function(feature) {
 *   // feature == feature
 * });
 */
function featureEach(layer, callback) {
    if (layer.type === 'Feature') {
        callback(layer, 0);
    } else if (layer.type === 'FeatureCollection') {
        for (var i = 0; i < layer.features.length; i++) {
            callback(layer.features[i], i);
        }
    }
}
module.exports.featureEach = featureEach;

/**
 * Get all coordinates from any GeoJSON object, returning an array of coordinate
 * arrays.
 *
 * @name coordAll
 * @param {Object} layer any GeoJSON object
 * @returns {Array<Array<number>>} coordinate position array
 */
function coordAll(layer) {
    var coords = [];
    coordEach(layer, function (coord) {
        coords.push(coord);
    });
    return coords;
}
module.exports.coordAll = coordAll;

/**
 * Iterate over each geometry in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name geomEach
 * @param {Object} layer any GeoJSON object
 * @param {Function} callback a method that takes (value)
 * @example
 * var point = {
 *   type: 'Feature',
 *   geometry: { type: 'Point', coordinates: [0, 0] },
 *   properties: {}
 * };
 * turfMeta.geomEach(point, function(geom) {
 *   // geom is the point geometry
 * });
 */
function geomEach(layer, callback) {
    var i, j, g, geometry, stopG,
        geometryMaybeCollection,
        isGeometryCollection,
        isFeatureCollection = layer.type === 'FeatureCollection',
        isFeature = layer.type === 'Feature',
        stop = isFeatureCollection ? layer.features.length : 1;

  // This logic may look a little weird. The reason why it is that way
  // is because it's trying to be fast. GeoJSON supports multiple kinds
  // of objects at its root: FeatureCollection, Features, Geometries.
  // This function has the responsibility of handling all of them, and that
  // means that some of the `for` loops you see below actually just don't apply
  // to certain inputs. For instance, if you give this just a
  // Point geometry, then both loops are short-circuited and all we do
  // is gradually rename the input until it's called 'geometry'.
  //
  // This also aims to allocate as few resources as possible: just a
  // few numbers and booleans, rather than any temporary arrays as would
  // be required with the normalization approach.
    for (i = 0; i < stop; i++) {

        geometryMaybeCollection = (isFeatureCollection ? layer.features[i].geometry :
        (isFeature ? layer.geometry : layer));
        isGeometryCollection = geometryMaybeCollection.type === 'GeometryCollection';
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (g = 0; g < stopG; g++) {
            geometry = isGeometryCollection ?
            geometryMaybeCollection.geometries[g] : geometryMaybeCollection;

            if (geometry.type === 'Point' ||
                geometry.type === 'LineString' ||
                geometry.type === 'MultiPoint' ||
                geometry.type === 'Polygon' ||
                geometry.type === 'MultiLineString' ||
                geometry.type === 'MultiPolygon') {
                callback(geometry);
            } else if (geometry.type === 'GeometryCollection') {
                for (j = 0; j < geometry.geometries.length; j++)
                    callback(geometry.geometries[j]);
            } else {
                throw new Error('Unknown Geometry Type');
            }
        }
    }
}
module.exports.geomEach = geomEach;

},{}],16:[function(require,module,exports){
var bearing = require('@turf/bearing');
var destination = require('@turf/destination');
var distance = require('@turf/distance');

/**
 * Takes two {@link Point|points} and returns a point midway between them.
 * The midpoint is calculated geodesically, meaning the curvature of the earth is taken into account.
 *
 * @name midpoint
 * @param {Feature<Point>} from first point
 * @param {Feature<Point>} to second point
 * @return {Feature<Point>} a point midway between `pt1` and `pt2`
 * @example
 * var pt1 = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [144.834823, -37.771257]
 *   }
 * };
 * var pt2 = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [145.14244, -37.830937]
 *   }
 * };
 *
 * var midpointed = turf.midpoint(pt1, pt2);
 * midpointed.properties['marker-color'] = '#f00';
 *
 *
 * var result = {
 *   "type": "FeatureCollection",
 *   "features": [pt1, pt2, midpointed]
 * };
 *
 * //=result
 */
module.exports = function (from, to) {
    var dist = distance(from, to, 'miles');
    var heading = bearing(from, to);
    var midpoint = destination(from, dist / 2, heading, 'miles');

    return midpoint;
};

},{"@turf/bearing":7,"@turf/destination":10,"@turf/distance":11}],17:[function(require,module,exports){
var distance = require('@turf/distance');

/**
 * Takes a bounding box and calculates the minimum square bounding box that
 * would contain the input.
 *
 * @name square
 * @param {Array<number>} bbox extent in [minX, minY, maxX, maxY] order
 * @return {Array<number>} a square surrounding `bbox`
 * @example
 * var bbox = [-20,-20,-15,0];
 *
 * var squared = turf.square(bbox);
 *
 * var features = {
 *   "type": "FeatureCollection",
 *   "features": [
 *     turf.bboxPolygon(bbox),
 *     turf.bboxPolygon(squared)
 *   ]
 * };
 *
 * //=features
 */
module.exports = function (bbox) {
    var horizontalDistance = distance(bbox.slice(0, 2), [bbox[2], bbox[1]], 'miles');
    var verticalDistance = distance(bbox.slice(0, 2), [bbox[0], bbox[3]], 'miles');
    if (horizontalDistance >= verticalDistance) {
        var verticalMidpoint = (bbox[1] + bbox[3]) / 2;
        return [
            bbox[0],
            verticalMidpoint - ((bbox[2] - bbox[0]) / 2),
            bbox[2],
            verticalMidpoint + ((bbox[2] - bbox[0]) / 2)
        ];
    } else {
        var horizontalMidpoint = (bbox[0] + bbox[2]) / 2;
        return [
            horizontalMidpoint - ((bbox[3] - bbox[1]) / 2),
            bbox[1],
            horizontalMidpoint + ((bbox[3] - bbox[1]) / 2),
            bbox[3]
        ];
    }
};

},{"@turf/distance":11}],18:[function(require,module,exports){
module.exports.RADIUS = 6378137;
module.exports.FLATTENING = 1/298.257223563;
module.exports.POLAR_RADIUS = 6356752.3142;

},{}]},{},[1]);
