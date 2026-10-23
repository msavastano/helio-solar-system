export const J2000 = 2451545.0;
export const AU_KM = 149597870.7;
export const DAY_MS = 86400000;

const DEG = Math.PI / 180;
const RAD = 180 / Math.PI;

// NASA/JPL Approximate Positions of the Planets, Table 1 (1800–2050 AD).
// Earth uses the table's EM Bary row. Pluto is the removed legacy Table 1 row.
export const ORBITAL_ELEMENTS = Object.freeze({
  Mercury: [0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749, 252.25032350, 149472.67411175, 77.45779628, 0.16047689, 48.33076593, -0.12534081],
  Venus: [0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890, 181.97909950, 58517.81538729, 131.60246718, 0.00268329, 76.67984255, -0.27769418],
  Earth: [1.00000261, 0.00000562, 0.01671123, -0.00004392, -0.00001531, -0.01294668, 100.46457166, 35999.37244981, 102.93768193, 0.32327364, 0.0, 0.0],
  Mars: [1.52371034, 0.00001847, 0.09339410, 0.00007882, 1.84969142, -0.00813131, -4.55343205, 19140.30268499, -23.94362959, 0.44441088, 49.55953891, -0.29257343],
  Jupiter: [5.20288700, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714, 34.39644051, 3034.74612775, 14.72847983, 0.21252668, 100.47390909, 0.20469106],
  Saturn: [9.53667594, -0.00125060, 0.05386179, -0.00050991, 2.48599187, 0.00193609, 49.95424423, 1222.49362201, 92.59887831, -0.41897216, 113.66242448, -0.28867794],
  Uranus: [19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939, 313.23810451, 428.48202785, 170.95427630, 0.40805281, 74.01692503, 0.04240589],
  Neptune: [30.06992276, 0.00026291, 0.00859048, 0.00005105, 1.77004347, 0.00035372, -55.12002969, 218.45945325, 44.96476227, -0.32241464, 131.78422574, -0.00508664],
  Pluto: [39.48211675, -0.00031596, 0.24882730, 0.00005170, 17.14001206, 0.00004818, 238.92903833, 145.20780515, 224.06891629, -0.04062942, 110.30393684, -0.01183482],
});

export const BODY_NAMES = Object.freeze(Object.keys(ORBITAL_ELEMENTS));

// NASA/JPL Planetary Satellite Mean Elements, epoch 2000-01-01.5 TDB.
// The Moon's inclination and sidereal period use NASA's Moon Fact Sheet.
export const SATELLITE_ELEMENTS = Object.freeze({
  Moon: { parent: 'Earth', aKm: 384400, e: 0.0554, omega: 318.15, M: 135.27, i: 5.145, node: 125.08, periodDays: 27.3217, frame: 'ecliptic' },
  Io: { parent: 'Jupiter', aKm: 421800, e: 0.004, omega: 49.1, M: 330.9, i: 0.0, node: 0.0, periodDays: 1.762732, frame: 'Laplace' },
  Europa: { parent: 'Jupiter', aKm: 671100, e: 0.009, omega: 45.0, M: 345.4, i: 0.5, node: 184.0, periodDays: 3.525463, frame: 'Laplace' },
  Ganymede: { parent: 'Jupiter', aKm: 1070400, e: 0.001, omega: 198.3, M: 324.8, i: 0.2, node: 58.5, periodDays: 7.155588, frame: 'Laplace' },
  Callisto: { parent: 'Jupiter', aKm: 1882700, e: 0.007, omega: 43.8, M: 87.4, i: 0.3, node: 309.1, periodDays: 16.690440, frame: 'Laplace' },
  Titan: { parent: 'Saturn', aKm: 1221900, e: 0.029, omega: 78.3, M: 11.7, i: 0.3, node: 78.6, periodDays: 15.945448, frame: 'Laplace' },
});

export const SATELLITE_NAMES = Object.freeze(Object.keys(SATELLITE_ELEMENTS));
const SATELLITE_EPOCH = J2000;

export function dateToJulian(date) {
  const milliseconds = date instanceof Date ? date.getTime() : Date.parse(date);
  if (!Number.isFinite(milliseconds)) throw new RangeError('Invalid date');
  return milliseconds / DAY_MS + 2440587.5;
}

export function julianToDate(jd) {
  if (!Number.isFinite(jd)) throw new RangeError('Invalid Julian date');
  return new Date((jd - 2440587.5) * DAY_MS);
}

export function julianCenturies(jd) {
  if (!Number.isFinite(jd)) throw new RangeError('Invalid Julian date');
  return (jd - J2000) / 36525;
}

export function wrap360(degrees) {
  return ((degrees % 360) + 360) % 360;
}

export function wrap180(degrees) {
  return wrap360(degrees + 180) - 180;
}

export function solveKepler(meanAnomaly, eccentricity, tolerance = 1e-6) {
  if (!Number.isFinite(meanAnomaly) || !Number.isFinite(eccentricity) || eccentricity < 0 || eccentricity >= 1) {
    throw new RangeError('Kepler solver requires finite angles and 0 ≤ e < 1');
  }
  if (!(tolerance > 0)) throw new RangeError('Tolerance must be positive');
  const M = wrap180(meanAnomaly);
  const eStar = RAD * eccentricity;
  let E = M + eStar * Math.sin(M * DEG);
  let deltaE = Infinity;
  let iterations = 0;
  while (Math.abs(deltaE) > tolerance && iterations < 50) {
    const deltaM = M - (E - eStar * Math.sin(E * DEG));
    deltaE = deltaM / (1 - eccentricity * Math.cos(E * DEG));
    E += deltaE;
    iterations += 1;
  }
  if (Math.abs(deltaE) > tolerance) throw new Error('Kepler iteration did not converge');
  return E;
}

export function elementsAt(body, jd) {
  const row = ORBITAL_ELEMENTS[body];
  if (!row) throw new RangeError(`Unknown body: ${body}`);
  const T = julianCenturies(jd);
  const [a0, aRate, e0, eRate, i0, iRate, L0, LRate, peri0, periRate, node0, nodeRate] = row;
  const a = a0 + aRate * T;
  const e = e0 + eRate * T;
  const i = i0 + iRate * T;
  const L = L0 + LRate * T;
  const peri = peri0 + periRate * T;
  const node = node0 + nodeRate * T;
  return { body, jd, T, a, e, i, L, peri, node, omega: peri - node, M: wrap180(L - peri) };
}

export function positionFromElements(elements, eccentricAnomaly) {
  const { a, e, i, peri, node } = elements;
  const omega = (peri - node) * DEG;
  const I = i * DEG;
  const O = node * DEG;
  const E = eccentricAnomaly * DEG;
  const xp = a * (Math.cos(E) - e);
  const yp = a * Math.sqrt(1 - e * e) * Math.sin(E);
  const cw = Math.cos(omega);
  const sw = Math.sin(omega);
  const cO = Math.cos(O);
  const sO = Math.sin(O);
  const ci = Math.cos(I);
  const si = Math.sin(I);
  return [
    (cw * cO - sw * sO * ci) * xp + (-sw * cO - cw * sO * ci) * yp,
    (cw * sO + sw * cO * ci) * xp + (-sw * sO + cw * cO * ci) * yp,
    (sw * si) * xp + (cw * si) * yp,
  ];
}

export function vectorLength(vector) {
  return Math.hypot(vector[0], vector[1], vector[2]);
}

export function vectorDistance(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

export function eclipticLongitude(position) {
  return wrap360(Math.atan2(position[1], position[0]) * RAD);
}

export function planetState(body, jd) {
  const elements = elementsAt(body, jd);
  const E = solveKepler(elements.M, elements.e);
  const position = positionFromElements(elements, E);
  return {
    ...elements,
    E,
    position,
    distance: vectorLength(position),
    longitude: eclipticLongitude(position),
    latitude: Math.atan2(position[2], Math.hypot(position[0], position[1])) * RAD,
  };
}

export function sampleOrbit(body, jd, segments = 360) {
  if (!Number.isInteger(segments) || segments < 3) throw new RangeError('Orbit needs at least three segments');
  const elements = elementsAt(body, jd);
  return Array.from({ length: segments }, (_, index) => positionFromElements(elements, index * 360 / segments));
}

export function apsides(body, jd) {
  const elements = elementsAt(body, jd);
  return {
    perihelion: positionFromElements(elements, 0),
    aphelion: positionFromElements(elements, 180),
  };
}

function satelliteElementsAt(name, jd) {
  const satellite = SATELLITE_ELEMENTS[name];
  if (!satellite) throw new RangeError(`Unknown satellite: ${name}`);
  if (!Number.isFinite(jd)) throw new RangeError('Invalid Julian date');
  return {
    ...satellite,
    a: satellite.aKm,
    peri: satellite.node + satellite.omega,
    M: wrap180(satellite.M + 360 * (jd - SATELLITE_EPOCH) / satellite.periodDays),
  };
}

export function satelliteState(name, jd) {
  const elements = satelliteElementsAt(name, jd);
  const E = solveKepler(elements.M, elements.e);
  const position = positionFromElements(elements, E);
  return { ...elements, E, position, distanceKm: vectorLength(position) };
}

export function sampleSatelliteOrbit(name, jd, segments = 180) {
  if (!Number.isInteger(segments) || segments < 3) throw new RangeError('Orbit needs at least three segments');
  const elements = satelliteElementsAt(name, jd);
  return Array.from({ length: segments }, (_, index) => positionFromElements(elements, index * 360 / segments));
}
