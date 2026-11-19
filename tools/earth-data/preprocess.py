"""Build the 0.5-degree Earth datasets used by tectonics.js's present-day Earth preset.

Inputs (in ./raw, see README.md for download commands):
  etopo1_ice_0p1.nc            ETOPO1 Ice Surface, 0.1 deg (NOAA NGDC / Amante & Eakins 2009), via ERDDAP etopo180
  age.2020.1.GTS2012.6m.grd    Seafloor age, 0.1 deg (Seton et al. 2020, EarthByte)
  PB2002_plates.dig.txt        Plate outlines (Bird 2003)
  PB2002_poles.dat.txt         Euler poles relative to the Pacific plate (Bird 2003)

Outputs (in ./out): earth_0p5.npz with every intermediate product, plus PNG previews.
Then run export_js.py to write noncompiled/generators/EarthData.js.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nc3  # noqa: E402

RAW = os.path.join(HERE, 'raw')
OUT = os.path.join(HERE, 'out')
os.makedirs(OUT, exist_ok=True)

RES = 0.5                      # output pixel size, degrees
W, H = int(360 / RES), int(180 / RES)

# pixel centres: row 0 is the northernmost band, column 0 starts at 180 W
pix_lat = 90.0 - (np.arange(H) + 0.5) * RES
pix_lon = -180.0 + (np.arange(W) + 0.5) * RES
LAT, LON = np.meshgrid(pix_lat, pix_lon, indexing='ij')

# exact pixel areas on the unit sphere (steradians)
band_top = np.radians(90.0 - np.arange(H) * RES)
band_bot = np.radians(90.0 - (np.arange(H) + 1) * RES)
band_area = np.radians(RES) * (np.sin(band_top) - np.sin(band_bot))
AREA = np.repeat(band_area[:, None], W, axis=1)


def xyz(lat_deg, lon_deg):
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


# --------------------------------------------------------------------------------------------
# 1. block-average the 0.1 degree grids into 0.5 degree pixels
# --------------------------------------------------------------------------------------------
def block_indices(lat_samples, lon_samples):
    """Map gridline-registered sample coordinates to output pixel rows/cols."""
    rows = np.clip(np.floor((90.0 - lat_samples) / RES).astype(int), 0, H - 1)
    cols = np.clip(np.floor((lon_samples + 180.0) / RES).astype(int), 0, W - 1)
    return rows, cols


def block_mean(values, lat_samples, lon_samples, valid=None):
    """Area-weighted mean of `values` (lat x lon) inside each output pixel."""
    rows, cols = block_indices(lat_samples, lon_samples)
    weight_row = np.cos(np.radians(lat_samples))
    weight_row = np.maximum(weight_row, 1e-6)          # keep polar samples (cos=0) from vanishing
    wgrid = np.repeat(weight_row[:, None], len(lon_samples), axis=1)
    if valid is None:
        valid = np.isfinite(values)
    wv = np.where(valid, wgrid, 0.0)
    vv = np.where(valid, values, 0.0)
    flat = (rows[:, None] * W + cols[None, :]).ravel()
    sum_w = np.bincount(flat, weights=wv.ravel(), minlength=W * H).reshape(H, W)
    sum_v = np.bincount(flat, weights=(vv * wv).ravel(), minlength=W * H).reshape(H, W)
    sum_all = np.bincount(flat, weights=wgrid.ravel(), minlength=W * H).reshape(H, W)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = sum_v / sum_w
    return mean, sum_w / sum_all


print('reading ETOPO1 ...')
_, _, v = nc3.read(os.path.join(RAW, 'etopo1_ice_0p1.nc'))
e_lat = v['latitude'][2].astype(np.float64)
e_lon = v['longitude'][2].astype(np.float64)
alt = v['altitude'][2].astype(np.float64)
alt[alt == 32767] = np.nan
keep = e_lon < 180.0                                     # drop duplicated 180 E column
elev, _ = block_mean(alt[:, keep], e_lat, e_lon[keep])
assert np.isfinite(elev).all(), 'ETOPO gap'
print('  elevation 0.5deg: min %.0f max %.0f, ocean fraction %.3f' %
      (elev.min(), elev.max(), (AREA * (elev < 0)).sum() / AREA.sum()))

print('reading Seton et al. 2020 age grid ...')
_, _, v = nc3.read(os.path.join(RAW, 'age.2020.1.GTS2012.6m.grd'))
a_lat = v['lat'][2].astype(np.float64)
a_lon = v['lon'][2].astype(np.float64)
age = v['z'][2].astype(np.float64)
keep = a_lon < 180.0
age_mean, ocean_frac = block_mean(age[:, keep], a_lat, a_lon[keep])
cont_frac = 1.0 - ocean_frac
print('  oceanic-crust area fraction %.3f' % ((AREA * ocean_frac).sum() / AREA.sum()))

# --------------------------------------------------------------------------------------------
# 2. PB2002 plates
# --------------------------------------------------------------------------------------------
def read_dig(path):
    plates = []
    name, pts = None, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('***'):
                plates.append((name, np.array(pts)))
                name, pts = None, []
            elif name is None:
                name = s
            else:
                lon, lat = (float(t) for t in s.split(','))
                pts.append((lat, lon))
    return plates


def read_poles(path):
    poles = {}
    with open(path) as f:
        for line in f:
            if len(line) < 30:
                continue
            code = line[0:2].strip()
            lat, lon, rate = float(line[2:11]), float(line[11:21]), float(line[21:30])
            poles[code] = (lat, lon, rate)
    return poles


NAMES = {
    'AF': 'Africa', 'AM': 'Amur', 'AN': 'Antarctica', 'AP': 'Altiplano', 'AR': 'Arabia',
    'AS': 'Aegean Sea', 'AT': 'Anatolia', 'AU': 'Australia', 'BH': 'Birds Head', 'BR': 'Balmoral Reef',
    'BS': 'Banda Sea', 'BU': 'Burma', 'CA': 'Caribbean', 'CL': 'Caroline', 'CO': 'Cocos',
    'CR': 'Conway Reef', 'EA': 'Easter', 'EU': 'Eurasia', 'FT': 'Futuna', 'GP': 'Galapagos',
    'IN': 'India', 'JF': 'Juan de Fuca', 'JZ': 'Juan Fernandez', 'KE': 'Kermadec', 'MA': 'Mariana',
    'MN': 'Manus', 'MO': 'Maoke', 'MS': 'Molucca Sea', 'NA': 'North America', 'NB': 'North Bismarck',
    'ND': 'North Andes', 'NH': 'New Hebrides', 'NI': "Niuafo'ou", 'NZ': 'Nazca', 'OK': 'Okhotsk',
    'ON': 'Okinawa', 'PA': 'Pacific', 'PM': 'Panama', 'PS': 'Philippine Sea', 'RI': 'Rivera',
    'SA': 'South America', 'SB': 'South Bismarck', 'SC': 'Scotia', 'SL': 'Shetland', 'SO': 'Somalia',
    'SS': 'Solomon Sea', 'SU': 'Sunda', 'SW': 'Sandwich', 'TI': 'Timor', 'TO': 'Tonga',
    'WL': 'Woodlark', 'YA': 'Yangtze',
}

plates = read_dig(os.path.join(RAW, 'PB2002_plates.dig.txt'))
poles = read_poles(os.path.join(RAW, 'PB2002_poles.dat.txt'))
codes = [p[0] for p in plates]
assert len(codes) == 52 and len(set(codes)) == 52, codes
assert set(codes) == set(poles) == set(NAMES), (set(codes) ^ set(poles), set(codes) ^ set(NAMES))

P = xyz(LAT, LON).reshape(-1, 3)


def turning_sum(V):
    """Sum of signed turning angles of a closed spherical polygon (vertices V, unclosed)."""
    a = V
    b = np.roll(V, -1, axis=0)
    c = np.roll(V, -2, axis=0)
    n1 = np.cross(a, b)
    n2 = np.cross(b, c)
    n1 /= np.linalg.norm(n1, axis=1)[:, None]
    n2 /= np.linalg.norm(n2, axis=1)[:, None]
    return np.arctan2(np.einsum('ij,ij->i', b, np.cross(n1, n2)), np.einsum('ij,ij->i', n1, n2)).sum()


def winding(V, pts):
    """Winding angle of polygon V (unclosed, n x 3) around each point of pts (m x 3)."""
    total = np.zeros(len(pts))
    n = len(V)
    for k in range(n):
        a = V[k]
        b = V[(k + 1) % n]
        axb = np.cross(a, b)
        s = pts @ axb
        c = a @ b - (pts @ a) * (pts @ b)
        total += np.arctan2(s, c)
    return total


# sanity check of the orientation convention on a small CCW ring around the north pole
ring = xyz(np.full(36, 80.0), np.arange(36) * 10.0)
assert turning_sum(ring) > 0, 'orientation test failed'
w_test = winding(ring, xyz(np.array([89.0, 0.0, -89.0]), np.array([0.0, 0.0, 0.0])))
assert abs(w_test[0] - 2 * np.pi) < 1e-6 and abs(w_test[1]) < 1e-6 and abs(w_test[2] + 2 * np.pi) < 1e-6, w_test

print('rasterizing %d plates ...' % len(plates))
inside_count = np.zeros(len(P), dtype=np.int16)
plate_id = np.full(len(P), 255, dtype=np.uint8)
plate_area_exact = {}
for idx, (code, latlon) in enumerate(plates):
    V = xyz(latlon[:, 0], latlon[:, 1])
    # drop closing duplicate and consecutive duplicates
    keep = np.ones(len(V), bool)
    keep[1:] = np.linalg.norm(V[1:] - V[:-1], axis=1) > 1e-12
    V = V[keep]
    if np.linalg.norm(V[0] - V[-1]) < 1e-12:
        V = V[:-1]
    tau = turning_sum(V)
    if tau < 0:                         # small region is on the right: reverse so it is on the left
        V = V[::-1]
        tau = turning_sum(V)
    plate_area_exact[code] = 2 * np.pi - tau
    # bounding cap prefilter
    centre = V.mean(axis=0)
    centre /= np.linalg.norm(centre)
    cos_r = (V @ centre).min()
    cand = np.nonzero(P @ centre >= min(cos_r, 1.0) - 0.02)[0] if cos_r > 0.05 else np.arange(len(P))
    wnd = winding(V, P[cand])
    inside = cand[wnd > np.pi]
    inside_count[inside] += 1
    plate_id[inside] = np.where(plate_id[inside] == 255, idx, plate_id[inside])

n_gap = int((inside_count == 0).sum())
n_overlap = int((inside_count > 1).sum())
print('  pixels with no plate: %d, with >1 plate: %d' % (n_gap, n_overlap))

# fill any gaps from the nearest assigned pixel (in 3D)
if n_gap:
    gaps = np.nonzero(plate_id == 255)[0]
    assigned = np.nonzero(plate_id != 255)[0]
    for g in gaps:
        d = P[assigned] @ P[g]
        plate_id[g] = plate_id[assigned[np.argmax(d)]]
plate_id = plate_id.reshape(H, W)

# area check against the polygons' own spherical areas
total_px_area = 0
print('  %-3s %-16s %9s %9s %6s' % ('id', 'name', 'poly sr', 'pixel sr', 'px'))
for idx, code in enumerate(codes):
    px_area = AREA[plate_id == idx].sum()
    total_px_area += px_area
    print('  %-3s %-16s %9.5f %9.5f %6d' % (code, NAMES[code], plate_area_exact[code], px_area, (plate_id == idx).sum()))
print('  sum of polygon areas %.5f sr (4pi = %.5f)' % (sum(plate_area_exact.values()), 4 * np.pi))

# --------------------------------------------------------------------------------------------
# 3. Euler poles: Pacific-fixed -> no-net-rotation (NNR) frame
# --------------------------------------------------------------------------------------------
def pole_vector(lat, lon, rate_deg_per_myr):
    return xyz(np.array(lat), np.array(lon)) * math.radians(rate_deg_per_myr)


omega_rel = {c: pole_vector(*poles[c]) for c in codes}
Pg = P.reshape(H, W, 3)
M_total = np.zeros((3, 3))
rhs = np.zeros(3)
for idx, code in enumerate(codes):
    sel = plate_id == idx
    p = Pg[sel]
    w = AREA[sel]
    M = np.eye(3) * w.sum() - np.einsum('i,ij,ik->jk', w, p, p)
    M_total += M
    rhs += M @ omega_rel[code]
omega_pa_nnr = -np.linalg.solve(M_total, rhs)


def to_pole(vec):
    rate = np.linalg.norm(vec)
    if rate == 0:
        return 0.0, 0.0, 0.0
    u = vec / rate
    return math.degrees(math.asin(u[2])), math.degrees(math.atan2(u[1], u[0])), math.degrees(rate)


pa = to_pole(omega_pa_nnr)
print('Pacific in NNR frame: lat %.2f lon %.2f rate %.4f deg/Myr' % pa)
print('  (published NNR-NUVEL-1A Pacific: lat -63.0 lon 107.4 rate 0.641 deg/Myr)')
omega_nnr = {c: omega_rel[c] + omega_pa_nnr for c in codes}
# check: net rotation of the whole lithosphere should now be ~0
net = sum((np.eye(3) * AREA[plate_id == i].sum() - np.einsum('i,ij,ik->jk', AREA[plate_id == i], Pg[plate_id == i], Pg[plate_id == i])) @ omega_nnr[c]
          for i, c in enumerate(codes))
print('  residual net rotation |sum M w| = %.3e (rad/Myr sr)' % np.linalg.norm(net))
pole_nnr = {c: to_pole(omega_nnr[c]) for c in codes}

# --------------------------------------------------------------------------------------------
# 4. fit of the crust constants to real depth-age pairs
# --------------------------------------------------------------------------------------------
ocean = (cont_frac < 0.01) & np.isfinite(age_mean)
ages = age_mean[ocean]
depths = elev[ocean]
areas = AREA[ocean]
# Fit only 0-105 Myr: that is where seafloor follows half-space cooling. Older floor flattens and is
# raised by plateaus/seamounts, which a model whose density rises linearly with age cannot express.
bins = np.arange(0, 110, 5)
centres, medians, weights = [], [], []
for lo in bins[:-1]:
    sel = (ages >= lo) & (ages < lo + 5)
    if sel.sum() < 20:
        continue
    order = np.argsort(depths[sel])
    cw = np.cumsum(areas[sel][order])
    med = depths[sel][order][np.searchsorted(cw, cw[-1] / 2)]
    centres.append(lo + 2.5)
    medians.append(med)
    weights.append(areas[sel].sum())
centres, medians, weights = map(np.array, (centres, medians, weights))
c1, c0 = np.polyfit(centres, medians, 1, w=np.sqrt(weights))
print('depth-age fit (median per 5 Myr bin, 0-105 Myr): depth = %.1f + %.2f * age' % (c0, c1))
for a_, m_ in zip(centres[::3], medians[::3]):
    print('   age %5.1f  median depth %7.1f  fit %7.1f' % (a_, m_, c0 + c1 * a_))

np.savez_compressed(os.path.join(OUT, 'earth_0p5.npz'),
                    elev=elev.astype(np.float32), cont_frac=cont_frac.astype(np.float32),
                    age_mean=age_mean.astype(np.float32), plate_id=plate_id, codes=np.array(codes),
                    pole_nnr=np.array([pole_nnr[c] for c in codes]),
                    pole_rel=np.array([poles[c] for c in codes]),
                    fit=np.array([c0, c1]), area=AREA.astype(np.float32))

# previews
try:
    from PIL import Image
    rng = np.random.default_rng(3)
    pal = rng.integers(40, 255, size=(256, 3)).astype(np.uint8)
    Image.fromarray(pal[plate_id]).resize((W * 2, H * 2), Image.NEAREST).save(os.path.join(OUT, 'preview_plates.png'))
    e8 = np.clip((elev + 8000) / 15000 * 255, 0, 255).astype(np.uint8)
    rgb = np.stack([e8, e8, e8], -1)
    rgb[elev < 0] = (rgb[elev < 0] * np.array([0.3, 0.5, 1.0])).astype(np.uint8)
    Image.fromarray(rgb).resize((W * 2, H * 2), Image.NEAREST).save(os.path.join(OUT, 'preview_elev.png'))
    cf = (cont_frac * 255).astype(np.uint8)
    Image.fromarray(cf).resize((W * 2, H * 2), Image.NEAREST).save(os.path.join(OUT, 'preview_continental.png'))
except ImportError:
    pass
print('done')
