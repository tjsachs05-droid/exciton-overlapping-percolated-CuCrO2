"""
cucro2_percolation.py
=====================
Blender script: build the CuCrO2 delafossite structure from a CIF and colour
the A-site sublattice to show a PERCOLATION THRESHOLD, in the style of the
Pd_x Cu_(1-x) CrO2 figure.

Run: Blender -> Scripting -> Open -> edit the USER PARAMETERS block -> Alt+P.

The physics being drawn
-----------------------
In delafossite ABO2 the A-site cations (Cu here) form a TRIANGULAR lattice in
the ab-plane, stacked along c and separated by layers of edge-sharing CrO6
octahedra. Substituting a fraction x of those A sites with a metallic species
(Pd) turns each substituted site into a node of a conducting network. Because
conduction is confined to the A-plane, this is 2D SITE percolation on a
triangular lattice -- for which the threshold is exactly

        p_c = 1/2

which is why the reference figure marks the insulator->metal crossover at 0.5.
This script assigns each A site as metallic with probability METAL_FRACTION,
finds the connected clusters of metallic sites (6 in-plane nearest neighbours
at distance a), reports whether a cluster spans the sample, and can draw the
"metallic channel" network as cylinders between adjacent metallic sites.

COMPOSITION_MODE decides how the substitution is laid out. By default it is
"islands": a mostly-host sheet with ISLAND_COUNT large metallic islands
dropped on it. "gradient" instead ramps x across the sheet, from no
substitution at one end through p_c = 0.5 in the middle to complete
replacement at the other, so one picture holds the insulator, the threshold
and the metal. "uniform" is the old behaviour, one composition everywhere,
and takes COMPOSITION_SERIES for several panels side by side.

The sheet is one A-plane with CRO2_SLABS complete CrO2 slabs around it -- one
above it by default, pulled LAYER_GAP Angstrom clear so the two layers read
separately.

Notes
-----
* 1 Blender unit = 1 Angstrom.
* Colours are given as sRGB hex strings (as sampled from the reference figure)
  and converted to linear internally, which is what Blender's Base Color wants.
* If CIF_PATH does not exist, a built-in CuCrO2 structure is used instead, so
  the script always runs.
* Transparency shows in Material Preview or Rendered viewport shading.
* SETUP_SCENE adds a camera, lights and render settings. BACKGROUND_MODE
  picks between the deep blue field of the reference figure, a flat colour,
  and a transparent RGBA film for layering the render over a slide. The
  camera is a perspective one by default, looking down at CAMERA_ELEVATION.
* LAYER_GAP pulls the Cu planes and the CrO2 slabs apart vertically, and
  SHOW_VERTICAL_BONDS drops the struts that run between them.
"""

import bpy
import colorsys
import math
import os
import random
import re
from math import atan2, cos, hypot, pi, radians, sin, sqrt, floor
from mathutils import Vector

# ============================================================================
#  USER PARAMETERS
# ============================================================================

CIF_PATH = r"C:\path\to\CuCrO2.cif"     # <-- EDIT (falls back to built-in data)

# --- supercell / orientation -------------------------------------------------
N_CELLS        = (16, 16, 1) # unit cells along a, b, c. Keep nz small: c = 17 A.
                             # A composition gradient needs a wide sheet to
                             # read -- with only a handful of sites per slice
                             # the statistics are noise, not a ramp.
VIEW_DIRECTION = "001"       # "001" looks down c (the figure's view); also
                             # accepts "100", "110", "111", etc.
LATTICE_SCALE  = 1.0         # multiplies the CIF lattice parameters

# --- percolation -------------------------------------------------------------
A_SITE_ELEMENT   = "Cu"      # element on the A site that gets substituted
METAL_FRACTION   = 0.50      # x: probability an A site is metallic (Pd-like).
                             # Used only by COMPOSITION_MODE = "uniform".
RANDOM_SEED      = 1         # change for a different random configuration

# --- how the composition is laid out ------------------------------------------
COMPOSITION_MODE      = "islands"
                              # "islands"  : a mostly-host sheet with a few
                              #   large metallic islands in it
                              # "gradient" : x ramps across the sheet, from no
                              #   substitution through p_c = 1/2 to complete
                              #   replacement
                              # "uniform"  : one composition everywhere, from
                              #   METAL_FRACTION or COMPOSITION_SERIES

# --- metallic islands ---------------------------------------------------------
ISLAND_COUNT          = 3     # how many metallic islands to drop on the sheet
ISLAND_RADIUS         = 0.22  # mean island radius, as a fraction of the
                              # sheet's shorter side. Three at 0.22 cover
                              # roughly a third of it, so the host stays the
                              # clear majority.
ISLAND_RADIUS_JITTER  = 0.22  # +/- this fraction on each island's radius
ISLAND_WOBBLE         = 0.35  # how far the outline departs from a circle.
                              # 0 gives discs, which read as drawn-on rather
                              # than grown.
ISLAND_EDGE           = 0.30  # width of the soft rim, as a fraction of the
                              # radius. The probability falls across it
                              # instead of switching, which frays the coast
                              # the way a real substituted alloy does.
ISLAND_FILL           = 0.97  # probability a site well inside an island is
                              # metallic
ISLAND_BACKGROUND     = 0.03  # probability a site out on the host sheet is
                              # metallic anyway -- a light sprinkle, so the
                              # sheet reads as an alloy rather than a mask
ISLAND_SPACING        = 0.85  # keep centres this many combined radii apart,
                              # so the islands stay distinct
ISLAND_INSET          = 0.55  # keep centres this many radii inside the edge

# --- composition gradient (COMPOSITION_MODE = "gradient") ---------------------
GRADIENT_AXIS         = "x"   # "x" or "y": which way the composition ramps
GRADIENT_MIN          = 0.00  # x at the low end (left, for GRADIENT_AXIS "x")
GRADIENT_MID          = 0.50  # x through the middle band
GRADIENT_MAX          = 1.00  # x at the high end (right)
GRADIENT_LEFT_BAND    = 0.22  # fraction of the width held flat at GRADIENT_MIN
GRADIENT_MID_BAND     = 0.20  # fraction held flat at GRADIENT_MID, centred
GRADIENT_RIGHT_BAND   = 0.22  # fraction held flat at GRADIENT_MAX
GRADIENT_SMOOTH       = True  # ease the ramps between bands instead of a
                              # straight line, so the bands blend in
GRADIENT_REPORT_BANDS = 7     # slices the console report breaks the sheet into
                             # (None = different every run)
COMPOSITION_SERIES = []      # e.g. [0.35, 0.50, 1.00] -> one panel per value,
                             # laid out along +X. Empty = single panel at
                             # METAL_FRACTION.
PANEL_GAP        = 10.0      # Angstrom of empty space between panels

CONNECT_INTERLAYER = False   # False: clusters are strictly in-plane (2D, the
                             # physically relevant case for delafossites).
                             # True: also link metallic sites in adjacent
                             # A-planes (3D percolation, p_c is much lower).

CRO2_SLABS = 1               # Complete CrO2 slabs to keep around the
                             # A-plane(s): 0 = the bare A-plane on its own,
                             # 1 = one slab above it (the two-layer figure),
                             # 2 = one either side (a CrO2 | Cu/Pd | CrO2
                             # sandwich).
N_A_PLANES = 1               # How many A-planes to keep, counted from the
                             # bottom. The R-3m hexagonal cell stacks THREE
                             # A-planes per unit cell (ABC), which overlap
                             # when viewed down c -- so 1 gives the clean
                             # single-plane view of the reference figure.
                             # Set to 0 or None to keep the whole slab.

# --- colouring ---------------------------------------------------------------
COLOR_MODE = "species"       # "species"  : metallic vs host colours (the figure)
                             # "spanning" : spanning cluster highlighted, other
                             #              metallic sites dimmed
                             # "clusters" : a distinct colour per cluster

# sRGB hex, sampled from the reference figure. Edit freely.
COLORS = {
    "metal_site":  "#A64DF6",   # metallic A site (Pd in the figure) - purple
    "host_site":   "#E8B04A",   # insulating A site (Cu)             - gold
    "Cr":          "#7474D3",   # B site                             - indigo
    "O":           "#E84C2D",   # oxygen                             - red
    "bond":        "#B9B9C2",   # chemical bonds (Cu-O, Cr-O)        - grey
    "channel":     "#A64DF6",   # metallic-channel network           - purple
    "dimmed":      "#C9A8E8",   # non-spanning metallic sites in "spanning" mode
    "guide":       "#BFBFBF",   # faint full A-lattice guide network
}

# Display radii in Angstrom (not physical radii -- these set the look).
RADII = {
    "metal_site": 0.62,
    "host_site":  0.62,
    "Cr":         0.50,
    "O":          0.26,
}
DEFAULT_RADIUS = 0.45

# --- layer stacking ----------------------------------------------------------
LAYER_GAP        = 7.0       # extra vertical distance (A) inserted between the
                             # Cu (A) planes and the CrO2 slabs on either side
                             # of them. 0 = the true crystal spacing; a couple
                             # of Angstrom pulls the layers apart so the
                             # sandwich is legible in a figure. It is a display
                             # exaggeration only: bonds, neighbours and the
                             # percolation analysis are all computed on the
                             # real geometry first, so the bonds simply stretch.
SLAB_TOL         = 1.4       # (A) height gap that separates one O-Cr-O slab
                             # from the next layer. Between 1.0 and 1.8 for
                             # CuCrO2; only used to group atoms into layers.

# --- what to draw ------------------------------------------------------------
SHOW_ATOMS        = True
SHOW_BONDS        = True     # chemical Cu-O / Cr-O bonds. With the vertical
                            # ones excluded below, what is left is the CrO6
                            # octahedral network inside the CrO2 slab.
SHOW_VERTICAL_BONDS = False  # False drops every bond that runs between layers
                             # (the vertical Cu-O struts, and the interlayer
                             # rungs of the channel network when
                             # CONNECT_INTERLAYER is on), leaving only the
                             # in-plane connectivity
VERTICAL_BOND_ANGLE = 40.0   # degrees: a bond counts as "vertical" when it
                             # sits within this angle of the stacking axis
SHOW_CHANNELS     = True     # purple network linking adjacent metallic sites
SHOW_LATTICE_GUIDE = False   # faint network over ALL A sites (figure's dashes)
SHOW_A_SITES      = True
SHOW_B_SITES      = True     # Cr
SHOW_OXYGEN       = True

# --- transparency (0.0 = opaque, 1.0 = invisible) ----------------------------
TRANSPARENCY = {
    "atoms":   0.0,
    "bonds":   0.0,
    "channel": 0.0,
    "guide":   0.55,
}

# --- geometry / quality ------------------------------------------------------
BOND_RADIUS      = 0.07      # chemical bond cylinder radius
CHANNEL_RADIUS   = 0.11      # metallic-channel cylinder radius
GUIDE_RADIUS     = 0.03
BOND_FACTOR      = 1.15      # bond if d < BOND_FACTOR * (r_cov_i + r_cov_j)
MAX_BOND_LENGTH  = 2.4       # hard cutoff for chemical bonds (A)
MIN_BOND_LENGTH  = 0.40
NEIGHBOR_TOL     = 0.25      # tolerance (A) for "in-plane nearest neighbour"
PLANE_TOL        = 0.30      # tolerance (A) for grouping atoms into a plane
SPHERE_SEGMENTS  = 24        # a 16x16 sheet is ~1800 atoms, so the mesh cost
SPHERE_RINGS     = 16        # per sphere matters; 24x16 is still smooth at
                             # the size these render
CYL_SEGMENTS     = 16
CLEAR_SCENE      = True

# --- scene / render ----------------------------------------------------------
SETUP_SCENE            = True   # camera, lights and render settings
BACKGROUND_MODE        = "gradient"
                                # "gradient"    -> the deep blue field of the
                                #   reference image, dark at the bottom and
                                #   brighter towards the top
                                # "flat"        -> a single BACKGROUND_COLOR
                                # "transparent" -> alpha channel, for layering
                                #   this render over something else. This is
                                #   the base layer of the pair, so it carries
                                #   the blue and exciton.py stays transparent
                                #   on top of it; flip both to "transparent"
                                #   instead if the slide itself is blue.
BACKGROUND_COLOR       = "#0A2A5E"  # deep blue, sampled from Reference 1
BACKGROUND_TOP_COLOR   = "#12539E"  # brighter blue at the top of the gradient
BACKGROUND_MARGIN      = 1.15   # how far the backdrop oversizes the frame
BACKGROUND_STEPS       = 48     # rows the gradient is built from
BACKGROUND_SMOOTH      = True   # ease the gradient instead of ramping it
ORTHOGRAPHIC           = False  # False gives a normal perspective camera,
                                # which is what makes the three-layer sandwich
                                # read as a solid rather than a flat pattern
CAMERA_LENS            = 50.0   # focal length in mm, for the perspective camera
CAMERA_ELEVATION       = 55.0   # degrees above the plane of the sheet. 90 is
                                # straight down (a plan view, where the
                                # sandwich is edge-on and invisible); lower
                                # angles show its thickness.
CAMERA_AZIMUTH         = -90.0  # degrees around the sheet, measured from +x.
                                # -90 puts the camera on the -y side, which
                                # lays the GRADIENT_AXIS left-to-right across
                                # the frame.
CAMERA_MARGIN          = 1.18   # >1 leaves air around the structure. The
                                # camera fits the geometry exactly, so this is
                                # pure breathing room -- at 1.0 the outermost
                                # atoms sit right on the frame edge.
AMBIENT_STRENGTH       = 0.35   # world light. It only lights the atoms; the
                                # world itself stays invisible on a
                                # transparent film.
ADD_LIGHTS             = True   # key/fill/rim suns, so the atoms are shaded
                                # spheres rather than flat silhouettes
KEY_LIGHT_ENERGY       = 4.0
RENDER_ENGINE          = "EEVEE"    # "EEVEE" or "CYCLES"
RESOLUTION             = (2600, 2000)
RENDER_SAMPLES         = 128
VIEW_TRANSFORM         = "Standard"   # "Standard" keeps the colours exactly
                                      # as sampled above, which is what a
                                      # structure figure wants; "AgX" or
                                      # "Filmic" tone-map them photographically
OUTPUT_PATH            = ""     # e.g. r"C:\figures\percolation.png"
RENDER_NOW             = False  # True = render straight to OUTPUT_PATH

# ============================================================================
#  Built-in CuCrO2 structure (used if CIF_PATH is missing)
# ============================================================================
#  Poienar et al., J. Solid State Chem. 185 (2012) 56: R-3m, a = 2.9747 A,
#  c = 17.1038 A. O 6c z = 0.1100 reproduces Cu-O = 1.881 A, Cr-O = 1.972 A.

BUILTIN_CELL  = (2.9747, 2.9747, 17.1038, 90.0, 90.0, 120.0)
BUILTIN_SITES = [("Cu", 0.0, 0.0, 0.0),
                 ("Cr", 0.0, 0.0, 0.5),
                 ("O",  0.0, 0.0, 0.11)]
BUILTIN_SG    = 166                      # R-3m

# covalent radii, only used for the chemical-bond distance criterion
COVALENT_RADII = {"Cu": 1.32, "Cr": 1.39, "O": 0.66, "Pd": 1.39, "Ag": 1.45,
                  "Pt": 1.36, "Al": 1.21, "Fe": 1.32, "Ga": 1.22, "Rh": 1.42,
                  "In": 1.42, "Sc": 1.70, "Y": 1.90, "La": 2.07}
DEFAULT_COVALENT = 1.40

CIF_UNKNOWN = ("?", ".")


def hex_to_linear(h):
    """sRGB hex string -> linear RGB triple (what Blender Base Color expects)."""
    h = h.strip().lstrip("#")
    if len(h) == 3:
        h = "".join(ch * 2 for ch in h)
    srgb = [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]

    def lin(u):
        return u / 12.92 if u <= 0.04045 else ((u + 0.055) / 1.055) ** 2.4
    return tuple(lin(u) for u in srgb)


def covalent_radius(elem):
    return COVALENT_RADII.get(elem, DEFAULT_COVALENT)


# ============================================================================
#  Space-group symmetry: generators -> full group by closure
# ============================================================================

_CENTERING = {
    "P": [(0.0, 0.0, 0.0)],
    "C": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)],
    "I": [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    "F": [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    "R": [(0.0, 0.0, 0.0), (2/3, 1/3, 1/3), (1/3, 2/3, 2/3)],
}

_SG_GENERATORS = {
    1:   ([], "P"),
    2:   (["-x,-y,-z"], "P"),
    12:  (["-x,y,-z", "-x,-y,-z"], "C"),
    14:  (["-x,y+1/2,-z+1/2", "-x,-y,-z"], "P"),
    15:  (["-x,y,-z+1/2", "-x,-y,-z"], "C"),
    62:  (["-x+1/2,-y,z+1/2", "-x,y+1/2,-z", "-x,-y,-z"], "P"),
    123: (["-y,x,z", "x,-y,z", "-x,-y,-z"], "P"),
    136: (["-y+1/2,x+1/2,z+1/2", "y,x,-z", "-x,-y,-z"], "P"),
    139: (["-y,x,z", "x,-y,z", "-x,-y,-z"], "I"),
    148: (["-y,x-y,z", "-x,-y,-z"], "R"),
    160: (["-y,x-y,z", "-y,-x,z"], "R"),
    166: (["-y,x-y,z", "y,x,-z", "-x,-y,-z"], "R"),      # R-3m : CuCrO2
    167: (["-y,x-y,z", "y,x,-z+1/2", "-x,-y,-z"], "R"),
    186: (["-y,x-y,z", "-x,-y,z+1/2", "-y,-x,z"], "P"),
    194: (["-y,x-y,z+1/2", "y,x,-z", "-x,-y,-z"], "P"),
    221: (["-y,x,z", "z,x,y", "-x,-y,-z"], "P"),
    225: (["-y,x,z", "z,x,y", "-x,-y,-z"], "F"),
    229: (["-y,x,z", "z,x,y", "-x,-y,-z"], "I"),
}

_HM_TO_NUMBER = {
    "p1": 1, "p-1": 2, "c2/m": 12, "p2_1/c": 14, "p21/c": 14, "c2/c": 15,
    "pnma": 62, "p4/mmm": 123, "p4_2/mnm": 136, "p42/mnm": 136, "i4/mmm": 139,
    "r-3": 148, "r3m": 160, "r-3m": 166, "r-3c": 167, "p6_3mc": 186,
    "p63mc": 186, "p6_3/mmc": 194, "p63/mmc": 194, "pm-3m": 221, "pm3m": 221,
    "fm-3m": 225, "fm3m": 225, "im-3m": 229, "im3m": 229,
}

_IDENTITY_OP = ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0])


def _op_from_string(op):
    comps = op.lower().replace(" ", "").split(",")
    if len(comps) != 3:
        raise ValueError(f"Bad symmetry operation: {op!r}")

    def ev(expr, x, y, z):
        return eval(expr, {"__builtins__": {}}, {"x": x, "y": y, "z": z})

    t = [float(ev(c, 0, 0, 0)) for c in comps]
    R = [[float(ev(c, *e)) - t[k] for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
         for k, c in enumerate(comps)]
    return R, t


def _compose(A, B):
    RA, tA = A
    RB, tB = B
    R = [[sum(RA[i][k] * RB[k][j] for k in range(3)) for j in range(3)]
         for i in range(3)]
    t = [(sum(RA[i][k] * tB[k] for k in range(3)) + tA[i]) % 1.0 for i in range(3)]
    return R, t


def _op_key(op, grid=24):
    R, t = op
    return (tuple(int(round(x)) for row in R for x in row),
            tuple(int(round((v % 1.0) * grid)) % grid for v in t))


def _generate_group(generator_strings, centering):
    gens = [_IDENTITY_OP] + [_op_from_string(s) for s in generator_strings]
    ops = {_op_key(g): g for g in gens}
    changed = True
    while changed and len(ops) < 400:
        changed = False
        for a in list(ops.values()):
            for g in gens:
                c = _compose(a, g)
                k = _op_key(c)
                if k not in ops:
                    ops[k] = c
                    changed = True
    full = {}
    for R, t in ops.values():
        for cv in _CENTERING[centering]:
            nt = [(t[i] + cv[i]) % 1.0 for i in range(3)]
            full[_op_key((R, nt))] = (R, nt)
    return list(full.values())


def apply_symmetry(sites, symops, tol=1e-3):
    out = []
    for elem, fx, fy, fz in sites:
        f = (fx, fy, fz)
        for R, t in symops:
            p = [(R[k][0] * f[0] + R[k][1] * f[1] + R[k][2] * f[2] + t[k]) % 1.0
                 for k in range(3)]
            p = [0.0 if abs(v - 1.0) < tol or abs(v) < tol else v for v in p]
            if not any(e2 == elem and
                       all(min(abs(p[k] - q[k]), 1.0 - abs(p[k] - q[k])) < tol
                           for k in range(3))
                       for e2, q in out):
                out.append((elem, p))
    return out


# ============================================================================
#  CIF parsing (multi-block aware, '?' tolerant)
# ============================================================================

_NUM_RE = re.compile(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_CELL_KEYS = ("_cell_length_a", "_cell_length_b", "_cell_length_c",
              "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")


def _num(s):
    m = _NUM_RE.match(s.strip())
    if not m:
        raise ValueError(f"Cannot parse number from CIF value: {s!r}")
    return float(m.group(0))


def _tokenize(line):
    tokens, i, n = [], 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c in "'\"":
            q, j = c, i + 1
            while j < n and not (line[j] == q and (j + 1 == n or line[j + 1] in " \t")):
                j += 1
            tokens.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and line[j] not in " \t":
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


def _new_block():
    return {"scalars": {}, "symops": [], "sites": []}


def _process_loop(headers, rows, block):
    for tag in ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"):
        if tag in headers:
            col = headers.index(tag)
            for r in rows:
                if r[col].strip() not in CIF_UNKNOWN:
                    block["symops"].append(r[col])
            return
    if "_atom_site_fract_x" in headers:
        ix = headers.index("_atom_site_fract_x")
        iy = headers.index("_atom_site_fract_y")
        iz = headers.index("_atom_site_fract_z")
        it = headers.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in headers else None
        il = headers.index("_atom_site_label") if "_atom_site_label" in headers else None
        for r in rows:
            coords = (r[ix].strip(), r[iy].strip(), r[iz].strip())
            if any(v in CIF_UNKNOWN for v in coords):
                continue
            raw = (r[it] if it is not None else (r[il] if il is not None else "X")).strip()
            if raw in CIF_UNKNOWN:
                continue
            m = re.match(r"([A-Za-z]{1,2})", raw)
            elem = m.group(1).capitalize() if m else "X"
            block["sites"].append((elem, _num(coords[0]), _num(coords[1]), _num(coords[2])))


def _resolve_symops(block, label):
    if block["symops"]:
        ops = []
        for s in block["symops"]:
            try:
                ops.append(_op_from_string(s))
            except Exception:
                print(f"[cucro2] WARNING: skipping unparsable symop {s!r}")
        if ops:
            return ops
    sc = block["scalars"]
    name = sc.get("_symmetry_space_group_name_h-m",
                  sc.get("_space_group_name_h-m_alt", "")).replace(" ", "").lower()
    num_str = sc.get("_symmetry_int_tables_number",
                     sc.get("_space_group_it_number", "")).strip()
    number = int(num_str) if num_str.isdigit() else _HM_TO_NUMBER.get(name)
    if number in _SG_GENERATORS:
        gens, centering = _SG_GENERATORS[number]
        ops = _generate_group(gens, centering)
        print(f"[cucro2] {label}: no symop loop; generated {len(ops)} operations "
              f"for space group {name or number} (No. {number}).")
        return ops
    print(f"[cucro2] WARNING: {label}: no symmetry operations available "
          f"(space group '{name or num_str or 'unknown'}'). Using P1 -- the "
          "structure will look incomplete. Re-save the CIF from VESTA to get "
          "a full symop loop.")
    return [_IDENTITY_OP]


def parse_cif(path):
    with open(path, "r", errors="replace") as f:
        lines = [ln.rstrip("\n") for ln in f]

    blocks = [_new_block()]
    cur = blocks[0]
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        low = line.lower()
        if low.startswith("data_"):
            cur = _new_block()
            blocks.append(cur)
            i += 1
            continue
        if line.startswith(";"):
            i += 1
            while i < n and not lines[i].strip().startswith(";"):
                i += 1
            i += 1
            continue
        if low.startswith("loop_"):
            i += 1
            headers = []
            while i < n and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip().split()[0].lower())
                i += 1
            values = []
            while i < n:
                s = lines[i].strip()
                if not s or s.startswith("#"):
                    i += 1
                    if not s and values:
                        break
                    continue
                if s.startswith("_") or s.lower().startswith(("loop_", "data_")):
                    break
                if s.startswith(";"):
                    buf = [s[1:]]
                    i += 1
                    while i < n and not lines[i].strip().startswith(";"):
                        buf.append(lines[i])
                        i += 1
                    i += 1
                    values.append(" ".join(buf).strip())
                    continue
                values.extend(_tokenize(s))
                i += 1
            if headers:
                ncol = len(headers)
                rows = [values[k:k + ncol] for k in range(0, len(values), ncol)
                        if len(values[k:k + ncol]) == ncol]
                _process_loop(headers, rows, cur)
            continue
        if line.startswith("_"):
            toks = _tokenize(line)
            if len(toks) >= 2:
                cur["scalars"][toks[0].lower()] = toks[1]
                i += 1
            else:
                key = toks[0].lower()
                i += 1
                if i < n and lines[i].strip().startswith(";"):
                    i += 1
                    while i < n and not lines[i].strip().startswith(";"):
                        i += 1
                    i += 1
                elif i < n:
                    cur["scalars"][key] = lines[i].strip()
                    i += 1
            continue
        i += 1

    def has_cell(b):
        return all(k in b["scalars"] and b["scalars"][k].strip() not in CIF_UNKNOWN
                   for k in _CELL_KEYS)

    chosen = next((b for b in blocks if has_cell(b) and b["sites"]), None)
    if chosen is None:
        raise ValueError(f"{path}: no data block with both a cell and atom sites")
    cell = tuple(_num(chosen["scalars"][k]) for k in _CELL_KEYS)
    return {"cell": cell,
            "symops": _resolve_symops(chosen, os.path.basename(path)),
            "sites": chosen["sites"]}


def load_structure():
    """CIF if available, else the built-in CuCrO2 data."""
    if os.path.isfile(CIF_PATH):
        print(f"[cucro2] Reading {CIF_PATH}")
        return parse_cif(CIF_PATH)
    print(f"[cucro2] CIF not found at {CIF_PATH} -- using built-in CuCrO2 "
          "(R-3m, a = 2.9747 A, c = 17.1038 A, O z = 0.1100).")
    gens, centering = _SG_GENERATORS[BUILTIN_SG]
    return {"cell": BUILTIN_CELL,
            "symops": _generate_group(gens, centering),
            "sites": BUILTIN_SITES}


# ============================================================================
#  Geometry
# ============================================================================

def lattice_vectors(a, b, c, alpha, beta, gamma):
    al, be, ga = radians(alpha), radians(beta), radians(gamma)
    va = Vector((a, 0.0, 0.0))
    vb = Vector((b * cos(ga), b * sin(ga), 0.0))
    cx = c * cos(be)
    cy = c * (cos(al) - cos(be) * cos(ga)) / sin(ga)
    cz = sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return va, vb, Vector((cx, cy, cz))


def view_quaternion(view, va, vb, vc):
    try:
        u, v, w = (int(ch) for ch in view.strip())
    except Exception:
        raise ValueError(f"VIEW_DIRECTION must be a 3-index string like "
                         f"'001', '100', '110', '111'; got {view!r}")
    d = u * va + v * vb + w * vc
    if d.length < 1e-9:
        raise ValueError(f"Degenerate view direction {view!r}")
    return d.normalized().rotation_difference(Vector((0.0, 0.0, 1.0)))


def stacking_axis(q, va, vb):
    """
    Unit normal of the A-site planes, in the rotated (view) frame.

    The layers are perpendicular to c*, i.e. to a x b -- which for a hexagonal
    cell is c itself. For the default VIEW_DIRECTION = "001" this comes out as
    +Z, which is what the rest of the script's plane bookkeeping assumes.
    """
    n = va.cross(vb)
    if n.length < 1e-9:
        return Vector((0.0, 0.0, 1.0))
    n = (q @ n.normalized())
    return n if n.z >= 0.0 else -n


def build_supercell(struct):
    """Return list of [elem, Vector pos], centred in XY and resting at z >= 0."""
    a, b, c, al, be, ga = struct["cell"]
    s = LATTICE_SCALE
    va, vb, vc = lattice_vectors(a * s, b * s, c * s, al, be, ga)
    frac = apply_symmetry(struct["sites"], struct["symops"])
    nx, ny, nz = N_CELLS
    q = view_quaternion(VIEW_DIRECTION, va, vb, vc)

    atoms = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                for elem, (fx, fy, fz) in frac:
                    cart = (ix + fx) * va + (iy + fy) * vb + (iz + fz) * vc
                    atoms.append([elem, q @ cart])
    atoms = _crop_to_layers(atoms)

    cx = sum(p[1].x for p in atoms) / len(atoms)
    cy = sum(p[1].y for p in atoms) / len(atoms)
    zmin = min(p[1].z for p in atoms)
    for p in atoms:
        p[1] = Vector((p[1].x - cx, p[1].y - cy, p[1].z - zmin))
    return atoms, (va, vb, vc), stacking_axis(q, va, vb)


def _ease(u):
    u = min(max(u, 0.0), 1.0)
    return u * u * (3.0 - 2.0 * u) if GRADIENT_SMOOTH else u


def gradient_fraction(t):
    """
    The metallic fraction at normalised position t across the sheet: 0 at the
    low end of GRADIENT_AXIS, 1 at the high end.

    Three flat bands -- GRADIENT_MIN on the left, GRADIENT_MID through the
    middle, GRADIENT_MAX on the right -- joined by ramps. The middle band is
    the interesting one: at 0.5 it sits exactly on the 2D triangular-lattice
    site-percolation threshold, so one sheet shows the insulator, the
    threshold and the metal side by side.
    """
    t = min(max(t, 0.0), 1.0)
    left = min(max(GRADIENT_LEFT_BAND, 0.0), 1.0)
    right = min(max(GRADIENT_RIGHT_BAND, 0.0), 1.0)
    mid = min(max(GRADIENT_MID_BAND, 0.0), 1.0)
    mid_lo, mid_hi = 0.5 - 0.5 * mid, 0.5 + 0.5 * mid

    if left > mid_lo or 1.0 - right < mid_hi:
        raise ValueError(
            "the gradient bands overlap: GRADIENT_LEFT_BAND and "
            "GRADIENT_RIGHT_BAND must each leave room for half of "
            f"GRADIENT_MID_BAND (left {left}, mid {mid}, right {right})")

    if t <= left:
        return GRADIENT_MIN
    if t >= 1.0 - right:
        return GRADIENT_MAX
    if mid_lo <= t <= mid_hi:
        return GRADIENT_MID
    if t < mid_lo:
        u = (t - left) / max(mid_lo - left, 1e-9)
        return GRADIENT_MIN + (GRADIENT_MID - GRADIENT_MIN) * _ease(u)
    u = (t - mid_hi) / max((1.0 - right) - mid_hi, 1e-9)
    return GRADIENT_MID + (GRADIENT_MAX - GRADIENT_MID) * _ease(u)


def _sheet_bounds(a_positions):
    xs = [p.x for p in a_positions] or [0.0]
    ys = [p.y for p in a_positions] or [0.0]
    return min(xs), max(xs), min(ys), max(ys)


def island_targets(a_positions, rng):
    """
    Scatter ISLAND_COUNT metallic islands over an otherwise host sheet.

    An island is a blob rather than a disc: its radius is modulated by a few
    harmonics of the polar angle, so the outline is irregular, and the rim is
    soft -- the probability falls from ISLAND_FILL to ISLAND_BACKGROUND across
    a band ISLAND_EDGE wide. Between them those two frays the coastline into
    something a substituted alloy might actually produce, instead of a circle
    someone drew on.

    Returns (per-site probability, [(cx, cy, radius, harmonics)]).
    """
    n = max(int(ISLAND_COUNT), 0)
    x0, x1, y0, y1 = _sheet_bounds(a_positions)
    base = ISLAND_RADIUS * min(x1 - x0, y1 - y0)

    islands = []
    for _ in range(n):
        r = max(base * (1.0 + ISLAND_RADIUS_JITTER * (rng.random() * 2.0 - 1.0)),
                1e-3)
        harmonics = [(rng.uniform(0.4, 1.0), rng.uniform(0.0, 2.0 * pi), k)
                     for k in (2, 3, 5)]
        inset = r * ISLAND_INSET
        lox, hix = x0 + inset, x1 - inset
        loy, hiy = y0 + inset, y1 - inset
        if lox > hix:
            lox = hix = 0.5 * (x0 + x1)
        if loy > hiy:
            loy = hiy = 0.5 * (y0 + y1)
        cx, cy = 0.5 * (lox + hix), 0.5 * (loy + hiy)
        for _attempt in range(200):
            cx, cy = rng.uniform(lox, hix), rng.uniform(loy, hiy)
            if all(hypot(cx - ox, cy - oy) >= (r + orad) * ISLAND_SPACING
                   for ox, oy, orad, _ in islands):
                break
        islands.append((cx, cy, r, harmonics))

    def blob_radius(island, angle):
        _, _, r, harmonics = island
        m = 1.0
        for amp, phase, k in harmonics:
            m += ISLAND_WOBBLE * amp * cos(k * angle + phase) / len(harmonics)
        return r * max(m, 0.25)

    lo_p = min(max(ISLAND_BACKGROUND, 0.0), 1.0)
    hi_p = min(max(ISLAND_FILL, 0.0), 1.0)
    targets = []
    for p in a_positions:
        best = lo_p
        for island in islands:
            cx, cy, _, _ = island
            dx, dy = p.x - cx, p.y - cy
            rim = blob_radius(island, atan2(dy, dx))
            edge = max(ISLAND_EDGE * rim, 1e-6)
            # 1 just inside the rim, 0 just outside it
            u = (rim + 0.5 * edge - hypot(dx, dy)) / edge
            best = max(best, lo_p + (hi_p - lo_p) * _ease_clamped(u))
        targets.append(best)
    return targets, islands


def _ease_clamped(u):
    u = min(max(u, 0.0), 1.0)
    return u * u * (3.0 - 2.0 * u)


def site_targets(a_positions, x_metal, rng):
    """
    Target metallic fraction for every A site, plus whatever the composition
    mode wants to report about itself.
    """
    if COMPOSITION_MODE == "islands":
        targets, islands = island_targets(a_positions, rng)
        return targets, {"islands": islands}

    if COMPOSITION_MODE == "gradient":
        if GRADIENT_AXIS not in ("x", "y"):
            raise ValueError("GRADIENT_AXIS must be 'x' or 'y'; got "
                             f"{GRADIENT_AXIS!r}")
        vals = [getattr(p, GRADIENT_AXIS) for p in a_positions]
        lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
        span = max(hi - lo, 1e-9)
        where = [(v - lo) / span for v in vals]
        return [gradient_fraction(t) for t in where], {"where": where}

    if COMPOSITION_MODE != "uniform":
        raise ValueError("COMPOSITION_MODE must be 'islands', 'gradient' or "
                         f"'uniform'; got {COMPOSITION_MODE!r}")
    frac = METAL_FRACTION if x_metal is None else x_metal
    return [frac] * len(a_positions), {}


def _crop_to_layers(atoms):
    """
    Keep N_A_PLANES A-planes plus CRO2_SLABS complete CrO2 slabs around them:
    none, one above, or one either side.

    The cuts run just inside the neighbouring A-planes rather than midway to
    them. An O-Cr-O slab sits centred between two A-planes, so a cut at the
    midpoint would slice the slab in half and leave one of its oxygen sheets
    behind.
    """
    if not N_A_PLANES:
        return atoms
    levels = sorted({round(p.z / PLANE_TOL) * PLANE_TOL
                     for e, p in atoms if e == A_SITE_ELEMENT})
    keep = max(int(N_A_PLANES), 1)
    if len(levels) < keep:
        return atoms
    slabs = min(max(int(CRO2_SLABS), 0), 2)

    # with a slab wanted on both sides the kept planes have to come from the
    # middle of the stack, so there is a neighbour to cut against either way
    start = max((len(levels) - keep) // 2, 0) if slabs >= 2 else 0
    window = levels[start:start + keep]
    gaps = [levels[i + 1] - levels[i] for i in range(len(levels) - 1)]
    spacing = min(gaps) if gaps else PLANE_TOL * 2

    has_below, has_above = start > 0, start + keep < len(levels)
    if (slabs >= 2 and not has_below) or (slabs >= 1 and not has_above):
        print("[cucro2] WARNING: not enough A-planes to take "
              f"{slabs} complete CrO2 slab(s) -- raise N_CELLS[2].")

    below = levels[start - 1] if has_below else window[0] - spacing
    above = levels[start + keep] if has_above else window[-1] + spacing
    lo = below + PLANE_TOL if slabs >= 2 else window[0] - PLANE_TOL
    hi = above - PLANE_TOL if slabs >= 1 else window[-1] + PLANE_TOL
    return [p for p in atoms if lo <= p[1].z <= hi]


def layer_groups(atoms, stack_axis):
    """
    Split the atoms into the alternating layers of the delafossite stack: flat
    A-site (Cu) planes, and the O-Cr-O slabs sandwiched between them.

    A-site atoms are grouped by height to within PLANE_TOL. Everything else is
    sorted by height and cut wherever the gap to the next atom exceeds
    SLAB_TOL, which is what separates one O-Cr-O slab from the next.

    Returns [(mean height, [atom indices])], ordered bottom to top.
    """
    s = [stack_axis.dot(p) for _, p in atoms]
    groups = {}
    rest = []
    for i, (elem, _) in enumerate(atoms):
        if elem == A_SITE_ELEMENT:
            groups.setdefault(("A", round(s[i] / PLANE_TOL)), []).append(i)
        else:
            rest.append(i)

    rest.sort(key=lambda i: s[i])
    slab, last, n = [], None, 0
    for i in rest:
        if last is not None and s[i] - last > SLAB_TOL:
            groups[("B", n)] = slab
            slab, n = [], n + 1
        slab.append(i)
        last = s[i]
    if slab:
        groups[("B", n)] = slab

    out = [(sum(s[i] for i in ids) / len(ids), ids) for ids in groups.values()]
    out.sort(key=lambda t: t[0])
    return out


def layer_offsets(atoms, stack_axis):
    """
    Per-atom displacement that pulls consecutive layers LAYER_GAP further
    apart along the stacking axis, keeping the stack centred on where it was.
    """
    zero = Vector((0.0, 0.0, 0.0))
    if abs(LAYER_GAP) < 1e-9:
        return [zero] * len(atoms)
    groups = layer_groups(atoms, stack_axis)
    n = len(groups)
    off = [zero] * len(atoms)
    for k, (_, ids) in enumerate(groups):
        d = stack_axis * (LAYER_GAP * (k - (n - 1) * 0.5))
        for i in ids:
            off[i] = d
    return off


def is_vertical(p1, p2, stack_axis):
    """True when the segment runs between layers rather than within one."""
    d = p2 - p1
    if d.length < 1e-9:
        return False
    return abs(d.dot(stack_axis)) / d.length > cos(radians(VERTICAL_BOND_ANGLE))


def drop_vertical(segments, stack_axis):
    """Filter out interlayer segments when SHOW_VERTICAL_BONDS is off."""
    if SHOW_VERTICAL_BONDS:
        return segments
    return [(p, q) for p, q in segments if not is_vertical(p, q, stack_axis)]


def chemical_bonds(atoms):
    """Distance-based bonds using a spatial hash grid."""
    cell = MAX_BOND_LENGTH
    grid = {}
    for i, (_, p) in enumerate(atoms):
        grid.setdefault((floor(p.x / cell), floor(p.y / cell), floor(p.z / cell)),
                        []).append(i)
    bonds = []
    for i, (ei, pi) in enumerate(atoms):
        ri = covalent_radius(ei)
        kx, ky, kz = floor(pi.x / cell), floor(pi.y / cell), floor(pi.z / cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((kx + dx, ky + dy, kz + dz), ()):
                        if j <= i:
                            continue
                        ej, pj = atoms[j]
                        d = (pj - pi).length
                        cutoff = min(BOND_FACTOR * (ri + covalent_radius(ej)),
                                     MAX_BOND_LENGTH)
                        if MIN_BOND_LENGTH < d <= cutoff:
                            bonds.append((i, j))
    return bonds


# ============================================================================
#  Percolation on the A-site sublattice
# ============================================================================

class _DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, i):
        while self.p[i] != i:
            self.p[i] = self.p[self.p[i]]
            i = self.p[i]
        return i

    def union(self, i, j):
        ri, rj = self.find(i), self.find(j)
        if ri != rj:
            self.p[ri] = rj


def _grid_candidates(positions, ids, cell):
    """Spatial hash of the given site ids; returns {cell key: [ids]}."""
    grid = {}
    for i in ids:
        p = positions[i]
        grid.setdefault((floor(p.x / cell), floor(p.y / cell)), []).append(i)
    return grid


def _plane_min_distance(positions, ids):
    """Shortest in-plane site-site distance, via a spatial hash (O(n))."""
    xs = [positions[i].x for i in ids]
    ys = [positions[i].y for i in ids]
    area = max((max(xs) - min(xs)) * (max(ys) - min(ys)), 1e-6)
    cell = max(2.0 * sqrt(area / len(ids)), 1e-3)
    for _ in range(8):                       # widen until every site has company
        grid = _grid_candidates(positions, ids, cell)
        d_min, lonely = None, False
        for i in ids:
            p = positions[i]
            kx, ky = floor(p.x / cell), floor(p.y / cell)
            best = None
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for j in grid.get((kx + dx, ky + dy), ()):
                        if j == i:
                            continue
                        d = (positions[j] - p).length
                        if best is None or d < best:
                            best = d
            if best is None:
                lonely = True
                break
            d_min = best if d_min is None else min(d_min, best)
        if not lonely:
            return d_min, cell
        cell *= 2.0
    # fall back to brute force if the grid never settled
    d_min = min((positions[i] - positions[j]).length
                for k, i in enumerate(ids) for j in ids[k + 1:])
    return d_min, max(d_min * 2.0, 1e-3)


def a_site_neighbors(a_positions):
    """
    Nearest-neighbour pairs on the A-site sublattice.

    In-plane: the shortest A-A distance within a plane (= the lattice
    parameter a for a triangular lattice), giving 6 neighbours per interior
    site. Optionally also links nearest sites in adjacent planes.
    Uses a spatial hash, so this stays fast for large supercells.
    """
    planes = {}
    for i, p in enumerate(a_positions):
        planes.setdefault(round(p.z / PLANE_TOL), []).append(i)

    pairs = []
    nn_in_plane = None
    for ids in planes.values():
        if len(ids) < 2:
            continue
        d_min, cell = _plane_min_distance(a_positions, ids)
        nn_in_plane = d_min if nn_in_plane is None else min(nn_in_plane, d_min)
        cutoff = d_min + NEIGHBOR_TOL
        cell = max(cell, cutoff)
        grid = _grid_candidates(a_positions, ids, cell)
        for i in ids:
            p = a_positions[i]
            kx, ky = floor(p.x / cell), floor(p.y / cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for j in grid.get((kx + dx, ky + dy), ()):
                        if j <= i:
                            continue
                        if (a_positions[j] - p).length <= cutoff:
                            pairs.append((i, j))

    if CONNECT_INTERLAYER:
        keys = sorted(planes)
        for k1, k2 in zip(keys, keys[1:]):
            g1, g2 = planes[k1], planes[k2]
            if not g1 or not g2:
                continue
            d_min = min((a_positions[i] - a_positions[j]).length
                        for i in g1 for j in g2)
            cutoff = d_min + NEIGHBOR_TOL
            grid = _grid_candidates(a_positions, g2, max(cutoff, 1e-3))
            for i in g1:
                p = a_positions[i]
                kx, ky = floor(p.x / cutoff), floor(p.y / cutoff)
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for j in grid.get((kx + dx, ky + dy), ()):
                            if (a_positions[j] - p).length <= cutoff:
                                pairs.append((i, j))
    return pairs, nn_in_plane


def percolation_analysis(a_positions, is_metal):
    """
    Cluster the metallic A sites. Returns:
      labels   : cluster id per A site (-1 for host sites)
      clusters : {cluster_id: [site indices]}
      channel_pairs : neighbour pairs where BOTH sites are metallic
      spanning : set of cluster ids that touch both opposite edges (x or y)
      nn       : in-plane nearest-neighbour distance
    """
    pairs, nn = a_site_neighbors(a_positions)
    dsu = _DSU(len(a_positions))
    channel_pairs = []
    for i, j in pairs:
        if is_metal[i] and is_metal[j]:
            dsu.union(i, j)
            channel_pairs.append((i, j))

    clusters = {}
    labels = [-1] * len(a_positions)
    for i in range(len(a_positions)):
        if is_metal[i]:
            r = dsu.find(i)
            clusters.setdefault(r, []).append(i)
            labels[i] = r

    xs = [p.x for p in a_positions]
    ys = [p.y for p in a_positions]
    if xs:
        x_lo, x_hi = min(xs), max(xs)
        y_lo, y_hi = min(ys), max(ys)
    edge = (nn or 1.0) * 0.75
    spanning = set()
    for cid, ids in clusters.items():
        cx = [a_positions[i].x for i in ids]
        cy = [a_positions[i].y for i in ids]
        if ((min(cx) <= x_lo + edge and max(cx) >= x_hi - edge) or
                (min(cy) <= y_lo + edge and max(cy) >= y_hi - edge)):
            spanning.add(cid)
    return labels, clusters, channel_pairs, spanning, nn


# ============================================================================
#  Mesh templates
# ============================================================================

def sphere_geometry(radius, segs, rings):
    verts = [(0.0, 0.0, radius)]
    for i in range(1, rings):
        phi = math.pi * i / rings
        z, r = radius * cos(phi), radius * sin(phi)
        for j in range(segs):
            th = 2.0 * math.pi * j / segs
            verts.append((r * cos(th), r * sin(th), z))
    verts.append((0.0, 0.0, -radius))
    last = len(verts) - 1

    def idx(i, j):
        return 1 + (i - 1) * segs + (j % segs)

    faces = [(0, idx(1, j), idx(1, j + 1)) for j in range(segs)]
    for i in range(1, rings - 1):
        for j in range(segs):
            faces.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)))
    faces += [(last, idx(rings - 1, j + 1), idx(rings - 1, j)) for j in range(segs)]
    return verts, faces


def cylinder_geometry(radius, segs):
    bot = [(radius * cos(2 * math.pi * j / segs),
            radius * sin(2 * math.pi * j / segs), 0.0) for j in range(segs)]
    top = [(x, y, 1.0) for x, y, _ in bot]
    verts = bot + top + [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    bc, tc = 2 * segs, 2 * segs + 1
    faces = []
    for j in range(segs):
        j2 = (j + 1) % segs
        faces.append((j, j2, segs + j2, segs + j))
        faces.append((bc, j2, j))
        faces.append((tc, segs + j, segs + j2))
    return verts, faces


# ============================================================================
#  Blender helpers
# ============================================================================

def clear_scene():
    for obj in list(bpy.data.objects):
        if obj.type in {"MESH", "CAMERA", "LIGHT"}:
            bpy.data.objects.remove(obj, do_unlink=True)
    for m in list(bpy.data.meshes):
        if m.users == 0:
            bpy.data.meshes.remove(m)
    for m in list(bpy.data.materials):
        if m.users == 0:
            bpy.data.materials.remove(m)


def get_collection(name, parent=None):
    col = bpy.data.collections.get(name)
    if col is None:
        col = bpy.data.collections.new(name)
        (parent or bpy.context.scene.collection).children.link(col)
    return col


def make_material(name, linear_rgb, alpha):
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = (*linear_rgb, 1.0)
        bsdf.inputs["Roughness"].default_value = 0.35
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = alpha
    mat.diffuse_color = (*linear_rgb, alpha)
    if alpha < 1.0:
        for attr, val in (("blend_method", "BLEND"),
                          ("shadow_method", "HASHED"),
                          ("surface_render_method", "BLENDED")):
            try:
                setattr(mat, attr, val)
            except Exception:
                pass
    return mat


def add_mesh_object(name, verts, faces, mat, collection, colors=None,
                    attr_name="Color", smooth=True):
    if not verts:
        return None
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.validate()
    if mat:
        mesh.materials.append(mat)
    for poly in mesh.polygons:
        poly.use_smooth = smooth
    if colors is not None:
        try:
            attr = mesh.color_attributes.new(name=attr_name,
                                             type="FLOAT_COLOR",
                                             domain="POINT")
            for i, c in enumerate(colors):
                attr.data[i].color = (c[0], c[1], c[2], 1.0)
        except Exception as e:
            print(f"[cucro2] could not write colour attribute: {e}")
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    return obj


def blob_spheres(positions, radius):
    tv, tf = sphere_geometry(radius, SPHERE_SEGMENTS, SPHERE_RINGS)
    V, F = [], []
    for pos in positions:
        base = len(V)
        V.extend((x + pos.x, y + pos.y, z + pos.z) for x, y, z in tv)
        F.extend(tuple(base + k for k in f) for f in tf)
    return V, F


def blob_cylinders(segments, radius):
    cv, cf = cylinder_geometry(radius, CYL_SEGMENTS)
    V, F = [], []
    for p1, p2 in segments:
        d = p2 - p1
        L = d.length
        if L < 1e-6:
            continue
        rot = Vector((0.0, 0.0, 1.0)).rotation_difference(d.normalized())
        base = len(V)
        for x, y, z in cv:
            w = p1 + rot @ Vector((x, y, z * L))
            V.append((w.x, w.y, w.z))
        F.extend(tuple(base + k for k in f) for f in cf)
    return V, F


# ============================================================================
#  Build one composition panel
# ============================================================================

def build_panel(struct, x_metal, panel_index, x_offset, rng):
    """Build one sheet. x_metal is None when the composition ramps across it."""
    tag = "gradient" if x_metal is None else f"x{x_metal:.2f}".replace(".", "p")
    col = get_collection(f"CuCrO2_{tag}")
    off = Vector((x_offset, 0.0, 0.0))

    atoms, _, stack_axis = build_supercell(struct)

    a_idx = [i for i, (e, _) in enumerate(atoms) if e == A_SITE_ELEMENT]
    if not a_idx:
        raise ValueError(f"No '{A_SITE_ELEMENT}' atoms found. Elements present: "
                         f"{sorted({e for e, _ in atoms})}. Set A_SITE_ELEMENT "
                         "to one of these.")
    a_pos = [atoms[i][1] for i in a_idx]

    target, info = site_targets(a_pos, x_metal, rng)
    is_metal = [rng.random() < t for t in target]
    labels, clusters, channel_pairs, spanning, nn = percolation_analysis(a_pos, is_metal)

    # ---------------- report ----------------
    n_a = len(a_idx)
    n_m = sum(is_metal)
    n_planes = len({round(p.z / PLANE_TOL) for p in a_pos})
    per_plane = n_a / n_planes
    biggest_id, biggest = None, 0
    for cid, ids in clusters.items():
        if len(ids) > biggest:
            biggest_id, biggest = cid, len(ids)

    if "islands" in info:
        print(f"\n--- panel {panel_index}: {len(info['islands'])} metallic "
              f"island(s) on a host sheet ---")
    elif "where" in info:
        print(f"\n--- panel {panel_index}: composition gradient along "
              f"{GRADIENT_AXIS}, x = {GRADIENT_MIN:.2f} -> {GRADIENT_MID:.2f} "
              f"-> {GRADIENT_MAX:.2f} ---")
    else:
        side = ("above threshold" if x_metal > 0.5 else
                "below threshold" if x_metal < 0.5 else "at threshold")
        print(f"\n--- panel {panel_index}: x = {x_metal:.3f} ({side}) ---")
    print(f"  A-planes: {n_planes}   A sites: {n_a} ({per_plane:.0f} per plane)")
    print(f"  metallic: {n_m} (actual x = {n_m / n_a:.3f})")
    print(f"  in-plane A-A nearest-neighbour distance: {nn:.3f} A "
          f"(= lattice parameter a)")
    denom = n_a if CONNECT_INTERLAYER else per_plane
    scope = "all A sites" if CONNECT_INTERLAYER else "the A sites in one plane"
    print(f"  metallic clusters: {len(clusters)}   largest: {biggest} sites "
          f"({biggest / denom * 100:.1f}% of {scope})")
    print(f"  spanning cluster: "
          f"{'YES -- percolating (metallic)' if spanning else 'no -- disconnected (insulating)'}")
    print("  2D triangular-lattice site percolation threshold: p_c = 0.5")

    if "islands" in info:
        _report_islands(a_pos, info["islands"], is_metal, labels, clusters)
    elif "where" in info:
        _report_gradient(a_pos, info["where"], target, is_metal, labels,
                         biggest_id, nn)

    # ---------------- layer gap ----------------
    # Everything above ran on the true crystal geometry. Now that the bonds
    # and the percolation analysis are settled, the layers can be pulled
    # apart for legibility: the bonds computed above simply stretch, which is
    # the point -- it is the vertical Cu-O struts that show the sandwich.
    bond_pairs = chemical_bonds(atoms) if SHOW_BONDS else []
    shift = layer_offsets(atoms, stack_axis)
    if abs(LAYER_GAP) > 1e-9:
        for i, d in enumerate(shift):
            atoms[i][1] = atoms[i][1] + d
        a_pos = [atoms[i][1] for i in a_idx]      # rebind after the shift
        print(f"  layer gap: {LAYER_GAP:+.2f} A between "
              f"{len(layer_groups(atoms, stack_axis))} layers")
    if not SHOW_VERTICAL_BONDS:
        print(f"  vertical bonds excluded (within "
              f"{VERTICAL_BOND_ANGLE:.0f} deg of the stacking axis)")

    # ---------------- materials ----------------
    a_atoms = TRANSPARENCY.get("atoms", 0.0)
    mat_metal = make_material(f"{tag}_A_metal", hex_to_linear(COLORS["metal_site"]), 1 - a_atoms)
    mat_host = make_material(f"{tag}_A_host", hex_to_linear(COLORS["host_site"]), 1 - a_atoms)
    mat_dim = make_material(f"{tag}_A_dimmed", hex_to_linear(COLORS["dimmed"]), 1 - a_atoms)
    mat_b = make_material(f"{tag}_B", hex_to_linear(COLORS["Cr"]), 1 - a_atoms)
    mat_o = make_material(f"{tag}_O", hex_to_linear(COLORS["O"]), 1 - a_atoms)
    mat_bond = make_material(f"{tag}_bond", hex_to_linear(COLORS["bond"]),
                             1 - TRANSPARENCY.get("bonds", 0.0))
    mat_chan = make_material(f"{tag}_channel", hex_to_linear(COLORS["channel"]),
                             1 - TRANSPARENCY.get("channel", 0.0))
    mat_guide = make_material(f"{tag}_guide", hex_to_linear(COLORS["guide"]),
                              1 - TRANSPARENCY.get("guide", 0.0))

    # ---------------- A sites ----------------
    if SHOW_ATOMS and SHOW_A_SITES:
        r_m = RADII.get("metal_site", DEFAULT_RADIUS)
        r_h = RADII.get("host_site", DEFAULT_RADIUS)
        host_pos = [a_pos[k] + off for k in range(len(a_idx)) if not is_metal[k]]
        V, F = blob_spheres(host_pos, r_h)
        add_mesh_object(f"{tag}_A_host", V, F, mat_host, col)

        if COLOR_MODE == "clusters":
            for n, (cid, ids) in enumerate(sorted(clusters.items(),
                                                  key=lambda t: -len(t[1]))):
                hue = (n * 0.6180339887) % 1.0
                rgb = colorsys.hsv_to_rgb(hue, 0.62, 0.95)
                m = make_material(f"{tag}_cluster{n:03d}",
                                  tuple(v ** 2.2 for v in rgb), 1 - a_atoms)
                V, F = blob_spheres([a_pos[i] + off for i in ids], r_m)
                add_mesh_object(f"{tag}_cluster{n:03d}_{len(ids)}sites", V, F, m, col)
        elif COLOR_MODE == "spanning":
            span_ids = [i for i in range(len(a_idx))
                        if labels[i] in spanning and is_metal[i]]
            other_ids = [i for i in range(len(a_idx))
                         if is_metal[i] and labels[i] not in spanning]
            V, F = blob_spheres([a_pos[i] + off for i in span_ids], r_m)
            add_mesh_object(f"{tag}_A_metal_spanning", V, F, mat_metal, col)
            V, F = blob_spheres([a_pos[i] + off for i in other_ids], r_m)
            add_mesh_object(f"{tag}_A_metal_isolated", V, F, mat_dim, col)
        else:                                    # "species"
            V, F = blob_spheres([a_pos[k] + off for k in range(len(a_idx))
                                 if is_metal[k]], r_m)
            add_mesh_object(f"{tag}_A_metal", V, F, mat_metal, col)

    # ---------------- B sites and oxygen ----------------
    if SHOW_ATOMS:
        by_elem = {}
        for i, (e, p) in enumerate(atoms):
            if i in set(a_idx):
                continue
            by_elem.setdefault(e, []).append(p + off)
        for e, ps in by_elem.items():
            if e == "O":
                if not SHOW_OXYGEN:
                    continue
                mat, r = mat_o, RADII.get("O", DEFAULT_RADIUS)
            else:
                if not SHOW_B_SITES:
                    continue
                mat = mat_b
                r = RADII.get(e, RADII.get("Cr", DEFAULT_RADIUS))
            V, F = blob_spheres(ps, r)
            add_mesh_object(f"{tag}_{e}", V, F, mat, col)

    # ---------------- chemical bonds ----------------
    if SHOW_BONDS:
        segs = drop_vertical([(atoms[i][1] + off, atoms[j][1] + off)
                              for i, j in bond_pairs], stack_axis)
        V, F = blob_cylinders(segs, BOND_RADIUS)
        add_mesh_object(f"{tag}_bonds", V, F, mat_bond, col)

    # ---------------- metallic channel network ----------------
    if SHOW_CHANNELS and channel_pairs:
        if COLOR_MODE == "spanning":
            span = drop_vertical([(a_pos[i] + off, a_pos[j] + off)
                                  for i, j in channel_pairs
                                  if labels[i] in spanning], stack_axis)
            rest = drop_vertical([(a_pos[i] + off, a_pos[j] + off)
                                  for i, j in channel_pairs
                                  if labels[i] not in spanning], stack_axis)
            V, F = blob_cylinders(span, CHANNEL_RADIUS)
            add_mesh_object(f"{tag}_channels_spanning", V, F, mat_chan, col)
            V, F = blob_cylinders(rest, CHANNEL_RADIUS)
            add_mesh_object(f"{tag}_channels_isolated", V, F, mat_dim, col)
        else:
            segs = drop_vertical([(a_pos[i] + off, a_pos[j] + off)
                                  for i, j in channel_pairs], stack_axis)
            V, F = blob_cylinders(segs, CHANNEL_RADIUS)
            add_mesh_object(f"{tag}_channels", V, F, mat_chan, col)

    # ---------------- faint guide over the whole A lattice ----------------
    if SHOW_LATTICE_GUIDE:
        pairs, _ = a_site_neighbors(a_pos)
        segs = drop_vertical([(a_pos[i] + off, a_pos[j] + off)
                              for i, j in pairs], stack_axis)
        V, F = blob_cylinders(segs, GUIDE_RADIUS)
        add_mesh_object(f"{tag}_A_lattice_guide", V, F, mat_guide, col)

    width = (max(p.x for p in a_pos) - min(p.x for p in a_pos)) if a_pos else 0.0

    # bounds of everything drawn in this panel, padded by the largest display
    # radius, so the camera can frame the whole series
    pad = max(list(RADII.values()) + [DEFAULT_RADIUS])
    pts = [p + off for _, p in atoms]
    bounds = (min(p.x for p in pts) - pad, max(p.x for p in pts) + pad,
              min(p.y for p in pts) - pad, max(p.y for p in pts) + pad,
              min(p.z for p in pts) - pad, max(p.z for p in pts) + pad)
    return width, bounds


def _report_islands(a_pos, islands, is_metal, labels, clusters):
    """
    What each island actually came out as. The interesting number is the last
    one: an island only conducts if its metallic sites join up into a single
    cluster, so a count above 1 means the island is internally broken.
    """
    n_metal = sum(1 for m in is_metal if m)
    print(f"  host sheet: {len(a_pos) - n_metal} sites "
          f"({(len(a_pos) - n_metal) / max(len(a_pos), 1) * 100:.0f}%), "
          f"metallic: {n_metal} ({n_metal / max(len(a_pos), 1) * 100:.0f}%)")
    print("     island   centre (A)         radius   sites   metallic   clusters")
    for n, (cx, cy, r, _) in enumerate(islands, 1):
        near = [k for k, p in enumerate(a_pos)
                if hypot(p.x - cx, p.y - cy) <= r * (1.0 + ISLAND_WOBBLE)]
        met = [k for k in near if is_metal[k]]
        cl = {labels[k] for k in met}
        biggest = max((sum(1 for k in met if labels[k] == c) for c in cl),
                      default=0)
        print(f"     {n:^6d}   ({cx:6.1f}, {cy:6.1f})   {r:6.2f}   "
              f"{len(near):5d}   {len(met):5d}      {len(cl)}"
              + (f" (largest {biggest})" if len(cl) > 1 else " -- connected"))
    stray = sum(1 for k, p in enumerate(a_pos) if is_metal[k] and
                all(hypot(p.x - cx, p.y - cy) > r * (1.0 + ISLAND_WOBBLE)
                    for cx, cy, r, _ in islands))
    print(f"  stray metallic sites out on the host sheet: {stray}")


def _report_gradient(a_pos, where, target, is_metal, labels, biggest_id, nn):
    """
    Break the sheet into slices along the gradient and report what each one
    came out as. This is where the threshold shows: the slices below p_c hold
    only small islands, and past it nearly every metallic site joins the one
    cluster that runs across the sheet.
    """
    nb = max(int(GRADIENT_REPORT_BANDS), 1)
    bands = [[] for _ in range(nb)]
    for k, t in enumerate(where):
        bands[min(int(t * nb), nb - 1)].append(k)

    print(f"  composition across {GRADIENT_AXIS}, in {nb} slices:")
    print("     slice     target x   actual x   metallic   in the largest cluster")
    for b, ids in enumerate(bands):
        if not ids:
            continue
        tgt = sum(target[k] for k in ids) / len(ids)
        met = [k for k in ids if is_metal[k]]
        act = len(met) / len(ids)
        big = sum(1 for k in met if labels[k] == biggest_id)
        share = (big / len(met) * 100.0) if met else 0.0
        lo, hi = b / nb, (b + 1) / nb
        print(f"     {lo:.2f}-{hi:.2f}   {tgt:8.3f}   {act:8.3f}   "
              f"{len(met):4d}/{len(ids):<4d}  {big:4d} ({share:5.1f}%)")

    # Which way the largest cluster reaches matters here in a way it does not
    # for a uniform sheet: a gradient is expected to percolate ACROSS the
    # ramp, over on the metallic side, while staying disconnected from the
    # insulating end, so end-to-end spanning is not the interesting question.
    ids = [k for k in range(len(a_pos)) if labels[k] == biggest_id]
    if not ids:
        return
    edge = (nn or 1.0) * 0.75
    for axis in ("x", "y"):
        vals = [getattr(p, axis) for p in a_pos]
        cvals = [getattr(a_pos[k], axis) for k in ids]
        full = max(vals) - min(vals)
        reach = max(cvals) - min(cvals)
        spans = (min(cvals) <= min(vals) + edge and
                 max(cvals) >= max(vals) - edge)
        role = "along the gradient" if axis == GRADIENT_AXIS else "across it"
        print(f"  largest cluster reaches {reach / max(full, 1e-9) * 100:5.1f}% "
              f"of the sheet in {axis} ({role}): "
              f"{'spans' if spans else 'does not span'}")


# ============================================================================
#  Camera, lights and render settings
# ============================================================================

def _set_render_engine(scene, name):
    """Set the engine by identifier; the EEVEE one is renamed between releases."""
    wanted = {"EEVEE": ("BLENDER_EEVEE_NEXT", "BLENDER_EEVEE"),
              "CYCLES": ("CYCLES",)}.get(str(name).upper(), (str(name),))
    try:
        available = {i.identifier for i in
                     scene.render.bl_rna.properties["engine"].enum_items}
    except Exception:
        available = set()
    for ident in wanted:
        if available and ident not in available:
            continue
        try:
            scene.render.engine = ident
            return ident
        except Exception:
            continue
    return scene.render.engine


def _mix(c1, c2, t):
    return tuple(a + (b - a) * t for a, b in zip(c1, c2))


def emission_vertex_material(name, attr_name):
    """Flat emission driven by a per-vertex colour attribute."""
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (300, 0)
    emis = nt.nodes.new("ShaderNodeEmission")
    emis.location = (100, 0)
    emis.inputs["Strength"].default_value = 1.0
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.location = (-150, 0)
    attr.attribute_name = attr_name
    nt.links.new(attr.outputs["Color"], emis.inputs["Color"])
    nt.links.new(emis.outputs["Emission"], out.inputs["Surface"])
    return mat


def add_background_plane(scene, cam, target, up_hint, bottom, top, reach):
    """
    The background is an emissive quad placed behind everything and sized to
    fill the frame, with the gradient baked into its vertex colours.

    This replaces compositing a generated image behind the render, which is
    where the background was going missing: bpy.data.images.new() makes an
    image whose source is GENERATED, and Blender rebuilds a generated image's
    buffer from its own settings whenever it re-evaluates it. The gradient
    poked in from Python was therefore not there at render time -- the Image
    node handed the compositor an empty frame, and the background came out
    transparent. Geometry and vertex colours have nothing to regenerate, and
    behave identically in EEVEE and Cycles under either projection.
    """
    loc = Vector(cam.location)
    forward = (target - loc)
    if forward.length < 1e-9:
        return None
    forward = forward.normalized()
    right = forward.cross(up_hint)
    right = right.normalized() if right.length > 1e-6 else Vector((1.0, 0.0, 0.0))
    up = right.cross(forward).normalized()

    centre = target + forward * (reach * 2.0 + 1.0)
    res_x, res_y = RESOLUTION
    aspect = res_x / float(res_y)
    if cam.data.type == "ORTHO":
        s = cam.data.ortho_scale
        half_w = 0.5 * (s if aspect >= 1.0 else s * aspect)
        half_h = 0.5 * (s / aspect if aspect >= 1.0 else s)
    else:
        tan_h = 18.0 / max(cam.data.lens, 1e-3)     # 36 mm sensor, half-angle
        half_w = (centre - loc).length * tan_h
        half_h = half_w / aspect
    half_w *= BACKGROUND_MARGIN
    half_h *= BACKGROUND_MARGIN

    rows = max(int(BACKGROUND_STEPS), 2)
    verts, faces, colors = [], [], []
    for i in range(rows + 1):
        t = i / rows
        c = _mix(bottom, top, _ease_clamped(t) if BACKGROUND_SMOOTH else t)
        y = (2.0 * t - 1.0) * half_h
        for side in (-1.0, 1.0):
            p = centre + right * (side * half_w) + up * y
            verts.append((p.x, p.y, p.z))
            colors.append(c)
    for i in range(rows):
        b = 2 * i
        faces.append((b, b + 1, b + 3, b + 2))

    mat = emission_vertex_material("Background", "BackgroundColor")
    obj = add_mesh_object("Background", verts, faces, mat,
                          get_collection("CuCrO2_Background"), colors=colors,
                          attr_name="BackgroundColor", smooth=False)
    if obj is None:
        return None
    # a backdrop, not a light: keep it out of everything but the camera ray
    for attr in ("visible_diffuse", "visible_glossy", "visible_transmission",
                 "visible_volume_scatter", "visible_shadow"):
        try:
            setattr(obj, attr, False)
        except Exception:
            pass
    need = (centre - loc).length + max(half_w, half_h) * 2.0
    try:
        if cam.data.clip_end < need:
            cam.data.clip_end = need * 1.2
    except Exception:
        pass
    return obj


def setup_render(scene):
    engine = _set_render_engine(scene, RENDER_ENGINE)
    r = scene.render
    r.resolution_x, r.resolution_y = RESOLUTION
    r.resolution_percentage = 100
    # The backdrop is an object in the scene, so the film can stay
    # transparent in every mode: where there is no backdrop the render simply
    # comes out with an alpha channel.
    r.film_transparent = True

    # RGBA is the part people miss: with the film transparent but the colour
    # mode left at RGB, the alpha channel is dropped on save and the
    # background comes back black.
    try:
        r.image_settings.file_format = "PNG"
        r.image_settings.color_mode = "RGBA"
        r.image_settings.color_depth = "16"
    except Exception:
        pass

    if engine == "CYCLES":
        try:
            scene.cycles.samples = RENDER_SAMPLES
            scene.cycles.use_denoising = True
        except Exception:
            pass
    else:
        for attr in ("taa_render_samples", "taa_samples"):
            if hasattr(scene.eevee, attr):
                try:
                    setattr(scene.eevee, attr, RENDER_SAMPLES)
                except Exception:
                    pass
    if VIEW_TRANSFORM:
        try:
            scene.view_settings.view_transform = VIEW_TRANSFORM
        except Exception as e:
            print(f"[cucro2] view transform {VIEW_TRANSFORM!r} not available: {e}")

    if OUTPUT_PATH:
        r.filepath = os.path.abspath(os.path.expanduser(OUTPUT_PATH))
    return engine


def setup_world(scene):
    """
    The world is only ever a neutral fill light here -- the visible background
    is the backdrop object, so the world never has to double as one.
    """
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    if not bg:
        return
    bg.inputs["Color"].default_value = (1.0, 1.0, 1.0, 1.0)
    bg.inputs["Strength"].default_value = AMBIENT_STRENGTH


def _add_sun(scene, name, direction, energy, softness=0.25):
    data = bpy.data.lights.new(name, "SUN")
    data.energy = energy
    try:
        data.angle = softness
    except Exception:
        pass
    obj = bpy.data.objects.new(name, data)
    scene.collection.objects.link(obj)
    d = Vector(direction).normalized()
    obj.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()
    return obj


def setup_lights(scene):
    """Key, fill and rim -- enough to make the spheres read as spheres."""
    _add_sun(scene, "Key", (-0.5, -0.7, -1.0), KEY_LIGHT_ENERGY)
    _add_sun(scene, "Fill", (0.8, 0.4, -0.6), KEY_LIGHT_ENERGY * 0.35)
    _add_sun(scene, "Rim", (0.1, 0.9, -0.25), KEY_LIGHT_ENERGY * 0.45)


def drawn_points():
    """
    Every vertex of every mesh built so far, in world space.

    The camera frames these rather than a bounding box worked out from atom
    positions and a radius. An estimate has to guess what the drawn geometry
    actually reaches -- sphere tessellation, bond caps, the widest display
    radius -- and anything it misses is something the frame clips off. The
    vertices cannot be wrong about it.

    Every mesh here is built in world coordinates with its object left at the
    origin, so the vertices can be read as they are.
    """
    pts = []
    for obj in bpy.data.objects:
        if obj.type != "MESH" or obj.data is None:
            continue
        for v in obj.data.vertices:
            co = getattr(v, "co", v)
            pts.append((co[0], co[1], co[2]))
    return pts


def setup_camera(scene, points):
    """
    Frame the structure from CAMERA_ELEVATION degrees above its plane, at
    CAMERA_AZIMUTH around it. build_supercell has already rotated the crystal
    so that VIEW_DIRECTION points along +Z, which is the axis the elevation is
    measured from.

    The framing is exact rather than a guess. Project every point onto the
    camera's own right and up axes; those two components do not change as the
    camera slides along its view direction, so each point sets a lower bound
    on the distance and the largest of them frames the lot.
    """
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    target = Vector((0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)),
                     0.5 * (min(zs) + max(zs))))
    corners = [Vector(p) for p in points]

    elev = radians(min(max(CAMERA_ELEVATION, -89.0), 89.0))
    azim = radians(CAMERA_AZIMUTH)
    # unit vector from the structure toward the camera
    eye = Vector((cos(elev) * cos(azim), cos(elev) * sin(azim), sin(elev)))
    forward = -eye
    right = forward.cross(Vector((0.0, 0.0, 1.0)))
    right = right.normalized() if right.length > 1e-6 else Vector((1.0, 0.0, 0.0))
    up = right.cross(forward).normalized()

    res_x, res_y = RESOLUTION
    aspect = res_x / float(res_y)

    cam_data = bpy.data.cameras.new("CuCrO2Camera")
    if ORTHOGRAPHIC:
        cam_data.type = "ORTHO"
        half_h = max(abs((c - target).dot(right)) for c in corners) * CAMERA_MARGIN
        half_v = max(abs((c - target).dot(up)) for c in corners) * CAMERA_MARGIN
        # ortho_scale spans the longer image axis, so convert the other one
        # through the aspect ratio before taking the maximum
        span = (max(2 * half_h, 2 * half_v * aspect) if aspect >= 1.0
                else max(2 * half_v, 2 * half_h / aspect))
        cam_data.ortho_scale = span
        depth = max(abs((c - target).dot(forward)) for c in corners)
        dist = depth * 2.0 + span
    else:
        cam_data.type = "PERSP"
        cam_data.lens = CAMERA_LENS
        try:
            cam_data.sensor_fit = "HORIZONTAL"     # makes the maths definite
        except Exception:
            pass
        tan_h = 18.0 / max(CAMERA_LENS, 1e-3)      # 36 mm sensor, half-angle
        tan_v = tan_h / aspect
        dist = 0.0
        for c in corners:
            d = c - target
            ahead = d.dot(forward)
            dist = max(dist,
                       abs(d.dot(right)) * CAMERA_MARGIN / tan_h - ahead,
                       abs(d.dot(up)) * CAMERA_MARGIN / tan_v - ahead)
        dist = max(dist, 1.0)

    reach = max((c - target).length for c in corners)
    cam_data.clip_start = max(0.01, (dist - reach) * 0.5)
    cam_data.clip_end = dist + reach * 4.0 + 100.0

    cam = bpy.data.objects.new("CuCrO2Camera", cam_data)
    scene.collection.objects.link(cam)
    cam.location = target + eye * dist
    cam.rotation_euler = forward.to_track_quat("-Z", "Y").to_euler()
    scene.camera = cam

    # how much of the frame the structure actually uses. A wide flat sheet
    # seen from a low elevation is much wider than it is tall, so it fills the
    # frame across and leaves air above and below -- raise CAMERA_ELEVATION or
    # widen RESOLUTION if that bothers you.
    if ORTHOGRAPHIC:
        fw = 2 * max(abs((c - target).dot(right)) for c in corners) / span
        fh = 2 * max(abs((c - target).dot(up)) for c in corners) / (span / aspect
                                                                   if aspect >= 1.0
                                                                   else span)
    else:
        fw = max(abs((c - target).dot(right)) /
                 max((c - target).dot(forward) + dist, 1e-9) for c in corners) / tan_h
        fh = max(abs((c - target).dot(up)) /
                 max((c - target).dot(forward) + dist, 1e-9) for c in corners) / tan_v
    print(f"  framing: the structure fills {fw * 100:.0f}% of the frame across "
          f"and {fh * 100:.0f}% of it up, from {len(points)} vertices")
    return cam


def setup_scene(bounds):
    if BACKGROUND_MODE not in ("gradient", "flat", "transparent"):
        raise ValueError("BACKGROUND_MODE must be 'gradient', 'flat' or "
                         f"'transparent'; got {BACKGROUND_MODE!r}")
    scene = bpy.context.scene
    engine = setup_render(scene)
    setup_world(scene)
    if ADD_LIGHTS:
        setup_lights(scene)
    # frame what was actually drawn; the estimated bounds are only a fallback
    # for the case where nothing was built at all
    points = drawn_points()
    if not points:
        x0, x1, y0, y1, z0, z1 = bounds
        points = [(x, y, z) for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)]
    cam = setup_camera(scene, points)

    if BACKGROUND_MODE in ("gradient", "flat"):
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        zs = [p[2] for p in points]
        target = Vector((0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)),
                         0.5 * (min(zs) + max(zs))))
        reach = max((Vector(p) - target).length for p in points)
        bottom = hex_to_linear(BACKGROUND_COLOR)
        top = (hex_to_linear(BACKGROUND_TOP_COLOR)
               if BACKGROUND_MODE == "gradient" else bottom)
        add_background_plane(scene, cam, target, Vector((0.0, 0.0, 1.0)),
                             bottom, top, reach)
    scene.use_nodes = False
    return engine


# ============================================================================
#  Main
# ============================================================================

def main():
    if CLEAR_SCENE:
        clear_scene()

    struct = load_structure()
    if COMPOSITION_MODE in ("islands", "gradient"):
        if COMPOSITION_SERIES:
            print(f"[cucro2] COMPOSITION_MODE is {COMPOSITION_MODE!r}, so "
                  "COMPOSITION_SERIES is ignored: this is one sheet, not a "
                  "series.")
        comps = [None]          # None means "the mode decides, site by site"
    else:
        comps = COMPOSITION_SERIES if COMPOSITION_SERIES else [METAL_FRACTION]
    rng = random.Random(RANDOM_SEED)

    x_off = 0.0
    box = None
    for n, x in enumerate(comps):
        if x is not None and not 0.0 <= x <= 1.0:
            raise ValueError(f"composition must be between 0 and 1; got {x}")
        width, bounds = build_panel(struct, x, n + 1, x_off, rng)
        box = bounds if box is None else (
            min(box[0], bounds[0]), max(box[1], bounds[1]),
            min(box[2], bounds[2]), max(box[3], bounds[3]),
            min(box[4], bounds[4]), max(box[5], bounds[5]))
        x_off += width + PANEL_GAP

    engine = None
    if SETUP_SCENE and box is not None:
        engine = setup_scene(box)

    print("\n[cucro2] Done. Metallic sites are coloured "
          f"{COLORS['metal_site']}, host sites {COLORS['host_site']}. "
          "Switch the viewport to Material Preview to see the colours.")
    if SETUP_SCENE:
        print(f"[cucro2] camera: "
              f"{'orthographic' if ORTHOGRAPHIC else f'{CAMERA_LENS:.0f} mm perspective'}"
              f", {CAMERA_ELEVATION:.0f} deg above the sheet at "
              f"{CAMERA_AZIMUTH:.0f} deg around it")
    if SETUP_SCENE:
        described = {"transparent": "transparent (RGBA PNG)",
                     "flat": f"flat {BACKGROUND_COLOR}",
                     "gradient": f"gradient {BACKGROUND_COLOR} -> "
                                 f"{BACKGROUND_TOP_COLOR}"}[BACKGROUND_MODE]
        print(f"[cucro2] engine {engine}, background {described}")

    if RENDER_NOW:
        if not OUTPUT_PATH:
            print("[cucro2] RENDER_NOW is on but OUTPUT_PATH is empty; "
                  "nothing was written.")
        else:
            path = os.path.abspath(os.path.expanduser(OUTPUT_PATH))
            folder = os.path.dirname(path)
            if folder and not os.path.isdir(folder):
                os.makedirs(folder, exist_ok=True)
            bpy.context.scene.render.filepath = path
            print(f"[cucro2] rendering to {path} ...")
            bpy.ops.render.render(write_still=True)
            print("[cucro2] render written.")


main()
