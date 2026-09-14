"""
exciton.py
==========
Blender script: a glowing exciton -- a bound electron-hole pair -- drawn as
a small electron in orbit around a larger hole. Behind the electron runs a
motion tail built as a band of concentric strands lying flat in the orbit
plane, like the rings of Saturn, brightest where it meets the electron and
fading away as it sweeps round. A volumetric aura sits directly between the
two particles. The real dipole field lines are available too, as an option.
Renders on a transparent background so the image can be dropped straight onto
a PowerPoint slide.

Run: Blender -> Scripting -> Open -> edit USER PARAMETERS -> Alt+P.

What is actually being drawn
----------------------------
The electron sits on a circle of radius ORBIT_RADIUS about the hole, at
ORBIT_PHASE degrees around it. The tail is TRAIL_STRANDS concentric strands
spread over TRAIL_WIDTH, each one a tube whose thickness and colour follow
how far behind the electron that point is -- so the band swells and brightens
into the electron and fades out behind it, while the rest of each circle
stays as a faint ring at TRAIL_RING_LEVEL. The outer strands trail slightly
longer than the inner ones (TRAIL_SHEAR) so the tail feathers out rather than
ending on a straight edge.

With SHOW_FIELD_LINES on, the field lines are not decorative arcs either:
they are streamlines of the electric field of two equal and opposite point
charges (+ at the hole, - at the electron), integrated with RK4 through

        E(r) = (r - r_h)/|r - r_h|^3  -  (r - r_e)/|r - r_e|^3

Every line therefore starts on the hole, ends on the electron, and bows
outward by exactly the amount the dipole geometry dictates. Lines launched
close to the axis stay in a tight bundle through the middle; lines launched at
larger angles swing out to the sides. FIELD_LINE_AXIS_BIAS controls how many
are launched near the axis, which is what makes the middle dense.

Colours default to the ones sampled from the reference image: a hot pink-red
hole at the centre and a bright cyan electron going round it. Swap them if
you prefer the other convention -- nothing else depends on the choice.

Shading
-------
Everything glows via Emission shaders, so switch the viewport to Material
Preview or Rendered to see it. The particles are not flat discs: their
emission colour and strength follow a half-Lambert gradient across the
surface plus a rim term at the silhouette, so they read as spheres without
needing any scene lights. The aura and the particle haloes are continuous
gradients -- a real volume for the aura, a view-angle alpha falloff for the
haloes -- rather than stacks of nested shells. CAMERA_ELEVATION decides how
far above the orbit plane the camera sits; at 0 the tail is seen exactly
edge-on and collapses to a straight line.

SETUP_SCENE builds the camera, the render settings and a compositor Glare
node for the bloom. BACKGROUND_MODE picks between a transparent RGBA film
(the default, so this figure lays over the cucro2perc render), the deep blue
field of Reference 1, and a flat colour. Whenever the film is transparent the
compositor pushes the bloom into the alpha channel too, so the glow survives
being composited instead of being clipped away.
"""

import bpy
import math
import os
from math import cos, sin, sqrt, pi
from mathutils import Vector, Matrix

# ============================================================================
#  USER PARAMETERS
# ============================================================================

# --- how tightly the pair is bound ------------------------------------------
# One dial for the whole model. Everything written below describes a STRONGLY
# bound exciton -- compact, hot, saturated. BINDING_STRENGTH blends each of the
# quantities listed in the WEAK_ block toward a loose, faint, washed-out one,
# so 1.0 is a tightly bound pair and 0.0 is a barely bound one. It drives
# several things at once on purpose: at these sizes a single quantity changing
# by 30% does not read, but the orbit, the aura, the tail and the colour all
# moving together do.
BINDING_STRENGTH     = 1.0

WEAK_ORBIT_RADIUS    = 6.4    # the pair drifts apart as the binding weakens
WEAK_AURA_DENSITY    = 0.10   # and the cloud between them thins out
WEAK_AURA_EMISSION   = 1.2
WEAK_AURA_MARGIN     = 2.2    # a loose halo rather than a tight bridge
WEAK_TRAIL_LENGTH    = 0.90   # a long wispy tail instead of a short hot one
WEAK_TRAIL_LEVEL     = 1.0
WEAK_GLOW_OUTER_SCALE = 2.3   # a bigger, fainter halo on each particle
WEAK_GLOW_ALPHA      = 0.26
WEAK_FUNNEL_LEVEL    = 0.35   # a dimmer funnel down to the sheet
WEAK_EMISSION_SCALE  = 0.65   # everything dims
WEAK_DESATURATION    = 0.55   # ... and washes toward grey, so weakly bound
                              # pairs recede into the background

# --- the pair: an electron in orbit around a hole ---------------------------
CENTER           = (0.0, 0.0, 0.0)   # the hole sits here
ORBIT_RADIUS     = 3.6        # how far the electron orbits from the hole
ORBIT_NORMAL     = (0.0, 0.0, 1.0)   # the axis the electron circles about
ORBIT_PHASE      = 40.0       # degrees around the orbit, measured from the
                              # right of frame: 0 puts the electron at the
                              # right-hand edge of the path, 90 nearest the
                              # camera, 270 furthest behind the hole
ORBIT_DIRECTION  = 1          # +1 or -1, which way round the electron
                              # travels. The trail streams behind it either
                              # way.

ELECTRON_RADIUS  = 1.10       # the electron is the smaller one, out on the
HOLE_RADIUS      = 1.25       # orbit; the hole is the larger one at the
                              # centre, which is what it is being orbited

# Colours as sRGB hex, sampled from Reference 1: a hot pink-red hole at the
# centre with a bright cyan electron going round it, which is the pairing
# that stays legible against the deep blue field of that figure.
ELECTRON_COLOR   = "#12C2FF"   # bright cyan core
ELECTRON_RIM     = "#8FE6FF"   # pale cyan highlight
HOLE_COLOR       = "#EF4E7F"   # hot pink-red core
HOLE_RIM         = "#FFA6C8"   # pale pink highlight
FIELD_LINE_COLOR = "#FBDB8C"   # gold field lines, as in the reference

ELECTRON_EMISSION = 1.5        # emission strength of each particle. Push
HOLE_EMISSION     = 1.3        # these much past ~2.5 and the whole sphere
                               # clips to white, which is what made the
                               # particles read as flat coloured circles in
                               # the first place -- the shading gradient has
                               # to stay inside the displayable range.

# --- the motion tail --------------------------------------------------------
# The tail is a DASHED band: several concentric strands lying flat in the orbit
# plane, like the rings of Saturn, each one solid where it leaves the electron
# and then breaking into dashes that shorten and space further apart behind it
# until they run out altogether.
SHOW_TRAIL           = True
TRAIL_STRANDS        = 5      # concentric strands making up the band. 1 gives
                              # a single dashed line instead of a band.
TRAIL_WIDTH          = 0.95   # width of the whole band, centred on the orbit
TRAIL_STRAND_RADIUS  = 0.055  # thickness of a strand where it meets the
                              # electron
TRAIL_LENGTH         = 0.50   # fraction of the loop the tail reaches back over
TRAIL_FALLOFF        = 1.6    # >1 dims the tail faster along its length
TRAIL_TAIL_WIDTH     = 0.30   # thickness at the far end, as a fraction of
                              # TRAIL_STRAND_RADIUS
TRAIL_SHEAR          = 0.35   # outer strands trail this much longer than
                              # inner ones, so the tail feathers out instead
                              # of ending square
TRAIL_EDGE_FADE      = 0.55   # how much dimmer the outermost strands are than
                              # the middle one: 0 = all equal, 1 = edges dark
TRAIL_LEVEL          = 2.6    # brightness where the band meets the electron
TRAIL_END_LEVEL      = 0.45   # brightness of the last dashes before they go

# How the dashes break up going backwards. The dash length and the gap are
# set separately rather than as a duty cycle of one period: tie them together
# and the growing period drags the early dashes longer before the falling duty
# gets on top of it, so the tail briefly strengthens before it fades. Set
# apart, each one only ever moves the way it should.
TRAIL_SOLID          = 0.18   # fraction of the tail nearest the electron that
                              # stays unbroken before the dashes start
TRAIL_DASH_LENGTH    = 0.026  # the first dash, as a fraction of the whole loop
TRAIL_DASH_GAP       = 0.008  # the first gap
TRAIL_DASH_FALLOFF   = 1.2    # >1 shrinks the dashes away faster
TRAIL_DASH_GROWTH    = 5.0    # the gap opens to this multiple by the end,
                              # which is what makes the dashes less frequent
TRAIL_DASH_MIN       = 0.12   # a dash shorter than this fraction of
                              # TRAIL_DASH_LENGTH is where the tail ends
TRAIL_DASH_TAPER     = 0.32   # taper on each dash, so they are lens shaped
                              # rather than blunt cylinders

TRAIL_COLOR          = "#8FE4FF"   # pale cyan, as in Reference 1
TRAIL_HEAD_COLOR     = ""     # "" = the electron's own highlight colour
TRAIL_EMISSION       = 2.0
TRAIL_RESOLUTION     = 320    # points per loop, sampled per dash
TRAIL_SEGMENTS       = 6      # cross-section resolution of a strand

# --- the light funnel down to the sheet -------------------------------------
# Borrowed from Reference 1: a cone of light dropping from the exciton to the
# layer below it, with a pool of light where it lands. It is what ties the
# exciton to the lattice rather than leaving it floating over the top.
SHOW_FUNNEL          = True
FUNNEL_DIRECTION     = (0.0, 0.0, -1.0)  # which way the funnel drops
FUNNEL_DROP          = 7.5    # how far down it reaches
# The funnel is an hourglass: it flares out of the exciton, pinches to a
# waist, then opens again where it lands. Each end is set on its own, so the
# mouth at the exciton and the spread on the sheet are independent.
FUNNEL_TOP_RADIUS    = 2.4    # radius where it leaves the exciton
FUNNEL_WAIST_RADIUS  = 0.50   # radius at the pinch
FUNNEL_BOTTOM_RADIUS = 4.2    # radius where it lands
FUNNEL_WAIST         = 0.42   # where the pinch sits along the drop: 0 at the
                              # exciton, 1 at the sheet. Below 0.5 the lower
                              # half is the longer one, which gives the wide
                              # slow flare down onto the layer.
FUNNEL_TOP_FLARE     = 1.8    # >1 holds each half near the waist radius and
FUNNEL_BOTTOM_FLARE  = 2.4    # then opens it near its end, which is what
                              # makes the curve an hourglass rather than two
                              # straight cones. 1 gives straight cones, below
                              # 1 flares immediately into a bell.
FUNNEL_TOP_LEVEL     = 1.5    # brightness at the exciton end
FUNNEL_WAIST_LEVEL   = 1.8    # brightness at the pinch, where the light is
                              # at its most concentrated
FUNNEL_BOTTOM_LEVEL  = 0.5    # brightness where it lands
FUNNEL_LEVEL         = 1.0    # overall multiplier, driven by the binding
                              # strength
FUNNEL_ALPHA         = 0.5
FUNNEL_FACING_FALLOFF = 1.6   # fades the cone out at its own silhouette, so
                              # it reads as a beam instead of a solid shell
FUNNEL_COLOR         = ""     # "" = the hole's highlight colour
FUNNEL_RINGS         = 40     # steps down the funnel
FUNNEL_SEGMENTS      = 56     # steps around it

SHOW_FUNNEL_POOL     = True   # the disc of light where the funnel lands
FUNNEL_POOL_RADIUS   = 5.0
FUNNEL_POOL_LEVEL    = 1.1
FUNNEL_POOL_ALPHA    = 0.55
FUNNEL_POOL_FALLOFF  = 2.2    # >1 pulls the pool in tighter around the centre
FUNNEL_POOL_STEPS    = 32

# --- field lines (optional) -------------------------------------------------
SHOW_FIELD_LINES     = False  # the full dipole streamline bundle. Off by
                              # default now that the tail and the aura carry
                              # the picture; turn it on for the
                              # physics-diagram version.
FIELD_LINE_RINGS     = 7      # launch angles (rings of lines)
FIELD_LINE_AZIMUTHS  = 18     # lines around the axis per ring
FIELD_LINE_AXIS_BIAS = 2.6    # >1 pushes launch angles toward the axis, which
                              # is what crowds the lines into the middle.
                              # 1.0 = even spread, 4.0 = very dense cord.
LAUNCH_ANGLE_MIN     = 4.0    # degrees from the axis
LAUNCH_ANGLE_MAX     = 78.0   # degrees; beyond ~90 lines swing far out

FIELD_LINE_RADIUS    = 0.045  # tube radius
FIELD_LINE_SEGMENTS  = 7      # cross-section resolution of each tube
FIELD_LINE_EMISSION  = 1.6    # scaled to the same range as the particles
MID_GLOW_BOOST       = 2.2    # extra brightness at the midpoint of each line
                              # (this is what makes the central cord white-hot)
MID_GLOW_WIDTH       = 0.30   # 0..1, how much of the line the boost covers

INTEGRATION_STEP     = 0.045  # RK4 step length; smaller = smoother, slower
MAX_STEPS            = 20000
FIELD_LINE_TAPER     = 0.12   # fraction of each end that tapers to a point

# --- particle shading (this is what makes them read as 3-D spheres) --------
PARTICLE_SHADING     = 0.75   # 0 = flat disc, 1 = fully shaded sphere
PARTICLE_LIGHT_DIR   = (-0.55, -1.0, 0.6)  # direction the fake key light
                              # comes from, in world space. It only drives the
                              # gradient -- no actual lamp is created.
PARTICLE_TERMINATOR  = 0.42   # 0..1, how bright the unlit side stays
PARTICLE_CORE_BOOST  = 1.35   # overall multiplier on the lit side
PARTICLE_RIM_POWER   = 3.0    # tightness of the bright silhouette rim
PARTICLE_RIM_BOOST   = 1.25   # how much the rim brightens the edge

# --- glow around each particle (one smooth gradient, no shells) ------------
SHOW_PARTICLE_GLOW   = True
GLOW_OUTER_SCALE     = 1.6    # halo radius / particle radius. Kept tight: a
                              # big soft halo swallows the particle it is
                              # meant to be lighting, and where two of them
                              # overlap the middle washes out to white. Let
                              # the compositor bloom do the soft work.
GLOW_ALPHA           = 0.60   # peak alpha of the halo, reached just outside
                              # the particle's edge
GLOW_FALLOFF         = 2.6    # >1 fades faster toward the halo's outer edge
GLOW_EMISSION        = 1.8

# --- aura enclosing the whole exciton (a single emissive volume) -----------
SHOW_AURA            = True
AURA_MARGIN          = 1.2        # how far the aura reaches past the two
                                  # particles. The aura is a spindle sitting
                                  # directly between the hole and the
                                  # electron, so this is what makes it a tight
                                  # bridge of light or a loose halo.
AURA_DENSITY         = 0.45       # density at the core of the aura. This
                                  # and AURA_EMISSION multiply; too much
                                  # density hazes over the particles instead
                                  # of glowing around them.
AURA_FALLOFF         = 2.4        # >1 = tighter, faster fade to nothing
AURA_EMISSION        = 2.5
AURA_RESOLUTION      = (64, 32)   # (azimuthal, polar) segments of the shell
                                  # that bounds the volume

# --- scene / render ---------------------------------------------------------
SETUP_SCENE            = True   # camera, render settings, and bloom
BACKGROUND_MODE        = "transparent"
                                # "transparent" -> alpha channel, so this
                                #   render drops straight onto a slide or over
                                #   the cucro2perc render. This is the overlay
                                #   of the pair, so it defaults to
                                #   transparent while cucro2perc carries the
                                #   blue underneath it.
                                # "gradient"    -> the same deep blue field as
                                #   Reference 1, for using this figure on its
                                #   own
                                # "flat"        -> a single BACKGROUND_COLOR
BACKGROUND_COLOR       = "#0A2A5E"  # deep blue, sampled from Reference 1
BACKGROUND_TOP_COLOR   = "#12539E"  # brighter blue at the top of the gradient
BACKGROUND_MARGIN      = 1.15   # how far the backdrop oversizes the frame
BACKGROUND_STEPS       = 48     # rows the gradient is built from
BACKGROUND_SMOOTH      = True   # ease the gradient instead of ramping it
CAMERA_ELEVATION       = 20.0   # degrees above the orbit plane. A camera
                                # level with that plane sees the tail exactly
                                # edge-on, as a straight line; this is what
                                # opens it out into an ellipse.
ADD_GLOW_COMPOSITOR    = True
GLOW_IN_ALPHA          = True   # also push the bloom into the alpha channel;
                                # without this the glow around the exciton is
                                # invisible once the PNG is composited, and in
                                # "gradient" mode it would not blend onto the
                                # blue either
GLOW_ALPHA_GAIN        = 1.6    # how strongly the bloom opens up the alpha
RENDER_ENGINE          = "EEVEE"   # "EEVEE" or "CYCLES". Both render the
                                   # volumetric aura; Cycles is cleaner and
                                   # slower.
CAMERA_LENS            = 50.0   # focal length in mm
CAMERA_MARGIN          = 1.18   # >1 leaves air around the exciton. The camera
                                # fits the geometry exactly, so this is pure
                                # breathing room.
RESOLUTION             = (2200, 2200)
RENDER_SAMPLES         = 128
VIEW_TRANSFORM         = "Standard"   # "Standard" keeps the colours exactly
                                      # as specified above, which is what a
                                      # figure wants. "AgX" or "Filmic" roll
                                      # the highlights off photographically
                                      # but desaturate the bright cores.
OUTPUT_PATH            = ""     # e.g. r"C:\figures\exciton.png". Empty just
                                # sets the scene up without rendering.
RENDER_NOW             = False  # True = render straight to OUTPUT_PATH
CLEAR_SCENE            = True

SPHERE_SEGMENTS      = 64     # particle mesh resolution
SPHERE_RINGS         = 48


# ============================================================================
#  Colour helpers
# ============================================================================

def hex_to_linear(h):
    """sRGB hex -> linear RGB, which is what Blender shader inputs expect."""
    h = h.strip().lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    srgb = [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]

    def lin(u):
        return u / 12.92 if u <= 0.04045 else ((u + 0.055) / 1.055) ** 2.4
    return tuple(lin(u) for u in srgb)


def mix(c1, c2, t):
    return tuple(a + (b - a) * t for a, b in zip(c1, c2))


# ============================================================================
#  Binding strength
# ============================================================================
#  BINDING_STRENGTH blends the model between a tightly bound pair and a barely
#  bound one. Each entry below names a parameter and the WEAK_ parameter that
#  holds its value at zero binding; the value written in USER PARAMETERS is
#  the strong end.
# ============================================================================

_STRENGTH_DRIVEN = {
    "ORBIT_RADIUS": "WEAK_ORBIT_RADIUS",
    "AURA_DENSITY": "WEAK_AURA_DENSITY",
    "AURA_EMISSION": "WEAK_AURA_EMISSION",
    "AURA_MARGIN": "WEAK_AURA_MARGIN",
    "TRAIL_LENGTH": "WEAK_TRAIL_LENGTH",
    "TRAIL_LEVEL": "WEAK_TRAIL_LEVEL",
    "GLOW_OUTER_SCALE": "WEAK_GLOW_OUTER_SCALE",
    "GLOW_ALPHA": "WEAK_GLOW_ALPHA",
    "FUNNEL_LEVEL": "WEAK_FUNNEL_LEVEL",
}

# emission strengths are scaled rather than replaced, so the balance between
# them survives
_STRENGTH_SCALED = ("ELECTRON_EMISSION", "HOLE_EMISSION", "TRAIL_EMISSION",
                    "GLOW_EMISSION", "FUNNEL_POOL_LEVEL")

_STRONG_REFERENCE = {}


def binding():
    """The binding strength, clamped to 0..1."""
    return min(max(BINDING_STRENGTH, 0.0), 1.0)


def apply_binding_strength():
    """
    Rewrite the strength-driven parameters for the current BINDING_STRENGTH.

    The first call records the values as written -- the strong end -- so that
    calling it again re-blends from those rather than from whatever the last
    call left behind.
    """
    g = globals()
    if not _STRONG_REFERENCE:
        for name in list(_STRENGTH_DRIVEN) + list(_STRENGTH_SCALED):
            _STRONG_REFERENCE[name] = g[name]

    t = binding()
    for name, weak_name in _STRENGTH_DRIVEN.items():
        weak = g[weak_name]
        g[name] = weak + (_STRONG_REFERENCE[name] - weak) * t
    scale = WEAK_EMISSION_SCALE + (1.0 - WEAK_EMISSION_SCALE) * t
    for name in _STRENGTH_SCALED:
        g[name] = _STRONG_REFERENCE[name] * scale
    return t


def pair_color(hex_string):
    """
    A colour from the palette, washed toward grey as the binding weakens.

    Desaturation is what actually makes a loosely bound pair recede: dimming
    alone just makes a small dark version of the same thing, while a washed
    out one reads as faint even next to a bright neighbour.
    """
    c = hex_to_linear(hex_string)
    fade = (1.0 - binding()) * WEAK_DESATURATION
    if fade <= 1e-6:
        return c
    grey = sum(c) / 3.0
    return mix(c, (grey, grey, grey), fade)


# ============================================================================
#  Dipole field
# ============================================================================

def dipole_field(p, r_pos, r_neg):
    """E of +1 at r_pos and -1 at r_neg (Gaussian units, k = 1)."""
    d1 = p - r_pos
    d2 = p - r_neg
    n1 = d1.length
    n2 = d2.length
    if n1 < 1e-9 or n2 < 1e-9:
        return Vector((0.0, 0.0, 0.0))
    return d1 / (n1 ** 3) - d2 / (n2 ** 3)


def _unit_field(p, r_pos, r_neg):
    e = dipole_field(p, r_pos, r_neg)
    L = e.length
    return e / L if L > 1e-12 else Vector((0.0, 0.0, 0.0))


def trace_field_line(start, r_pos, r_neg, r_stop, step=None, max_steps=None):
    """
    RK4-integrate a streamline of the dipole field from `start` until it
    reaches within r_stop of the negative charge. Returns a list of points,
    or None if the line failed to terminate.
    """
    step = INTEGRATION_STEP if step is None else step
    max_steps = MAX_STEPS if max_steps is None else max_steps
    bound = (r_pos - r_neg).length * 12.0 + 10.0
    center = (r_pos + r_neg) * 0.5

    pts = [start.copy()]
    p = start.copy()
    for _ in range(max_steps):
        k1 = _unit_field(p, r_pos, r_neg)
        if k1.length < 0.5:
            return None
        k2 = _unit_field(p + k1 * (step * 0.5), r_pos, r_neg)
        k3 = _unit_field(p + k2 * (step * 0.5), r_pos, r_neg)
        k4 = _unit_field(p + k3 * step, r_pos, r_neg)
        d = (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if d.length < 1e-9:
            return None
        p = p + d.normalized() * step
        pts.append(p.copy())

        if (p - r_neg).length <= r_stop:
            return pts
        if (p - center).length > bound:
            return None
    return None


def launch_directions(axis, n_rings, n_azim, a_min_deg, a_max_deg, bias):
    """
    Launch directions on the sphere around the positive charge.

    The polar angle is measured from `axis` (which points at the other
    charge). Angles are distributed as a_min + (a_max-a_min) * u**bias with u
    evenly spaced, so bias > 1 crowds the launches toward the axis -- those
    are exactly the lines that run through the middle of the pair.
    """
    axis = Vector(axis).normalized()
    ref = Vector((0.0, 0.0, 1.0))
    if abs(axis.dot(ref)) > 0.95:
        ref = Vector((1.0, 0.0, 0.0))
    u1 = axis.cross(ref).normalized()
    u2 = axis.cross(u1).normalized()

    a_min, a_max = math.radians(a_min_deg), math.radians(a_max_deg)
    dirs = []
    for i in range(n_rings):
        u = (i / (n_rings - 1)) if n_rings > 1 else 0.0
        alpha = a_min + (a_max - a_min) * (u ** bias)
        # stagger every other ring so the lines don't stack into sheets
        phase = (pi / n_azim) if (i % 2) else 0.0
        for j in range(n_azim):
            phi = 2.0 * pi * j / n_azim + phase
            d = (axis * cos(alpha)
                 + (u1 * cos(phi) + u2 * sin(phi)) * sin(alpha))
            dirs.append((d.normalized(), u))
    return dirs


def build_field_lines(r_electron, r_hole):
    """Trace every field line from the hole (+) to the electron (-)."""
    axis_to_electron = (r_electron - r_hole).normalized()
    dirs = launch_directions(axis_to_electron, FIELD_LINE_RINGS,
                             FIELD_LINE_AZIMUTHS, LAUNCH_ANGLE_MIN,
                             LAUNCH_ANGLE_MAX, FIELD_LINE_AXIS_BIAS)
    start_r = HOLE_RADIUS * 1.02
    stop_r = ELECTRON_RADIUS * 1.02

    lines, failed = [], 0
    for d, ring_u in dirs:
        pts = trace_field_line(r_hole + d * start_r, r_hole, r_electron, stop_r)
        if pts and len(pts) > 3:
            lines.append(pts)
        else:
            failed += 1
    return lines, failed


# ============================================================================
#  Mesh construction
# ============================================================================

def sphere_mesh_data(radius, segs, rings, scale=(1.0, 1.0, 1.0)):
    """UV sphere as (verts, faces), optionally scaled into an ellipsoid."""
    sx, sy, sz = scale
    verts = [(0.0, 0.0, radius * sz)]
    for i in range(1, rings):
        phi = pi * i / rings
        z, r = radius * cos(phi), radius * sin(phi)
        for j in range(segs):
            th = 2.0 * pi * j / segs
            verts.append((r * cos(th) * sx, r * sin(th) * sy, z * sz))
    verts.append((0.0, 0.0, -radius * sz))
    last = len(verts) - 1

    def idx(i, j):
        return 1 + (i - 1) * segs + (j % segs)

    faces = [(0, idx(1, j), idx(1, j + 1)) for j in range(segs)]
    for i in range(1, rings - 1):
        for j in range(segs):
            faces.append((idx(i, j), idx(i + 1, j),
                          idx(i + 1, j + 1), idx(i, j + 1)))
    faces += [(last, idx(rings - 1, j + 1), idx(rings - 1, j))
              for j in range(segs)]
    return verts, faces


def _initial_normal(t):
    """Any unit vector perpendicular to tangent t."""
    ref = Vector((0.0, 0.0, 1.0))
    if abs(t.dot(ref)) > 0.9:
        ref = Vector((1.0, 0.0, 0.0))
    return t.cross(ref).normalized()


def tube_from_polyline(points, radius, segs, taper=0.0):
    """
    Sweep a circular cross-section along a polyline using parallel transport,
    so the tube does not twist. Returns (verts, faces, param) where param is
    the 0..1 position along the line for each vertex (used for colouring).

    `radius` is either one number for the whole tube or a radius per point,
    which is what lets the orbit trail swell at the electron and thin out
    behind it.
    """
    n = len(points)
    if n < 2:
        return [], [], []

    tangents = []
    for i in range(n):
        if i == 0:
            t = points[1] - points[0]
        elif i == n - 1:
            t = points[-1] - points[-2]
        else:
            t = points[i + 1] - points[i - 1]
        if t.length < 1e-9:
            t = Vector((0.0, 0.0, 1.0))
        tangents.append(t.normalized())

    # parallel transport the frame along the curve
    normals = [_initial_normal(tangents[0])]
    for i in range(1, n):
        prev_t, cur_t = tangents[i - 1], tangents[i]
        v = prev_t.cross(cur_t)
        nv = normals[-1]
        if v.length < 1e-9:
            normals.append(nv)
            continue
        angle = math.atan2(v.length, prev_t.dot(cur_t))
        rot = Matrix.Rotation(angle, 3, v.normalized())
        nrm = (rot @ nv)
        nrm = (nrm - cur_t * nrm.dot(cur_t))
        normals.append(nrm.normalized() if nrm.length > 1e-9 else nv)

    # arclength parameterisation, for colour and taper
    lengths = [0.0]
    for i in range(1, n):
        lengths.append(lengths[-1] + (points[i] - points[i - 1]).length)
    total = lengths[-1] if lengths[-1] > 1e-9 else 1.0
    params = [L / total for L in lengths]

    def taper_scale(s):
        if taper <= 0.0:
            return 1.0
        e = min(s, 1.0 - s) / taper
        if e >= 1.0:
            return 1.0
        return e * e * (3.0 - 2.0 * e)          # smoothstep

    varying = not isinstance(radius, (int, float))

    verts, vparam = [], []
    for i in range(n):
        t = tangents[i]
        u = normals[i]
        w = t.cross(u).normalized()
        r = (radius[min(i, len(radius) - 1)] if varying else radius)
        r *= taper_scale(params[i])
        for j in range(segs):
            a = 2.0 * pi * j / segs
            p = points[i] + u * (r * cos(a)) + w * (r * sin(a))
            verts.append((p.x, p.y, p.z))
            vparam.append(params[i])

    faces = []
    for i in range(n - 1):
        b0, b1 = i * segs, (i + 1) * segs
        for j in range(segs):
            j2 = (j + 1) % segs
            faces.append((b0 + j, b0 + j2, b1 + j2, b1 + j))
    return verts, faces, vparam


# ============================================================================
#  Blender scene helpers
# ============================================================================

def clear_scene():
    for obj in list(bpy.data.objects):
        if obj.type in {"MESH", "CAMERA", "LIGHT"}:
            bpy.data.objects.remove(obj, do_unlink=True)
    for block in (bpy.data.meshes, bpy.data.materials):
        for b in list(block):
            if b.users == 0:
                block.remove(b)


def get_collection(name):
    col = bpy.data.collections.get(name)
    if col is None:
        col = bpy.data.collections.new(name)
        bpy.context.scene.collection.children.link(col)
    return col


def _set_blend(mat, alpha, single_layer=True):
    """Enable alpha blending in a way that works across Blender versions."""
    if alpha >= 1.0:
        return
    for attr, val in (("blend_method", "BLEND"),
                      ("shadow_method", "NONE"),
                      ("surface_render_method", "BLENDED")):
        try:
            setattr(mat, attr, val)
        except Exception:
            pass
    if single_layer:
        # only the camera-facing half of a transparent sphere contributes, so
        # the gradient is not doubled up by the far side of the same surface
        for attr in ("show_transparent_back", ):
            try:
                setattr(mat, attr, False)
            except Exception:
                pass
        try:
            mat.use_backface_culling = True
        except Exception:
            pass


def emission_material(name, color, strength, alpha=1.0):
    """Emission shader, optionally alpha-blended."""
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
    emis.location = (0, 0)
    emis.inputs["Color"].default_value = (*color, 1.0)
    emis.inputs["Strength"].default_value = strength

    if alpha >= 1.0:
        nt.links.new(emis.outputs["Emission"], out.inputs["Surface"])
    else:
        trans = nt.nodes.new("ShaderNodeBsdfTransparent")
        trans.location = (0, -150)
        mixer = nt.nodes.new("ShaderNodeMixShader")
        mixer.location = (150, 0)
        mixer.inputs["Fac"].default_value = alpha
        nt.links.new(trans.outputs["BSDF"], mixer.inputs[1])
        nt.links.new(emis.outputs["Emission"], mixer.inputs[2])
        nt.links.new(mixer.outputs["Shader"], out.inputs["Surface"])
    mat.diffuse_color = (*color, alpha)
    _set_blend(mat, alpha)
    return mat


def vertex_color_emission_material(name, attr_name, strength, alpha=1.0):
    """Emission driven by a per-vertex colour attribute."""
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (450, 0)
    emis = nt.nodes.new("ShaderNodeEmission")
    emis.location = (100, 0)
    emis.inputs["Strength"].default_value = strength
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.location = (-150, 0)
    attr.attribute_name = attr_name
    nt.links.new(attr.outputs["Color"], emis.inputs["Color"])
    if alpha >= 1.0:
        nt.links.new(emis.outputs["Emission"], out.inputs["Surface"])
    else:
        trans = nt.nodes.new("ShaderNodeBsdfTransparent")
        trans.location = (100, -160)
        mixer = nt.nodes.new("ShaderNodeMixShader")
        mixer.location = (280, 0)
        mixer.inputs["Fac"].default_value = alpha
        nt.links.new(trans.outputs["BSDF"], mixer.inputs[1])
        nt.links.new(emis.outputs["Emission"], mixer.inputs[2])
        nt.links.new(mixer.outputs["Shader"], out.inputs["Surface"])
        _set_blend(mat, alpha, single_layer=False)
    return mat


# ============================================================================
#  Gradient shader helpers
#
#  Everything below builds node graphs by hand. Two small ideas do all the
#  work:
#
#    * |dot(N, I)| is the cosine of the angle between the surface normal and
#      the view ray: 1 where a surface faces the camera head-on and 0 exactly
#      at its silhouette. Driving alpha or brightness with it gives a smooth
#      radial gradient across a sphere with no visible edge -- which is what
#      replaces the old stacks of nested shells.
#
#    * dot(N, L) against a fixed direction is a lambert term. Wrapped to
#      0..1 ("half lambert") it shades a sphere from a lit side to a shaded
#      side without any lamp in the scene, so an emissive particle stops
#      looking like a flat coloured circle.
# ============================================================================

def _node(nt, kind, x, y):
    n = nt.nodes.new(kind)
    n.location = (x, y)
    return n


def _set_in(node, name, value):
    """Set an input by name, ignoring names this Blender version does not have."""
    try:
        node.inputs[name].default_value = value
        return True
    except Exception:
        return False


def _link_in(nt, socket, node, name):
    try:
        nt.links.new(socket, node.inputs[name])
        return True
    except Exception:
        return False


def _math(nt, op, x, y, value1=None, value2=None, value3=None, clamp=False):
    n = _node(nt, "ShaderNodeMath", x, y)
    n.operation = op
    n.use_clamp = clamp
    for k, v in ((0, value1), (1, value2), (2, value3)):
        if v is not None and len(n.inputs) > k:
            n.inputs[k].default_value = v
    return n


def _facing(nt, x, y):
    """(geometry node, facing value) with facing = 1 head-on, 0 at the edge."""
    geo = _node(nt, "ShaderNodeNewGeometry", x, y)
    dot = _node(nt, "ShaderNodeVectorMath", x + 190, y)
    dot.operation = "DOT_PRODUCT"
    nt.links.new(geo.outputs["Normal"], dot.inputs[0])
    nt.links.new(geo.outputs["Incoming"], dot.inputs[1])
    fac = _math(nt, "ABSOLUTE", x + 360, y, clamp=True)
    nt.links.new(dot.outputs["Value"], fac.inputs[0])
    return geo, fac


def _two_stop_ramp(nt, x, y, c0, c1, c_mid=None, mid_pos=0.55):
    ramp = _node(nt, "ShaderNodeValToRGB", x, y)
    e0, e1 = ramp.color_ramp.elements[0], ramp.color_ramp.elements[1]
    e0.position = 0.0
    e0.color = (*c0, 1.0)
    e1.position = 1.0
    e1.color = (*c1, 1.0)
    if c_mid is not None:
        em = ramp.color_ramp.elements.new(mid_pos)
        em.color = (*c_mid, 1.0)
    return ramp


def particle_material(name, core_color, rim_color, strength):
    """
    A glowing sphere that actually looks like a sphere.

    The emission colour is read off a ramp driven by a half-lambert gradient,
    and the emission strength is that same gradient times a rim term that
    brightens the silhouette. No scene light is involved, so it renders
    identically in Material Preview, EEVEE and Cycles.
    """
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    geo, facing = _facing(nt, -1180, -280)

    # --- half-lambert shading from a fixed, scene-independent direction ----
    L = Vector(PARTICLE_LIGHT_DIR)
    L = L.normalized() if L.length > 1e-9 else Vector((0.0, 0.0, 1.0))
    ldot = _node(nt, "ShaderNodeVectorMath", -820, 200)
    ldot.operation = "DOT_PRODUCT"
    ldot.inputs[1].default_value = (L.x, L.y, L.z)
    nt.links.new(geo.outputs["Normal"], ldot.inputs[0])

    half = _math(nt, "MULTIPLY_ADD", -640, 200,
                 value2=0.5, value3=0.5, clamp=True)     # (dot + 1) / 2
    nt.links.new(ldot.outputs["Value"], half.inputs[0])

    t = min(max(PARTICLE_TERMINATOR, 0.0), 1.0)
    lit = _math(nt, "MULTIPLY_ADD", -460, 200,
                value2=1.0 - t, value3=t, clamp=True)    # never fully black
    nt.links.new(half.outputs["Value"], lit.inputs[0])

    s = min(max(PARTICLE_SHADING, 0.0), 1.0)
    grad = _math(nt, "MULTIPLY_ADD", -280, 200,
                 value2=s, value3=1.0 - s, clamp=True)   # lerp(1, lit, s)
    nt.links.new(lit.outputs["Value"], grad.inputs[0])

    ramp = _two_stop_ramp(nt, -80, 220,
                          tuple(c * 0.4 for c in core_color),
                          rim_color, c_mid=core_color, mid_pos=0.6)
    nt.links.new(grad.outputs["Value"], ramp.inputs["Fac"])

    # --- rim light at the silhouette --------------------------------------
    inv = _math(nt, "SUBTRACT", -280, -280, value1=1.0, clamp=True)
    nt.links.new(facing.outputs["Value"], inv.inputs[1])
    rim = _math(nt, "POWER", -100, -280, value2=max(PARTICLE_RIM_POWER, 0.01))
    nt.links.new(inv.outputs["Value"], rim.inputs[0])
    rim_mul = _math(nt, "MULTIPLY_ADD", 80, -280,
                    value2=PARTICLE_RIM_BOOST, value3=1.0)
    nt.links.new(rim.outputs["Value"], rim_mul.inputs[0])

    base = _math(nt, "MULTIPLY", 80, 20, value2=strength * PARTICLE_CORE_BOOST)
    nt.links.new(grad.outputs["Value"], base.inputs[0])
    total = _math(nt, "MULTIPLY", 280, -120)
    nt.links.new(base.outputs["Value"], total.inputs[0])
    nt.links.new(rim_mul.outputs["Value"], total.inputs[1])

    emis = _node(nt, "ShaderNodeEmission", 480, 80)
    nt.links.new(ramp.outputs["Color"], emis.inputs["Color"])
    nt.links.new(total.outputs["Value"], emis.inputs["Strength"])
    out = _node(nt, "ShaderNodeOutputMaterial", 680, 80)
    nt.links.new(emis.outputs["Emission"], out.inputs["Surface"])

    mat.diffuse_color = (*core_color, 1.0)
    return mat


def halo_material(name, inner_color, outer_color, strength, alpha, falloff):
    """
    One sphere whose alpha follows |dot(N, I)| ** falloff, so it reaches
    exactly zero at its own silhouette -- a continuous halo gradient in place
    of the old nested shells, with no outline and none of the faint concentric
    rings a stack of shells leaves behind.

    Only the half of the sphere facing *away* from the camera is shaded (the
    Backfacing term). The particle itself is opaque and sits at the centre, so
    that far half is hidden everywhere the particle covers it and visible
    everywhere outside it: the halo glows around the particle instead of
    washing a flat veil across it and flattening its shading back out.
    """
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    geo, facing = _facing(nt, -900, 0)
    prof = _math(nt, "POWER", -460, 0, value2=max(falloff, 0.01))
    nt.links.new(facing.outputs["Value"], prof.inputs[0])
    scaled = _math(nt, "MULTIPLY", -280, 0, value2=max(min(alpha, 1.0), 0.0),
                   clamp=True)
    nt.links.new(prof.outputs["Value"], scaled.inputs[0])
    a = _math(nt, "MULTIPLY", -100, -60, clamp=True)
    nt.links.new(scaled.outputs["Value"], a.inputs[0])
    nt.links.new(geo.outputs["Backfacing"], a.inputs[1])

    ramp = _two_stop_ramp(nt, -280, 220, outer_color, inner_color)
    nt.links.new(prof.outputs["Value"], ramp.inputs["Fac"])

    emis = _node(nt, "ShaderNodeEmission", 20, 200)
    emis.inputs["Strength"].default_value = strength
    nt.links.new(ramp.outputs["Color"], emis.inputs["Color"])
    trans = _node(nt, "ShaderNodeBsdfTransparent", 20, 40)
    mixer = _node(nt, "ShaderNodeMixShader", 220, 120)
    nt.links.new(a.outputs["Value"], mixer.inputs["Fac"])
    nt.links.new(trans.outputs["BSDF"], mixer.inputs[1])
    nt.links.new(emis.outputs["Emission"], mixer.inputs[2])
    out = _node(nt, "ShaderNodeOutputMaterial", 420, 120)
    nt.links.new(mixer.outputs["Shader"], out.inputs["Surface"])

    mat.diffuse_color = (*inner_color, alpha)
    # the far side of the halo is the side that renders, so back faces must
    # not be culled
    _set_blend(mat, alpha, single_layer=False)
    return mat


def beam_material(name, attr_name, strength, alpha, facing_falloff=0.0):
    """
    Emission from a vertex colour, faded out by that colour's own alpha.

    With `facing_falloff` above zero the alpha is also multiplied by
    |dot(N, I)|, so the surface disappears at its own silhouette. That is what
    turns the funnel from a visible cone-shaped shell into something that
    reads as a shaft of light: no outline anywhere, brightest where you are
    looking straight through the most of it.
    """
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    attr = _node(nt, "ShaderNodeAttribute", -700, 0)
    attr.attribute_name = attr_name

    a = _math(nt, "MULTIPLY", -300, -140, value2=max(min(alpha, 1.0), 0.0),
              clamp=True)
    nt.links.new(attr.outputs["Alpha"], a.inputs[0])
    if facing_falloff > 0.0:
        _, facing = _facing(nt, -900, -340)
        prof = _math(nt, "POWER", -480, -340, value2=facing_falloff)
        nt.links.new(facing.outputs["Value"], prof.inputs[0])
        shaped = _math(nt, "MULTIPLY", -120, -140, clamp=True)
        nt.links.new(a.outputs["Value"], shaped.inputs[0])
        nt.links.new(prof.outputs["Value"], shaped.inputs[1])
        a = shaped

    emis = _node(nt, "ShaderNodeEmission", -120, 120)
    emis.inputs["Strength"].default_value = strength
    nt.links.new(attr.outputs["Color"], emis.inputs["Color"])
    trans = _node(nt, "ShaderNodeBsdfTransparent", -120, -20)
    mixer = _node(nt, "ShaderNodeMixShader", 120, 60)
    nt.links.new(a.outputs["Value"], mixer.inputs["Fac"])
    nt.links.new(trans.outputs["BSDF"], mixer.inputs[1])
    nt.links.new(emis.outputs["Emission"], mixer.inputs[2])
    out = _node(nt, "ShaderNodeOutputMaterial", 320, 60)
    nt.links.new(mixer.outputs["Shader"], out.inputs["Surface"])

    # both walls of the funnel should show through each other
    _set_blend(mat, min(alpha, 0.999), single_layer=False)
    return mat


def aura_volume_material(name, hole_color, electron_color, density, falloff,
                         strength):
    """
    Emissive volume whose density falls off from the centre of the aura to its
    surface, tinted along its length from the hole's colour at the hole end to
    the electron's at the electron end.

    Why the old version rendered as nothing: it fed the *world-space* distance
    from the origin straight into a colour ramp. A ramp's Fac is clamped to
    0..1, and the aura is several Blender units across, so every point more
    than one unit from the centre came out at density zero. Here the aura mesh
    is a plain unit sphere and its size lives in the object's scale instead,
    so Texture Coordinate -> Object hands back a normalised radius that really
    does run 0 at the core to 1 at the surface. The falloff is then explicit
    math -- clamp(1 - r) ** AURA_FALLOFF -- instead of a two-stop ramp.
    """
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    tex = _node(nt, "ShaderNodeTexCoord", -1000, 0)
    length = _node(nt, "ShaderNodeVectorMath", -800, -80)
    length.operation = "LENGTH"
    nt.links.new(tex.outputs["Object"], length.inputs[0])

    inv = _math(nt, "SUBTRACT", -620, -80, value1=1.0, clamp=True)  # 1 - r
    nt.links.new(length.outputs["Value"], inv.inputs[1])
    prof = _math(nt, "POWER", -440, -80, value2=max(falloff, 0.01))
    nt.links.new(inv.outputs["Value"], prof.inputs[0])

    dens = _math(nt, "MULTIPLY", -260, -160, value2=max(density, 0.0))
    nt.links.new(prof.outputs["Value"], dens.inputs[0])
    emit = _math(nt, "MULTIPLY", -260, -320, value2=max(strength, 0.0))
    nt.links.new(prof.outputs["Value"], emit.inputs[0])

    # colour gradient along the spindle: object Z runs -1 at the hole end to
    # +1 at the electron end, because the object is rotated to put its local
    # +Z along the hole -> electron direction
    sep = _node(nt, "ShaderNodeSeparateXYZ", -800, 220)
    nt.links.new(tex.outputs["Object"], sep.inputs[0])
    zt = _math(nt, "MULTIPLY_ADD", -620, 220, value2=0.5, value3=0.5,
               clamp=True)
    nt.links.new(sep.outputs["Z"], zt.inputs[0])
    ramp = _two_stop_ramp(nt, -440, 220, hole_color, electron_color)
    nt.links.new(zt.outputs["Value"], ramp.inputs["Fac"])

    vol = _node(nt, "ShaderNodeVolumePrincipled", 0, 0)
    _set_in(vol, "Anisotropy", 0.0)
    _link_in(nt, ramp.outputs["Color"], vol, "Color")
    _link_in(nt, ramp.outputs["Color"], vol, "Emission Color")
    _link_in(nt, dens.outputs["Value"], vol, "Density")
    _link_in(nt, emit.outputs["Value"], vol, "Emission Strength")

    out = _node(nt, "ShaderNodeOutputMaterial", 260, 0)
    nt.links.new(vol.outputs["Volume"], out.inputs["Volume"])
    return mat


def add_object(name, verts, faces, mat, collection, colors=None,
               attr_name="LineColor"):
    if not verts or not faces:
        return None
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.validate()
    for poly in mesh.polygons:
        poly.use_smooth = True
    if mat:
        mesh.materials.append(mat)
    if colors is not None:
        try:
            attr = mesh.color_attributes.new(name=attr_name,
                                             type="FLOAT_COLOR",
                                             domain="POINT")
            for i, c in enumerate(colors):
                # a fourth component is per-vertex alpha, which is what the
                # funnel and its pool fade themselves out with
                attr.data[i].color = (c[0], c[1], c[2],
                                      c[3] if len(c) > 3 else 1.0)
        except Exception as e:
            print(f"[exciton] could not write colour attribute: {e}")
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    return obj


# ============================================================================
#  Builders
# ============================================================================

def build_particles(col, r_electron, r_hole):
    """The electron and the hole, shaded so they read as spheres."""
    for label, centre, radius, core_hex, rim_hex, strength in (
            ("Electron", r_electron, ELECTRON_RADIUS, ELECTRON_COLOR,
             ELECTRON_RIM, ELECTRON_EMISSION),
            ("Hole", r_hole, HOLE_RADIUS, HOLE_COLOR, HOLE_RIM,
             HOLE_EMISSION)):
        mat = particle_material(label, pair_color(core_hex),
                                pair_color(rim_hex), strength)
        v, f = sphere_mesh_data(radius, SPHERE_SEGMENTS, SPHERE_RINGS)
        obj = add_object(f"Exciton_{label}", v, f, mat, col)
        if obj:
            obj.location = centre


def build_particle_glow(col, r_electron, r_hole):
    """
    One halo sphere per particle, fading out as a continuous gradient.

    This replaces the old nested GLOW_SHELLS: a single surface whose alpha
    follows |dot(N, I)| ** GLOW_FALLOFF is transparent exactly at its own
    silhouette, so the halo has no outline at all -- whereas a stack of
    shells always leaves faint concentric rings where they overlap.
    """
    for label, centre, radius, core_hex, rim_hex in (
            ("Electron", r_electron, ELECTRON_RADIUS, ELECTRON_COLOR,
             ELECTRON_RIM),
            ("Hole", r_hole, HOLE_RADIUS, HOLE_COLOR, HOLE_RIM)):
        mat = halo_material(f"{label}_Halo", pair_color(rim_hex),
                            pair_color(core_hex), GLOW_EMISSION,
                            GLOW_ALPHA, GLOW_FALLOFF)
        v, f = sphere_mesh_data(1.0, 48, 32)
        obj = add_object(f"Exciton_{label}_Halo", v, f, mat, col)
        if obj:
            obj.location = centre
            r = radius * GLOW_OUTER_SCALE
            obj.scale = (r, r, r)


def build_field_line_objects(col, lines):
    """One merged mesh for all field lines, coloured per vertex."""
    e_col = hex_to_linear(ELECTRON_COLOR)
    h_col = hex_to_linear(FIELD_LINE_COLOR)
    V, F, C = [], [], []
    for pts in lines:
        v, f, params = tube_from_polyline(pts, FIELD_LINE_RADIUS,
                                          FIELD_LINE_SEGMENTS,
                                          FIELD_LINE_TAPER)
        if not v:
            continue
        base = len(V)
        V.extend(v)
        F.extend(tuple(base + i for i in face) for face in f)
        for s in params:
            c = mix(h_col, e_col, s)
            # brighten the middle of the line into a white-hot cord
            g = math.exp(-((s - 0.5) / max(MID_GLOW_WIDTH, 1e-6)) ** 2)
            boost = 1.0 + (MID_GLOW_BOOST - 1.0) * g
            C.append(tuple(min(v_ * boost, 8.0) for v_ in c))
    mat = vertex_color_emission_material("FieldLines", "LineColor",
                                         FIELD_LINE_EMISSION)
    add_object("Exciton_FieldLines", V, F, mat, col, colors=C)


def camera_basis(axis):
    """
    (toward-camera, screen-horizontal) unit vectors for the default camera.

    setup_camera places the camera along `side`, so `side` points from the
    exciton at the viewer and `right` runs across the frame. The orbit is
    built on this basis, which is what lets ORBIT_PHASE be quoted relative to
    the view -- and what keeps the orbit from being seen exactly edge-on,
    where it would render as a straight line.
    """
    ref = Vector((1.0, 0.0, 0.0))
    if abs(axis.dot(ref)) > 0.9:
        ref = Vector((0.0, 1.0, 0.0))
    side = axis.cross(ref).normalized()
    right = axis.cross(side).normalized()
    return side, right


def scene_extent():
    """Radius of everything drawn, which is what the camera has to frame."""
    return (ORBIT_RADIUS + max(ELECTRON_RADIUS, HOLE_RADIUS)
            + 0.5 * max(TRAIL_WIDTH, 0.0))


def aura_shape(r_electron, r_hole):
    """
    The aura is a spindle sitting directly between the two particles: a
    prolate spheroid centred on the midpoint, its long axis running from the
    hole to the electron and both ends reaching AURA_MARGIN past them.

    Returns (centre, axis, along, across).
    """
    mid = (r_electron + r_hole) * 0.5
    d = r_electron - r_hole
    axis = d.normalized() if d.length > 1e-9 else Vector((0.0, 0.0, 1.0))
    half = d.length * 0.5
    biggest = max(ELECTRON_RADIUS, HOLE_RADIUS)
    return mid, axis, half + biggest + AURA_MARGIN, biggest + AURA_MARGIN


def orbit_frame(normal):
    """
    (e1, e2) spanning the orbit plane, with e1 to the right of frame and e2
    toward the camera. Tying the frame to the view is what makes ORBIT_PHASE
    mean the same thing whatever ORBIT_NORMAL is set to.
    """
    side, right = camera_basis(normal)
    return right, side


def orbit_position(centre, normal, degrees, radius=None):
    e1, e2 = orbit_frame(normal)
    th = math.radians(degrees)
    r = ORBIT_RADIUS if radius is None else radius
    return centre + (e1 * cos(th) + e2 * sin(th)) * r


def orbit_points(centre, normal, n, radius=None):
    """
    A full circle, starting where the electron is and running backwards along
    its direction of travel. Point i is therefore i/n of the way *behind* the
    electron, which is exactly the parameter the tail needs.
    """
    e1, e2 = orbit_frame(normal)
    phase = math.radians(ORBIT_PHASE)
    sweep = -2.0 * pi * (1.0 if ORBIT_DIRECTION >= 0 else -1.0)
    r = ORBIT_RADIUS if radius is None else radius
    pts = []
    for i in range(n + 1):
        th = phase + sweep * i / n
        pts.append(centre + (e1 * cos(th) + e2 * sin(th)) * r)
    return pts


def build_aura(col, r_electron, r_hole):
    """
    The aura is a single emissive volume sitting directly between the two
    particles: density AURA_DENSITY at the midpoint, falling smoothly to
    nothing at the surface of a spindle that reaches past both of them.

    The mesh is a plain unit sphere and the shape lives entirely in the
    object's scale and rotation. That is not cosmetic -- it is what makes the
    volume shader work, because the shader reads Texture Coordinate -> Object,
    which is normalised by the object transform and so runs 0 at the core to 1
    at the surface however large the aura is.
    """
    e_col = pair_color(ELECTRON_RIM)
    h_col = pair_color(HOLE_RIM)
    mid, axis, along, across = aura_shape(r_electron, r_hole)
    segs, rings = AURA_RESOLUTION

    mat = aura_volume_material("Aura_Volume", h_col, e_col, AURA_DENSITY,
                               AURA_FALLOFF, AURA_EMISSION)
    v, f = sphere_mesh_data(1.0, segs, rings)
    obj = add_object("Exciton_Aura", v, f, mat, col)
    if obj is None:
        return
    obj.location = mid
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = Vector((0.0, 0.0, 1.0)).rotation_difference(axis)
    obj.scale = (across, across, along)
    _enable_volumes()


def dash_spans(length):
    """
    Where the tail is drawn and where it is not, as (start, end) pairs measured
    in fractions of the loop behind the electron.

    The first span is unbroken: the tail leaves the electron as a solid
    stroke. After TRAIL_SOLID it breaks up, and from there the dashes shorten
    while the gaps between them open out, so they arrive less and less often
    on the way back. The tail ends where the dash length falls below
    TRAIL_DASH_MIN -- the point at which it has visibly run out, rather than
    an arbitrary cut.
    """
    length = min(max(length, 1e-3), 1.0)
    solid_end = min(max(TRAIL_SOLID, 0.0), 1.0) * length
    spans = []
    if solid_end > 1e-6:
        spans.append((0.0, solid_end))

    if TRAIL_DASH_LENGTH <= 1e-6:
        return spans              # no dashes wanted: the solid head is the tail
    dash0 = TRAIL_DASH_LENGTH
    gap0 = max(TRAIL_DASH_GAP, 1e-5)
    dashed = max(length - solid_end, 1e-9)
    floor = TRAIL_DASH_MIN * dash0
    s = solid_end
    for _ in range(4000):                      # a bound, not a limit
        if s >= length:
            break
        u = (s - solid_end) / dashed           # 0 where the dashes start, 1 at the end
        dash = dash0 * max(1.0 - u, 0.0) ** max(TRAIL_DASH_FALLOFF, 0.01)
        if dash <= floor:
            break                              # the dashes have shrunk to nothing
        gap = gap0 * (1.0 + (TRAIL_DASH_GROWTH - 1.0) * u)
        spans.append((s, min(s + dash, length)))
        s += dash + gap
    return spans


def arc_points(centre, normal, radius, start, end, n):
    """Points along the orbit between two positions behind the electron."""
    e1, e2 = orbit_frame(normal)
    phase = math.radians(ORBIT_PHASE)
    sweep = -2.0 * pi * (1.0 if ORBIT_DIRECTION >= 0 else -1.0)
    pts = []
    for i in range(n + 1):
        behind = start + (end - start) * (i / n)
        th = phase + sweep * behind
        pts.append(centre + (e1 * cos(th) + e2 * sin(th)) * radius)
    return pts


def build_motion_trail(col, centre, normal):
    """
    The motion tail: a band of concentric strands lying flat in the orbit
    plane, each one solid where it leaves the electron and then breaking into
    dashes that shorten and space out behind it until they run out.

    Three things fade together along a strand, which is what sells it as
    motion rather than as a drawn dotted line: the dashes get shorter and
    further apart, each dash gets thinner, and the colour cools from the
    electron's own highlight to the dim tail colour. Across the band, the
    outer strands trail TRAIL_SHEAR longer than the inner ones so the tail
    feathers instead of ending on a straight edge, and they are dimmer by
    TRAIL_EDGE_FADE so the band keeps a bright spine.
    """
    strands = max(int(TRAIL_STRANDS), 1)
    head = pair_color(TRAIL_HEAD_COLOR or ELECTRON_RIM)
    tail = pair_color(TRAIL_COLOR)
    falloff = max(TRAIL_FALLOFF, 0.01)
    resolution = max(int(TRAIL_RESOLUTION), 24)

    V, F, C = [], [], []
    dashes = 0
    for j in range(strands):
        # -1 at the inner edge of the band, +1 at the outer edge
        u = 0.0 if strands == 1 else (2.0 * j / (strands - 1) - 1.0)
        radius = ORBIT_RADIUS + 0.5 * TRAIL_WIDTH * u
        if radius <= 1e-3:
            continue
        level = max(1.0 - TRAIL_EDGE_FADE * abs(u), 0.0)
        length = min(max(TRAIL_LENGTH * (1.0 + TRAIL_SHEAR * u), 1e-3), 1.0)
        hot = tuple(v * TRAIL_LEVEL * level for v in head)
        cold = tuple(v * TRAIL_END_LEVEL * level for v in tail)

        for start, end in dash_spans(length):
            steps = max(int(resolution * (end - start)), 3)
            pts = arc_points(centre, normal, radius, start, end, steps)
            radii, colors = [], []
            for i in range(steps + 1):
                behind = start + (end - start) * (i / steps)
                w = max(0.0, 1.0 - behind / length) ** falloff
                radii.append(TRAIL_STRAND_RADIUS *
                             (TRAIL_TAIL_WIDTH + (1.0 - TRAIL_TAIL_WIDTH) * w))
                colors.append(mix(cold, hot, w))

            v, f, _ = tube_from_polyline(pts, radii, TRAIL_SEGMENTS,
                                         taper=TRAIL_DASH_TAPER)
            if not v:
                continue
            base = len(V)
            V.extend(v)
            F.extend(tuple(base + i for i in face) for face in f)
            for c in colors:
                C.extend([c] * TRAIL_SEGMENTS)
            dashes += 1

    if not V:
        return 0
    mat = vertex_color_emission_material("MotionTrail", "TrailColor",
                                         TRAIL_EMISSION)
    add_object("Exciton_MotionTrail", V, F, mat, col, colors=C,
               attr_name="TrailColor")
    return dashes


def funnel_profile(t, top, waist, bottom):
    """
    The hourglass, as a function of how far down the funnel you are: `top` at
    the exciton, pinching to `waist`, opening out to `bottom` at the sheet.

    Each half is shaped on its own so the two ends are independent -- the
    mouth at the exciton and the spread where it lands rarely want the same
    curve. The exponents act away from the waist, so above 1 each half hugs
    the pinch and then opens near its own end; that is what makes the
    silhouette an hourglass rather than two straight cones.

    Used for the radius and, with the three brightness levels, for the colour
    down the funnel as well.
    """
    w = min(max(FUNNEL_WAIST, 0.0), 1.0)
    if t <= w:
        u = (w - t) / w if w > 1e-6 else 0.0        # 1 at the top, 0 at the waist
        return waist + (top - waist) * u ** max(FUNNEL_TOP_FLARE, 0.01)
    u = (t - w) / (1.0 - w) if w < 1.0 - 1e-6 else 0.0   # 0 at the waist, 1 at the end
    return waist + (bottom - waist) * u ** max(FUNNEL_BOTTOM_FLARE, 0.01)


def build_light_funnel(col, centre):
    """
    A cone of light dropping from the exciton to the layer below it, as in
    Reference 1.

    It is an hourglass: a mouth at the exciton, a pinch at FUNNEL_WAIST, and a
    wide landing on the sheet, with each end shaped independently. Its alpha
    is written per vertex so it fades along its length, and the material fades
    it again at its own silhouette -- without that second term the cone has a
    visible outline and reads as a solid shell instead of a shaft of light.
    """
    drop = Vector(FUNNEL_DIRECTION)
    if drop.length < 1e-9:
        return None
    drop = drop.normalized()
    ref = Vector((1.0, 0.0, 0.0))
    if abs(drop.dot(ref)) > 0.9:
        ref = Vector((0.0, 1.0, 0.0))
    e1 = drop.cross(ref).normalized()
    e2 = drop.cross(e1).normalized()

    colour = pair_color(FUNNEL_COLOR or HOLE_RIM)
    rings = max(int(FUNNEL_RINGS), 2)
    segs = max(int(FUNNEL_SEGMENTS), 6)

    verts, faces, colors = [], [], []
    for i in range(rings + 1):
        t = i / rings                          # 0 at the exciton, 1 at the sheet
        r = max(funnel_profile(t, FUNNEL_TOP_RADIUS, FUNNEL_WAIST_RADIUS,
                               FUNNEL_BOTTOM_RADIUS), 0.0)
        level = funnel_profile(t, FUNNEL_TOP_LEVEL, FUNNEL_WAIST_LEVEL,
                               FUNNEL_BOTTOM_LEVEL) * FUNNEL_LEVEL
        ring_centre = centre + drop * (FUNNEL_DROP * t)
        # brightness may run past 1 for the glow; the alpha may not
        c = (tuple(v * max(level, 0.0) for v in colour)
             + (min(max(level, 0.0), 1.0),))
        for k in range(segs):
            a = 2.0 * pi * k / segs
            p = ring_centre + (e1 * cos(a) + e2 * sin(a)) * r
            verts.append((p.x, p.y, p.z))
            colors.append(c)
    for i in range(rings):
        for k in range(segs):
            k2 = (k + 1) % segs
            b0, b1 = i * segs, (i + 1) * segs
            faces.append((b0 + k, b0 + k2, b1 + k2, b1 + k))

    mat = beam_material("Funnel", "FunnelColor", 1.0, FUNNEL_ALPHA,
                        FUNNEL_FACING_FALLOFF)
    obj = add_object("Exciton_Funnel", verts, faces, mat, col, colors=colors,
                     attr_name="FunnelColor")

    if SHOW_FUNNEL_POOL:
        build_funnel_pool(col, centre + drop * FUNNEL_DROP, e1, e2, colour)
    return obj


def build_funnel_pool(col, centre, e1, e2, colour):
    """
    The disc of light where the funnel lands: bright at the middle, fading to
    nothing at its rim.

    A flat disc is the one place the silhouette trick cannot help -- every
    point on it faces the camera the same way -- so the falloff is written
    into the vertex alpha instead.
    """
    steps = max(int(FUNNEL_POOL_STEPS), 8)
    rings = 16
    verts, faces, colors = [], [], []

    centre_level = FUNNEL_POOL_LEVEL * FUNNEL_LEVEL
    verts.append((centre.x, centre.y, centre.z))
    colors.append(tuple(v * centre_level for v in colour) + (FUNNEL_POOL_ALPHA,))
    for i in range(1, rings + 1):
        t = i / rings
        level = centre_level * max(1.0 - t, 0.0) ** max(FUNNEL_POOL_FALLOFF, 0.01)
        alpha = FUNNEL_POOL_ALPHA * max(1.0 - t, 0.0) ** max(FUNNEL_POOL_FALLOFF, 0.01)
        c = tuple(v * level for v in colour) + (alpha,)
        for k in range(steps):
            a = 2.0 * pi * k / steps
            p = centre + (e1 * cos(a) + e2 * sin(a)) * (FUNNEL_POOL_RADIUS * t)
            verts.append((p.x, p.y, p.z))
            colors.append(c)

    def idx(ring, k):
        return 1 + (ring - 1) * steps + (k % steps)

    faces.extend((0, idx(1, k), idx(1, k + 1)) for k in range(steps))
    for i in range(1, rings):
        for k in range(steps):
            faces.append((idx(i, k), idx(i + 1, k),
                          idx(i + 1, k + 1), idx(i, k + 1)))

    mat = beam_material("FunnelPool", "PoolColor", 1.0, 1.0, 0.0)
    return add_object("Exciton_FunnelPool", verts, faces, mat, col,
                      colors=colors, attr_name="PoolColor")


def _enable_volumes():
    """Make sure the render engine will actually show the volumetric aura."""
    scene = bpy.context.scene
    try:
        ev = scene.eevee
        for attr, val in (("use_volumetric_lights", True),
                          ("volumetric_tile_size", "2"),
                          ("volumetric_samples", 128),
                          ("volumetric_start", 0.05),
                          ("volumetric_end", 1000.0)):
            if hasattr(ev, attr):
                try:
                    setattr(ev, attr, val)
                except Exception:
                    pass
    except Exception:
        pass


def _set_render_engine(scene, name):
    """Set the engine by its identifier, which is renamed between releases."""
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


def setup_render(scene):
    """Resolution, samples, and the transparent RGBA film."""
    engine = _set_render_engine(scene, RENDER_ENGINE)
    r = scene.render
    r.resolution_x, r.resolution_y = RESOLUTION
    r.resolution_percentage = 100
    # The backdrop is an object in the scene, so the film can stay
    # transparent in every mode: where there is no backdrop the render simply
    # comes out with an alpha channel.
    r.film_transparent = True

    # PNG + RGBA: without the RGBA colour mode the alpha channel is thrown
    # away on save and the background comes out black
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
            scene.cycles.film_transparent_glass = True
        except Exception:
            pass
    else:
        for attr in ("taa_render_samples", "taa_samples"):
            if hasattr(scene.eevee, attr):
                try:
                    setattr(scene.eevee, attr, RENDER_SAMPLES)
                except Exception:
                    pass
    _enable_volumes()

    if VIEW_TRANSFORM:
        try:
            scene.view_settings.view_transform = VIEW_TRANSFORM
        except Exception as e:
            print(f"[exciton] view transform {VIEW_TRANSFORM!r} not available: {e}")

    if OUTPUT_PATH:
        r.filepath = os.path.abspath(os.path.expanduser(OUTPUT_PATH))
    return engine


def setup_world(scene):
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    if bg:
        # the backdrop object carries the background, so the world only ever
        # has to stay out of the way -- at any strength it would tint the aura
        bg.inputs["Color"].default_value = (0.0, 0.0, 0.0, 1.0)
        bg.inputs["Strength"].default_value = 0.0


def drawn_points():
    """
    Every vertex of every mesh built so far, in world space.

    The camera frames these rather than an estimated radius. An estimate has
    to guess how far the drawn geometry actually reaches -- the width of the
    tail band, the aura's spindle, the halo spheres -- and anything it misses
    is something the frame clips off. The vertices cannot be wrong about it.

    Meshes here are built in world coordinates, but the aura and the haloes
    carry their size in the object transform, so those are expanded by hand.
    """
    pts = []
    for obj in bpy.data.objects:
        if obj.type != "MESH" or obj.data is None:
            continue
        loc = Vector(getattr(obj, "location", (0.0, 0.0, 0.0)))
        sx, sy, sz = tuple(getattr(obj, "scale", (1.0, 1.0, 1.0)))
        scaled = abs(sx - 1.0) + abs(sy - 1.0) + abs(sz - 1.0) > 1e-9
        radius = max(abs(sx), abs(sy), abs(sz))
        for v in obj.data.vertices:
            co = getattr(v, "co", v)
            if scaled:
                # a rotated, non-uniformly scaled object: its own rotation is
                # unknown here, so take the enclosing sphere of the scale
                pts.append((loc.x + co[0] * radius, loc.y + co[1] * radius,
                            loc.z + co[2] * radius))
            else:
                pts.append((loc.x + co[0], loc.y + co[1], loc.z + co[2]))
    return pts


def setup_camera(scene, centre, axis, points=None):
    """
    Look at the exciton from side-on, raised CAMERA_ELEVATION degrees above
    the orbit plane. The elevation is what opens the orbit out into an
    ellipse: from dead level the electron's path is seen edge-on and renders
    as a straight line.

    The distance is fitted rather than guessed: every point is projected onto
    the camera's right and up axes, and since those components do not change
    as the camera slides along its view direction, each point sets a lower
    bound on the distance and the largest of them frames the lot.
    """
    side, _ = camera_basis(axis)
    elev = math.radians(CAMERA_ELEVATION)
    eye = (side * cos(elev) + axis * sin(elev)).normalized()
    forward = -eye
    right = forward.cross(Vector((0.0, 0.0, 1.0)))
    right = right.normalized() if right.length > 1e-6 else Vector((1.0, 0.0, 0.0))
    up = right.cross(forward).normalized()

    if not points:
        r = scene_extent()
        points = [(centre.x + x * r, centre.y + y * r, centre.z + z * r)
                  for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    target = Vector((0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)),
                     0.5 * (min(zs) + max(zs))))

    res_x, res_y = RESOLUTION
    aspect = res_x / float(res_y)
    tan_h = 18.0 / max(CAMERA_LENS, 1e-3)          # 36 mm sensor, half-angle
    tan_v = tan_h / aspect

    dist = 0.0
    reach = 0.0
    for p in points:
        d = Vector(p) - target
        ahead = d.dot(forward)
        dist = max(dist,
                   abs(d.dot(right)) * CAMERA_MARGIN / tan_h - ahead,
                   abs(d.dot(up)) * CAMERA_MARGIN / tan_v - ahead)
        reach = max(reach, d.length)
    dist = max(dist, 1.0)

    cam_data = bpy.data.cameras.new("ExcitonCamera")
    cam_data.lens = CAMERA_LENS
    try:
        cam_data.sensor_fit = "HORIZONTAL"         # makes the maths definite
    except Exception:
        pass
    cam_data.clip_start = max(0.01, (dist - reach) * 0.5)
    cam_data.clip_end = dist + reach * 4.0 + 100.0

    cam = bpy.data.objects.new("ExcitonCamera", cam_data)
    scene.collection.objects.link(cam)
    cam.location = target + eye * dist
    cam.rotation_euler = forward.to_track_quat("-Z", "Y").to_euler()
    scene.camera = cam

    fw = max(abs((Vector(p) - target).dot(right)) /
             max((Vector(p) - target).dot(forward) + dist, 1e-9)
             for p in points) / tan_h
    fh = max(abs((Vector(p) - target).dot(up)) /
             max((Vector(p) - target).dot(forward) + dist, 1e-9)
             for p in points) / tan_v
    print(f"[exciton] framing: fills {fw * 100:.0f}% of the frame across and "
          f"{fh * 100:.0f}% up, from {len(points)} vertices")
    return cam
    return cam


def _ease_clamped(u):
    u = min(max(u, 0.0), 1.0)
    return u * u * (3.0 - 2.0 * u)


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
    behave identically in EEVEE and Cycles.
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
    tan_h = 18.0 / max(getattr(cam.data, "lens", 50.0), 1e-3)
    half_w = (centre - loc).length * tan_h * BACKGROUND_MARGIN
    half_h = half_w / aspect

    rows = max(int(BACKGROUND_STEPS), 2)
    verts, faces, colors = [], [], []
    for i in range(rows + 1):
        t = i / rows
        c = mix(bottom, top, _ease_clamped(t) if BACKGROUND_SMOOTH else t)
        y = (2.0 * t - 1.0) * half_h
        for side in (-1.0, 1.0):
            p = centre + right * (side * half_w) + up * y
            verts.append((p.x, p.y, p.z))
            colors.append(c)
    for i in range(rows):
        b = 2 * i
        faces.append((b, b + 1, b + 3, b + 2))

    mat = vertex_color_emission_material("Background", "BackgroundColor", 1.0)
    obj = add_object("Exciton_Background", verts, faces, mat,
                     get_collection("Exciton"), colors=colors,
                     attr_name="BackgroundColor")
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


def setup_compositor(scene):
    """
    Glare for the bloom, and on a transparent film a second pass that folds
    the bloom into the alpha channel.

    Glare only adds colour; it does not touch alpha. On a transparent film
    every pixel of glow outside the geometry therefore keeps alpha = 0 and
    vanishes the moment the PNG is composited onto a slide. Taking the
    luminance of the glare result and maxing it into the existing alpha is
    what keeps the halo visible there.
    """
    scene.use_nodes = True
    nt = scene.node_tree
    nt.nodes.clear()

    rl = nt.nodes.new("CompositorNodeRLayers")
    rl.location = (-500, 0)
    comp = nt.nodes.new("CompositorNodeComposite")
    comp.location = (500, 0)

    image_out = rl.outputs["Image"]
    if ADD_GLOW_COMPOSITOR:
        glare = nt.nodes.new("CompositorNodeGlare")
        glare.location = (-250, 0)
        for gt in ("BLOOM", "FOG_GLOW"):
            try:
                glare.glare_type = gt
                break
            except Exception:
                continue
        for attr, val in (("quality", "HIGH"), ("mix", 0.1),
                          ("threshold", 0.6), ("size", 8)):
            if hasattr(glare, attr):
                try:
                    setattr(glare, attr, val)
                except Exception:
                    pass
        nt.links.new(rl.outputs["Image"], glare.inputs["Image"])
        image_out = glare.outputs["Image"]

    if scene.render.film_transparent and GLOW_IN_ALPHA:
        try:
            bw = nt.nodes.new("CompositorNodeRGBToBW")
            bw.location = (-60, -220)
            nt.links.new(image_out, bw.inputs["Image"])

            gain = nt.nodes.new("CompositorNodeMath")
            gain.location = (110, -220)
            gain.operation = "MULTIPLY"
            gain.inputs[1].default_value = GLOW_ALPHA_GAIN
            nt.links.new(bw.outputs["Val"], gain.inputs[0])

            biggest = nt.nodes.new("CompositorNodeMath")
            biggest.location = (280, -220)
            biggest.operation = "MAXIMUM"
            biggest.use_clamp = True
            nt.links.new(rl.outputs["Alpha"], biggest.inputs[0])
            nt.links.new(gain.outputs["Value"], biggest.inputs[1])

            setalpha = nt.nodes.new("CompositorNodeSetAlpha")
            setalpha.location = (330, 0)
            try:
                setalpha.mode = "REPLACE_ALPHA"
            except Exception:
                pass
            nt.links.new(image_out, setalpha.inputs["Image"])
            nt.links.new(biggest.outputs["Value"], setalpha.inputs["Alpha"])
            image_out = setalpha.outputs["Image"]
        except Exception as e:
            print(f"[exciton] glow-in-alpha skipped: {e}")

    nt.links.new(image_out, comp.inputs["Image"])


def setup_scene(centre, axis):
    if BACKGROUND_MODE not in ("gradient", "flat", "transparent"):
        raise ValueError("BACKGROUND_MODE must be 'gradient', 'flat' or "
                         f"'transparent'; got {BACKGROUND_MODE!r}")
    scene = bpy.context.scene
    engine = setup_render(scene)
    setup_world(scene)
    points = drawn_points()
    cam = setup_camera(scene, centre, axis, points)

    if BACKGROUND_MODE in ("gradient", "flat"):
        bottom = hex_to_linear(BACKGROUND_COLOR)
        top = (hex_to_linear(BACKGROUND_TOP_COLOR)
               if BACKGROUND_MODE == "gradient" else bottom)
        reach = max((Vector(p) - centre).length for p in points) if points \
            else scene_extent()
        add_background_plane(scene, cam, centre, Vector((0.0, 0.0, 1.0)),
                             bottom, top, reach)

    if ADD_GLOW_COMPOSITOR:
        try:
            setup_compositor(scene)
        except Exception as e:
            print(f"[exciton] compositor skipped: {e}")
    else:
        scene.use_nodes = False
    return engine


# ============================================================================
#  Main
# ============================================================================

def main():
    if CLEAR_SCENE:
        clear_scene()

    t = apply_binding_strength()
    print(f"[exciton] binding strength {t:.2f} "
          f"({'strongly' if t >= 0.66 else 'weakly' if t <= 0.33 else 'loosely'}"
          f" bound): orbit {ORBIT_RADIUS:.2f}, aura density "
          f"{AURA_DENSITY:.2f}, tail over {TRAIL_LENGTH * 100:.0f}% of the loop")

    normal = Vector(ORBIT_NORMAL)
    if normal.length < 1e-9:
        raise ValueError("ORBIT_NORMAL must be a non-zero vector")
    normal = normal.normalized()
    centre = Vector(CENTER)
    r_hole = centre                                   # the hole is orbited
    r_electron = orbit_position(centre, normal, ORBIT_PHASE)

    if ORBIT_RADIUS <= ELECTRON_RADIUS + HOLE_RADIUS:
        print("[exciton] WARNING: ORBIT_RADIUS is smaller than the two radii "
              "combined -- the electron overlaps the hole.")

    col = get_collection("Exciton")

    if SHOW_FIELD_LINES:
        lines, failed = build_field_lines(r_electron, r_hole)
        print(f"[exciton] traced {len(lines)} field lines "
              f"({failed} discarded), "
              f"{sum(len(l) for l in lines)} integration points total")
        if lines:
            build_field_line_objects(col, lines)

    if SHOW_TRAIL:
        dashes = build_motion_trail(col, centre, normal)
        print(f"[exciton] electron at {ORBIT_PHASE:.0f} deg; tail of "
              f"{TRAIL_STRANDS} strands over {TRAIL_WIDTH:.2f} units, "
              f"{dashes} dashes, solid for the first "
              f"{TRAIL_SOLID * 100:.0f}% of it")

    if SHOW_FUNNEL:
        build_light_funnel(col, centre)
        print(f"[exciton] light funnel dropping {FUNNEL_DROP:.1f} units: "
              f"{FUNNEL_TOP_RADIUS:.1f} at the exciton, pinched to "
              f"{FUNNEL_WAIST_RADIUS:.1f} at {FUNNEL_WAIST * 100:.0f}% down, "
              f"{FUNNEL_BOTTOM_RADIUS:.1f} where it lands"
              + (f", pool of {FUNNEL_POOL_RADIUS:.1f}" if SHOW_FUNNEL_POOL
                 else ""))

    build_particles(col, r_electron, r_hole)

    if SHOW_PARTICLE_GLOW:
        build_particle_glow(col, r_electron, r_hole)

    if SHOW_AURA:
        build_aura(col, r_electron, r_hole)

    engine = None
    if SETUP_SCENE:
        engine = setup_scene(centre, normal)

    print("[exciton] Done. Electron at "
          f"({r_electron.x:.2f}, {r_electron.y:.2f}, {r_electron.z:.2f}), "
          f"hole at ({r_hole.x:.2f}, {r_hole.y:.2f}, {r_hole.z:.2f}).")
    if SETUP_SCENE:
        described = {"transparent": "transparent (RGBA PNG)",
                     "flat": f"flat {BACKGROUND_COLOR}",
                     "gradient": f"gradient {BACKGROUND_COLOR} -> "
                                 f"{BACKGROUND_TOP_COLOR}"}[BACKGROUND_MODE]
        print(f"[exciton] engine {engine}, background {described}")
    print("[exciton] The aura is a volume: use Rendered shading (not Material "
          "Preview) to see it.")

    if RENDER_NOW:
        if not OUTPUT_PATH:
            print("[exciton] RENDER_NOW is on but OUTPUT_PATH is empty; "
                  "nothing was written.")
        else:
            path = os.path.abspath(os.path.expanduser(OUTPUT_PATH))
            folder = os.path.dirname(path)
            if folder and not os.path.isdir(folder):
                os.makedirs(folder, exist_ok=True)
            bpy.context.scene.render.filepath = path
            print(f"[exciton] rendering to {path} ...")
            bpy.ops.render.render(write_still=True)
            print("[exciton] render written.")


main()
