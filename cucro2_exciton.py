"""
cucro2_exciton.py
=================
Blender script: excitons over a percolating CuCrO2 sheet, in one scene.

The lattice is a Cu/Pd plane under a CrO2 slab, with the A-site substitution
laid out as large metallic islands on an otherwise host sheet (or as a
gradient, or uniformly). Over it float excitons -- an electron orbiting a
hole, with a dashed motion tail and a funnel of light dropping to the plane
below. Strongly bound ones gather on the metallic islands; weakly bound ones
are scattered out over the host sheet, smaller and fainter.

Run: Blender -> Scripting -> Open -> edit USER PARAMETERS -> Alt+P.

Why one file
------------
The two halves were separate scripts, overlaid afterwards. They are together
because they have to be:

  * Placement needs the lattice. "Strong excitons on the islands" can only be
    answered by whatever decided where the islands are.
  * One camera, or the perspective does not match. An exciton over the far
    edge of the sheet has to be smaller and higher in frame than one over the
    near edge, by exactly what the lens and camera position dictate. Two
    renders overlaid by hand cannot track that.
  * Occlusion. An exciton behind the CrO2 slab, or passing behind an atom,
    needs real depth sorting.

DRAW still lets the two come out as separate images when that is wanted --
but framed by the same camera, so they line up by construction rather than by
eye.

1 Blender unit = 1 Angstrom, except for the excitons: they are drawn far
larger than life, and EXCITON_SCALE is where that is admitted.
"""

import bpy
import colorsys
import math
import os
import random
import re
from math import atan2, cos, hypot, pi, radians, sin, sqrt, floor
from mathutils import Matrix, Vector

# ============================================================================
#  USER PARAMETERS
# ============================================================================

# --- what to build ----------------------------------------------------------
DRAW = "both"                # "both"     -> the whole figure
                             # "lattice"  -> the sheet on its own
                             # "excitons" -> the excitons on their own, on a
                             #    transparent film, for layering over the
                             #    lattice pass later
                             # Either single pass is framed by the camera that
                             # frames the WHOLE scene, so the two line up
                             # exactly when composited.


# ============================================================================
#  THE LATTICE
# ============================================================================

CIF_PATH = r"C:\path\to\CuCrO2.cif"     # <-- EDIT (falls back to built-in data)

# --- supercell / orientation -------------------------------------------------
N_CELLS        = (36, 36, 1) # unit cells along a, b, c. Keep nz small: c = 17 A.
                             # A composition gradient needs a wide sheet to
                             # read -- with only a handful of sites per slice
                             # the statistics are noise, not a ramp.
VIEW_DIRECTION = "001"       # "001" looks down c (the figure's view); also
                             # accepts "100", "110", "111", etc.
LATTICE_SCALE  = 1.2         # multiplies the CIF lattice parameters

# --- percolation -------------------------------------------------------------
A_SITE_ELEMENT   = "Cu"      # element on the A site that gets substituted
METAL_FRACTION   = 0.50      # x: probability an A site is metallic (Pd-like).
                             # Used only by COMPOSITION_MODE = "uniform".
RANDOM_SEED      = 6         # change for a different random configuration

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
ISLAND_RADIUS         = 0.18  # mean island radius, as a fraction of the
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
SHOW_BONDS        = False     # chemical Cu-O / Cr-O bonds. With the vertical
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
SHOW_LATTICE_GUIDE = True   # faint network over ALL A sites (figure's dashes)
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
ATOM_SEGMENTS    = 24        # a 16x16 sheet is ~1800 atoms, so the mesh cost
ATOM_RINGS       = 16        # per sphere matters; 24x16 is still smooth at
                             # the size these render
CYL_SEGMENTS     = 16
CLEAR_SCENE      = True


# ============================================================================
#  EXCITONS: how many, how tightly bound, and where they sit
# ============================================================================

STRONG_EXCITONS      = 7      # how many strongly bound pairs to place
WEAK_EXCITONS        = 5      # how many weakly bound ones
STRONG_BINDING       = 1.00   # binding strength of the strong ones, 0..1
WEAK_BINDING         = 0.15   # ... and of the weak ones

# The strength dial drives the whole model at once -- orbit radius, particle
# size, field lines, tail, funnel, brightness and colour -- at these sizes
# one quantity changing by 30% does not read, but everything loosening and
# fading together does. Each WEAK_ value below is the far end of the quantity
# named after it; the value in the exciton model block is the strong end.
WEAK_PARTICLE_SCALE  = 0.55   # a weakly bound pair's electron and hole are
                              # this fraction of the size, keeping the ratio
                              # between the two
WEAK_ORBIT_RADIUS    = 6.4    # the pair drifts apart as the binding weakens
WEAK_FIELD_LINE_RINGS = 3     # the loose end of the field line bundle: this
WEAK_FIELD_LINE_AZIMUTHS = 8  # many rings of this many lines each, against
                              # FIELD_LINE_RINGS x FIELD_LINE_AZIMUTHS in the
                              # model block. Both counts are rounded after the
                              # blend, so at WEAK_BINDING = 0.15 a weak pair
                              # comes out at 4 x 10 = 40 lines, not 3 x 8
WEAK_TRAIL_LENGTH    = 0.90   # a long wispy tail instead of a short hot one
WEAK_TRAIL_LEVEL     = 1.0
WEAK_GLOW_OUTER_SCALE = 2.3   # a bigger, fainter halo on each particle
WEAK_GLOW_ALPHA      = 0.26
WEAK_FUNNEL_LEVEL    = 0.35   # a dimmer funnel down to the sheet
WEAK_EMISSION_SCALE  = 0.65   # everything dims
WEAK_DESATURATION    = 0.55   # ... and washes toward grey, so weakly bound
                              # pairs recede into the background

# --- where they go ----------------------------------------------------------
EXCITON_SCALE        = 3.0    # size of the whole exciton model relative to
                              # the numbers in the model block below. The
                              # lattice is drawn to scale and an exciton is
                              # not: at 1.0 a pair is a few Angstrom across
                              # and vanishes against the sheet.
EXCITON_HEIGHT       = 0.62   # how high they float above the Cu plane, as a
                              # fraction of the gap up to the CrO2 slab. The
                              # funnel drops from there to the plane, so this
                              # also sets how long the funnel is.
EXCITON_HEIGHT_JITTER = 0.12  # +/- this much of the gap, so they do not sit
                              # in a flat row
EXCITON_TILT         = 14.0   # degrees of random tilt on each orbit plane.
                              # Without it every ellipse is identical and the
                              # field reads as wallpaper.
EXCITON_SPACING      = 1.5    # keep pairs this many orbit diameters apart
EXCITON_SEED         = 11     # change for a different arrangement

STRONG_SPREAD        = 0.85   # how far out into its island a strong exciton
                              # may sit, as a fraction of the island radius.
                              # Small values cluster them on the middle of
                              # each island, near 1 pushes them to the coast.
WEAK_CLEARANCE       = 1.25   # weak excitons keep this many island radii
                              # away from every island, so the two
                              # populations do not mix


# ============================================================================
#  THE EXCITON MODEL
#  These describe a STRONGLY bound pair; see the WEAK_ block above.
# ============================================================================

# --- the pair: an electron in orbit around a hole ---------------------------
ORBIT_RADIUS     = 3.6        # how far the electron orbits from the hole,
                              # before EXCITON_SCALE
# The orbit plane is parallel to the sheet, tilted per exciton by
# EXCITON_TILT, and the phase is dealt at random -- so ORBIT_NORMAL and
# ORBIT_PHASE are gone; the placement decides both.

ELECTRON_RADIUS  = 1.10       # the electron is the smaller one, out on the
HOLE_RADIUS      = 1.25       # orbit; the hole is the larger one at the
                              # centre, which is what it is being orbited.
                              # Both shrink with the binding strength, by
                              # WEAK_PARTICLE_SCALE.

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
SHOW_FIELD_LINES     = True   # the dipole streamline bundle binding the
                              # pair. How many lines each exciton gets follows
                              # its binding strength: the counts here are the
                              # strong end and WEAK_FIELD_LINE_* the other, so
                              # a loosely bound pair is held by visibly fewer.
FIELD_LINE_RINGS     = 7      # launch angles (rings of lines), for a pair at
FIELD_LINE_AZIMUTHS  = 18     # full binding: 7 x 18 = 126 lines, against 4 x
                              # 10 = 40 for the weak population at
                              # WEAK_BINDING 0.15
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

INTEGRATION_STEP     = 0.045  # RK4 step length, before EXCITON_SCALE. It
                              # scales with the exciton, so a larger pair
                              # costs no more steps per line than a small one
                              # -- otherwise tracing slows by the cube of the
                              # scale for no extra detail.
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

# --- mesh resolution ---------------------------------------------------------
PARTICLE_SEGMENTS    = 48     # mesh resolution of the electron and the hole
PARTICLE_RINGS       = 32

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
CAMERA_ELEVATION       = 0.0   # degrees above the plane of the sheet. 90 is
                                # straight down (a plan view, where the
                                # sandwich is edge-on and invisible); lower
                                # angles show its thickness.
CAMERA_AZIMUTH         = 60.0  # degrees around the sheet, measured from +x.
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


ADD_GLOW_COMPOSITOR    = True   # Glare on the whole frame, which is what makes
                                # the excitons bloom
GLOW_IN_ALPHA          = True   # push the bloom into the alpha channel too,
                                # so it survives being composited; without it
                                # every pixel of glow outside the geometry
                                # keeps alpha 0 and is clipped away
GLOW_ALPHA_GAIN        = 1.6


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


# ============================================================================
#  Shared helpers
# ============================================================================

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


def mix(c1, c2, t):
    return tuple(a + (b - a) * t for a, b in zip(c1, c2))


def smoothstep(u):
    u = min(max(u, 0.0), 1.0)
    return u * u * (3.0 - 2.0 * u)


def clear_scene():
    for obj in list(bpy.data.objects):
        if obj.type in {"MESH", "CAMERA", "LIGHT"}:
            bpy.data.objects.remove(obj, do_unlink=True)
    for block in (bpy.data.meshes, bpy.data.materials):
        for b in list(block):
            if b.users == 0:
                block.remove(b)


def get_collection(name, parent=None):
    col = bpy.data.collections.get(name)
    if col is None:
        col = bpy.data.collections.new(name)
        (parent or bpy.context.scene.collection).children.link(col)
    return col


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


def add_object(name, verts, faces, mat, collection, colors=None,
               attr_name="LineColor", smooth=True):
    if not verts or not faces:
        return None
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.validate()
    for poly in mesh.polygons:
        poly.use_smooth = smooth
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


def blob_spheres(positions, radius):
    tv, tf = sphere_mesh_data(radius, ATOM_SEGMENTS, ATOM_RINGS)
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
#  Space-group symmetry: generators -> full group by closure
# ============================================================================

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

# ============================================================================
#  CIF parsing (multi-block aware, '?' tolerant)
# ============================================================================

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
#  Lattice geometry
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
#  How the substitution is laid out
# ============================================================================

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
        return GRADIENT_MIN + (GRADIENT_MID - GRADIENT_MIN) * smoothstep(u)
    u = (t - mid_hi) / max((1.0 - right) - mid_hi, 1e-9)
    return GRADIENT_MID + (GRADIENT_MAX - GRADIENT_MID) * smoothstep(u)


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
            best = max(best, lo_p + (hi_p - lo_p) * smoothstep(u))
        targets.append(best)
    return targets, islands


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
#  Building the sheet
# ============================================================================

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
        add_object(f"{tag}_A_host", V, F, mat_host, col)

        if COLOR_MODE == "clusters":
            for n, (cid, ids) in enumerate(sorted(clusters.items(),
                                                  key=lambda t: -len(t[1]))):
                hue = (n * 0.6180339887) % 1.0
                rgb = colorsys.hsv_to_rgb(hue, 0.62, 0.95)
                m = make_material(f"{tag}_cluster{n:03d}",
                                  tuple(v ** 2.2 for v in rgb), 1 - a_atoms)
                V, F = blob_spheres([a_pos[i] + off for i in ids], r_m)
                add_object(f"{tag}_cluster{n:03d}_{len(ids)}sites", V, F, m, col)
        elif COLOR_MODE == "spanning":
            span_ids = [i for i in range(len(a_idx))
                        if labels[i] in spanning and is_metal[i]]
            other_ids = [i for i in range(len(a_idx))
                         if is_metal[i] and labels[i] not in spanning]
            V, F = blob_spheres([a_pos[i] + off for i in span_ids], r_m)
            add_object(f"{tag}_A_metal_spanning", V, F, mat_metal, col)
            V, F = blob_spheres([a_pos[i] + off for i in other_ids], r_m)
            add_object(f"{tag}_A_metal_isolated", V, F, mat_dim, col)
        else:                                    # "species"
            V, F = blob_spheres([a_pos[k] + off for k in range(len(a_idx))
                                 if is_metal[k]], r_m)
            add_object(f"{tag}_A_metal", V, F, mat_metal, col)

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
            add_object(f"{tag}_{e}", V, F, mat, col)

    # ---------------- chemical bonds ----------------
    if SHOW_BONDS:
        segs = drop_vertical([(atoms[i][1] + off, atoms[j][1] + off)
                              for i, j in bond_pairs], stack_axis)
        V, F = blob_cylinders(segs, BOND_RADIUS)
        add_object(f"{tag}_bonds", V, F, mat_bond, col)

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
            add_object(f"{tag}_channels_spanning", V, F, mat_chan, col)
            V, F = blob_cylinders(rest, CHANNEL_RADIUS)
            add_object(f"{tag}_channels_isolated", V, F, mat_dim, col)
        else:
            segs = drop_vertical([(a_pos[i] + off, a_pos[j] + off)
                                  for i, j in channel_pairs], stack_axis)
            V, F = blob_cylinders(segs, CHANNEL_RADIUS)
            add_object(f"{tag}_channels", V, F, mat_chan, col)

    # ---------------- faint guide over the whole A lattice ----------------
    if SHOW_LATTICE_GUIDE:
        pairs, _ = a_site_neighbors(a_pos)
        segs = drop_vertical([(a_pos[i] + off, a_pos[j] + off)
                              for i, j in pairs], stack_axis)
        V, F = blob_cylinders(segs, GUIDE_RADIUS)
        add_object(f"{tag}_A_lattice_guide", V, F, mat_guide, col)

    width = (max(p.x for p in a_pos) - min(p.x for p in a_pos)) if a_pos else 0.0

    # bounds of everything drawn in this panel, padded by the largest display
    # radius, so the camera can frame the whole series
    pad = max(list(RADII.values()) + [DEFAULT_RADIUS])
    pts = [p + off for _, p in atoms]
    bounds = (min(p.x for p in pts) - pad, max(p.x for p in pts) + pad,
              min(p.y for p in pts) - pad, max(p.y for p in pts) + pad,
              min(p.z for p in pts) - pad, max(p.z for p in pts) + pad)

    # What the excitons need to know about the sheet they sit over: where the
    # islands are, how high the A-plane is, and how much room there is between
    # it and the slab above. The islands come out in panel-local coordinates,
    # so they take the panel offset here.
    islands = [(cx + off.x, cy + off.y, r, h)
               for cx, cy, r, h in info.get("islands", [])]
    plane_z = sum(p.z for p in a_pos) / len(a_pos) if a_pos else 0.0
    heights = [h for h, _ in layer_groups(atoms, stack_axis)]
    layer_gap = (max(heights) - min(heights)) if len(heights) > 1 else 6.0
    return {
        "width": width,
        "bounds": bounds,
        "bounds_xy": (bounds[0], bounds[1], bounds[2], bounds[3]),
        "islands": islands,
        "plane_z": plane_z,
        "layer_gap": layer_gap,
        "a_positions": [p + off for p in a_pos],
        "is_metal": is_metal,
    }


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
#  Gradient shaders
#  Two small ideas do most of the work here:
#  
#    * |dot(N, I)| is the cosine between the surface normal and the view
#      ray: 1 where a surface faces the camera head-on and 0 exactly at
#      its silhouette. Driving alpha or brightness with it gives a smooth
#      radial gradient with no visible edge anywhere.
#  
#    * dot(N, L) against a fixed direction is a lambert term. Wrapped to
#      0..1 it shades a sphere from a lit side to a shaded side with no
#      lamp in the scene, so an emissive particle stops looking like a
#      flat coloured circle.
# ============================================================================

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


# ============================================================================
#  The dipole field, for the optional field lines
#  Streamlines of the field of two equal and opposite point charges,
#  integrated with RK4 through
#  
#        E(r) = (r - r_h)/|r - r_h|^3  -  (r - r_e)/|r - r_e|^3
#  
#  so every line starts on the hole, ends on the electron, and bows out
#  by exactly what the dipole geometry dictates.
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


def build_field_lines(ex):
    """Trace every field line from this exciton's hole (+) to its electron (-)."""
    r_electron, r_hole = ex.electron(), ex.centre
    axis_to_electron = (r_electron - r_hole).normalized()
    dirs = launch_directions(axis_to_electron, ex.field_rings,
                             ex.field_azimuths, LAUNCH_ANGLE_MIN,
                             LAUNCH_ANGLE_MAX, FIELD_LINE_AXIS_BIAS)
    start_r = ex.hole_radius * 1.02
    stop_r = ex.electron_radius * 1.02

    # the step is a length, so it scales with the pair: without this a pair
    # drawn three times larger costs three times the steps per line
    step = INTEGRATION_STEP * ex.scale

    lines, failed = [], 0
    for d, _ring_u in dirs:
        pts = trace_field_line(r_hole + d * start_r, r_hole, r_electron,
                               stop_r, step=step)
        if pts and len(pts) > 3:
            lines.append(pts)
        else:
            failed += 1
    return lines, failed


def build_field_line_objects(col, ex, lines):
    """One merged mesh for all of this exciton's field lines, coloured per vertex."""
    e_col = ex.color(ELECTRON_COLOR)
    h_col = ex.color(FIELD_LINE_COLOR)
    V, F, C = [], [], []
    for pts in lines:
        v, f, params = tube_from_polyline(pts, FIELD_LINE_RADIUS * ex.scale,
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
    mat = vertex_color_emission_material(f"FieldLines_{ex.tag}", "LineColor",
                                         FIELD_LINE_EMISSION * ex.emission)
    add_object("Exciton_FieldLines", V, F, mat, col, colors=C)

# ============================================================================
#  The motion tail and the light funnel, in profile
# ============================================================================

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
    return waist + (bottom - waist) * u ** max(FUNNEL_BOTTOM_FLARE, 0.01)# ============================================================================
#  One exciton
# ============================================================================

class Exciton:
    """
    One exciton, with every quantity that depends on its binding strength
    already worked out.

    The model parameters describe a strongly bound pair and each WEAK_ value
    is the far end of the quantity named after it, so `strength` is a straight
    blend between the two. Resolving it here rather than by rewriting module
    globals is what lets strong and weak pairs stand in the same scene.

    `scale` multiplies every length. The lattice is drawn to scale and the
    excitons are not -- at 1.0 a pair would be a few Angstrom across and
    invisible against the sheet.
    """

    def __init__(self, centre, strength, scale, normal, phase, direction):
        t = min(max(strength, 0.0), 1.0)

        def blend(weak, strong):
            return weak + (strong - weak) * t

        self.centre = centre
        self.strength = t
        self.normal = normal.normalized()
        self.phase = phase
        self.direction = 1.0 if direction >= 0 else -1.0
        self.scale = scale
        # materials are shared by every pair of the same strength, so two
        # populations cost two sets rather than one set each
        self.tag = f"{int(round(t * 100)):03d}"

        particle = blend(WEAK_PARTICLE_SCALE, 1.0) * scale
        self.electron_radius = ELECTRON_RADIUS * particle
        self.hole_radius = HOLE_RADIUS * particle
        self.orbit_radius = blend(WEAK_ORBIT_RADIUS, ORBIT_RADIUS) * scale

        self.emission = blend(WEAK_EMISSION_SCALE, 1.0)
        self.fade = (1.0 - t) * WEAK_DESATURATION

        # how much of the dipole bundle actually gets drawn. Rounded, because
        # rings and azimuths are counts, and floored at 1 so a pair is never
        # left with no lines at all
        self.field_rings = max(int(round(
            blend(WEAK_FIELD_LINE_RINGS, FIELD_LINE_RINGS))), 1)
        self.field_azimuths = max(int(round(
            blend(WEAK_FIELD_LINE_AZIMUTHS, FIELD_LINE_AZIMUTHS))), 1)

        self.glow_scale = blend(WEAK_GLOW_OUTER_SCALE, GLOW_OUTER_SCALE)
        self.glow_alpha = blend(WEAK_GLOW_ALPHA, GLOW_ALPHA)

        self.trail_length = blend(WEAK_TRAIL_LENGTH, TRAIL_LENGTH)
        self.trail_level = blend(WEAK_TRAIL_LEVEL, TRAIL_LEVEL)
        self.trail_width = TRAIL_WIDTH * scale
        self.strand_radius = TRAIL_STRAND_RADIUS * scale

        self.funnel_level = blend(WEAK_FUNNEL_LEVEL, FUNNEL_LEVEL)
        self.funnel_top = FUNNEL_TOP_RADIUS * scale
        self.funnel_waist = FUNNEL_WAIST_RADIUS * scale
        self.funnel_bottom = FUNNEL_BOTTOM_RADIUS * scale
        self.pool_radius = FUNNEL_POOL_RADIUS * scale

    def color(self, hex_string):
        """A palette colour, washed toward grey as the binding weakens."""
        c = hex_to_linear(hex_string)
        if self.fade <= 1e-6:
            return c
        grey = sum(c) / 3.0
        return mix(c, (grey, grey, grey), self.fade)

    def frame(self):
        """
        (e1, e2) spanning the orbit plane, e1 to the right of frame and e2
        toward the camera. Tying the frame to the view is what makes the phase
        mean the same thing for every pair however its plane is tilted.
        """
        eye = camera_eye()
        e1 = self.normal.cross(eye)
        if e1.length < 1e-6:
            e1 = self.normal.cross(Vector((1.0, 0.0, 0.0)))
        e1 = e1.normalized()
        return e1, self.normal.cross(e1).normalized()

    def at(self, behind, radius=None):
        """A point on the orbit, `behind` fractions of a loop behind the electron."""
        e1, e2 = self.frame()
        r = self.orbit_radius if radius is None else radius
        th = radians(self.phase) - 2.0 * pi * self.direction * behind
        return self.centre + (e1 * cos(th) + e2 * sin(th)) * r

    def electron(self):
        return self.at(0.0)

    def reach(self):
        """Radius of everything this exciton draws, for spacing them out."""
        return self.orbit_radius + 0.5 * self.trail_width + self.electron_radius
# ============================================================================
#  Building one exciton
# ============================================================================

def build_particles(col, ex):
    """The electron and the hole, shaded so they read as spheres."""
    for label, centre, radius, core_hex, rim_hex, strength in (
            ("Electron", ex.electron(), ex.electron_radius, ELECTRON_COLOR,
             ELECTRON_RIM, ELECTRON_EMISSION),
            ("Hole", ex.centre, ex.hole_radius, HOLE_COLOR, HOLE_RIM,
             HOLE_EMISSION)):
        mat = particle_material(f"{label}_{ex.tag}", ex.color(core_hex),
                                ex.color(rim_hex), strength * ex.emission)
        v, f = sphere_mesh_data(radius, PARTICLE_SEGMENTS, PARTICLE_RINGS)
        obj = add_object(f"Exciton_{label}", v, f, mat, col)
        if obj:
            obj.location = centre


def build_particle_glow(col, ex):
    """
    One halo sphere per particle, fading out as a continuous gradient.

    Only the half of the sphere facing away from the camera is shaded. The
    particle is opaque and sits at the centre, so that far half is hidden
    wherever the particle covers it and visible everywhere outside it: the
    halo glows around the particle instead of washing a veil across it.
    """
    for label, centre, radius, core_hex, rim_hex in (
            ("Electron", ex.electron(), ex.electron_radius, ELECTRON_COLOR,
             ELECTRON_RIM),
            ("Hole", ex.centre, ex.hole_radius, HOLE_COLOR, HOLE_RIM)):
        mat = halo_material(f"{label}_Halo_{ex.tag}", ex.color(rim_hex),
                            ex.color(core_hex), GLOW_EMISSION * ex.emission,
                            ex.glow_alpha, GLOW_FALLOFF)
        v, f = sphere_mesh_data(1.0, 32, 22)
        obj = add_object(f"Exciton_{label}_Halo", v, f, mat, col)
        if obj:
            obj.location = centre
            r = radius * ex.glow_scale
            obj.scale = (r, r, r)


def build_motion_trail(col, ex):
    """
    A band of concentric strands lying flat in the orbit plane, each solid
    where it leaves the electron and then breaking into dashes that shorten
    and space out behind it until they run out.

    Three things fade together along a strand, which is what sells it as
    motion rather than a drawn dotted line: the dashes get shorter and further
    apart, each dash gets thinner, and the colour cools from the electron's
    own highlight to the dim tail colour. Across the band, the outer strands
    trail TRAIL_SHEAR longer so the tail feathers instead of ending on a
    straight edge, and they are dimmer by TRAIL_EDGE_FADE so it keeps a
    bright spine.
    """
    strands = max(int(TRAIL_STRANDS), 1)
    head = ex.color(TRAIL_HEAD_COLOR or ELECTRON_RIM)
    tail = ex.color(TRAIL_COLOR)
    falloff = max(TRAIL_FALLOFF, 0.01)
    resolution = max(int(TRAIL_RESOLUTION), 24)

    V, F, C = [], [], []
    dashes = 0
    for j in range(strands):
        # -1 at the inner edge of the band, +1 at the outer edge
        u = 0.0 if strands == 1 else (2.0 * j / (strands - 1) - 1.0)
        radius = ex.orbit_radius + 0.5 * ex.trail_width * u
        if radius <= 1e-3:
            continue
        level = max(1.0 - TRAIL_EDGE_FADE * abs(u), 0.0)
        length = min(max(ex.trail_length * (1.0 + TRAIL_SHEAR * u), 1e-3), 1.0)
        hot = tuple(v * ex.trail_level * level for v in head)
        cold = tuple(v * TRAIL_END_LEVEL * level for v in tail)

        for start, end in dash_spans(length):
            steps = max(int(resolution * (end - start)), 3)
            pts = [ex.at(start + (end - start) * (i / steps), radius)
                   for i in range(steps + 1)]
            radii, colors = [], []
            for i in range(steps + 1):
                behind = start + (end - start) * (i / steps)
                w = max(0.0, 1.0 - behind / length) ** falloff
                radii.append(ex.strand_radius *
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
    mat = vertex_color_emission_material(f"MotionTrail_{ex.tag}", "TrailColor",
                                         TRAIL_EMISSION * ex.emission)
    add_object("Exciton_MotionTrail", V, F, mat, col, colors=C,
               attr_name="TrailColor")
    return dashes


def build_light_funnel(col, ex, drop):
    """
    A funnel of light from the exciton down to the Cu plane, `drop` below it.

    It is an hourglass -- a mouth at the exciton, a pinch at FUNNEL_WAIST, a
    wide landing on the sheet -- with each end shaped on its own. Its alpha is
    written per vertex so it fades along its length, and the material fades it
    again at its own silhouette; without that second term the cone has a
    visible outline and reads as a solid shell rather than a shaft of light.
    """
    if drop <= 1e-6:
        return None
    fall = Vector(FUNNEL_DIRECTION)
    if fall.length < 1e-9:
        return None
    fall = fall.normalized()
    ref = Vector((1.0, 0.0, 0.0))
    if abs(fall.dot(ref)) > 0.9:
        ref = Vector((0.0, 1.0, 0.0))
    e1 = fall.cross(ref).normalized()
    e2 = fall.cross(e1).normalized()

    colour = ex.color(FUNNEL_COLOR or HOLE_RIM)
    rings = max(int(FUNNEL_RINGS), 2)
    segs = max(int(FUNNEL_SEGMENTS), 6)

    verts, faces, colors = [], [], []
    for i in range(rings + 1):
        t = i / rings                          # 0 at the exciton, 1 at the sheet
        r = max(funnel_profile(t, ex.funnel_top, ex.funnel_waist,
                               ex.funnel_bottom), 0.0)
        level = funnel_profile(t, FUNNEL_TOP_LEVEL, FUNNEL_WAIST_LEVEL,
                               FUNNEL_BOTTOM_LEVEL) * ex.funnel_level
        ring_centre = ex.centre + fall * (drop * t)
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
        build_funnel_pool(col, ex, ex.centre + fall * drop, e1, e2, colour)
    return obj


def build_funnel_pool(col, ex, centre, e1, e2, colour):
    """
    The disc of light where the funnel lands on the plane: bright at the
    middle, fading to nothing at its rim.

    A flat disc is the one place the silhouette trick cannot help -- every
    point on it faces the camera the same way -- so the falloff goes into the
    vertex alpha instead.
    """
    steps = max(int(FUNNEL_POOL_STEPS), 8)
    rings = 16
    verts, faces, colors = [], [], []

    centre_level = FUNNEL_POOL_LEVEL * ex.funnel_level
    verts.append((centre.x, centre.y, centre.z))
    colors.append(tuple(v * centre_level for v in colour) + (FUNNEL_POOL_ALPHA,))
    for i in range(1, rings + 1):
        t = i / rings
        edge = max(1.0 - t, 0.0) ** max(FUNNEL_POOL_FALLOFF, 0.01)
        c = tuple(v * centre_level * edge for v in colour) + (FUNNEL_POOL_ALPHA * edge,)
        for k in range(steps):
            a = 2.0 * pi * k / steps
            p = centre + (e1 * cos(a) + e2 * sin(a)) * (ex.pool_radius * t)
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


def build_exciton(parent, ex, index, drop):
    """Everything one exciton draws, in a collection of its own."""
    col = get_collection(f"Exciton_{index:02d}_s{ex.tag}", parent)
    lines_drawn = 0
    if SHOW_FIELD_LINES:
        lines, _failed = build_field_lines(ex)
        lines_drawn = len(lines)
        if lines:
            build_field_line_objects(col, ex, lines)
    dashes = build_motion_trail(col, ex) if SHOW_TRAIL else 0
    if SHOW_FUNNEL:
        build_light_funnel(col, ex, drop)
    build_particles(col, ex)
    if SHOW_PARTICLE_GLOW:
        build_particle_glow(col, ex)
    return dashes, lines_drawn
# ============================================================================
#  Placing the excitons over the sheet
# ============================================================================

def _too_close(x, y, placed, reach):
    """True if this spot crowds one already taken."""
    for px, py, preach in placed:
        if hypot(x - px, y - py) < (reach + preach) * EXCITON_SPACING:
            return True
    return False


def place_excitons(lattice, rng):
    """
    Work out where every exciton goes: the strongly bound ones over the
    metallic islands, the weakly bound ones scattered out on the host sheet.

    Strong pairs are dealt round the islands rather than sampled at random
    over them, so every island gets one before any island gets two -- with a
    handful of pairs and a handful of islands, random sampling leaves some
    islands bare often enough to matter.

    Weak pairs are rejected within WEAK_CLEARANCE island radii of every
    island, which is what keeps the two populations from mixing and reading
    as one scattered field.
    """
    islands = lattice["islands"]
    x0, x1, y0, y1 = lattice["bounds_xy"]
    plane_z = lattice["plane_z"]
    gap = lattice["layer_gap"]

    height = gap * EXCITON_HEIGHT
    jitter = gap * EXCITON_HEIGHT_JITTER
    placed = []
    out = []

    def spec(x, y, strength):
        z = plane_z + height + rng.uniform(-jitter, jitter)
        tilt = radians(EXCITON_TILT)
        normal = Vector((rng.uniform(-1.0, 1.0) * sin(tilt),
                         rng.uniform(-1.0, 1.0) * sin(tilt), 1.0))
        return Exciton(Vector((x, y, z)), strength, EXCITON_SCALE, normal,
                       rng.uniform(0.0, 360.0),
                       1 if rng.random() < 0.5 else -1)

    # --- strong: on the islands ---------------------------------------------
    want_strong = max(int(STRONG_EXCITONS), 0)
    if want_strong and not islands:
        print("[figure] no islands to gather the strong excitons on -- "
              "they will be spread over the sheet instead.")
    for n in range(want_strong):
        trial = None
        for _attempt in range(200):
            if islands:
                cx, cy, r, _ = islands[n % len(islands)]
                a = rng.uniform(0.0, 2.0 * pi)
                d = r * STRONG_SPREAD * sqrt(rng.random())
                x, y = cx + d * cos(a), cy + d * sin(a)
            else:
                x, y = rng.uniform(x0, x1), rng.uniform(y0, y1)
            cand = spec(x, y, STRONG_BINDING)
            if not _too_close(x, y, placed, cand.reach()):
                trial = (x, y, cand)
                break
            trial = (x, y, cand)
        x, y, cand = trial
        placed.append((x, y, cand.reach()))
        out.append(cand)

    # --- weak: out on the host sheet ----------------------------------------
    for _n in range(max(int(WEAK_EXCITONS), 0)):
        trial = None
        for _attempt in range(400):
            x, y = rng.uniform(x0, x1), rng.uniform(y0, y1)
            if any(hypot(x - cx, y - cy) < r * WEAK_CLEARANCE
                   for cx, cy, r, _ in islands):
                continue
            cand = spec(x, y, WEAK_BINDING)
            if not _too_close(x, y, placed, cand.reach()):
                trial = (x, y, cand)
                break
            trial = trial or (x, y, cand)
        if trial is None:
            print("[figure] nowhere left on the host sheet for a weak "
                  "exciton -- lower WEAK_CLEARANCE or EXCITON_SPACING.")
            break
        x, y, cand = trial
        placed.append((x, y, cand.reach()))
        out.append(cand)

    out.sort(key=lambda e: -e.strength)
    return out


def build_excitons(lattice, rng):
    """Place them, build them, and say what went where."""
    parent = get_collection("Excitons")
    excitons = place_excitons(lattice, rng)
    if not excitons:
        return []

    plane_z = lattice["plane_z"]
    total = 0
    traced = 0
    for n, ex in enumerate(excitons):
        drop = max(ex.centre.z - plane_z, 0.0)
        dashes, lines = build_exciton(parent, ex, n, drop)
        total += dashes
        traced += lines

    strong = [e for e in excitons if e.strength >= 0.5]
    weak = [e for e in excitons if e.strength < 0.5]

    where = "on the islands" if lattice["islands"] else "over the sheet"
    print(f"\n[figure] {len(excitons)} excitons: {len(strong)} strongly bound "
          f"{where}, {len(weak)} weakly bound out on the host sheet")
    for label, group in (("strong", strong), ("weak  ", weak)):
        if not group:
            continue
        e = group[0]
        lines = (f"{e.field_rings}x{e.field_azimuths} = "
                 f"{e.field_rings * e.field_azimuths} field lines"
                 if SHOW_FIELD_LINES else "field lines off")
        print(f"  {label}: binding {e.strength:.2f}, orbit "
              f"{e.orbit_radius:.1f} A, particles "
              f"{e.hole_radius:.1f}/{e.electron_radius:.1f} A, {lines}")
    print(f"  {total} tail dashes"
          + (f", {traced} field lines traced" if SHOW_FIELD_LINES else ""))
    print(f"  floating {excitons[0].centre.z - lattice['plane_z']:.1f} A over "
          f"the Cu plane, in a gap of {lattice['layer_gap']:.1f} A")
    return excitons


# ============================================================================
#  Camera, lights, backdrop and render settings
# ============================================================================

def camera_eye():
    """
    Unit vector from the scene toward the camera, from CAMERA_ELEVATION and
    CAMERA_AZIMUTH alone.

    The excitons need the view direction before the camera exists -- their
    orbit phase is quoted relative to it -- and only the camera's DISTANCE
    depends on the geometry, so the direction can be settled up front.
    """
    elev = radians(min(max(CAMERA_ELEVATION, -89.0), 89.0))
    azim = radians(CAMERA_AZIMUTH)
    return Vector((cos(elev) * cos(azim), cos(elev) * sin(azim), sin(elev)))


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
        c = mix(bottom, top, smoothstep(t) if BACKGROUND_SMOOTH else t)
        y = (2.0 * t - 1.0) * half_h
        for side in (-1.0, 1.0):
            p = centre + right * (side * half_w) + up * y
            verts.append((p.x, p.y, p.z))
            colors.append(c)
    for i in range(rows):
        b = 2 * i
        faces.append((b, b + 1, b + 3, b + 2))

    mat = vertex_color_emission_material("Background", "BackgroundColor", 1.0)
    obj = add_object("Background", verts, faces, mat,
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

    eye = camera_eye()
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

def apply_draw_mode():
    """
    Hide whichever half DRAW does not want.

    Both halves are always BUILT, and the camera is always fitted to all of
    it. That is the whole point of one file: the lattice pass and the exciton
    pass come out of the same camera, so they line up exactly when they are
    composited, instead of being matched by eye afterwards.
    """
    if DRAW not in ("both", "lattice", "excitons"):
        raise ValueError("DRAW must be 'both', 'lattice' or 'excitons'; got "
                         f"{DRAW!r}")
    if DRAW == "both":
        return
    hide_excitons = DRAW == "lattice"
    for name, obj in bpy.data.objects.items():
        if obj.type != "MESH":
            continue
        is_exciton = name.startswith("Exciton_")
        if is_exciton == hide_excitons:
            for attr in ("hide_render", "hide_viewport"):
                try:
                    setattr(obj, attr, True)
                except Exception:
                    pass
    print(f"[figure] DRAW = {DRAW!r}: the other half is built but hidden, so "
          "this pass is framed exactly as the full scene is")


def setup_scene(points, bounds):
    if BACKGROUND_MODE not in ("gradient", "flat", "transparent"):
        raise ValueError("BACKGROUND_MODE must be 'gradient', 'flat' or "
                         f"'transparent'; got {BACKGROUND_MODE!r}")
    scene = bpy.context.scene
    engine = setup_render(scene)
    setup_world(scene)
    if ADD_LIGHTS:
        setup_lights(scene)
    cam = setup_camera(scene, points)

    # the exciton pass has to come out on a transparent film, whatever the
    # background is set to, or it cannot be laid over the lattice pass
    want_backdrop = BACKGROUND_MODE in ("gradient", "flat") and DRAW != "excitons"
    if want_backdrop:
        x0, x1, y0, y1, z0, z1 = bounds
        target = Vector((0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1)))
        reach = max((Vector(p) - target).length for p in points)
        bottom = hex_to_linear(BACKGROUND_COLOR)
        top = (hex_to_linear(BACKGROUND_TOP_COLOR)
               if BACKGROUND_MODE == "gradient" else bottom)
        add_background_plane(scene, cam, target, Vector((0.0, 0.0, 1.0)),
                             bottom, top, reach)
    elif DRAW == "excitons" and BACKGROUND_MODE != "transparent":
        print("[figure] DRAW = 'excitons', so the backdrop is left off and "
              "the pass comes out with an alpha channel")

    if ADD_GLOW_COMPOSITOR:
        try:
            setup_compositor(scene)
        except Exception as e:
            print(f"[figure] compositor skipped: {e}")
    else:
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
            print(f"[figure] COMPOSITION_MODE is {COMPOSITION_MODE!r}, so "
                  "COMPOSITION_SERIES is ignored: this is one sheet, not a "
                  "series.")
        comps = [None]          # None means "the mode decides, site by site"
    else:
        comps = COMPOSITION_SERIES if COMPOSITION_SERIES else [METAL_FRACTION]
    rng = random.Random(RANDOM_SEED)

    x_off = 0.0
    box = None
    panels = []
    for n, x in enumerate(comps):
        if x is not None and not 0.0 <= x <= 1.0:
            raise ValueError(f"composition must be between 0 and 1; got {x}")
        panel = build_panel(struct, x, n + 1, x_off, rng)
        panels.append(panel)
        b = panel["bounds"]
        box = b if box is None else (min(box[0], b[0]), max(box[1], b[1]),
                                     min(box[2], b[2]), max(box[3], b[3]),
                                     min(box[4], b[4]), max(box[5], b[5]))
        x_off += panel["width"] + PANEL_GAP

    # the excitons go over the first panel, which is the one the islands the
    # report talks about belong to
    excitons = []
    if STRONG_EXCITONS or WEAK_EXCITONS:
        lattice = dict(panels[0])
        lattice["bounds_xy"] = (box[0], box[1], box[2], box[3]) if len(panels) == 1 \
            else panels[0]["bounds_xy"]
        excitons = build_excitons(lattice, random.Random(EXCITON_SEED))

    engine = None
    if SETUP_SCENE:
        points = drawn_points()
        if not points:
            points = [(x, y, z) for x in box[0:2] for y in box[2:4]
                      for z in box[4:6]]
        engine = setup_scene(points, box)
    apply_draw_mode()

    print(f"\n[figure] Done. Metallic sites are {COLORS['metal_site']}, host "
          f"sites {COLORS['host_site']}.")
    if SETUP_SCENE:
        described = {"transparent": "transparent (RGBA PNG)",
                     "flat": f"flat {BACKGROUND_COLOR}",
                     "gradient": f"gradient {BACKGROUND_COLOR} -> "
                                 f"{BACKGROUND_TOP_COLOR}"}[BACKGROUND_MODE]
        print(f"[figure] engine {engine}, background {described}, camera "
              f"{'orthographic' if ORTHOGRAPHIC else f'{CAMERA_LENS:.0f} mm'} "
              f"at {CAMERA_ELEVATION:.0f} deg above the sheet")

    if RENDER_NOW:
        if not OUTPUT_PATH:
            print("[figure] RENDER_NOW is on but OUTPUT_PATH is empty; "
                  "nothing was written.")
        else:
            path = os.path.abspath(os.path.expanduser(OUTPUT_PATH))
            folder = os.path.dirname(path)
            if folder and not os.path.isdir(folder):
                os.makedirs(folder, exist_ok=True)
            bpy.context.scene.render.filepath = path
            print(f"[figure] rendering to {path} ...")
            bpy.ops.render.render(write_still=True)
            print("[figure] render written.")


main()
