"""
exciton.py
==========
Blender script: a glowing exciton -- a bound electron-hole pair -- with real
dipole field lines running between the two particles, shaded 3-D particles,
a volumetric aura that fades smoothly from the core outward, and optional
rings around the pair. Renders on a transparent background so the image can
be dropped straight onto a PowerPoint slide.

Run: Blender -> Scripting -> Open -> edit USER PARAMETERS -> Alt+P.

What is actually being drawn
----------------------------
The field lines are not decorative arcs: they are streamlines of the electric
field of two equal and opposite point charges (+ at the hole, - at the
electron), integrated with RK4 through

        E(r) = (r - r_h)/|r - r_h|^3  -  (r - r_e)/|r - r_e|^3

Every line therefore starts on the hole, ends on the electron, and bows
outward by exactly the amount the dipole geometry dictates. Lines launched
close to the axis stay in a tight bundle through the middle; lines launched at
larger angles swing out to the sides. FIELD_LINE_AXIS_BIAS controls how many
are launched near the axis, which is what makes the middle dense.

Colours default to the ones sampled from the reference image: a warm
orange-red electron and a cool blue hole. Swap them if you prefer the more
common convention -- nothing else depends on the choice.

Shading
-------
Everything glows via Emission shaders, so switch the viewport to Material
Preview or Rendered to see it. The particles are not flat discs: their
emission colour and strength follow a half-Lambert gradient across the
surface plus a rim term at the silhouette, so they read as spheres without
needing any scene lights. The aura and the particle haloes are continuous
gradients -- a real volume for the aura, a view-angle alpha falloff for the
haloes -- rather than stacks of nested shells.

SETUP_SCENE builds the camera, the render settings and a compositor Glare
node for the bloom. With TRANSPARENT_BACKGROUND the film is transparent and
the output is RGBA, and the compositor pushes the bloom into the alpha
channel too, so the glow survives being composited over a slide.
"""

import bpy
import math
import os
from math import cos, sin, sqrt, pi
from mathutils import Vector, Matrix

# ============================================================================
#  USER PARAMETERS
# ============================================================================

# --- the pair ---------------------------------------------------------------
SEPARATION       = 6.0        # centre-to-centre distance (Blender units)
AXIS_DIRECTION   = (0.0, 0.0, 1.0)   # electron sits at +axis, hole at -axis
CENTER           = (0.0, 0.0, 0.0)

ELECTRON_RADIUS  = 1.25
HOLE_RADIUS      = 1.10

# Colours as sRGB hex (sampled from the reference image).
ELECTRON_COLOR   = "#F3724F"   # hot orange-red core
ELECTRON_RIM     = "#FEC45A"   # warmer rim / glow
HOLE_COLOR       = "#5FABF5"   # blue core
HOLE_RIM         = "#7CC3F6"   # lighter blue glow
FIELD_LINE_COLOR = "#FBDB8C"   # gold field lines
BACKGROUND_COLOR = "#0A1430"   # world colour, used only when
                               # TRANSPARENT_BACKGROUND is False

ELECTRON_EMISSION = 1.5        # emission strength of each particle. Push
HOLE_EMISSION     = 1.3        # these much past ~2.5 and the whole sphere
                               # clips to white, which is what made the
                               # particles read as flat coloured circles in
                               # the first place -- the shading gradient has
                               # to stay inside the displayable range.

# --- field lines ------------------------------------------------------------
SHOW_FIELD_LINES     = True
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
GLOW_OUTER_SCALE     = 2.4    # halo radius / particle radius
GLOW_ALPHA           = 0.60   # peak alpha of the halo, reached just outside
                              # the particle's edge
GLOW_FALLOFF         = 2.6    # >1 fades faster toward the halo's outer edge
GLOW_EMISSION        = 1.8

# --- aura enclosing the whole exciton (a single emissive volume) -----------
SHOW_AURA            = True
AURA_MARGIN          = 2.6        # how far the aura extends past the particles
AURA_DENSITY         = 0.25       # density at the core of the aura. This
                                  # and AURA_EMISSION multiply; too much
                                  # density hazes over the particles instead
                                  # of glowing around them.
AURA_FALLOFF         = 2.4        # >1 = tighter, faster fade to nothing
AURA_EMISSION        = 2.5
AURA_RESOLUTION      = (64, 32)   # (azimuthal, polar) segments of the shell
                                  # that bounds the volume

# --- rings around the exciton ----------------------------------------------
SHOW_RINGS           = True
RING_MODE            = "equator"  # "equator" -> rings stacked along the pair
                                  #   axis, lying on the aura's surface like
                                  #   lines of latitude.
                                  # "orbit"   -> great-circle rings all
                                  #   through the centre, fanned around the
                                  #   axis like a textbook atom diagram.
RING_COUNT           = 3
RING_SPREAD          = 0.55       # "equator": outermost ring position, as a
                                  # fraction of the aura's half-length
RING_SCALE           = 1.04       # ring size / aura size (>1 sits outside it)
RING_TILT            = 0.0        # degrees. "equator": tips the rings off the
                                  # axis. "orbit": angle between each ring's
                                  # normal and the axis (0 -> uses 72 deg, so
                                  # the rings do not all coincide).
RING_TUBE_RADIUS     = 0.07
RING_SEGMENTS        = 10         # cross-section resolution of the ring tube
RING_RESOLUTION      = 192        # points around each ring
RING_COLOR           = ""         # "" = gradient from the hole colour to the
                                  # electron colour along the axis; or a hex
                                  # string such as "#FBDB8C" for a flat colour
RING_EMISSION        = 1.8
RING_ALPHA           = 1.0        # <1 makes the rings translucent

# --- scene / render ---------------------------------------------------------
SETUP_SCENE            = True   # camera, render settings, and bloom
TRANSPARENT_BACKGROUND = True   # render with an alpha channel, so the image
                                # drops straight onto a slide
ADD_GLOW_COMPOSITOR    = True
GLOW_IN_ALPHA          = True   # also push the bloom into the alpha channel;
                                # without this the glow around the exciton is
                                # invisible once the PNG is composited
GLOW_ALPHA_GAIN        = 1.6    # how strongly the bloom opens up the alpha
RENDER_ENGINE          = "EEVEE"   # "EEVEE" or "CYCLES". Both render the
                                   # volumetric aura; Cycles is cleaner and
                                   # slower.
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

    verts, vparam = [], []
    for i in range(n):
        t = tangents[i]
        u = normals[i]
        w = t.cross(u).normalized()
        r = radius * taper_scale(params[i])
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


def aura_volume_material(name, near_color, far_color, density, falloff,
                         strength):
    """
    Emissive volume whose density falls off from the centre of the aura to its
    surface, tinted from the hole's colour at one end of the pair axis to the
    electron's at the other.

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

    # colour gradient along the pair axis (object Z runs -1 .. +1)
    sep = _node(nt, "ShaderNodeSeparateXYZ", -800, 220)
    nt.links.new(tex.outputs["Object"], sep.inputs[0])
    zt = _math(nt, "MULTIPLY_ADD", -620, 220, value2=0.5, value3=0.5,
               clamp=True)
    nt.links.new(sep.outputs["Z"], zt.inputs[0])
    ramp = _two_stop_ramp(nt, -440, 220, near_color, far_color)
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
                attr.data[i].color = (c[0], c[1], c[2], 1.0)
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
        mat = particle_material(label, hex_to_linear(core_hex),
                                hex_to_linear(rim_hex), strength)
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
        mat = halo_material(f"{label}_Halo", hex_to_linear(rim_hex),
                            hex_to_linear(core_hex), GLOW_EMISSION,
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


def _aura_shape():
    """Prolate spheroid enclosing both particles."""
    half = SEPARATION * 0.5
    long_axis = half + max(ELECTRON_RADIUS, HOLE_RADIUS) + AURA_MARGIN
    short_axis = max(ELECTRON_RADIUS, HOLE_RADIUS) + AURA_MARGIN
    return long_axis, short_axis


def build_aura(col, centre, axis):
    """
    The aura is a single emissive volume: density AURA_DENSITY at the centre
    of the pair, falling smoothly to nothing at the surface of a prolate
    spheroid that encloses both particles.

    The mesh is a unit sphere and the shape lives entirely in the object's
    scale and rotation. That is not cosmetic -- it is what makes the volume
    shader work, because the shader reads Texture Coordinate -> Object, which
    is normalised by the object transform and therefore runs 0 at the core to
    1 at the surface no matter how large the aura is.
    """
    e_col = hex_to_linear(ELECTRON_RIM)
    h_col = hex_to_linear(HOLE_RIM)
    long_axis, short_axis = _aura_shape()
    segs, rings = AURA_RESOLUTION

    mat = aura_volume_material("Aura_Volume", h_col, e_col, AURA_DENSITY,
                               AURA_FALLOFF, AURA_EMISSION)
    v, f = sphere_mesh_data(1.0, segs, rings)
    obj = add_object("Exciton_Aura", v, f, mat, col)
    if obj is None:
        return
    obj.location = centre
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = Vector((0.0, 0.0, 1.0)).rotation_difference(axis)
    obj.scale = (short_axis, short_axis, long_axis)
    _enable_volumes()


def _ring_points(k, long_axis, short_axis):
    """
    Points of ring k in the local frame where the pair axis is +Z.

    "equator" puts ring k at height z on the aura's surface, so its radius is
    the spheroid's own radius there and the ring hugs the aura like a line of
    latitude. "orbit" builds a great circle of the unit sphere in a plane
    tilted off the axis and then scales it by the spheroid, which lands the
    ring on that same surface.
    """
    n = max(int(RING_RESOLUTION), 12)
    count = max(int(RING_COUNT), 1)
    L = long_axis * RING_SCALE
    S = short_axis * RING_SCALE
    tilt = math.radians(RING_TILT)
    pts = []

    if RING_MODE == "orbit":
        psi = 2.0 * pi * k / count
        beta = tilt if abs(RING_TILT) > 1e-6 else math.radians(72.0)
        nrm = Vector((sin(beta) * cos(psi), sin(beta) * sin(psi), cos(beta)))
        e1 = nrm.cross(Vector((0.0, 0.0, 1.0)))
        if e1.length < 1e-6:
            e1 = nrm.cross(Vector((1.0, 0.0, 0.0)))
        e1.normalize()
        e2 = nrm.cross(e1).normalized()
        for i in range(n + 1):
            th = 2.0 * pi * i / n
            u = e1 * cos(th) + e2 * sin(th)
            pts.append(Vector((u.x * S, u.y * S, u.z * L)))
    else:
        if count == 1:
            frac = 0.0
        else:
            frac = (-1.0 + 2.0 * k / (count - 1)) * RING_SPREAD
        z = frac * L
        r = S * sqrt(max(1.0 - min(abs(frac), 1.0) ** 2, 1e-4))
        ca, sa = cos(tilt), sin(tilt)
        for i in range(n + 1):
            th = 2.0 * pi * i / n
            x, y0, z0 = r * cos(th), r * sin(th), z
            pts.append(Vector((x, y0 * ca - z0 * sa, y0 * sa + z0 * ca)))
    return pts


def build_rings(col, centre, axis):
    """Rings encircling the exciton, as tubes swept around each ring path."""
    long_axis, short_axis = _aura_shape()
    rot = Vector((0.0, 0.0, 1.0)).rotation_difference(axis).to_matrix()
    e_col = hex_to_linear(ELECTRON_RIM)
    h_col = hex_to_linear(HOLE_RIM)
    flat = hex_to_linear(RING_COLOR) if RING_COLOR else None

    V, F, C = [], [], []
    for k in range(max(int(RING_COUNT), 1)):
        local = _ring_points(k, long_axis, short_axis)
        world = [(rot @ p) + centre for p in local]
        v, f, _ = tube_from_polyline(world, RING_TUBE_RADIUS, RING_SEGMENTS)
        if not v:
            continue
        base = len(V)
        V.extend(v)
        F.extend(tuple(base + i for i in face) for face in f)
        if flat is None:
            # colour every ring vertex by where it sits along the pair axis:
            # the hole's colour at one end, the electron's at the other
            span = long_axis * RING_SCALE
            for p in local:
                t = min(max(0.5 + 0.5 * p.z / span, 0.0), 1.0)
                C.extend([mix(h_col, e_col, t)] * RING_SEGMENTS)

    if not V:
        return
    if flat is None:
        mat = vertex_color_emission_material("Rings", "RingColor",
                                             RING_EMISSION, RING_ALPHA)
        add_object("Exciton_Rings", V, F, mat, col, colors=C,
                   attr_name="RingColor")
    else:
        mat = emission_material("Rings", flat, RING_EMISSION, RING_ALPHA)
        add_object("Exciton_Rings", V, F, mat, col)


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
    r.film_transparent = bool(TRANSPARENT_BACKGROUND)

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
        if TRANSPARENT_BACKGROUND:
            # the film hides the world anyway; killing its strength keeps it
            # from tinting the volumetric aura
            bg.inputs["Color"].default_value = (0.0, 0.0, 0.0, 1.0)
            bg.inputs["Strength"].default_value = 0.0
        else:
            bg.inputs["Color"].default_value = (
                *hex_to_linear(BACKGROUND_COLOR), 1.0)
            bg.inputs["Strength"].default_value = 1.0


def setup_camera(scene, centre, axis):
    long_axis, _ = _aura_shape()
    ref = Vector((1.0, 0.0, 0.0))
    if abs(axis.dot(ref)) > 0.9:
        ref = Vector((0.0, 1.0, 0.0))
    side = axis.cross(ref).normalized()
    dist = long_axis * 3.2
    loc = centre + side * dist + axis * (long_axis * 0.15)

    cam_data = bpy.data.cameras.new("ExcitonCamera")
    cam = bpy.data.objects.new("ExcitonCamera", cam_data)
    scene.collection.objects.link(cam)
    cam.location = loc
    direction = (centre - loc).normalized()
    cam.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    scene.camera = cam
    return cam


def setup_compositor(scene):
    """
    Glare for the bloom -- and, on a transparent film, a second pass that
    folds the bloom into the alpha channel.

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

    if TRANSPARENT_BACKGROUND and GLOW_IN_ALPHA:
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
    scene = bpy.context.scene
    engine = setup_render(scene)
    setup_world(scene)
    setup_camera(scene, centre, axis)
    if ADD_GLOW_COMPOSITOR:
        try:
            setup_compositor(scene)
        except Exception as e:
            print(f"[exciton] compositor skipped: {e}")
    return engine


# ============================================================================
#  Main
# ============================================================================

def main():
    if CLEAR_SCENE:
        clear_scene()

    axis = Vector(AXIS_DIRECTION)
    if axis.length < 1e-9:
        raise ValueError("AXIS_DIRECTION must be a non-zero vector")
    axis = axis.normalized()
    centre = Vector(CENTER)
    r_electron = centre + axis * (SEPARATION * 0.5)
    r_hole = centre - axis * (SEPARATION * 0.5)

    if SEPARATION <= ELECTRON_RADIUS + HOLE_RADIUS:
        print("[exciton] WARNING: SEPARATION is smaller than the two radii "
              "combined -- the particles overlap.")

    col = get_collection("Exciton")

    if SHOW_FIELD_LINES:
        lines, failed = build_field_lines(r_electron, r_hole)
        print(f"[exciton] traced {len(lines)} field lines "
              f"({failed} discarded), "
              f"{sum(len(l) for l in lines)} integration points total")
        if lines:
            build_field_line_objects(col, lines)

    build_particles(col, r_electron, r_hole)

    if SHOW_PARTICLE_GLOW:
        build_particle_glow(col, r_electron, r_hole)

    if SHOW_AURA:
        build_aura(col, centre, axis)

    if SHOW_RINGS:
        build_rings(col, centre, axis)
        print(f"[exciton] {RING_COUNT} ring(s), mode {RING_MODE!r}")

    engine = None
    if SETUP_SCENE:
        engine = setup_scene(centre, axis)

    print("[exciton] Done. Electron at "
          f"({r_electron.x:.2f}, {r_electron.y:.2f}, {r_electron.z:.2f}), "
          f"hole at ({r_hole.x:.2f}, {r_hole.y:.2f}, {r_hole.z:.2f}).")
    if SETUP_SCENE:
        print(f"[exciton] engine {engine}, background "
              f"{'transparent (RGBA)' if TRANSPARENT_BACKGROUND else BACKGROUND_COLOR}")
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
