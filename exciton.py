"""
exciton.py
==========
Blender script: a glowing exciton -- a bound electron-hole pair -- with real
dipole field lines running between the two particles and a soft aura that
fades from opaque at the core to transparent at the edge.

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

Everything glows via Emission shaders, so switch the viewport to Material
Preview or Rendered to see it. SETUP_SCENE also builds a dark background,
a camera framing the pair, and a compositor Glare node for the bloom.
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
BACKGROUND_COLOR = "#0A1430"   # deep blue world background

ELECTRON_EMISSION = 6.0        # emission strength of each particle
HOLE_EMISSION     = 5.0

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
FIELD_LINE_EMISSION  = 4.0
MID_GLOW_BOOST       = 2.2    # extra brightness at the midpoint of each line
                              # (this is what makes the central cord white-hot)
MID_GLOW_WIDTH       = 0.30   # 0..1, how much of the line the boost covers

INTEGRATION_STEP     = 0.045  # RK4 step length; smaller = smoother, slower
MAX_STEPS            = 20000
FIELD_LINE_TAPER     = 0.12   # fraction of each end that tapers to a point

# --- glow around each particle ---------------------------------------------
SHOW_PARTICLE_GLOW   = True
GLOW_SHELLS          = 5      # nested shells per particle
GLOW_OUTER_SCALE     = 2.1    # outermost shell radius / particle radius
GLOW_INNER_ALPHA     = 0.34   # alpha of the innermost shell
GLOW_FALLOFF         = 2.0    # >1 fades faster toward the outside

# --- aura enclosing the whole exciton ---------------------------------------
SHOW_AURA            = True
AURA_MODE            = "shells"   # "shells" -> nested transparent surfaces.
                                  #   Reliable everywhere, incl. Material
                                  #   Preview and Cycles.
                                  # "volume" -> a single volumetric object.
                                  #   Smoother and physically nicer, needs
                                  #   Rendered shading to be visible.
AURA_SHELLS          = 9          # shells, for AURA_MODE = "shells"
AURA_MARGIN          = 2.6        # how far the aura extends past the particles
AURA_INNER_ALPHA     = 0.20       # alpha of the innermost aura shell
AURA_FALLOFF         = 2.4        # >1 = tighter, faster fade to transparent
AURA_EMISSION        = 1.1
AURA_DENSITY         = 0.35       # for AURA_MODE = "volume"
AURA_RESOLUTION      = (48, 24)   # (azimuthal, polar) segments

# --- scene ------------------------------------------------------------------
SETUP_SCENE          = True   # world colour, camera, and bloom in compositor
ADD_GLOW_COMPOSITOR  = True
CLEAR_SCENE          = True

SPHERE_SEGMENTS      = 48     # particle mesh resolution
SPHERE_RINGS         = 32


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


def _set_blend(mat, alpha):
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
    try:
        mat.show_transparent_back = False     # cleaner stacked shells
    except Exception:
        pass


def emission_material(name, color, strength, alpha=1.0):
    """Emission shader, optionally alpha-blended (for glow shells)."""
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


def vertex_color_emission_material(name, attr_name, strength):
    """Emission driven by a per-vertex colour attribute."""
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
    emis.inputs["Strength"].default_value = strength
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.location = (-150, 0)
    attr.attribute_name = attr_name
    nt.links.new(attr.outputs["Color"], emis.inputs["Color"])
    nt.links.new(emis.outputs["Emission"], out.inputs["Surface"])
    return mat


def volume_material(name, color, density, strength):
    """Emissive volume with a distance-driven density falloff."""
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (500, 0)
    vol = nt.nodes.new("ShaderNodeVolumePrincipled")
    vol.location = (250, 0)
    vol.inputs["Color"].default_value = (*color, 1.0)
    vol.inputs["Emission Color"].default_value = (*color, 1.0)
    vol.inputs["Emission Strength"].default_value = strength

    # density = ramp( distance from object centre ), so it thins out toward
    # the edge and is thickest at the core
    tex = nt.nodes.new("ShaderNodeTexCoord")
    tex.location = (-400, 0)
    length = nt.nodes.new("ShaderNodeVectorMath")
    length.operation = "LENGTH"
    length.location = (-200, 0)
    ramp = nt.nodes.new("ShaderNodeValToRGB")
    ramp.location = (0, 0)
    ramp.color_ramp.elements[0].position = 0.0
    ramp.color_ramp.elements[0].color = (density, density, density, 1.0)
    ramp.color_ramp.elements[1].position = 1.0
    ramp.color_ramp.elements[1].color = (0.0, 0.0, 0.0, 1.0)
    nt.links.new(tex.outputs["Object"], length.inputs[0])
    nt.links.new(length.outputs["Value"], ramp.inputs["Fac"])
    nt.links.new(ramp.outputs["Color"], vol.inputs["Density"])
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
    e_col = hex_to_linear(ELECTRON_COLOR)
    h_col = hex_to_linear(HOLE_COLOR)

    v, f = sphere_mesh_data(ELECTRON_RADIUS, SPHERE_SEGMENTS, SPHERE_RINGS)
    v = [(x + r_electron.x, y + r_electron.y, z + r_electron.z) for x, y, z in v]
    add_object("Exciton_Electron", v, f,
               emission_material("Electron", e_col, ELECTRON_EMISSION), col)

    v, f = sphere_mesh_data(HOLE_RADIUS, SPHERE_SEGMENTS, SPHERE_RINGS)
    v = [(x + r_hole.x, y + r_hole.y, z + r_hole.z) for x, y, z in v]
    add_object("Exciton_Hole", v, f,
               emission_material("Hole", h_col, HOLE_EMISSION), col)


def build_particle_glow(col, r_electron, r_hole):
    """Nested emissive shells around each particle: opaque core -> clear edge."""
    for label, centre, radius, hexcol in (
            ("Electron", r_electron, ELECTRON_RADIUS, ELECTRON_RIM),
            ("Hole", r_hole, HOLE_RADIUS, HOLE_RIM)):
        base = hex_to_linear(hexcol)
        for k in range(GLOW_SHELLS):
            t = (k + 1) / GLOW_SHELLS                    # 0..1 outward
            r = radius * (1.0 + (GLOW_OUTER_SCALE - 1.0) * t)
            # sample the falloff at the shell's mid-band so the outermost
            # shell keeps a small non-zero alpha instead of vanishing
            ta = (k + 0.5) / GLOW_SHELLS
            alpha = GLOW_INNER_ALPHA * (1.0 - ta) ** GLOW_FALLOFF
            if alpha < 0.004:
                continue
            mat = emission_material(f"{label}_Glow{k}", base,
                                    ELECTRON_EMISSION * 0.35, alpha)
            v, f = sphere_mesh_data(r, 32, 20)
            v = [(x + centre.x, y + centre.y, z + centre.z) for x, y, z in v]
            add_object(f"Exciton_{label}_Glow{k:02d}", v, f, mat, col)


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
    e_col = hex_to_linear(ELECTRON_RIM)
    h_col = hex_to_linear(HOLE_RIM)
    aura_col = mix(h_col, e_col, 0.5)
    long_axis, short_axis = _aura_shape()
    segs, rings = AURA_RESOLUTION

    # rotation taking +Z to the pair axis
    rot = Vector((0.0, 0.0, 1.0)).rotation_difference(axis).to_matrix()

    if AURA_MODE == "volume":
        mat = volume_material("Aura_Volume", aura_col, AURA_DENSITY,
                              AURA_EMISSION)
        v, f = sphere_mesh_data(1.0, segs, rings,
                                scale=(short_axis, short_axis, long_axis))
        v = [tuple((rot @ Vector(p)) + centre) for p in v]
        add_object("Exciton_Aura_Volume", v, f, mat, col)
        _enable_eevee_volumes()
        return

    for k in range(AURA_SHELLS):
        t = (k + 1) / AURA_SHELLS                        # 0..1 outward
        s = 0.45 + 0.55 * t                              # innermost not tiny
        # sample the falloff mid-band so the outermost shell still renders
        ta = (k + 0.5) / AURA_SHELLS
        alpha = AURA_INNER_ALPHA * (1.0 - ta) ** AURA_FALLOFF
        if alpha < 0.003:
            continue
        # warm toward the electron end, cool toward the hole end
        shell_col = mix(aura_col, e_col if k % 2 else h_col, 0.25)
        mat = emission_material(f"Aura_Shell{k}", shell_col,
                                AURA_EMISSION, alpha)
        v, f = sphere_mesh_data(1.0, segs, rings,
                                scale=(short_axis * s, short_axis * s,
                                       long_axis * s))
        v = [tuple((rot @ Vector(p)) + centre) for p in v]
        add_object(f"Exciton_Aura_Shell{k:02d}", v, f, mat, col)


def _enable_eevee_volumes():
    try:
        ev = bpy.context.scene.eevee
        for attr, val in (("use_volumetric_lights", True),
                          ("volumetric_tile_size", "2"),
                          ("volumetric_samples", 128)):
            if hasattr(ev, attr):
                setattr(ev, attr, val)
    except Exception:
        pass


def setup_scene(centre, axis):
    scene = bpy.context.scene

    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    if bg:
        bg.inputs["Color"].default_value = (*hex_to_linear(BACKGROUND_COLOR), 1.0)
        bg.inputs["Strength"].default_value = 1.0

    # camera looking at the pair, perpendicular to the axis
    long_axis, short_axis = _aura_shape()
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

    if ADD_GLOW_COMPOSITOR:
        try:
            scene.use_nodes = True
            nt = scene.node_tree
            nt.nodes.clear()
            rl = nt.nodes.new("CompositorNodeRLayers")
            rl.location = (-300, 0)
            glare = nt.nodes.new("CompositorNodeGlare")
            glare.location = (0, 0)
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
            comp = nt.nodes.new("CompositorNodeComposite")
            comp.location = (300, 0)
            nt.links.new(rl.outputs["Image"], glare.inputs["Image"])
            nt.links.new(glare.outputs["Image"], comp.inputs["Image"])
        except Exception as e:
            print(f"[exciton] compositor glare skipped: {e}")


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

    if SETUP_SCENE:
        setup_scene(centre, axis)

    print("[exciton] Done. Electron at "
          f"({r_electron.x:.2f}, {r_electron.y:.2f}, {r_electron.z:.2f}), "
          f"hole at ({r_hole.x:.2f}, {r_hole.y:.2f}, {r_hole.z:.2f}).")
    print("[exciton] Switch the viewport to Material Preview or Rendered "
          "to see the emission and the aura.")


main()
