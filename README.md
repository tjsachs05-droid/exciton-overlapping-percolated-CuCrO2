# exciton-overlapping-percolated-CuCrO2

Diagram of an exciton model overlapping percolated CuCrO2 (hexagonal crystal
structure). Consists of two Python Blender scripts that are overlaid with
transparent backgrounds in PowerPoint.

| script | draws |
| --- | --- |
| `cucro2perc.py` | a Cu/Pd plane under a CrO2 slab, with large metallic islands scattered over an otherwise gold host sheet (or a composition gradient, or a uniform alloy) |
| `exciton.py` | a small electron orbiting a larger hole, with a Saturn-ring motion tail behind it and a volumetric aura between the two |

Both take their palette from `Reference 1` in this repo: a deep blue field, a
hot pink-red hole at the centre, a bright cyan electron going round it, pale
cyan rings in the tail, and gold field lines.

## Running them

Blender → Scripting → Open → edit the `USER PARAMETERS` block at the top →
Alt+P. Each script builds its own scene, camera and render settings, so
`F12` renders straight away.

The aura in `exciton.py` is a volume, so the viewport has to be in
**Rendered** shading to see it — Material Preview does not draw volumes.

## Backgrounds

`BACKGROUND_MODE` in each script takes one of three values:

| mode | what you get |
| --- | --- |
| `"gradient"` | the deep blue field of `Reference 1`, dark at the bottom and brighter towards the top |
| `"flat"` | a single `BACKGROUND_COLOR` |
| `"transparent"` | a transparent film and an RGBA PNG, for layering |

The two scripts default to opposite ends of that: **`cucro2perc.py` is the
base layer and carries the blue**, **`exciton.py` is the overlay and stays
transparent** so it drops on top of it. If the slide itself is blue, set both
to `"transparent"` instead.

Whenever the film is transparent, `image_settings.color_mode = "RGBA"` keeps
the alpha channel through the save — leave it at RGB and the background comes
back black. `exciton.py` also folds its bloom into the alpha channel
(`GLOW_IN_ALPHA`), because the compositor's Glare node only adds colour;
without that step every pixel of glow outside the geometry keeps alpha 0 and
disappears the moment the image is composited.

The backdrop is an emissive quad placed behind the scene with the gradient
baked into its vertex colours — not a world shader and not a composited
image. A world-space gradient has nothing to vary over under an orthographic
camera, and a generated image is worse: `bpy.data.images.new()` makes an
image whose source is GENERATED, and Blender rebuilds a generated image's
buffer from its own settings whenever it re-evaluates it, so a gradient poked
in from Python is simply not there at render time and the background comes
out transparent. Geometry and vertex colours have nothing to regenerate.

The film stays transparent in every mode; the backdrop is what makes the
render opaque where it is present.

Set `OUTPUT_PATH` and `RENDER_NOW = True` in either script to render to a
file as soon as it finishes building.

## Framing

Both cameras fit the **vertices actually built** rather than an estimated
bounding box. An estimate has to guess how far the drawn geometry reaches —
sphere tessellation, bond caps, the widest display radius, the width of the
exciton's tail band — and whatever it misses is what the frame clips off.

The fit itself is exact: every point is projected onto the camera's own right
and up axes, and because those two components do not change as the camera
slides along its view direction, each point sets a lower bound on the
distance and the largest of them frames the lot. Each script prints how much
of the frame the result uses.

A wide flat sheet seen from a low `CAMERA_ELEVATION` is far wider than it is
tall, so it fills the frame across and leaves air above and below — at 5
degrees it uses about a quarter of the frame's height. Raise the elevation or
widen `RESOLUTION` if that bothers you.

## Knobs worth knowing

`exciton.py`

* `PARTICLE_SHADING`, `PARTICLE_LIGHT_DIR`, `PARTICLE_RIM_*` — how strongly
  the electron and hole read as 3-D spheres rather than flat discs.
* `ELECTRON_EMISSION` / `HOLE_EMISSION` — keep these below ~2.5; past that
  the sphere clips to white and the shading gradient disappears.
* `AURA_DENSITY`, `AURA_FALLOFF`, `AURA_EMISSION` — the volumetric aura.
* `ORBIT_RADIUS`, `ORBIT_NORMAL`, `ORBIT_PHASE`, `ORBIT_DIRECTION` — where
  the electron is on its circle around the hole and which way it travels.
  `ORBIT_PHASE` is measured from the right of frame, so 0 puts the electron
  at the right-hand edge of the path and 90 nearest the camera.
* The motion tail is a **band**, not a line: `TRAIL_STRANDS` concentric
  strands lying flat in the orbit plane over `TRAIL_WIDTH`, like the rings of
  Saturn. `TRAIL_LENGTH` and `TRAIL_FALLOFF` set how far back it reaches and
  how fast it fades; `TRAIL_SHEAR` makes the outer strands trail longer so
  the tail feathers out instead of ending on a straight edge;
  `TRAIL_EDGE_FADE` dims the outer strands so the band has a bright spine.
* `TRAIL_RING_LEVEL` sets how brightly the rest of each circle shows. At 0
  the circles are cut short and only the tail is drawn.
* `SHOW_FIELD_LINES` — the full dipole streamline bundle, off by default.
  The physics is unchanged; it is just no longer what carries the picture.
* `AURA_MARGIN` — the aura is a spindle sitting directly between the two
  particles, so this is what makes it a tight bridge of light or a loose halo.
* `CAMERA_ELEVATION` — how far above the orbit plane the camera sits. At 0
  the tail is seen edge-on and renders as a straight line, so this is what
  opens it out into an ellipse.

`cucro2perc.py`

* `COMPOSITION_MODE` picks how the substitution is laid out:
  * `"islands"` (default) — a mostly-host sheet with `ISLAND_COUNT` large
    metallic islands on it. `ISLAND_RADIUS` sizes them as a fraction of the
    sheet's shorter side; `ISLAND_WOBBLE` and `ISLAND_EDGE` keep the outlines
    irregular and the coastlines ragged rather than drawn-on;
    `ISLAND_BACKGROUND` sprinkles a few stray metallic sites over the host.
    Three islands at the defaults cover about a third of the sheet, so gold
    stays the clear majority. The console reports each island's size and
    whether its metallic sites join into a single cluster.
  * `"gradient"` — x ramps across the sheet, `GRADIENT_MIN` on one side,
    `GRADIENT_MID` through the middle, `GRADIENT_MAX` on the other, with
    `GRADIENT_*_BAND` setting how much of the width each holds flat.
  * `"uniform"` — one composition everywhere, from `METAL_FRACTION` or
    `COMPOSITION_SERIES` for several panels side by side.
* `N_CELLS` — the islands and the gradient both need a wide sheet to read;
  below about 200 A sites the statistics are noise rather than structure.
* `CRO2_SLABS` / `N_A_PLANES` — how many complete CrO2 slabs sit around the
  A-plane: 0 for the bare plane, 1 (the default) for one above it, 2 for a
  CrO2 | Cu/Pd | CrO2 sandwich. The cuts run just inside the neighbouring
  A-planes, not midway to them: an O-Cr-O slab sits centred between two
  A-planes, so a midpoint cut would slice it in half.
* `ORTHOGRAPHIC`, `CAMERA_ELEVATION`, `CAMERA_AZIMUTH`, `CAMERA_LENS` — the
  camera is a perspective one by default, looking down at the sheet from
  `CAMERA_ELEVATION` degrees above its plane, which is what makes the layers
  read as solids. The azimuth decides which way the sheet runs across the
  frame; -90 lays the x axis left to right.
* `CAMERA_MARGIN` — air around the structure. The camera fits the geometry
  exactly, so at 1.0 the outermost atoms sit precisely on the frame edge and
  this is pure breathing room.
* `COLOR_MODE` — `"species"`, `"spanning"` or `"clusters"`.
* `LAYER_GAP` — extra vertical distance between the Cu planes and the CrO2
  slabs. Bonds, neighbours and the percolation analysis are all computed on
  the true geometry first, so raising it stretches the vertical Cu-O struts
  instead of breaking them.
* `SHOW_VERTICAL_BONDS` / `VERTICAL_BOND_ANGLE` — drop the bonds that run
  between layers. On CuCrO2 this removes the Cu-O struts and leaves the CrO6
  octahedra intact, since only the Cu-O bonds lie along the stacking axis.
* `CAMERA_MARGIN`, `ADD_LIGHTS` — framing and shading.

`LAYER_GAP` defaults to 7 A, which lifts the CrO2 slab clear of the Cu/Pd
plane so the two read separately. With `SHOW_VERTICAL_BONDS = False` the
vertical Cu-O struts between them are dropped, leaving the CrO6 octahedral
network in the slab and the metallic channel network in the plane.
