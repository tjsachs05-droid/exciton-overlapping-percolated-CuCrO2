# exciton-overlapping-percolated-CuCrO2

Diagram of an exciton model overlapping percolated CuCrO2 (hexagonal crystal
structure). Consists of two Python Blender scripts that are overlaid with
transparent backgrounds in PowerPoint.

| script | draws |
| --- | --- |
| `cucro2perc.py` | the CuCrO2 delafossite lattice with the A-site sublattice coloured to show 2D site percolation (p_c = 1/2) |
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

The gradient is painted into a small image and composited behind the render
rather than being a world shader, which is the only approach that survives an
orthographic camera — under one, every view ray points the same way and a
world-space gradient has nothing left to vary over.

Set `OUTPUT_PATH` and `RENDER_NOW = True` in either script to render to a
file as soon as it finishes building.

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

* `METAL_FRACTION` / `COMPOSITION_SERIES` — one panel, or several side by
  side across the percolation threshold.
* `COLOR_MODE` — `"species"`, `"spanning"` or `"clusters"`.
* `LAYER_GAP` — extra vertical distance between the Cu planes and the CrO2
  slabs. Bonds, neighbours and the percolation analysis are all computed on
  the true geometry first, so raising it stretches the vertical Cu-O struts
  instead of breaking them.
* `SHOW_VERTICAL_BONDS` / `VERTICAL_BOND_ANGLE` — drop the bonds that run
  between layers. On CuCrO2 this removes the Cu-O struts and leaves the CrO6
  octahedra intact, since only the Cu-O bonds lie along the stacking axis.
* `ORTHOGRAPHIC`, `CAMERA_MARGIN`, `ADD_LIGHTS` — framing and shading.

`LAYER_GAP` shows best with `SHOW_BONDS = True` and `N_A_PLANES` above 1;
with the defaults (`SHOW_BONDS = False`, one A-plane kept) there is a single
gap in the model and no bonds drawn across it.
