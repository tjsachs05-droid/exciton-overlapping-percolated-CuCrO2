# exciton-overlapping-percolated-CuCrO2

Diagram of an exciton model overlapping percolated CuCrO2 (hexagonal crystal
structure). Consists of two Python Blender scripts that are overlaid with
transparent backgrounds in PowerPoint.

| script | draws |
| --- | --- |
| `cucro2perc.py` | the CuCrO2 delafossite lattice with the A-site sublattice coloured to show 2D site percolation (p_c = 1/2) |
| `exciton.py` | a bound electron-hole pair with its real dipole field lines, a volumetric aura, and optional rings |

## Running them

Blender → Scripting → Open → edit the `USER PARAMETERS` block at the top →
Alt+P. Each script builds its own scene, camera and render settings, so
`F12` renders straight away.

The aura in `exciton.py` is a volume, so the viewport has to be in
**Rendered** shading to see it — Material Preview does not draw volumes.

## Transparent backgrounds

Both scripts default to `TRANSPARENT_BACKGROUND = True`, which sets

* `render.film_transparent` — nothing is drawn behind the geometry, and
* `image_settings.color_mode = "RGBA"` — the alpha channel survives the save.

Save as PNG (the default) and drop the file onto a slide; the two renders
line up as separate layers. `exciton.py` additionally folds its bloom into
the alpha channel (`GLOW_IN_ALPHA`), because the compositor's Glare node only
adds colour — without that step every pixel of glow outside the geometry
keeps alpha 0 and disappears the moment the image is composited.

Set `OUTPUT_PATH` and `RENDER_NOW = True` in either script to render to a
file as soon as it finishes building.

## Knobs worth knowing

`exciton.py`

* `PARTICLE_SHADING`, `PARTICLE_LIGHT_DIR`, `PARTICLE_RIM_*` — how strongly
  the electron and hole read as 3-D spheres rather than flat discs.
* `ELECTRON_EMISSION` / `HOLE_EMISSION` — keep these below ~2.5; past that
  the sphere clips to white and the shading gradient disappears.
* `AURA_DENSITY`, `AURA_FALLOFF`, `AURA_EMISSION` — the volumetric aura.
* `SHOW_RINGS`, `RING_MODE` (`"equator"` or `"orbit"`), `RING_COUNT`,
  `RING_TILT` — rings encircling the pair.

`cucro2perc.py`

* `METAL_FRACTION` / `COMPOSITION_SERIES` — one panel, or several side by
  side across the percolation threshold.
* `COLOR_MODE` — `"species"`, `"spanning"` or `"clusters"`.
* `ORTHOGRAPHIC`, `CAMERA_MARGIN`, `ADD_LIGHTS` — framing and shading.
