# exciton-overlapping-percolated-CuCrO2

Excitons over a percolating CuCrO2 sheet, drawn in Blender as one figure:
`cucro2_exciton.py`.

The lattice is a Cu/Pd plane under a CrO2 slab, with the A-site substitution
laid out as large metallic islands on an otherwise gold host sheet. Over it
float excitons — an electron orbiting a hole, bound by dipole field lines,
with a dashed motion tail and a funnel of light dropping to the plane below. Strongly bound ones gather on the
islands; weakly bound ones are scattered out over the host sheet, smaller and
fainter.

The palette comes from `Reference 1` in this repo: a deep blue field, gold
host sites and purple metallic ones, a hot pink-red hole with a bright cyan
electron going round it, and pale cyan tail dashes.

## Running it

Blender → Scripting → Open → edit the `USER PARAMETERS` block at the top →
Alt+P. It builds the scene, camera, lights and render settings, so `F12`
renders straight away.

## Why one file

The two halves used to be separate scripts, overlaid afterwards. They are
together because they have to be:

* **Placement needs the lattice.** "Strong excitons on the islands" can only
  be answered by whatever decided where the islands are.
* **One camera, or the perspective does not match.** An exciton over the far
  edge of the sheet has to be smaller and higher in frame than one over the
  near edge, by exactly what the lens and camera position dictate.
* **Occlusion.** An exciton behind the CrO2 slab, or passing behind an atom,
  needs real depth sorting.

`DRAW` still gives separate images when they are wanted — `"lattice"` or
`"excitons"` instead of `"both"`. Both halves are always built and the camera
is always fitted to all of it, so the passes line up by construction rather
than by eye. The exciton pass always comes out on a transparent film.

## The excitons

| | |
| --- | --- |
| `STRONG_EXCITONS` / `WEAK_EXCITONS` | how many of each |
| `STRONG_BINDING` / `WEAK_BINDING` | how tightly each population is bound, 0 to 1 |
| `EXCITON_SCALE` | how much larger than life they are drawn |
| `EXCITON_HEIGHT` | how high they float, as a fraction of the gap to the CrO2 slab |
| `STRONG_SPREAD` | how far out into its island a strong exciton may sit |
| `WEAK_CLEARANCE` | how far weak ones stay clear of every island |

The binding strength drives the whole model at once — orbit radius, particle
size, field lines, tail, funnel, brightness and colour — because at these
sizes one
quantity changing by 30% does not read, but everything loosening and fading
together does. Each `WEAK_` parameter is the far end of the quantity named
after it, and the value in the exciton model block is the strong end.

Weak pairs get smaller particles through `WEAK_PARTICLE_SCALE`, which keeps
the electron/hole size ratio while shrinking both.

The field lines are the clearest reading of binding strength. `SHOW_FIELD_LINES`
turns the bundle on; a pair at full binding is held by `FIELD_LINE_RINGS` x
`FIELD_LINE_AZIMUTHS` streamlines and one at zero binding by
`WEAK_FIELD_LINE_RINGS` x `WEAK_FIELD_LINE_AZIMUTHS`. Both counts blend with
the strength like everything else and are then rounded, so at the defaults the
strong pairs carry 7 x 18 = 126 lines and the weak ones, sitting at
`WEAK_BINDING` 0.15 rather than 0, come out at 4 x 10 = 40. They are real
streamlines of the dipole field, not decorative arcs, integrated with RK4 from
the hole to the electron.

## Other knobs worth knowing

* `COMPOSITION_MODE` — `"islands"` (the default), `"gradient"`, or
  `"uniform"`. Islands are blobs with ragged coasts rather than discs:
  `ISLAND_COUNT`, `ISLAND_RADIUS`, `ISLAND_WOBBLE`, `ISLAND_EDGE`.
* `CRO2_SLABS` / `N_A_PLANES` / `LAYER_GAP` — how many CrO2 slabs sit around
  the A-plane and how far apart they are pulled. Bonds and the percolation
  analysis are computed on the true geometry first, so the gap stretches the
  vertical Cu-O struts rather than breaking them.
* `SHOW_VERTICAL_BONDS` — off by default, which drops those struts and leaves
  the CrO6 octahedral network in the slab and the channel network in the plane.
* `BACKGROUND_MODE` — `"gradient"`, `"flat"` or `"transparent"`. The backdrop
  is an emissive quad behind the scene with the gradient in its vertex
  colours, not a world shader and not a composited image: a world gradient has
  nothing to vary over under an orthographic camera, and a generated image is
  rebuilt by Blender from its own settings whenever it is re-evaluated, so a
  gradient written from Python is simply not there at render time.
* `ORTHOGRAPHIC`, `CAMERA_ELEVATION`, `CAMERA_AZIMUTH`, `CAMERA_LENS`,
  `CAMERA_MARGIN` — the camera fits the vertices actually built rather than an
  estimated bounding box, so nothing is clipped. At `CAMERA_MARGIN = 1.0` the
  outermost atom sits exactly on the frame edge.

A wide flat sheet seen from a low `CAMERA_ELEVATION` is far wider than it is
tall, so it fills the frame across and leaves air above and below. The script
prints how much of the frame it uses on every run.
