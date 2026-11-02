# Real-world constraints

How TreeSim's trees compare with real ones, the rules added to close the gap,
and what is still off. Parameters named here are in `params/default.json`.

## Why the trees were skinny

Measured from recorded frames (`tools/render_run.py` frame format):

| Run | Trees 2–5 m tall | Median trunk | Height / diameter |
| --- | --- | --- | --- |
| 20 m plot, year 150, before | 28 | 1.9 cm | 137 |
| 20 m plot, year 150, with secondary growth | 7 | 3.6 cm | 63 |
| 120 m plot, year 642, before | 197 | 1.6 cm | 155 |

Real trees rarely stand above a height/diameter (H/D) ratio of 100.
Open-grown trees run 20–40 and managed stands 60–80. Wind damage climbs
steeply past 80 and peaks near 100 ([Watt & Moore 2023][slender];
[Sheehan et al.][sheehan]; [COFORD windthrow guide][coford]).

The old model had no rule that tied a stem's girth to what stood above it.
Each limb thickened only when its growth brain asked to, and the founders'
brain asks only when a storm nearly breaks it. The trees also carry almost no
crown (see [Still off](#still-off)), so no rule forced them to thicken.

## What real trees do, and what the model now does

### Secondary growth (cambium)

Every woody axis of a real tree lays down a new ring each year. The tree does
not decide this limb by limb. Three constraints set how much wood each stem
needs. All three are now enforced, behind `wood.EnableSecondaryGrowth`,
and the cambium is paid before any other growth.

1. **Pipe model.** Each unit of leaf area is supplied by a fixed amount of
   conducting sapwood all the way to the roots (Shinozaki 1964). Its inverse,
   the Huber value, falls between about 1e-4 and 1e-3 m² of sapwood per m² of
   leaf. It is higher in dry climates and in tall trees
   ([McDowell et al. 2002][huber-height]; [Huber value review][huber]).
   Wood never shrinks, so pipes left by shed branches stay on as heartwood.
   `HuberValue = 2.5e-4`.
   Each root carries a share of the crown's demand in proportion to the
   absorbing surface beyond it.
2. **Continuity.** A limb is at least as thick as the thickest limb it bears.
   Leonardo da Vinci's rule, that the area of a stem equals the summed area
   of its branches, comes out of designing branches against wind loads
   ([Eloy 2011, APS abstract][eloy]).
3. **Elastic stability.** A column of diameter D buckles under its own
   weight above Greenhill's height `0.79 (E / ρg)^(1/3) D^(2/3)`. Trees stay
   about a factor of 4 below it in the open ([McMahon 1973][mcmahon];
   [review][mcmahon-review]). Shaded trees put height before girth and run
   closer to a factor of 2. The model scales the factor between
   `OpenBucklingSafetyFactor = 4` and `ShadedBucklingSafetyFactor = 1.5`
   with how much sun the tree's green parts get. E follows wood density, as
   it does in real wood ([Niklas & Spatz][niklas]). The default is 9 GPa for
   green wood at 506 kg/m³; green wood weighs about twice its dry mass.
   A stem's load against this limit is also fed to the growth brain as part
   of its torque margin.

The cambium has no fixed pace: it lays on what the tree can afford, keeping
each limb's shape. Open-grown trees put on 0.5–2.5 cm of diameter a year
([diameter growth studies][dbh]), which is what an energy-limited tree in
this model should land near.

### Only sapwood respires

Heartwood is dead tissue. Maintenance respiration scales with sapwood
volume, not stem volume ([Ryan et al. 1995][ryan]). The old model charged
upkeep on the whole volume, which made a thick trunk ruinous to keep. Now
only the outer `SapwoodDepthM = 0.04` m respires (real sapwood is 1–10 cm
deep).

### What it costs

Girth slows height growth. On the same seed of the 20 m run, the tallest
tree at year 150 reaches 3.0 m without the cambium, 2.5 m with it under the
old fixed growth rates, and 2.2 m with it now that energy alone limits
growth. Without rate caps the plot fills in (about 450 trees) by year 30
instead of 120, and small trees carry 13 limbs instead of 5. Trees 2–5 m
tall stand at a height/diameter ratio of 64 either way.

## Terrain

`terrain.TerrainEnabled` (off by default) raises the plot into a fractal
heightfield that is identical on every rank. Trees grow upright from wherever
they root. Each limb treats the ground as level at its tree's base, so
height, wind profile and what counts as underground are all measured from
there. Light is traced in 3D, so uphill trees shade downhill ones.
The ground itself casts no shadow.

| Effect | Real-world basis | Model |
| --- | --- | --- |
| Water collects in hollows, drains off ridges | Topographic wetness index `ln(a / tan β)`, from the area draining through a point and its slope, predicts soil moisture and site quality ([Sørensen et al. 2006][twi]; [Appalachian oak study][oak]) | D8 flow routing on the terrain grid. Water supply × (1 + `TerrainWetnessEffect` × standardized TWI) |
| Thin soil on steep slopes | Erosion strips slopes, and soil gathers below | Soil capacity × (1 − `TerrainSlopeSoilLoss` × slope) |
| Wind speeds up over crests | Jackson–Hunt theory. Near-surface speed-up of about 1.6× on low hills, up to 2× on steep ones ([Taylor & Teunissen 1987][hill-wind]; [Lubitz 2007][lubitz]) | Storm and design wind × (1 + `TerrainWindSpeedup` × standardized height above the surrounding 20 m), kept within 0.6–1.6 |
| Sunward slopes run dry | Aspect sets how much radiation reaches the soil, and equator-facing slopes are warmer and drier | Water demand × (1 + `TerrainSunSlopeDryness` on a 45° slope facing the noon sun) |
| Shorter season higher up | About 3–5 days lost per 100 m ([Moser et al. 2010][larch]; [Barnard et al. 2017][colorado]) | Growing season × (1 − `TerrainSeasonLossPer100M` × height / 100 m). This only matters for relief of hundreds of metres |

`params/hills_120m.json` runs the 120 m forest on hills with 25 m of relief.
Runs with terrain record `terrain.bin`. `render_run.py` draws it as a
hill-shaded surface in the corner view, as contours in the top view, and as
the ground profile in the side view. The live viewer draws it as its floor.

## Still off

The model's largest departure from real trees is no longer girth:

- **No crowns.** At year 642 of the 120 m run, a 2–5 m tree holds about
  0.02 m² of light-catching area. A real sapling that size holds several m²,
  and a closed forest holds 3–8 m² of leaf per m² of ground (leaf area index).
  This one holds about 0.0003. The founders' growth brain barely branches
  and never makes flat leaves, and 640 years of evolution have not found
  them. A stem only needs girth in proportion to the crown and wind it
  carries, so these trees will stay thin-looking whips until they grow
  crowns. Seeding founders with a leafy habit (side buds that start as thin
  blades and expand within a season) would change this. It would also change
  what the evolution experiment starts from.
- **Leaves earn too much.** A sunlit limb fixes CO₂ at the light-saturated
  rate `MaxLeafAssimilationMolCo2PerM2Second = 12e-6` for 45% of the year.
  That is about 170 mol CO₂ per m² of leaf per year. Real forests average
  about 30–40 (gross production of 1.5–2.5 kg C per m² of ground spread
  over a leaf area index of 4–5), because leaves are light-limited most of
  the day and stomata close in heat and drought. This matters once trees
  have crowns.
- **Leaf expansion.** Growth now has no fixed rate, so a blade can widen as
  fast as its tree can pay. Founders still widen only about 1/200 as fast as
  they lengthen, so blades must evolve.
- **Leaves are evergreen woody limbs.** There is no leaf fall or annual
  replacement.

[slender]: https://www.frontiersin.org/journals/forests-and-global-change/articles/10.3389/ffgc.2023.1188094/pdf
[sheehan]: https://www.scionresearch.com/__data/assets/pdf_file/0004/59395/NZJFS1221982SHEEHAN269-280.pdf
[coford]: https://www.coford.ie/media/coford/content/publications/2018/SM23ManagingWindthrow160418.pdf
[huber-height]: https://www.fsl.orst.edu/rna/Documents/publications/The%20relationship%20between%20tree%20height%20and%20leaf%20area%20sapwood%20area%20ratio%20pub3571.pdf
[huber]: https://prometheusprotocols.net/structure/architecture/water-balance-traits-leaf-to-sapwood-area-ratio/
[eloy]: https://meetings-archive.aps.org/dfd/2011/l24/6
[mcmahon]: https://www.sccs.swarthmore.edu/users/08/bblonder/phys120/docs/mcmahon.pdf
[mcmahon-review]: https://pmc.ncbi.nlm.nih.gov/articles/PMC524850
[niklas]: https://repository.si.edu/handle/10088/18591
[dbh]: https://www.afs-journal.org/articles/forest/full_html/2010/08/f09393/F2.html
[ryan]: https://www7.nau.edu/mpcer/direnet/publications/publications_r/files/Ryan_et_al_1995.pdf
[twi]: https://hess.copernicus.org/articles/10/101/2006/
[oak]: https://vtechworks.lib.vt.edu/items/27dfcb1e-3f75-4f84-8612-4b2b89ad3fbd
[hill-wind]: https://data-ww3.ifremer.fr/BIB/Taylor_Teunissen_BLM1987.pdf
[lubitz]: https://research.engineering.ucdavis.edu/wind/wp-content/uploads/sites/17/2014/03/lubitz-2007.pdf
[larch]: https://www.climatology.uni-mainz.de/files/2012/03/Moser_2010_TreePhys.pdf
[colorado]: https://scholar.colorado.edu/concern/articles/c821gk39d?locale=it
