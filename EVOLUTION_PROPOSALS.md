# Proposals: a more diverse evolutionary drive

The goal is a world where several different kinds of tree persist side by side, without adding species, fitness functions or hand-written niches. Each proposal below is one rule, reuses structures that already exist, and adds at most a couple of constants. Three of them remove genes.

The predicted effects are hypotheses. Proposals 0, 1, 2, 4 and 7 have been implemented, along with three root factors that were not on the original list. Two trials are reported below; proposals 3, 5 and 6 are untested. The "what the sim does today" notes are from reading the code and from the headless runs done while tuning branching.

## Why the sim converges on one form today

Diversity needs at least one of four things: trade-offs (no single best design), variation in space or time (different best designs in different places), an advantage to being rare, and limited mixing (so local winners are not swamped). The sim currently has almost none of them.

| What the sim does today | Consequence |
|---|---|
| Seed output is capped by the attempt rate (0.35 per year) and a fixed clutch of 3. One reproduction event costs about 0.19 MJ; the median adult holds about 0.6 MJ. | A 1 m stick already reproduces at the maximum rate. Extra size, leaf area or branches buy almost nothing. |
| Shade comes only from wood. A limb's shadow is 1.25× its radius wide and lets 35% through. With 4,000 trees on the plot, average occlusion was about 0.01. | There is effectively no competition for light, so height and crown shape are not under selection. |
| Soil starts uniform, recovers to the same value everywhere, and the sun never moves. | One environment, so one optimum. |
| Seeds land uniformly 3–40 m away on a 100 m plot. | The population is well mixed; a local variant cannot hold a patch. |
| Several costs are themselves genes: branch cost, seed cost, reproduction overhead, minimum reserve, establishment probability (mutable up to 0.95). | These only race to their clamps. That is directional change, not diversity, and it erodes the costs that would create trade-offs. |
| Mutation is a small Gaussian step on every gene, every generation. Integer genes are rounded, so a clutch of 3 needs a 17% jump to change. | Populations cannot cross a valley to reach a different design, as the split-signal cliff showed. |
| A limb's brain sees its own length, light and load, but not where it sits in the tree. | A trunk and a twig in the same state must behave the same, so many architectures cannot be expressed at all. |

## Proposal 0: measure diversity first

Give every founding tree a lineage id that seeds inherit unchanged, and record it in the snapshot.

- **Why:** without a number, "more diverse" is a matter of looking at the screen. Surviving lineage count and the spread of a few traits (height at maturity, limbs per tree, seed energy, wood density) make every later proposal testable.
- **Cost:** one integer in the genome, one more colour mode in the renderer, two or three lines on the metrics panel.
- **Success test for everything below:** after the plot fills, more than one lineage persists, and trait histograms have more than one peak.

## Proposal 1: make seed output follow energy

At a reproduction event, the clutch is whatever the reproductive budget affords, instead of a fixed count.

- **Why:** this is the missing link between a tree's body and its fitness. Once energy buys offspring, leaf area, branching, height and root investment all come under selection, and different ways of earning energy can compete.
- **Cost:** one line where the clutch is chosen; `seedsPerReproduction` becomes a minimum, or is removed.
- **Risk:** the population already grows without limit (4,000+ trees by year 140). This needs proposal 2 or 3 alongside it so that establishment, not seed supply, is the bottleneck.

## Proposal 2: one seed-size gene instead of four

Keep `seedInitialEnergy` and derive the rest from it:

- seed cost = seed energy × the existing provisioning multiplier;
- establishment probability rises with seed energy and saturates;
- dispersal distance falls with seed energy.

Remove `seedEnergyCost`, `propaguleEstablishmentProbability` and `seedsPerReproduction` as independent genes.

- **Why:** this is the classic many-small versus few-large trade-off. Small seeds travel and colonise gaps; large seeds establish under competition. Neither wins everywhere, and shorter dispersal for heavy seeds gives local lineages a chance to hold ground.
- **Cost:** two short formulas with one constant each; three genes fewer.

## Proposal 3: shade from leaves, on the existing grid

Reuse the 2 m soil grid as a light grid. Each tree adds its leaf area to the cells it covers, tagged with height. Light reaching a limb is reduced by the leaf area above it in that cell, with exponential falloff. This replaces the per-limb ray test against other trees.

- **Why:** it makes light a contested resource. Tall trees gain, but only by paying for trunk, transport and breakage risk, and the shaded floor becomes a second habitat for trees that can live on little light. It is also the main source of an advantage to being rare: in a tall forest, a cheap understory tree faces little competition of its own kind.
- **Pair with:** a small maintenance cost per unit of leaf area. Today more leaf is always better up to what the roots supply. With a cost, dense foliage pays in sun and loses in shade, which separates sun and shade strategies.
- **Cost:** one accumulation pass per substep and one constant for the falloff. It should be cheaper than the current occluder queries, which dominate run time.

## Proposal 4: uneven ground

Let each soil cell recover toward its own capacity instead of toward 1.0, with capacity set by a fixed smooth gradient across the plot.

- **Why:** rich and poor ground favour different root investment, seed size and growth speed. A gradient also makes the result visible as bands of different trees.
- **Cost:** one function of position and one changed line in soil recovery.

## Proposal 5: storms

Once a year draw a gust strength, and for that step add a sideways load proportional to each limb's height and exposed area before the existing breakage check.

- **Why:** disturbance is what keeps fast colonisers in the game against slow dominant trees. Storms would fell tall thin trees, and the gaps and detritus they leave are openings for small-seeded pioneers. It also gives wood density, thickness and bending strength a reason to vary.
- **Cost:** one random number per year and one term in the torque sum. No new death rule: the existing mechanics do the rest.

## Proposal 6: clean up the genome and the mutation rule

- **Remove the cost genes:** branch cost, reproduction overhead and minimum reserve become constants. They are prices set by the world, and letting them evolve only erases them.
- **Rare large mutations:** with a small probability per gene (around 2%), draw the step from a distribution ten times wider. Small steps refine a design; rare large ones let a lineage reach a different one.
- **Cost:** fewer genes, and one extra branch in the mutation lambda.

## Proposal 7: let a limb know where it is

Add two brain inputs: the limb's branch depth, and its tip height relative to the tree's height.

- **Why:** this is what lets a trunk keep climbing while side branches spread, or the reverse. Conifer-like and broadleaf-like shapes need limbs in the same local state to behave differently by position.
- **Cost:** two inputs and 14 more weights. No hidden layer.
- **Related, optional:** allow a limb to split again after it has grown a further set length, instead of only once. Today every branch on a limb comes from a single whorl.

## First trial: proposals 0, 4 and 7

Proposals 0, 4 and 7 are implemented. To separate their effects, four builds were compared, two runs each: neither change (baseline), the soil gradient only, the position inputs only, and both. Every run started from 600 founders spread over the whole plot, each its own lineage, and ran for 30 minutes of wall time. All eight reached year 100, so that is where they are compared.

| At year 100 | Baseline | Gradient only | Position inputs only | Both |
|---|---|---|---|---|
| Trees | 3,053 / 3,940 | 866 / 829 | 6,922 / 3,796 | 741 / 637 |
| Surviving lineages (of 600) | 104 / 92 | 16 / 19 | 86 / 96 | 20 / 25 |
| Effective lineages | 31.2 / 18.9 | 3.2 / 2.6 | 1.8 / 14.6 | 10.3 / 10.7 |
| Adults that have branched | 7.1% / 5.6% | 6.6% / 4.0% | 3.5% / 6.0% | 2.5% / 2.1% |
| Seed energy, sd / mean | 0.24 / 0.21 | 0.21 / 0.22 | 0.22 / 0.21 | 0.32 / 0.19 |
| Mean spread per brain weight | 0.22 / 0.22 | 0.22 / 0.20 | 0.21 / 0.22 | 0.23 / 0.22 |

What the runs show:

- **Neither change raised diversity in this window.** Trait spread grew at the same rate in all four builds, which is what accumulating mutation alone would give. Wood density, leaf area, growth rate, root uptake and the split bias look the same as the two rows shown.
- **The gradient cut lineage count.** Founders on poor ground failed early (600 trees fell to about 180–260 by year 40), and the population stayed four to five times smaller. Fewer trees means faster loss of lineages; these runs cannot separate that from any effect of the gradient itself.
- **No local adaptation along the gradient.** The correlation between each genetic trait and position along the gradient stayed within ±0.2, the same range as in the builds with uniform soil. Adult height did follow fertility (correlation 0.07–0.26), but that is the same genome growing better on better ground.
- **No consistent direction of selection.** Trait means moved by less than 0.1 standard deviations per decade over years 50–100, with the sign differing between replicates.
- **One possible sweep with the position inputs.** In one of the two position-input runs a few lineages took over (effective lineages fell to 1.8 while the population grew to 6,922), and the weights from branch depth to elongation and to splitting rose by about 0.15 standard deviations per decade. The other replicate showed nothing similar, so this is a lead, not a result.

Limits of the trial: two runs per build, about ten generations, and no run had filled the plot (populations were still growing at year 100). The test set out in proposal 0, "after the plot fills", has not been reached.

This fits the diagnosis above. With seed output capped, richer ground does not change which design is best, and with seeds travelling 3–40 m across a 100 m gradient any local variant is swamped. Proposals 1 and 2 are what would give the gradient something to select on.

## Second trial: proposals 1 and 2, and roots that matter

Implemented for this trial:

- **Proposal 1:** the clutch is whatever the reproductive budget affords, up to 60 seeds per event.
- **Proposal 2:** seed size is one gene. A seed costs four times its reserve, establishes with probability 1 / (1 + (16 kJ / (soil × reserve))²), and travels up to 20 m × √(default size / size), kept between 3 and 40 m. The seedling starts with the reserve its parent packed. The separate seed-cost, establishment and clutch-size genes are gone.
- **Roots forage where they are:** each root limb takes up from, and depletes, its own soil cell, and roots branch by the same rule as shoots.
- **Root upkeep:** absorbing root surface costs 10 MJ per m² per year.
- **Anchorage:** a tree uproots when shoot weight × height exceeds 10 × root weight × root reach.
- **Leaf area supported per root area** no longer mutates; it was a free gain.

The steep gradient plus the first upkeep setting (30 MJ) made the default 15-tree start die out in about a third of runs, so the poor edge was raised from 30% to 50% fertility and upkeep lowered to 10 MJ. With those values 24 of 24 sixty-year runs survived.

Six runs from 300 founders, 25 minutes each: four with everything above, and two controls with upkeep and anchorage switched off. All reached year 80. "Drive" is the change in a trait's mean per decade, in units of its spread within the population.

| At year 80 | Full (four runs) | Control (two runs) |
|---|---|---|
| Trees | 3,721 – 4,264 | 4,937 / 5,184 |
| Surviving lineages (of 300) | 39 – 53 | 62 / 62 |
| Root uptake gene, mean (starts at 1.00) | 1.07 – 1.23 | 1.13 / 1.15 |
| Root uptake gene, drive | +0.07 to +0.23 | +0.14 / +0.15 |
| Wood density, mean (starts at 506) | 390 – 447 | 400 / 450 |
| Wood density, drive | −0.13 to −0.28 | −0.12 / −0.25 |
| Seed energy, mean (starts at 11.7 kJ) | 11.8 – 13.2 kJ | 11.8 / 12.1 kJ |
| Seed energy, drive | +0.01 to +0.12 | +0.01 / +0.04 |
| Root share of biomass | 23.0% – 23.6% | 22.8% / 23.7% |
| Adults with branched roots | 0.8% – 2.6% | 2.6% / 2.0% |
| Adults with branched shoots | 0.7% – 2.1% | 3.2% / 2.2% |

What the runs show:

- **There is now consistent selection.** Root uptake rose and wood density fell in all six runs; leaf area rose in all six by a smaller amount. In the first trial no trait moved in the same direction across replicates.
- **Roots respond through the cheapest route.** The uptake-per-area gene climbs, while root share of biomass, root growth rate and root branching have barely moved. It climbs just as fast in the controls, so the pressure comes from soil running short as the plot fills (mean resource supply fell from 0.69 to about 0.43), and the 10 MJ upkeep is too small to hold it back.
- **Seed size is drifting up**, toward the 16 kJ and above that the establishment rule favours, clearly in three of the four full runs. There is still no sign of larger seeds on poorer ground: the correlation with position is slightly negative in all six runs but never beyond −0.09.
- **Anchorage and upkeep may hold branching back a little**: shoot branching is 0.7–2.1% with them and 2.2–3.2% without. That is two control runs, so it is weak evidence.
- **Lineage loss is as fast as before.** Nothing here was expected to slow it.

Limits: eight or so generations, populations still growing at year 80, and the controls only isolate upkeep and anchorage. Per-root uptake, root branching and the seed rules were on in every run.

The obvious next step for roots is to close the cheap route: make uptake per area a constant, or raise its upkeep until it pays to build root instead. Then the pressure that is now visible would have to act on root size and shape.

One loophole to know about: a seedling's first shoot and root are free. Their size is set by genes and not paid for out of the seed's reserve, so seed size does not yet limit how big a seedling starts.

## Suggested order

1. Proposal 0, so there is a baseline.
2. Proposals 1 and 2 together. They make size matter and set up the seed trade-off.
3. Proposal 3. This is the largest single change and the one most likely to produce coexisting forms.
4. Proposal 4 or 5, whichever is more interesting to watch. Both add variation; doing both at once makes it hard to tell which one worked.
5. Proposals 6 and 7 at any point; they are independent of the rest.

## Deliberately left out

- **Sexual recombination.** It adds variation, but needs mate choice and a crossover rule. Worth revisiting only if lineages stay too uniform after the changes above.
- **A hidden layer in the brain.** Position inputs give most of the expressive gain for far less.
- **Root architecture and soil depth layers.** Soil is charged only at the cell under the trunk, so root shape has nothing to act on yet. This is a second niche axis for later.
- **Seasons or a moving sun.** Variation in time is useful, but storms give it with one number instead of a calendar.
