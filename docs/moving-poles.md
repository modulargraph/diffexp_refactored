# Moving poles of the master bases

A denominator factor of a level's connection A(x, eps) that mixes the merge
parameter with epsilon, such as `11*x*eps-5*x-eps+7`, is a singular point
whose position depends on the dimension. This note records how many there
are, what removes them, and what that buys. Measurements are for 2.2.0 on the
4-core container of [performance.md](performance.md),
one thread (`OMP_NUM_THREADS=1`), built-in IBP solver with its default seeds.
The machine was shared, so times within about 10 % of each other are equal.

The short answer: the default bases have few moving poles (none in the kite
and the four-loop banana, one in the planar double box, three in the double
pentagon). Each goes away when one master is exchanged, in three of the four
cases for another dotted integral that the level's reduction already
contains, and then the connection gets smaller with it. That is worth a
third of the double box's numerical time and nothing measurable on a
double-pentagon sample. The epsilon windows do not change at all: they are
not caused by the moving poles.

## Counts

`diffexp_master_basis_report CONFIGURATION --cache DIRECTORY` factors every
denominator of every level's connection. Distinct irreducible factors by what
they depend on; "order" is the sum of the highest powers of the mixed factors;
"eps pole" is the highest power of 1/eps in an entry.

| Family | Level | Masters | x only | eps only | Mixed (order) | eps pole |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| kite | 1 / 2 / 3 | 16 / 10 / 3 | 9 / 11 / 6 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| box_triangle | 1 / 2 / 3 / 4 | 4 / 4 / 3 / 2 | 3 / 2 / 5 / 4 | 1 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 1 / 0 / 0 / 0 |
| double_box_planar | 1 | 12 | 4 | 2 | 0 | 2 |
| | 2 | 11 | 5 | 2 | 1 (1) | 0 |
| | 3 / 4 / 5 | 7 / 3 / 2 | 6 / 7 / 4 | 1 / 0 / 0 | 0 / 0 / 0 | 1 / 0 / 0 |
| banana4 | 1 / 2 / 3 | 15 / 7 / 3 | 6 / 6 / 5 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| henn_double_pentagon_x0 | 1 | 69 | 8 | 5 | 2 (2) | 2 |
| | 2 | 33 | 14 | 6 | 1 (1) | 2 |
| | 3 / 4 | 23 / 11 | 27 / 20 | 3 / 2 | 0 / 0 | 2 / 1 |
| | 5 / 6 | 3 / 2 | 7 / 4 | 0 / 0 | 0 / 0 | 0 / 0 |

The mixed factors, all linear in x and all simple poles:

| Level | Factor | Position at eps = 0 | Rows |
| --- | --- | --- | --- |
| double box 2 | `11*x*eps-5*x-eps+7` | 7/5 | the two top-sector masters |
| pentagon 1 | `3*x*eps+x+6*eps+4` | -4 | one six-line sector and the top sector |
| pentagon 1 | `9*x*eps+x-33*eps-17` | 17 | four sectors and the top sector |
| pentagon 2 | `12*x*eps-x-6*eps` | 0 | eleven masters |

Only the last one comes near the integration interval: it sits at
6 eps/(12 eps - 1), which is -0.068 at eps = 1/100 and merges with the
endpoint x = 0 at eps = 0.

The internal basis scan of the pipeline ranks by connection size within 256
trials. On the double box it installs a level-2 basis of 7197 characters
(from 8002) that still has the mixed factor. Its ranking is not the
obstacle: allowed 4000 trials (it uses 1336, 1.8 s) the same scan arrives at
2311 characters and no mixed factor. The exchanges that matter are among its
candidates but beyond its trial bound.

## What removes them

A moving pole of a Feynman integral is apparent: the integrals are regular
there, and the pole says that the masters become linearly dependent at
f(x, eps) = 0. `include/diffexp/moving_poles.hpp` states this as a lattice
closure. The basis without the pole is unique up to changes that are regular
at f, and it is reached by replacing one master per step with a row that has
a simple pole at f. Three kinds of change were tried.

**Rescaling masters (diagonal T).** Not possible for any of the four
factors. The pole sits in two rows of one sector (or more) and in the entries
between them, so the difference constraints on the exponents contain a cycle
(`--rescale 1` of the report tool).

**Another integral of the same sector.** Works for every factor. On the
double box all 125 single exchanges (29 integrals with a cached reduction, at
every position where their reduction is nonzero) were tried: five lower the
mixed order and all five remove it. Each of the five uses an integral whose
own reduction contains the mixed factor (three of the five such integrals).
This is what the lattice says: only an integral outside the span of the
masters at f = 0 can help.

| Level | Master | Exchanged for | Connection characters |
| --- | --- | --- | --- |
| double box 2 | `[1,1,1,1,2]` | `[2,1,1,1,1]` | 8002 -> 3321 |
| pentagon 2 | `[1,1,2,0,1,0]` | `[2,1,1,0,1,0]` | 257309 -> 222802 |
| pentagon 1 | `[1,1,0,1,1,1,2]`, `[1,1,1,2,0,1,1]` | `[2,1,0,1,1,1,1]` with one numerator, and the closure row below | 514677 -> 620301 |

An exchange is admitted only if the other denominators of the integral's
reduction add no singular point in x to the connection. For the second
pentagon factor the twelve integrals tried all add one (without that
condition level 1 ends at 469094 characters, with one more factor in x), so
the closure row is used there.

**The closure row itself.** If row i of A has a pole of order p at f, the
row r/f with r = (f^p A_i) reduced modulo f is the missing generator. It needs
no integral and removes all four factors in one step each, but the result is
larger (double box level 2: 14226 characters) because r carries the epsilon
dependence of the position x0(eps). It is the fallback. For a factor that is
not linear in x the reduction is arithmetic in the residue field at f
(polynomials in x modulo f over the rational functions of epsilon).

**Seed preference.** The solver's default ordering already prefers dotted
masters at equal degree. With the other ordering it offers (numerators
kept, `--dots-first 1` in the report tool, a separate cache) the bases are
worse, which confirms the rule of thumb of Smirnov and Smirnov
(arXiv:2002.08042) and Usovitsch (arXiv:2002.08173):

| Family | Mixed factors, default ordering | Mixed factors, numerators kept |
| --- | --- | --- |
| kite | 0 | 0 |
| double_box_planar | 1 (level 2) | 2 (levels 1, 2; one quadratic in x) |
| banana4 | 0 | 3 (levels 1, 1, 2; two not linear in x) |

On these worse bases an exchange for a dotted integral removes three of the
five factors (connection sizes 3425 -> 1486, 8417 -> 4279, 25394 -> 5726).
The two factors of banana4 level 1 need closure rows: the linear one is a
relation of order three (three steps), the quadratic one takes one step in
its residue field. All nine mixed factors met in either ordering are removed.

## What it buys

Two options install the result. `diffexp ft` uses `remove` by default and
`keep` for an `--epsilon-shift` sample, where the removal costs exact work
and buys nothing measurable:

- `--moving-poles remove` (`ft`): an internal basis per level, the smallest
  connection among cached-integral exchanges and closure rows. The level's
  requests to its child do not change, and the automatic basis scan of the
  level starts from this basis (`scan_physical_bases`).
- `--moving-poles exchange` (`prepare`, `ft`): the level's masters
  themselves, cached integrals of the same sector only. The levels below are
  prepared for the new masters and every later stage sees an ordinary level.

Commands (the double box twice back to back, the others once):

```sh
export OMP_NUM_THREADS=1
diffexp ft examples/feynman/double_box_planar.json --cache C --no-numerical-cache --json
diffexp ft examples/feynman/double_box_planar.json --cache C --no-numerical-cache --json --moving-poles remove
diffexp ft examples/feynman/double_box_planar.json --cache C --no-numerical-cache --json --moving-poles exchange
# sampled: add --epsilon-shift 1/100 --epsilon-order 0
```

| Run | keep | remove | exchange |
| --- | ---: | ---: | ---: |
| Double box, numerical seconds | 19.6, 20.0 | 13.7, 14.3 | 13.9, 13.4 |
| ... transport | 12.6, 12.6 | 8.8, 9.4 | 9.0, 8.7 |
| ... level 2 lower / upper leg | 2.1 / 6.1 | 1.6 / 3.0 | 1.7 / 3.5 |
| ... level 2 connection after the scan, characters | 7197 | 2800 | 2800 |
| ... charts, level 2 lower / upper | 29 / 24 | 29 / 24 | 29 / 24 |
| ... child asked through epsilon, levels 1-5 | 3, 5, 7, 9, 9 | 3, 5, 7, 9, 9 | 3, 5, 7, 9, 9 |
| Double box sample at 1/100, numerical seconds | 3.56, 3.59 | 3.12, 3.27 | 3.30, 3.26 |
| ... level 2 lower / upper leg | 0.35 / 0.55 | 0.26 / 0.33 | 0.28 / 0.41 |
| ... charts, level 2 lower / upper | 27 / 19 | 26 / 19 | 26 / 19 |
| Kite, numerical seconds | 5.2 | 5.1 | |
| Four-loop banana, numerical seconds | 49.7, 52.0 | 53.6, 48.6 | 50.3 |
| ... transport | 27.6 | 27.6 | 26.8 |

The kite and the banana have no moving pole; the option leaves their
computation identical (coefficients equal to every digit, same transport
time). Their totals scatter because the endpoint series run on every hardware
thread and feel the other jobs. On the double box the coefficients of the
three modes agree to 1.5e-39 relative to the largest, expanded, and to 7e-42
sampled.

The gain on the double box, close to a third of the numerical time, is the
smaller level-2 connection of the better master (and of the scan that
continues from it). The moving pole itself costs nothing there: it is at
x = 7/5, the chart counts do not change, and one chart fewer in a sampled leg
is the whole effect on the path.

One sample of the double pentagon at eps = 1/100 (the same commands with
`--epsilon-shift 1/100 --epsilon-order 0`; keep was run before and after):

| | keep | remove | keep again |
| --- | ---: | ---: | ---: |
| Numerical seconds | 292.2 | 277.2 | 262.5 |
| ... transport | 256.0 | 240.9 | 228.3 |
| ... level 2 lower / upper leg | 38.1 / 124.0 | 36.5 / 109.1 | 34.7 / 103.6 |
| Charts, level 2 lower / upper | 35 / 58 (24 halved) | 35 / 58 (24 halved) | 35 / 58 (24 halved) |
| Charts, level 1 lower / upper | 23 / 33 | 24 / 33 | 23 / 33 |
| Endpoint overlap, levels 1 / 2 / 3 | 1/128, 1/2048, 1/256 | the same | the same |
| Connection after the scan, level 1 / 2, characters | 148126 / 65846 | 157495 / 61263 | |
| Exact removal before the numerical phase, seconds | | 22.1 | |

The two keep runs are one computation and differ by 11 % (the load of the
machine fell during the series); the remove run lies between them. Removing
the three moving poles changes no chart count, no halving and no endpoint
overlap, and the value agrees to 2.0e-35 relative (its radius is 5e-35). The
pole of level 2 at x = -0.068 is not what limits that level: its overlap
stays at 1/2048 without it. Presumably the factors in x alone do (one of
them, `112*x^2-120*x+9`, vanishes at x = 0.990). `exchange` was not run on
this family, because it prepares levels 2 to 6 again.

## The epsilon windows have another cause

The structural losses of a level (`--profile 1` of the report tool: the
lowest epsilon order of the two endpoint operators and of the forcing, the
loss of the anchor gauge, and the child window they add up to):

| Level | Lower / upper endpoint low | Forcing low | Anchor loss | Child loss |
| --- | --- | ---: | ---: | ---: |
| double box 1 | -2 / -1 | 0 | 2 | 3 |
| double box 2, with its moving pole | -2 / -3 | -1 | 0 | 2 |
| double box 2, moving pole removed | -2 / -3 | -1 | 0 | 2 |
| double box 3 | 0 / -1 | 0 | 1 | 2 |
| double box 4 | -2 / -1 | -1 | 0 | 2 |
| pentagon 2, with its moving pole | -1 / -1 | 0 | 2 | 3 |
| pentagon 2, moving pole removed | -1 / -1 | 0 | 2 | 3 |

Levels without a moving pole lose as much as the two that have one, and
removing the pole changes no entry, also for the pentagon factor that merges
with the endpoint. The kite (no mixed factor, no epsilon in any denominator)
still asks its levels for epsilon through 1, 3, 4, the banana through 3, 6,
8. What a level loses is set by the poles in epsilon of its endpoint
operators, of its observable rows and of its anchor gauge, i.e. by how the
child's masters at the anchor combine into the level's outputs.

## What is not done

- **Narrower windows.** They need a change of what a level asks of its
  child, so that the combinations which enter with 1/eps^k are requested as
  such instead of being rebuilt from separately truncated masters. The
  `exchange_masters` hook of `recursion::Options` is the place where a level's
  public masters can change; a criterion in terms of epsilon poles is the
  missing part. An internal basis cannot do it (the windows are invariant
  under it, see [epsilon-demand-design.md](epsilon-demand-design.md)).
- **Singular points in x alone.** They are many (27 factors at pentagon
  level 3) and they, not the moving poles, set the endpoint overlaps and the
  detours of the paths. `moving_poles::Options::fixed_points` (`--fixed 1` of
  the report tool) applies the same closure to every factor in x that
  vanishes at neither end of the interval. On the double box none of the 18
  such factors closes within eight steps. That is how a genuine singular
  point behaves (the step bound is no proof), and it is what one expects: a
  level's integrals have thresholds in x once the anchors below it are fixed
  numbers. No choice of masters removes those.
- **Certificates.** An exchanged level is the certified closure transformed
  by an exact, verified change of basis; the IBP witnesses of the old masters
  are dropped rather than rebuilt. A level prepared with the exchanged masters
  from the start would keep them.
- **Trial bounds.** Twelve integrals at two positions are tried per factor
  and step. On pentagon level 1 this takes 20 s; an integral without a new
  denominator for its second factor may exist beyond the twelve.
