# ACC 2027: Certified Branch Recovery

Research code for **Sharp Distinguishability Bounds and Certified Branch
Recovery on Self-Intersecting Paths**, by Stanislav Kim and Anton Pyrkin.

The code checks which directed branches can explain a finite position history
when the speed profile is unknown and position errors have known bounds.
It includes a straight-line/Euclidean-ball solver, a cubic-path decoder,
quartic/quintic comparisons, and fixed checks of the geometric threshold,
curvature bound, informative intermediate measurements, and computational cost.

## Quick start

Use Python 3.10 or later. No third-party packages are needed. Run commands
from the repository root, without `-O`: validation uses assertions.

```bash
python3 --version
python3 demo.py
python3 -m unittest discover -s . -p 'test_*.py'
```

The demo returns `UNIQUE`, `AMBIGUOUS`, and `INCOMPATIBLE` for three short
cubic-path histories. On Windows, use `python` instead of `python3` if needed.
The polynomial benchmark requires Linux, macOS, or WSL for per-call timeouts.

## Reproduce the numerical evidence

Every run requires a new output directory and refuses to overwrite an existing
one. Generated files under `results/` are ignored by Git. These are fixed,
deterministic checks; there is no training, random search, or external service
call. Ground truth is stored separately and used for validation, not supplied
to a decoder.

### Geometric boundary and informative histories

```bash
python3 review_checks.py --output results/geometric-01 --protocol README.md
```

This command covers the straight-branch threshold, the curved sufficient rule,
the intermediate-measurement example, and the no-pruning cost comparison.
`--protocol README.md` archives this description with the source snapshot.

| Check | Fixed setup |
| --- | --- |
| Straight branches | Euclidean error radius `b=0.01`, arc-speed band `[1, 1.2]`, unit directions `(+h,k,0)` and `(-h,k,0)` with `(h,k)` in `(5/13,12/13)`, `(3/5,4/5)`, `(4/5,3/5)`, `(1,0)`. Vary `alpha = T*h/(2*b)` over `0.8, 0.99, 1, 1.01, 1.2`, using 2, 3, 9, and 33 samples. |
| Informative histories | Quartic path, scale `1/2`, coordinate-error radius `0.012`, endpoint times `-0.005, 0.005`, interior times `-0.003, -0.0025, -0.002`, and parameter shifts `-0.0002, 0, 0.0002`. Endpoint observations are shared midpoints; the interior positive-branch position has an x-error of `-0.0108`. |
| Curvature | Symmetric circular arcs with positive-branch crossing tangent `(3/5,4/5)`, curvatures `0.5, 1, 2`, anchor radii `0, 0.05, 0.1`, arc charts `[-0.3,0.3]`, and `T=j/200` for `j=1,...,50`. |
| Cost without pruning | Retain the informative observations and densify to `n=3,9,33` samples at fixed window length, with `B=18,24,32` bisections. Compare DC/MI/DMI coordinate counts and time MI/DMI on the nine informative records in three alternating-order repetitions. |

The straight solver uses rational square-root enclosures for line/ball
preimages, with 40-bit bounds and exact handling of perfect squares. It does
not use directional pruning. The curvature calculation uses rational Taylor
enclosures to compare the general lower bound with the exact minimum projected
chord; it does not compute the exact curved-history distinguishability boundary.

Expected outcomes:

| Check | Result |
| --- | --- |
| 80 midpoint evaluations | 48 ambiguous at/below `alpha=1`; 32 incompatible above |
| 640 model-valid evaluations | 628 unique, 12 ambiguous; all 256 above `alpha=1` correctly unique |
| 80 out-of-plane controls | All incompatible |
| Nine informative family members | Endpoints ambiguous; full history uniquely positive; no pruning |
| 450 circular grid points | 369 satisfy strict chart coverage; the general bound certifies 296, the exact projected-chord bound 305 |
| `n=3,9,33`; `B=18,24,32` | All unique with no pruning; respectively `168n`, `216n`, `280n` coordinate evaluations for each method |

Counts are evaluations, not independent trajectories or random trials. With
two samples the constant/alternating speed schedules coincide. The exact
decimal example is checked separately and also belongs to the nine-member
family: the output reports 10 informative records, not 10 distinct cases.
For `alpha>1`, the symmetric midpoint record is incompatible with both branches;
only model-valid inputs inherit the uniqueness guarantee. Below the threshold,
some records are still uniquely identifiable.

The run saves `source/`, `manifest.json`, `inputs.json`, `truth.json`,
`results.jsonl`, `curvature.json`, `timings.json`, and `summary.json`.
These contain source snapshots and hashes, decoder results and witnesses,
current-parameter bounds, curvature enclosures, and timing samples.

### Polynomial comparison and precision check

DC checks pairwise temporal difference constraints. MI propagates feasible
parameter intervals. DMI adds directional exclusion before interval propagation.
All three share the coordinate inverse routine. Their agreement checks temporal
feasibility; it is not an independent validation of that shared inverse routine.
Emitted witnesses are also checked in the original equations.

First check all 60 histories once, without timing repetitions:

```bash
python3 benchmark.py --check-only --output results/polynomial-check-01
python3 check_precision.py results/polynomial-check-01
```

For the timing comparison, run one correctness pass followed by five timed
passes. Method order rotates, and methods do not share computed caches.

```bash
python3 benchmark.py --output results/polynomial-timing-01
```

The limits are 5 seconds per call and 180 seconds for the suite. A successful
full run writes 1080 call records; a check-only run writes 180. Partial files
after an exception are not a completed run.

The fixed examples use quartic and quintic paths with scales `1/2, 1, 2`,
nominal durations `0.02, 0.08`, and coordinate-error bounds `0.0005, 0.01`.
There are 48 five-sample nominal histories, six ambiguous histories, and six
incompatible histories. All three methods should agree on these decisions.
DC and MI use 57024 coordinate evaluations per batch; DMI uses 28512.

All 48 nominal histories already satisfy the endpoint directional condition
(minimum margin `0.007`). This batch measures inverse-computation savings.
The nine informative histories above demonstrate the additional information
from intermediate samples; they allow no pruning and give equal MI/DMI
coordinate counts.

The benchmark uses 24 inverse bisections. `check_precision.py` verifies the
derivative, error, and chart-margin premises for all 48 nominal histories:
the parameter margin is `mu=e/230`, so 18 bisections already suffice for
their inner witnesses. This check does not tune the decoder's precision.

The two-zero polynomial controls use timestamps `-0.005, 0.005` and error
bounds `0.06, 0.0005`. The cubic demo uses `-0.01, 0.01` for its two-zero
controls; these are different examples.

Outputs are `inputs.json`, `truth.json`, `calls.jsonl`, `manifest.json`, and
`summary.json`. They include source hashes, checked decisions and witnesses,
parameter bounds, coordinate counts, and elapsed times.

### Reference timings and reproducibility

The manuscript reports complete-batch medians, not worst-case per-call latency:

| Batch | DC (s) | MI (s) | DMI (s) |
| --- | ---: | ---: | ---: |
| 60 polynomial histories, five timed repetitions | 1.539531 | 1.502498 | 0.768793 |
| Nine informative histories, three timed repetitions | Not timed | 0.188216 | 0.188633 |

These reference runs used Python 3.10.12 on Linux x86-64 (kernel 5.4.210-39.1,
glibc 2.35), an Intel Xeon Platinum 8468 CPU, one Python process, and
standard-library rational arithmetic. Batch timings exclude setup and file
output. Elapsed times depend on the machine and Python version; new runs
record their own environment in `manifest.json`.

The geometric and informative-history implementation was introduced in
[`dd8dfa0`](https://github.com/stalkim/ACC_2027/commit/dd8dfa04e7c3a222059232116658a802f20a2a62).
The polynomial reference comparison used
[`a193dcd`](https://github.com/stalkim/ACC_2027/commit/a193dcd7d104c82f0e5321e70755f2fcef638547).
Source hashes in each run identify the implementation actually executed.
The reference geometric `summary.json` SHA-256 is
`0db661f932c864a0882d175b064e41a468b3e366ccd98002b12ff882bf04e9ec`.
New summaries have different timing values and therefore different hashes.

## Models and decisions

The two models use different error sets and progress variables:

| Model | Position-error bounds | Progress bounds |
| --- | --- | --- |
| Straight lines and circular-arc checks | Euclidean radius `b=0.01` | Arc-length rate in `[1, 1.2]` |
| Cubic, quartic, and quintic paths | A common bound for all three coordinates at each timestamp; it may vary between timestamps | Parameter rate in `[9/10, 11/10]` |

The polynomial charts are `[-5/4, -3/4]` and `[3/4, 5/4]`. Their histories
must satisfy the decoder's coverage checks. No separate arc-speed bound is
imposed in this model: changing the curve's shape at fixed parameter rate
need not preserve physical speed. The error model imposes no additional
constraints across samples or coordinates; statistical independence is not
assumed, and correlations within the bounds are not used.

Inputs accept rational strings such as `"1/50"` and `"0.0005"`. Prefer them to
binary floating-point values when specifying exact decimal bounds.

| Decision | Meaning |
| --- | --- |
| `UNIQUE` | One branch has a checked feasible history and the other is excluded |
| `AMBIGUOUS` | Both branches have checked feasible histories |
| `INCOMPATIBLE` | Neither branch can explain the data under the model |
| `UNRESOLVED` | Finite-precision bounds did not settle the decision |
| `OUTSIDE_COVERAGE` | The polynomial decoder's required chart-coverage checks did not pass |

A nonempty outer interval alone does not prove that a feasible history exists.
Every emitted witness is checked, and valid current parameters are checked
for inclusion in their branch's outer bound.

## Code map

| File | Purpose |
| --- | --- |
| [straight_history.py](straight_history.py) | Line/Euclidean-ball preimages, rational square-root bounds, and history certificates |
| [monotone_history.py](monotone_history.py) | Rational inverse bounds, interval propagation, and cubic-path decoder |
| [extension_history.py](extension_history.py) | Quartic/quintic paths and the DC, MI, DMI methods |
| [review_checks.py](review_checks.py) | Straight boundary, circular-arc bounds, informative histories, scaling, and no-pruning timings |
| [benchmark.py](benchmark.py) | Deterministic 60-history polynomial generation, validation, and timing |
| [check_precision.py](check_precision.py) | Positive-margin precision check on generated nominal histories |
| [demo.py](demo.py) | Three small cubic-path examples |
| `test_*.py` | Interval, root-enclosure, decision, temporal-feasibility, and generator tests |

## Scope

This repository contains the numerical checks described above. It does not
include the manuscript, archived closest-point/continuation and branch-and-bound
experiments, or a closed-loop controller. The code assumes known geometry,
valid error/rate bounds, positive progress, and covered polynomial charts with
connected measurement preimages. It does not establish arbitrary-curve,
large-map, unknown-geometry, global-lap, or hard-real-time performance.

The temporal constraint formulation follows R. Dechter, I. Meiri, and J. Pearl,
["Temporal constraint networks"](https://www.sciencedirect.com/science/article/pii/0004370291900066),
*Artificial Intelligence*, 49 (1991), 61-95. Interval reachability and
set-membership estimation are established methods; this repository implements
specializations for the stated path models. DC is not a COMMA reimplementation,
and no speedup over COMMA or an optimized external implementation is claimed.
