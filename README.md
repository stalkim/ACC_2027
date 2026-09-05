Python code for branch recovery from finite position histories on
self-intersecting polynomial curves. The speed profile is unknown; its
parameter rate and the coordinate measurement errors have known bounds.

The repository contains the cubic-path decoder and the quartic/quintic
comparison used in Experiment 3. It does not include the manuscript, the
earlier full set-membership baselines, or the complete Experiments 1/2 pipeline.

## Quick start

Python 3.10 or later. No third-party packages are needed.

```bash
python3 demo.py
python3 -m unittest discover -s . -p 'test_*.py'
```

The demo returns `UNIQUE`, `AMBIGUOUS`, and `INCOMPATIBLE` for three short
histories. Use normal Python execution, without `-O`: the checks use assertions.
On Windows, `python` can be used instead of `python3`.

## Comparing the methods

The comparison script needs Linux, macOS, or WSL for per-call timeouts.
First check all 60 histories once, without timing repetitions:

```bash
python3 benchmark.py --check-only --output results/check-01
python3 check_precision.py results/check-01
```

For the timing comparison, run one correctness pass followed by five timed
passes. Method order rotates, and the methods do not share computed caches.

```bash
python3 benchmark.py --output results/benchmark-01
```

The limits are 5 seconds per call and 180 seconds for the suite. A successful
full run writes 1080 call records. The check-only run writes 180. Each output
directory must be new; an existing directory is never overwritten.

Outputs are `inputs.json`, `truth.json`, `calls.jsonl`, `manifest.json`, and
`summary.json`. They include source hashes, decisions, feasible histories,
parameter bounds, coordinate evaluation counts, and elapsed times. Ground
truth is used by the checker, not passed to the decoder. Generated outputs
are ignored by Git. Partial files after an exception are not a completed run.

The fixed examples contain 48 unique-branch histories, six ambiguous histories,
and six incompatible histories. All three methods should agree on these
decisions. DC and MI use 57024 coordinate evaluations per batch; DMI uses 28512.
Elapsed times depend on the machine and are not expected to match exactly.

## Code

| File | Purpose |
| --- | --- |
| `monotone_history.py` | Rational inverse bounds, interval propagation, and cubic-path decoder |
| `extension_history.py` | Quartic/quintic paths and the DC, MI, DMI methods |
| `benchmark.py` | Deterministic example generation, validation, and timing |
| `check_precision.py` | Check the positive-margin precision bound on the generated nominal examples |
| `demo.py` | Three small cubic-path examples |
| `test_*.py` | Tests for the interval operations, decisions, and example generator |

DC checks pairwise temporal difference constraints. MI propagates feasible
parameter intervals. DMI adds directional exclusion before interval propagation.
They share the coordinate inverse routine; DC is not an external implementation
of a map-matching package.

## Model and limits

The charts are `[-5/4, -3/4]` and `[3/4, 5/4]`; the parameter rate is in
`[9/10, 11/10]`. Measurements have three coordinates and an independent common
coordinate-error bound at each timestamp. The histories must satisfy the
coverage checks in the decoder. The code is for these polynomial families,
not arbitrary curves, flight dynamics, or closed-loop control.

Inputs accept rational strings such as `"1/50"` and `"0.0005"`. Prefer them to
binary floating-point values when specifying exact decimal bounds. The
benchmark uses 24 inverse bisections. It does not tune precision per case.

`UNIQUE` requires a feasible history on one branch and exclusion of the other.
`AMBIGUOUS` means both branches have checked feasible histories.
`INCOMPATIBLE` means neither branch can explain the data under the model.
`UNRESOLVED` means the finite-precision bounds did not settle the decision;
`OUTSIDE_COVERAGE` means the required chart-coverage conditions did not pass.
A nonempty outer interval alone does not prove a feasible history exists.

The quartic/quintic benchmark uses scales `1/2, 1, 2`, nominal durations
`0.02, 0.08`, and error bounds `0.0005, 0.01`. Its two-zero controls use
timestamps `-0.005, 0.005` and bounds `0.06, 0.0005`. The cubic demo uses
`-0.01, 0.01` for its two-zero controls; these are different examples.

The temporal constraint formulation follows R. Dechter, I. Meiri, and J. Pearl,
["Temporal constraint networks"](https://www.sciencedirect.com/science/article/pii/0004370291900066),
*Artificial Intelligence*, 49 (1991), 61–95.
Interval reachability and set-membership estimation are established methods;
this repository implements their polynomial-path specialization.
