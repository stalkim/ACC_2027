# Geometric boundary and informative-history checks

These checks supplement the polynomial inverse-cost benchmark. They require
only Python 3.10 or newer and the existing `monotone_history.py` and
`extension_history.py` in this repository. No dependency installation is needed.
The existing benchmark, decoder and README are unchanged.

From the repository root:

```bash
python3 -m unittest discover -p 'test_*.py'
python3 review_checks.py --output results/review-verification-001
```

Use `python` instead of `python3` if that is your Windows interpreter name.
Each execution requires a new output directory and refuses to overwrite one.
There is no training, random search or external service call.

## Fixed comparisons

- Straight lines: Euclidean radius 0.01, arc-speed band [1, 1.2], rational
  unit directions (+h,k), (-h,k), with (h,k) = (5/13,12/13), (3/5,4/5),
  (4/5,3/5), (1,0). Set alpha = T*h/0.02 to 0.8, 0.99, 1, 1.01, 1.2;
  use 2, 3, 9 and 33 samples. A separate line/ball preimage solver uses
  rational square-root enclosures, not directional pruning.
- Informative histories: quartic path, lambda=1/2, coordinate error 0.012,
  endpoint times +/-0.005, interior times -0.003, -0.0025, -0.002 and
  parameter shifts -0.0002, 0, 0.0002. Endpoints are shared midpoints;
  the interior positive position has x-error -0.0108. Endpoints must be
  ambiguous, the full record uniquely positive, with zero pruning.
- Curvature: symmetric circular arcs, tangent (3/5,4/5), curvature
  0.5, 1, 2, anchor radius 0, 0.05, 0.1, arc charts [-0.3,0.3]. The grid
  T=j/200, j=1,...,50 is checked against chart coverage. Rational Taylor
  enclosures compare the general lower bound with the exact minimum
  projected chord; this is not an exact curved-history separation threshold.
- Cost: retain the informative samples and densify to 3, 9 and 33 samples
  at fixed window length, with 18, 24 and 32 bisections. Compare DC/MI/DMI
  coordinate counts. Measure MI/DMI on the nine informative records in
  three rotating-order repetitions; there is no expected pruning speedup.

## Reference outcomes

The reference executions used Python 3.10.12 on Linux x86-64
(kernel 5.4.210-39.1, glibc 2.35), an Intel Xeon Platinum 8468 CPU,
one Python process and standard-library rational arithmetic. No third-party
numerical libraries, GPU or shared decoder caches were used. Batch timings
exclude setup and source/input output. The additional checks completed in
7.62 seconds including that output; this is not the timed batch duration.
Actual times vary by machine and Python version. Each new run records its
own Python version and platform in `manifest.json`.

The 60-history benchmark uses five rotating-order timing repetitions;
the nine informative histories use three. The manuscript reports complete-batch
medians, not worst-case per-call latency. Its reference medians are:

| Batch | DC (s) | MI (s) | DMI (s) |
|---|---:|---:|---:|
| 60 polynomial histories | 1.539531 | 1.502498 | 0.768793 |
| Nine informative histories | Not timed | 0.188216 | 0.188633 |

The public baseline used for the polynomial comparison is commit
`a193dcd7d104c82f0e5321e70755f2fcef638547`. The article links to the repository
root; source hashes and run manifests identify the exact implementation.

| Check | Reference outcome |
|---|---|
| 80 midpoint evaluations | 48 ambiguous at/below alpha=1; 32 incompatible above |
| 640 valid-history evaluations | 628 unique, 12 ambiguous; all 256 above alpha=1 correctly unique |
| 80 out-of-plane controls | All incompatible |
| Nine informative family members | Endpoints ambiguous; full history uniquely positive; no pruning |
| 369 chart-covered circular grid points | General bound sufficient at 296; exact projected-chord bound at 305 |
| n=3,9,33; B=18,24,32 | All unique with no pruning; respectively 168n,216n,280n coordinate calls |
| Nine-record timings | MI median 0.188216 s; DMI median 0.188633 s |

Counts are evaluations, not independent trajectories or random trials.
For two samples the constant/alternating schedule settings coincide.
The exact decimal example is checked separately and is also one of the
nine family members; it is not a tenth independent informative case.
For alpha>1 the symmetric midpoint record is incompatible, not uniquely
identified: only model-valid inputs inherit the uniqueness guarantee.
Below the threshold some records are still unique.

Every emitted witness is checked, and each declared valid current parameter
must remain in its branch's outer. The public polynomial decoder is unchanged.
No speedup is claimed over COMMA or an optimized external implementation.

## Saved evidence

The run saves source snapshots and hashes, input data, separate truth,
full decoder results, curvature enclosures, per-call timings and a summary.
Truth never enters a decoder. Reference `summary.json` SHA-256:
`0db661f932c864a0882d175b064e41a468b3e366ccd98002b12ff882bf04e9ec`.
New runs will have different timing values and therefore different summary hashes.
