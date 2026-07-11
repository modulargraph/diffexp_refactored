# Retained parity-oracle logs

These two logs are historical DiffExp 1 Feynman-trick ladder records retained
as transparent parity targets for the curated values in `Docs/Results.md`.
They are not independent pySecDec validations and they are not timing claims
for DiffExp 2.

| File | Scope | Recorded finite coefficient |
| --- | --- | ---: |
| `l2_bubsun.log` | equal-mass bubble and sunrise, `D=2-2 eps`, `p^2=-1`, working precision 200 | bubble `0.8608178819280080777765623653643473221`; sunrise `2.2367927002126465108229117827758723767` |
| `l2_banana.log` | equal-mass three-loop banana, `D=2-2 eps`, `p^2=-1` | `8.2681045358689687593219901952256396` |

The raw logs predate the release facade and therefore contain old loader
banners and development diagnostics. They are shipped only as immutable
provenance for the displayed numbers; current workflows use `DiffExp2.m`,
`FeynmanTrick.m`, and `Scripts/run_ft_stepwise2.m`.
