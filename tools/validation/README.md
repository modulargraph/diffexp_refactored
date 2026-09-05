# Completed transport basis conversion

`rebase_completed_transport.cpp` validates the reuse of an ordinary completed
transport after adding explicit algebraic basis prefactors. Build it with the
same compiler, include and library flags as `diffexp`, then run:

```sh
rebase_completed_transport continuous-request.json prefactor-request.json continuous-response.json principal-response.json
```

The utility requires identical requests after removing `basis_prefactors` from
the second request. It recompiles the exact path, continues its roots, and calls
the production endpoint conversion helper on the saved value and error balls.
It does not integrate the differential equation again or read reference values.
The response preserves the original native timings and records conversion time
separately in `basis_rebase.seconds`.

This was used for the [128-digit original planar comparison](../../docs/validation/zzz-high-basis.json).
That report records input/output hashes, published basis data, the initial
convention mismatch and the independent endpoint comparison. To reproduce the
whole calculation with endpoint conversion directly in the normal native
process, use the `zzz-high` case in `scripts/check_original_mathematica.py`.
