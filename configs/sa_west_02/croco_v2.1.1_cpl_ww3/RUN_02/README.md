# RUN_02 — conservative wave effects only

CROCO↔WW3 coupled run for `sa_west_02`, identical to `RUN_01` except that **every
non-conservative wave effect is switched off**. The purpose is to isolate the conservative
wave forcing (Stokes drift, vortex force, Stokes–Coriolis, Bernoulli head) from the
wave-driven mixing and bottom-boundary-layer effects that are bundled with it in RUN_01.

## The comparison this run supports

| Run | Wave channels active |
|---|---|
| `croco_v2.1.1/C10_I01_GLORYS_ERA5` (no-wave) | none |
| **RUN_02** (this run) | conservative only (**A**) |
| `RUN_01` | A + surface mixing (B) + bottom streaming (C) + wave bottom friction (D) |

**RUN_02 − no-wave = a pure conservative-wave signal.** Bottom friction in RUN_02 is the
plain log-law with `Zob=1e-3`, physically identical to the no-wave run, so nothing else
contaminates the difference.

Note the limitation: **`RUN_01 − RUN_02` is a lump of B+C+D** and cannot be decomposed
further with these three runs. Separating B would need a fourth run with `MRL_WCI` mixing
intact but the bottom keys off.

## How it differs from RUN_01

Only `CROCO_IN/COMP/cppdefs.h` differs. Everything else — grid, forcing, boundary
conditions, OASIS `namcouple`, all WW3 namelists, `croco_inter.in`, MPI decomposition,
time step — is byte-identical to RUN_01.

```
#  define WCI_NO_WAVE_MIXING     (new; user block, ~line 93)
#  define WCI_NO_WAVE_STREAMING  (new; user block, ~line 93)
#  undef  BBL                    (was define)
#  undef  WAVE_STREAMING         (was define)
#  undef  WAVE_FRICTION          (was define)
```

Why each one:

- **`WCI_NO_WAVE_MIXING`** — a *new* CPP key requiring the patched `gls_mixing.F` (see
  `../HPC_SOURCE_PATCH.md`). CROCO hard-wires the surface breaking-TKE injection and the
  wave-height surface roughness inside `#ifdef MRL_WCI`, so they cannot be switched off
  from `cppdefs.h` alone without also killing the conservative terms. The patch guards
  those three blocks with `#if defined MRL_WCI && !defined WCI_NO_WAVE_MIXING`, so they
  fall back to the standard non-wave branches (`flux_top = 0`, Charnock roughness).
- **`WCI_NO_WAVE_STREAMING`** — a second *new* key, requiring a patched `cppdefs_dev.h`.
  **The `#undef WAVE_STREAMING` below is not sufficient on its own**: `cppdefs_dev.h:81`
  force-defines `WAVE_STREAMING` whenever `OW_COUPLING && MRL_WCI`, and it is included at
  the *end* of `cppdefs.h` (via `set_global_definitions.h`), so it overrides the undef.
  Without this key the run would still carry the bottom-streaming body force
  (`rhs3d.F:1807`, `step2d.F:1416`, `mrl_wci.F:1137`) — and would compile and run without
  complaint.
- **`BBL`** off — removes the wave–current bottom stress, and with it any dependence on
  the patched `bbl.F` (the zero-bottom-friction bug fix).
- **`WAVE_STREAMING` + `WAVE_FRICTION`** off — removes near-bed streaming and
  wave-enhanced roughness. With both off, the `#if (defined WAVE_FRICTION && !defined
  WKB_WWAVE) || defined WAVE_STREAMING` block at `mrl_wci.F:156` is not compiled at all,
  so the `#define Zob Zob_my` override never happens and `get_vbc` uses the same
  `Zob=1e-3` log-law as the no-wave run.

`MRL_WCI` itself stays **on** — it carries all the conservative machinery we want to keep.

## Source-code prerequisite

This run **will not behave correctly against an unpatched source tree**. Both files in
`$CROCO_SRC` (`/home/gfearon/croco-v2.1.1/OCEAN/`, set in `CROCO_IN/myenv_inter.sh`) must
carry the patches described in `../HPC_SOURCE_PATCH.md`.

If either patch is missing the corresponding key is **silently ignored** — the run
compiles and executes cleanly, it just quietly retains the physics you meant to remove.
So check before compiling:

```bash
grep -c "WCI_NO_WAVE_MIXING"    $CROCO_SRC/gls_mixing.F    # must be 3
grep -c "WCI_NO_WAVE_STREAMING" $CROCO_SRC/cppdefs_dev.h   # must be 1
```

Both patches are inert without their keys, so RUN_01 is unaffected and needs no recompile.

## Verified channel state

Checked locally by running `cpp` over these exact `cppdefs.h` files against a patched
source tree:

| Channel | RUN_01 | RUN_02 |
|---|---|---|
| Conservative — Stokes, vortex force, Stokes–Coriolis, BHD | ON | **ON** |
| Surface breaking-TKE mixing | ON | **OFF** |
| Bottom streaming | ON | **OFF** |
| Wave bottom friction / roughness | ON | **OFF** |
| `Zob → Zob_my` override block compiled | ON | **OFF** |

## Coupling fields

`namcouple` is unchanged from RUN_01 — all 13 fields are still exchanged. `FOC` (the
wave-to-ocean energy flux, `wepb`) is still received and still written to the history and
average files by `wrt_his`/`wrt_avg`; it simply no longer forces anything. That makes it a
useful diagnostic for confirming the run really did ignore it. `UBRX`/`UBRY` are likewise
received but unused with the bottom keys off.

## Running it

```bash
# 1. compile (script is now path-relative; run it from this run's CROCO_IN)
cd $CONFIG_DIR/RUN_02/CROCO_IN && ./jobcomp_lengau

# 2. submit the chained monthly jobs
cd $CONFIG_DIR && RUN_NAME=RUN_02 ./run_inter_chpc.bash
```

`run_inter_chpc.bash` defaults to `RUN_NAME=RUN_02`; `run_inter.pbs` falls back to
`RUN_01` if submitted by hand without the variable. Per-run `stdout`/`stderr` and the PBS
job name are set by `qsub -o/-e/-N` in `run_inter_chpc.bash`, because the `#PBS` directives
in `run_inter.pbs` are static and cannot interpolate `$RUN_NAME`.

Period and restart behaviour are unchanged from RUN_01: 2008-01-01 to 2013-12-31, cold
start (`RSTFLAG=0`) from the same GLORYS initial condition, in 4-month PBS chunks.

## Files not in git

`OASIS_IN/oce.nc` and `OASIS_IN/wav.nc` (the cold-start OASIS coupling restarts) are
matched by the repo's `*.nc*` ignore rule, so **they will not reach the HPC via git**.
They are byte-identical to RUN_01's, so either rsync them across with the rest of the run
directory or copy them on the HPC:

```bash
cp $CONFIG_DIR/RUN_01/OASIS_IN/{oce.nc,wav.nc} $CONFIG_DIR/RUN_02/OASIS_IN/
```

`CROCO_IN/COMP/` likewise only tracks `cppdefs.h` and `param_.h`; `param.h`, `Compile/`
and the `croco` executable are all regenerated by `jobcomp_lengau`.
