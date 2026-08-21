# HPC task: patch CROCO source to disable wave-driven mixing while keeping conservative wave effects

**Read this whole file before editing anything.**

## Context

We run CROCO v2.1.1 coupled to WW3 via OASIS (`OW_COUPLING` + `MRL_WCI`) for the
`sa_west_02` domain. We want a new experiment that isolates the **conservative** wave
effects (Stokes drift, vortex force, Bernoulli head, Stokes–Coriolis) from the
**wave-driven surface mixing** (breaking-induced TKE injection and wave-height surface
roughness).

In CROCO these are *not* separately switchable: the surface wave mixing is hard-wired
inside `#ifdef MRL_WCI` blocks in the GLS turbulence closure. Turning off `MRL_WCI` would
also kill all the conservative terms, which is not what we want. Hence this source patch.

There are **two independent edits**, introducing two new CPP keys:

| Key | File | Disables |
|---|---|---|
| `WCI_NO_WAVE_MIXING` | `gls_mixing.F` | surface breaking-TKE mixing and wave-height surface roughness |
| `WCI_NO_WAVE_STREAMING` | `cppdefs_dev.h` | the forced `#define WAVE_STREAMING` (near-bed streaming body force) |

**Both are required.** The second one is easy to miss: `cppdefs_dev.h` force-defines
`WAVE_STREAMING` whenever `OW_COUPLING && MRL_WCI`, and it is included at the *end* of
`cppdefs.h`, so it silently overrides any `#undef WAVE_STREAMING` a config sets. Without
edit 2 the run still carries bottom streaming — and compiles and runs without complaint.

## Target files

```
$HOME/croco-v2.1.1/OCEAN/gls_mixing.F
$HOME/croco-v2.1.1/OCEAN/cppdefs_dev.h
```

Note this is `$HOME/croco-v2.1.1`, **not** `$HOME/code/croco-v2.1.1`. Confirm the path
matches `CROCO_SRC` in
`somisana-croco/configs/sa_west_02/croco_v2.1.1_cpl_ww3/RUN_01/CROCO_IN/myenv_inter.sh`
(it should read `export CROCO_SRC=/home/gfearon/croco-v2.1.1/OCEAN/`). If it does not
match, **stop and report** rather than guessing.

## Before you start

1. Back up the pristine files:
   ```bash
   cd $HOME/croco-v2.1.1/OCEAN
   cp gls_mixing.F  gls_mixing.F.orig
   cp cppdefs_dev.h cppdefs_dev.h.orig
   ```
   If either `.orig` already exists, the patch may already have been applied — check with
   the verification step below before overwriting anything.

2. This source tree is shared by other configs, including the existing RUN_01. Both
   patches are **inert unless the new keys are defined** in a config's `cppdefs.h`, so
   existing runs are unaffected and RUN_01 needs no recompile. Do not define either key in
   the source tree itself.

## Edit 1 — `gls_mixing.F`

There are exactly **three** `# ifdef MRL_WCI` lines to change, all inside the
`! Boundary conditions` section of the GLS solve (around lines 673, 693 and 705 in the
unmodified v2.1.1 file). Locate them by surrounding context, not by line number alone.

Change each of these three occurrences of:

```fortran
# ifdef MRL_WCI
```

to:

```fortran
# if defined MRL_WCI && !defined WCI_NO_WAVE_MIXING
```

### Site 1 — TKE surface boundary condition (breaking TKE injection)

Identify it by this exact block:

```fortran
            ! Boundary conditions
            IF( ig == itke ) THEN
               DO i=imin,imax
                  ! surface
# ifdef MRL_WCI
                  trb_sfc = MAX( tke_min, cmu_fac*
     &                                       (zalp*wepb(i,j))**(2./3.) )
                  flux_top=zalp*wepb(i,j)
# else
                  trb_sfc = MAX( tke_min, cm0inv2*ustar_sfc_sq(i,j) )
                  flux_top = 0.
# endif
```

Only the `# ifdef MRL_WCI` line changes. With the key defined, the `#else` branch runs:
standard law-of-the-wall surface TKE and zero TKE flux — identical to an uncoupled run.

### Site 2 — surface roughness length

```fortran
            ELSE
               DO i=imin,imax
                  ! surface
# ifdef MRL_WCI
                  z0_s=MAX( Zosmin , 0.5*whrm(i,j) )
# else
                  z0_s=MAX( Zosmin , chk*ustar_sfc_sq(i,j) )   !<-- Charnock
# endif
```

With the key defined, `z0_s` reverts to the Charnock relation instead of scaling with
significant wave height.

### Site 3 — extra `psi`/`gls` surface flux term

```fortran
                  flux_top = -rn*cm0**(rp+1.)
     &                             *vonKar*OneOverSig(igls)
     &                             *(cff**(rm+0.5))*(lgthsc**rn)
# ifdef MRL_WCI
                  flux_top=flux_top - OneOverSig(itke)*OneOverSig(igls)
     &                                  *   cm0**(rp+1.)
     &                                  *   rm
     &                                  *   cff**(rm-0.5)
     &                                  *   lgthsc** (rn+1.)
     &                                  *   zalp*   wepb(i,j)
# endif
```

This block has no `#else`; with the key defined the wave term simply drops out of
`flux_top`.

## Do NOT change

- The variable **declarations** at roughly lines 140 and 165:
  ```fortran
  # ifdef MRL_WCI
        REAL            ::  cmu_fac, cb_wallE, gls_sigp_cb
  # endif
  ```
  and
  ```fortran
  # ifdef MRL_WCI
        REAL, PARAMETER ::  zalp    =  0.05    ! fraction of wave energy
  # endif
  ```
  Leave both under plain `MRL_WCI` — otherwise the compile fails on undeclared names.

- The `cmu_fac` / `gls_sigp_cb` **initialisation** block (around line 328, also
  `# ifdef MRL_WCI`). `cmu_fac` is only consumed at site 1, so leaving it computed and
  unused is harmless.

- The `# ifdef BBL` block for `z0_b` (bottom roughness, around line 714). Bottom-side
  behaviour is controlled by CPP keys in `cppdefs.h`, not by this patch.

- Any other `MRL_WCI` guard anywhere else in the source tree. `mrl_wci.F` must stay fully
  active — that is where the conservative terms we want to keep are computed.

## Edit 2 — `cppdefs_dev.h`

Around **line 78–83**, inside the `#ifdef OW_COUPLING` block, find:

```c
# ifdef MRL_WCI
#  undef  WAVE_ROLLER
#  define WAVE_STREAMING
#  define WAVE_RAMP
# endif
```

Change it to:

```c
# ifdef MRL_WCI
#  undef  WAVE_ROLLER
#  ifndef WCI_NO_WAVE_STREAMING
#   define WAVE_STREAMING
#  endif
#  define WAVE_RAMP
# endif
```

Only the `#  define WAVE_STREAMING` line is wrapped. **Do not** touch `WAVE_ROLLER` or
`WAVE_RAMP`, and do not alter the surrounding `#ifdef OW_COUPLING` block (it sets
`MPI_COMM_WORLD` and other coupling essentials).

Note the indentation: the new `#  ifndef` is at two spaces, and the `#   define` it wraps
moves to three.

## Verification

1. Confirm the `gls_mixing.F` changes are exactly three single-line hunks:
   ```bash
   cd $HOME/croco-v2.1.1/OCEAN
   diff gls_mixing.F.orig gls_mixing.F
   ```
   Expect three hunks, each changing `# ifdef MRL_WCI` to
   `# if defined MRL_WCI && !defined WCI_NO_WAVE_MIXING`.

2. Confirm the `cppdefs_dev.h` change is the single expected hunk:
   ```bash
   diff cppdefs_dev.h.orig cppdefs_dev.h
   ```

3. Confirm the key counts:
   ```bash
   grep -c "WCI_NO_WAVE_MIXING"    gls_mixing.F    # must be 3
   grep -c "WCI_NO_WAVE_STREAMING" cppdefs_dev.h   # must be 1
   ```

4. Confirm no stray edits to wave variables:
   ```bash
   grep -n "wepb\|whrm\|zalp\|cmu_fac" gls_mixing.F
   ```

5. **Prove RUN_01 is unaffected.** RUN_01's `cppdefs.h` defines neither new key, so both
   guards must evaluate exactly as the originals did. Verify by reverting the
   `gls_mixing.F` guard in a pipe and checking it reproduces the backup byte-for-byte —
   this should print nothing:
   ```bash
   cd $HOME/croco-v2.1.1/OCEAN
   sed 's/^# if defined MRL_WCI && !defined WCI_NO_WAVE_MIXING$/# ifdef MRL_WCI/' \
       gls_mixing.F | diff - gls_mixing.F.orig && echo "gls_mixing.F: guards only"
   ```
   For `cppdefs_dev.h` the diff from step 2 is the proof: it must show one hunk that only
   adds the `#  ifndef` / `#  endif` pair around the existing `define`, with no other line
   added, removed or reordered.

6. Report back:
   - both diffs,
   - the resolved `CROCO_SRC` path you patched,
   - confirmation that both `.orig` backups exist.

## What happens next (not your job unless asked)

The config side is already done in the `somisana-croco` repo: `RUN_02/CROCO_IN/COMP/cppdefs.h`
defines both `WCI_NO_WAVE_MIXING` and `WCI_NO_WAVE_STREAMING`, and turns off `BBL`,
`WAVE_STREAMING` and `WAVE_FRICTION`. See `RUN_02/README.md`. **Do not edit any
`cppdefs.h` as part of this task** — only the two source files above.

Once both patches are in, RUN_02 can be compiled with
`configs/sa_west_02/croco_v2.1.1_cpl_ww3/RUN_02/CROCO_IN/jobcomp_lengau`.

Expected channel state after patching (verified locally with `cpp` against these exact
`cppdefs.h` files):

| Channel | RUN_01 | RUN_02 |
|---|---|---|
| Conservative (`MRL_WCI`: Stokes, vortex force, BHD) | ON | **ON** |
| Surface breaking-TKE mixing | ON | **OFF** |
| Bottom streaming (`WAVE_STREAMING`) | ON | **OFF** |
| Wave bottom friction (`BBL`, `WAVE_FRICTION`) | ON | **OFF** |
