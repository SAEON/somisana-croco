Here is where we prepare all the files to be used in our CROCO configuration and where we run the model

This configuration is the `sa_eez_01` domain run with the **croco-v2.1.1** source code. It is a direct
port of the `croco_v2.0.1` C12 setup in the sibling directory - same grid, same CPP keys, same runtime
options - with the compile files rebuilt on top of the v2.1.1 source templates.

For inter-annual runs, there is a `jobcomp_lengau` file for compiling the croco code and a
`run_croco_inter.pbs` file for running the simulation. Both of these scripts use environment variables
which are set up in `myenv_inter.sh`

None of the netcdf files are copied to the remote repo, only the scripts used to generate them
(the grid file is the exception - it is small enough and static enough to keep under version control).

<pre>
Grid input
----------
`GRID` - this is not intended to be a configurable directory. i.e. if you want a new grid, create a new
domain e.g. sa\_eez\_02. Copied unchanged from `croco_v2.0.1/GRID`, so the grid is byte-identical between
the two source versions. `AGRIF_FixedGrids.in` is kept because `run_croco_inter.pbs` copies it
unconditionally, even though AGRIF is not used here.

Compile options
---------------
`C**`
Numbering is inherited from `croco_v2.0.1` so the two directories can be cross-referenced. Only C12 has
been carried over.

12 - as per the croco\_v2.0.1 C12: ONLINE ERA5 bulk forcing (ERA\_ECMWF, BULK\_MONTH\_1DIGIT), READ\_PATM,
     TIDES, CLIMATOLOGY (FRC\_BRY undefined), UV\_VIS2 + UV\_VIS\_SMAGO, VADV\_ADAPT\_IMP, GLS\_MIXING,
     no TIDERAMP, open boundaries E/W/S (OBC\_NORTH undefined).

`cppdefs.h` and `param_.h` were rebuilt from the croco-v2.1.1 templates rather than copied from
croco\_v2.0.1, because the REGIONAL block was reordered between the two releases (Lateral Forcing now
precedes Surface Forcing, ABL1D moved inside BULK\_FLUX, the RVTK\_DEBUG block was removed). The
effective key settings are identical to croco\_v2.0.1/C12 - verified key by key. Differences that remain
are all upstream v2.1.1 changes sitting behind keys we do not define (OA\_COUPLING block,
DIAGNOSTICS\_TS\_MLD\_CRIT, key\_pisces\_light renamed key\_pisces\_npzd).

Note there is no C13 here. C13 in croco\_v2.0.1 was C12 plus a trailing `#undef SPONGE_GRID` to take the
sponge width from x\_sponge in the .in file. This configuration uses C12, so SPONGE\_GRID stays defined
(cppdefs\_dev.h force-defines it) and set\_nudgcof.F uses a fixed 10 grid point sponge. The `XXX`
placeholder on the sponge line of the .in file is therefore never read - but it would be a landmine if
anyone switches to a C13-style build, so check it before doing so.

Runtime input (\*.in files)
---------------------------
`I**`
01 - baseline runtime options, for an inter-annual hindcast, writing daily averaged outputs for the full
     domain, and hourly outputs for the surface

Carried over from croco\_v2.0.1/I01 with three changes:
  - the ERA5 ONLINE path now points at /mnt/lustre/groups/ERTH1103/data/ERA5/eez\_for\_croco/
  - the ONLINE record end is 2025 12 (was 2019 12), to cover the run window
  - NLEVEL is 1, not 4 - I01 was the AGRIF 4-grid baseline, and C12 has AGRIF undefined

No new .in blocks are needed for v2.1.1. The release adds the `diag_mld_crit_threshold`,
`diag_mld_depth_ref` and `wave_maker` keywords, but all three sit behind CPP keys this configuration
does not define, and no keywords were removed.

Surface and boundary forcing
----------------------------
you'll need to have directory names corresponding to what you have specified in `myenv_inter.sh` i.e.
for `ATMOS_BULK` and `OGCM`.

`GLORYS` - initial and climatology files. C12 defines CLIMATOLOGY and undefines FRC_BRY, so we generate
           clm files (`make_clm_inter.bash`), not bry files. `myenv_inter.sh` has CLIMATOLOGY_FILES=1
           and BOUNDARY_FILES=0 to match.
`TPXO10` - tidal forcing.
`ERA5`   - not a directory here; surface forcing is read ONLINE straight from the raw ERA5 files at the
           path given at the bottom of `I01/croco_inter.in`. BULK_FILES=0 accordingly.

The make\_\*.bash scripts run on lengau and reference the raw input data under
/mnt/lustre/groups/ERTH1103/data/.

Run window
----------
1993-03 to 2025-11, cold start (RSTFLAG=0), submitted by `run_croco_inter_chpc.bash` as 99 dependent
4-month jobs. Note that window is not a whole number of 4-month jobs, so that script clamps the final
job to END_DATE - without the clamp the chain would run on to 2026-02 and off the end of the forcing.

Yorig is 1993 throughout.
</pre>
