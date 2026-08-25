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
     no TIDERAMP, open boundaries E/W/S (OBC\_NORTH undefined),
     PLUS OUTPUTS\_SURFACE defined - see "Surface output" below. This is the one deliberate CPP
     divergence from croco\_v2.0.1/C12; every other key matches it exactly.

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
     domain, and hourly averaged outputs for the surface

Carried over from croco\_v2.0.1/I01 with four changes:
  - the ERA5 ONLINE path now points at /mnt/lustre/groups/ERTH1103/data/ERA5/eez\_for\_croco/
  - the ONLINE record end is 2025 12 (was 2019 12), to cover the run window
  - NLEVEL is 1, not 4 - I01 was the AGRIF 4-grid baseline, and C12 has AGRIF undefined
  - surface output is actually enabled and uses the v2.1.1 filenames (see below)

Surface output
--------------
Hourly *averaged* surface fields, written to `croco_avg_surf.nc`. The instantaneous surface stream is
off.

  surf:                F   NUMHISSURF   ->  croco_his_surf.nc   (not written)
  surf_avg:            T   NUMAVGSURF   ->  croco_avg_surf.nc
  surf_history_fields: T
  surf_average_fields: T

  NH_AVGSURF=1  -> nwrtsurf_avg = 3600/240 = 15 timesteps exactly
  NH_HISSURF=24 -> inert while ldefsurf is F

IMPORTANT - this needs OUTPUTS\_SURFACE defined in cppdefs.h, and it is NOT defined in any of the
croco\_v2.0.1 sa\_eez\_01 compile options (C10, C12 and C13 all have it undef, as do both stock CROCO
trees). It is defined here, matching sa\_west\_02/croco\_v2.1.1/C10.

Without that key the whole `surf:`/`surf_avg:` machinery is compiled out: read\_inp.F only recognises
those keywords under `#if defined OUTPUTS_SURFACE && ! defined XIOS`, and an unrecognised keyword is a
warning that is silently disregarded (read\_inp.F ~3105), not an error. So croco\_v2.0.1 has been
carrying surf blocks in its .in files that never produced a file. If you compare output against a
croco\_v2.0.1 run, expect surface files here and none there.

For the surf streams specifically, the .in flags are real on/off switches rather than the usual
"create a new file" flags: output.F gates the writes on `.and. ldefsurf` / `.and. ldefsurf_avg`, and
step.F gates the averaging accumulation (`set_surf_avg`) on ldefsurf\_avg too, so nothing is computed
or written for the disabled stream. Neither flag is reassigned anywhere in the source.

Do not generalise that to the other streams - it is specific to surf. The history write
(output.F ~277) is gated only on nwrt and wrthis, NOT on ldefhis, so LDEFHIS=F there keeps its
conventional "do not create a new file" meaning and does not switch history output off.

Archiving
---------
`run_croco_inter.pbs` moves croco\_avg\_surf.nc out of the scratch dir into
$RUNNAME/output/croco\_surf\_avg\_Y<yyyy>M<mm>.nc at the end of each month, in the Archive block after
the "month finished properly" check. Unlike the other mv lines it is wrapped in an existence test that
prints a greppable WARNING if the file is absent, because a silent mv failure would otherwise cost a
month of surface output somewhere in a 99 job chain without anyone noticing.

The averaging window divides the month exactly - nwrtsurf\_avg = 3600/240 = 15 timesteps against
NUMTIMES of 10080/10440/10800/11160 for 28/29/30/31 day months, giving 672/696/720/744 hourly records
with no remainder - so the last average of the month is always written before the file is moved.

CROCO writes a single joined file here (PARALLEL\_FILES, NC4PAR and XIOS are all undefined, so no
per-node suffix; nrpfsurf\_avg=0, so no time-index suffix), which is why the plain croco\_avg\_surf.nc
name in the mv line is correct.

Filenames are croco\_his\_surf.nc / croco\_avg\_surf.nc, the convention used by croco\_v2.0.1, sa\_eez\_02,
sa\_west\_02 v1.3.1 and v2.0.1, and oceanmotion-models gulf\_01. In that last repo the name is a
load-bearing interface - tidal\_analysis.py globs croco\_avg\_surf*, and the WW3 hindcast reads
croco\_avg\_surf\_<month>\_current.nc / \_level.nc as forcing.

Deliberately NOT following sa\_west\_02/croco\_v2.1.1 here: that config names its files croco\_surf.nc /
croco\_surf\_avg.nc, but it runs surf\_avg=F, so croco\_surf\_avg.nc is never actually written anywhere.
Note also that sa\_west\_02 uses the opposite streams to this config - daily *instantaneous* surf, no
averages.

The closest match to this configuration is oceanmotion-models gulf\_01/croco\_v1.3.1/hindcast:
OUTPUTS\_SURFACE defined, surf=F, surf\_avg=T, NH\_AVGSURF=1, NH\_HISSURF=24, single joined file, and the
same 672/696/720/744 records per month. It differs only in DT0 (180 s there, 240 s here) and in having
no existence guard on the archive mv.

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

Monthly looping
---------------
`run_croco_inter.pbs` takes START_DATE and END_DATE (`qsub -v`), not the NY\_START/NY\_END/NM\_START/
NM\_END quadruple the older scripts used, and steps months with a single date-driven loop:

    CURRENT_DATE=$START_DATE
    while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do
        TIME=$(time_stamp "$CURRENT_DATE")
        ...
        CURRENT_DATE=$(date -d "$CURRENT_DATE +1 month" +%Y-%m-01)
    done

This follows sa\_west\_02/croco\_v2.1.1\_cpl\_ww3/run_inter.pbs. It replaces the nested year/month loops
with their NY\_END+1 / MONTH\_END=13 bookkeeping, and replaces ~36 lines of hand-rolled month-length and
leap-year arithmetic with `date -d "<1st of month> +1 month -1 day" +%d`. Restart file names come from
`date -d "$START_DATE -1 month"` rather than a manual month decrement. 325 lines -> 258.

Behaviour is unchanged. The old and new loop logic were both transcribed into a harness and driven
through the real 99-job submission chain: the sequence of TIME stamps, NY, NM, NDAYS and ini/restart
filenames is identical for all 393 months, and also for century-leap edge cases (1900, 2000, 2100),
year-boundary spans and single-month jobs.

TIME\_SCHED and NY\_SPIN are gone entirely - removed from both the run script and I01/myenv\_in.sh:

- TIME\_SCHED selected yearly (0) vs monthly (1) output files. This config only ever runs monthly, and
  the yearly branch was partly nonsensical anyway (it set NM\_START=1979 on a variable holding a month
  number). One file per month is now simply what the script does.
- NY\_SPIN forced spin-up years to 365 days by dropping the leap day. It was set to 0, so it never did
  anything, and it was measured from the *job's* start year rather than the run's - meaning it would
  have reset every 4 months across a chained run even if it had been switched on. Leap years are now
  just leap years.

`time_stamp()` honours MTH\_FORMAT from myenv\_in.sh rather than hardcoding %m, so the Y<yyyy>M<mm>
naming stays configurable and keeps matching the GLORYS/TPXO10 filenames.
</pre>
