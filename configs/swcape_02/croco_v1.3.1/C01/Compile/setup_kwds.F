#include "cppdefs.h"

      subroutine setup_kwds (ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new keywords were !!!!!!
!!!!!! introduced into "read_inp.F". To create or refresh this !!!!!!
!!!!!! file use compile and execute "checkkwds.F" as an        !!!!!!
!!!!!! independent program, or use commands "make checkkwds"   !!!!!!
!!!!!! or "make depend".                                       !!!!!!

      implicit none
      integer ierr, is,ie
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      is=1
#ifdef SEDIMENT
#endif
#ifdef FLOATS
#endif
#ifdef STATIONS
#endif
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
#endif
#ifdef WKB_WWAVE
#endif
#ifdef ONLINE
#endif
#ifdef MPI
#endif
#ifdef NBQ
#endif
#ifdef SOLVE3D
#endif
#ifdef DIAGNOSTICS_BIO
#endif
#ifdef USE_CALENDAR
#endif
#ifdef XIOS
#endif
#if defined BASIN
#elif defined CANYON
#elif defined EQUATOR
#elif defined GRAV_ADJ
# ifdef NBQ
# else
# endif
#elif defined ISOLITON
#elif defined INNERSHELF
#elif defined INTERNAL
#elif defined OVERFLOW
#elif defined RIVER
#elif defined SEAMOUNT
#elif defined SHELFRONT
#elif defined SOLITON
#elif defined UPWELLING
#elif defined VORTEX
#elif defined JET
#elif defined SHOREFACE
#elif defined RIP
#elif defined SANDBAR
#elif defined THACKER
#elif defined SINGLE_COLUMN && defined KATO_PHILIPS
#elif defined SINGLE_COLUMN && defined WILLIS_DEARDORFF
#elif defined SINGLE_COLUMN && defined DIURNAL_CYCLE
#elif defined SINGLE_COLUMN && defined FORCED_EKBBL
#elif defined SINGLE_COLUMN && defined FORCED_DBLEEK
#elif defined SINGLE_COLUMN && defined FORCED_NONROTBBL
#elif defined SINGLE_COLUMN && defined FORCED_OSCNONROTBBL
#elif defined TANK
#elif defined MOVING_BATHY
#elif defined ACOUSTIC
#elif defined SWASH
#elif defined KH_INST
#elif defined IGW
# ifdef NBQ
# else
# endif
#elif defined TS_HADV_TEST
# ifdef SOLID_BODY_ROT
# elif defined DIAGONAL_ADV
# elif defined SOLID_BODY_PER
# endif
#elif defined DUNE
# ifdef DUNE3D
# elif defined ANA_DUNE
# else
# endif
#elif defined SED_TOY
# ifdef SED_TOY_RESUSP
# elif defined SED_TOY_CONSOLID
# elif defined SED_TOY_FLOC
# elif defined SED_TOY_ROUSE
# endif
#elif defined TIDAL_FLAT
#elif defined ESTUARY
#else
#endif
#ifdef MPI
#else
#endif
#ifdef AGRIF
#ifdef AGRIF_ADAPTIVE
#else
#endif
#endif
#ifdef LOGFILE
# ifndef AGRIF
# else
# endif
# ifdef AGRIF
# endif
#endif
#ifdef AVERAGES
#endif
#ifdef TRACERS
# ifdef DIAGNOSTICS_TS
#  ifdef AVERAGES
#  endif
#  ifdef DIAGNOSTICS_TS_MLD
#   ifdef AVERAGES
#   endif
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_UV
#ifdef AVERAGES
#endif
#endif
#ifdef DIAGNOSTICS_VRT
#ifdef AVERAGES
#endif
#endif
#ifdef DIAGNOSTICS_EK
#ifdef AVERAGES
#endif
#endif
#if defined DIAGNOSTICS_PV && defined TRACERS
#ifdef AVERAGES
#endif
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
# ifdef AVERAGES
# endif
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
# ifdef AVERAGES
# endif
#endif
#ifdef DIAGNOSTICS_BIO
#ifdef AVERAGES
#endif
#endif
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='title'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef USE_CALENDAR
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='start_date'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef ANA_INITIAL
# endif
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='end_date'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='output_time_steps'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef XIOS
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='xios_origin_date'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='time_stepping'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef NBQ
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='time_stepping_nbq'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SOLVE3D
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='S-coord'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='initial'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef ANA_INITIAL
#endif
#if defined MPI && defined PARALLEL_FILES
#endif
#ifdef ANA_INITIAL
#endif
#ifndef ANA_GRID
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='grid'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
#if !defined NO_FRCFILE  && \
        ( defined TIDES \
          || (defined MRL_WCI     && !defined ANA_WWAVE && \
                                     !defined WKB_WWAVE && \
                                     !defined OW_COUPLING) \
          || (defined SFLX_CORR   && !defined ANA_SSS   && \
                                   !defined BULK_FLUX   && \
                                   !defined OA_COUPLING) \
          || (defined QCORRECTION && !defined ANA_SST   && \
                                   !defined OA_COUPLING) \
          || (defined BBL         && !defined ANA_BSEDIM  && \
                                   !defined SEDIMENT) \
          || (!defined ANA_STFLUX && !defined BULK_FLUX && \
                                   !defined OA_COUPLING && \
                                   defined SOLVE3D) \
          || (defined SALINITY    && !defined ANA_SSFLUX && \
                                   !defined BULK_FLUX && \
                                   !defined OA_COUPLING && \
                                    defined SOLVE3D) \
          || (!defined ANA_SRFLUX && !defined BULK_FLUX && \
                                   !defined OA_COUPLING && \
                                   defined SOLVE3D) \
          || (!defined ANA_SMFLUX && !defined BULK_FLUX && \
                                   !defined OA_COUPLING) \
          || ( defined BULK_FLUX    && !defined SALINITY ) \
          || ( defined WAVE_OFFLINE && !defined MUSTANG ) )
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='forcing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif  /* !defined NO_FRCFILE */
# if defined WAVE_OFFLINE && defined MUSTANG
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wave_offline'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif  /* !defined WAVE_OFFLINE */
#if defined BIOLOGY && defined PISCES
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='biology'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#endif
#if defined BULK_FLUX
# ifndef ONLINE /* ONLINE FORCING */
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bulk_forcing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  if defined MPI && defined PARALLEL_FILES
#  endif
# endif  /* ONLINE FORCING */
#endif
#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY  && !defined ANA_SSH) || \
      (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
      (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='climatology'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AGRIF
# endif
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AGRIF
# endif
#endif
#if !defined ANA_BRY && defined FRC_BRY
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='boundary'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef AGRIF
#  endif
#  if defined MPI && defined PARALLEL_FILES
#  endif
#  ifdef AGRIF
#  endif
#endif
#if defined WKB_WWAVE && !defined ANA_BRY_WKB
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wkb_boundary'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef AGRIF
#  endif
#  if defined MPI && defined PARALLEL_FILES
#  endif
#  ifdef AGRIF
#  endif
#endif
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='restart'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#ifdef USE_CALENDAR
#endif
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='history'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#ifdef USE_CALENDAR
#endif
#ifdef AVERAGES
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='averages'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
#ifdef USE_CALENDAR
#endif
#endif
#if defined DIAGNOSTICS_TS
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagnostics'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diag_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /* DIAGOSTICS_TS */
#if defined DIAGNOSTICS_UV
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagnosticsM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagM_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /*DIAGNOSTICS_UV */
#if defined DIAGNOSTICS_VRT
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_vrt'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_vrt_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /*DIAGNOSTICS_VRT */
#if defined DIAGNOSTICS_EK
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_ek'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_ek_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /*DIAGNOSTICS_EK */
#if defined DIAGNOSTICS_PV
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_pv'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_pv_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /*DIAGNOSTICS_PV */
# if defined DIAGNOSTICS_EDDY && ! defined XIOS
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_eddy'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_eddy_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /*DIAGNOSTICS_EDDY */
#if defined OUTPUTS_SURFACE && ! defined XIOS
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='surf'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='surf_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /*OUTPUTS_SURFACE */
#ifdef DIAGNOSTICS_BIO
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagnostics_bio'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AVERAGES
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbio_avg'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES
# endif
# endif
#endif /* DIAGNOSTICS_BIO */
#ifdef FLOATS
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='floats'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef AGRIF
#endif
# if defined MPI && defined PARALLEL_FILES
# endif
#ifdef AGRIF
#endif
#endif /* FLOATS */
#ifdef STATIONS
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='stations'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AGRIF
# endif
# if defined MPI && defined PARALLEL_FILES
# endif
# ifdef AGRIF
# endif
#endif /* STATIONS */
#ifdef ASSIMILATION
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='assimilation'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='primary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
# ifdef TRACERS
# endif
#endif
#ifdef SOLVE3D
#endif
#ifdef SOLVE3D
# ifdef TRACERS
# endif
#endif
#if !defined SOLVE3D && defined RIP
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='auxiliary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#elif defined SOLVE3D
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='auxiliary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef TEMPERATURE
# else
# endif
# ifdef SALINITY
# else
# endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
# else
# endif
# ifdef VIS_COEF_3D
# else
# endif
# ifdef DIF_COEF_3D
# else
# endif
# if defined LMD_SKPP || defined GLS_MIXING
# else
# endif
# ifdef LMD_BKPP
# else
# endif
# ifdef TEMPERATURE
# else
# endif
# ifdef SALINITY
# else
# endif
# ifdef TEMPERATURE
# else
# endif
# ifdef BULK_FLUX
# else
# endif
#ifdef BHFLUX
#else
#endif
#if defined BWFLUX && defined SALINITY
#else
#endif
# ifdef MORPHODYN
# else
# endif
# if defined BIOLOGY && !defined PISCES
#  ifdef BIO_NChlPZD
#   ifdef OXYGEN
#   endif
#  elif defined BIO_BioEBUS
#  endif
# endif
# ifdef TEMPERATURE
# endif
# ifdef SALINITY
# endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
# endif
# ifdef VIS_COEF_3D
# endif
# ifdef DIF_COEF_3D
# endif
# if defined LMD_SKPP || defined GLS_MIXING
# endif
# ifdef LMD_BKPP
# endif
# ifdef TEMPERATURE
# endif
# ifdef SALINITY
# endif
# ifdef TEMPERATURE
# endif
# if defined BULK_FLUX
# endif
# if defined BHFLUX
# endif
# if defined BWFLUX && defined SALINITY
# endif
# if defined BIOLOGY && !defined PISCES
#  ifdef BIO_NChlPZD
#   ifdef OXYGEN
#   endif
#  elif defined BIO_BioEBUS
#  endif
# endif
# ifdef TEMPERATURE
# endif
# ifdef SALINITY
# endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
#  endif
# ifdef VIS_COEF_3D
# endif
# ifdef DIF_COEF_3D
# endif
# if defined LMD_SKPP || defined GLS_MIXING
# endif
# ifdef LMD_BKPP
# endif
# ifdef BULK_FLUX
# endif
# ifdef MORPHODYN
# endif
# if defined BIOLOGY && !defined PISCES
#  ifdef BIO_NChlPZD
#   ifdef OXYGEN
#   endif
#  elif defined BIO_BioEBUS
#  endif
# endif /* BIOLOGY */
# ifdef TEMPERATURE
# endif
# ifdef SALINITY
# endif
# ifdef TEMPERATURE
# endif
# ifdef BHFLUX
# endif
# if defined BWFLUX && defined SALINITY
# endif
# ifdef GLS_MIXING
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='gls_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif /* SOLVE3D */
#ifdef AVERAGES
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='primary_averages'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
#  ifdef TRACERS
#  endif
# endif
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
#  ifdef TRACERS
#  endif
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='auxiliary_averages'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef TEMPERATURE
# else
# endif
# ifdef SALINITY
# else
# endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
# else
# endif
# ifdef VIS_COEF_3D
# else
# endif
# ifdef DIF_COEF_3D
# else
# endif
# if defined LMD_SKPP || defined GLS_MIXING
# else
# endif
# ifdef LMD_BKPP
# else
# endif
# ifdef TEMPERATURE
# else
# endif
# ifdef SALINITY
# else
# endif
# ifdef TEMPERATURE
# else
# endif
# ifdef BULK_FLUX
# else
# endif
#ifdef BHFLUX
#else
#endif
#if defined BWFLUX && defined SALINITY
#else
#endif
# ifdef MORPHODYN
# else
# endif
# if defined BIOLOGY && !defined PISCES
#  ifdef BIO_NChlPZD
#   ifdef OXYGEN
#   endif
#  elif defined BIO_BioEBUS
#  endif
# endif
#  ifdef TEMPERATURE
#  endif
#  ifdef SALINITY
#  endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
# endif
# ifdef VIS_COEF_3D
# endif
# ifdef DIF_COEF_3D
# endif
#  if defined LMD_SKPP || defined GLS_MIXING
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef TEMPERATURE
#  endif
#  ifdef SALINITY
#  endif
#  ifdef TEMPERATURE
#  endif
#  ifdef BULK_FLUX
#  endif
# if defined BHFLUX
#endif
# if defined BWFLUX && defined SALINITY
#endif
#  if defined BIOLOGY && !defined PISCES
#   ifdef BIO_NChlPZD
#    ifdef OXYGEN
#    endif
#   elif defined BIO_BioEBUS
#   endif
#  endif  /* BIOLOGY */
#  ifdef TEMPERATURE
#  endif
#  ifdef SALINITY
#  endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
#  endif
#  ifdef VIS_COEF_3D
#  endif
#  ifdef DIF_COEF_3D
#  endif
#  if defined LMD_SKPP || defined GLS_MIXING
#  endif
#  ifdef LMD_BKPP
#  endif
#  ifdef BULK_FLUX
#  endif
# ifdef MORPHODYN
# endif
#  if defined BIOLOGY && !defined PISCES
#   ifdef BIO_NChlPZD
#    ifdef OXYGEN
#    endif
#   elif defined BIO_BioEBUS
#   endif
#  endif   /* BIOLOGY */
#  ifdef TEMPERATURE
#  endif
#  ifdef SALINITY
#  endif
#  ifdef TEMPERATURE
# endif
# ifdef BHFLUX
# endif
# if defined BWFLUX && defined SALINITY
# endif
#  ifdef GLS_MIXING
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='gls_averages'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  endif
# endif /* SOLVE3D */
#endif /* AVERAGES */
#ifdef FLOATS
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='float_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef AGRIF
#endif
#ifdef AGRIF
#endif
#endif /* FLOATS */
#if defined DIAGNOSTICS_TS
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diag3D_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef TRACERS
# endif  /* TRACERS */
# ifdef DIAGNOSTICS_TS_MLD
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diag2D_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef TRACERS
#  endif  /* TRACERS */
# endif  /* DIAGNOSTICS_TS_MLD */
# ifdef AVERAGES
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diag3D_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef TRACERS
#  endif  /* TRACERS */
#  ifdef DIAGNOSTICS_TS_MLD
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diag2D_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#   ifdef TRACERS
#   endif  /* TRACERS */
#  endif  /* DIAGNOSTICS_TS_MLD */
# endif  /* AVERAGES */
#endif  /* DIAGNOSTICS_TS */
#if defined DIAGNOSTICS_UV
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagM_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagM_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_UV */
#if defined DIAGNOSTICS_VRT
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_vrt_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_vrt_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_VRT */
#if defined DIAGNOSTICS_EK
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_ek_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_ek_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_EK */
#if defined DIAGNOSTICS_PV
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_pv_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef TRACERS
# endif /* TRACERS  */
# ifdef AVERAGES
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_pv_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef TRACERS
#  endif /* TRACERS  */
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_PV */
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
      ie=is +25
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_eddy_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +25
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diags_eddy_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_EDDY */
#if defined OUTPUTS_SURFACE && ! defined XIOS
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='surf_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='surf_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif /* AVERAGES  */
#endif /*OUTPUTS_SURFACE*/
#ifdef DIAGNOSTICS_BIO
      ie=is +26
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbioFlux_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +27
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbioVSink_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      ie=is +28
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbioGasExc_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  endif
# ifdef AVERAGES
      ie=is +26
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbioFlux_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +27
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbioVSink_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      ie=is +28
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='diagbioGasExc_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  endif
# endif  /* AVERAGES */
#endif /* DIAGNOSTICS_BIO */
#ifdef STATIONS
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='station_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AGRIF
# endif
# ifdef AGRIF
# endif
#endif /* STATIONS */
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='rho0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#if defined UV_VIS2 || defined UV_VIS4 || defined SPONGE_VIS2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='lateral_visc'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#if defined UV_VIS2 || defined SPONGE_VIS2
#endif
#ifdef UV_VIS4
#endif
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bottom_drag'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='gamma2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#ifdef SOLVE3D
# if defined TS_DIF2 || defined SPONGE_DIF2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='tracer_diff2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef TRACERS
#  endif /* TRACERS */
# endif
# ifdef TS_DIF4
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='tracer_diff4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef TRACERS
#  endif /* TRACERS */
# endif
# if !defined LMD_MIXING && !defined BVF_MIXING && \
     !defined GLS_MIXING
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='vertical_mixing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef TRACERS
#  endif /* TRACERS */
#  ifdef TRACERS
#  endif /* TRACERS */
# endif
# ifdef BODYFORCE
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bodyforce'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif
#if  (defined SPONGE && !defined SPONGE_GRID)
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sponge'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#elif  (defined SPONGE && defined SPONGE_GRID)
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sponge'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#if  defined T_FRC_BRY     || defined M2_FRC_BRY    || \
     defined M3_FRC_BRY    || defined Z_FRC_BRY     || \
     defined TCLIMATOLOGY  || defined M2CLIMATOLOGY || \
     defined M3CLIMATOLOGY || defined ZCLIMATOLOGY
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nudg_cof'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if  defined AGRIF && !defined AGRIF_OBC_M2ORLANSKI && \
     !defined AGRIF_OBC_M3ORLANSKI && !defined AGRIF_OBC_TORLANSKI
# endif
# if defined AGRIF && !defined AGRIF_OBC_M2ORLANSKI && \
     !defined AGRIF_OBC_M3ORLANSKI && !defined AGRIF_OBC_TORLANSKI
# endif
#endif
#ifdef SOLVE3D
# if !defined NONLIN_EOS
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='lin_EOS_cff'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif
# ifdef SEDIMENT
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sediments'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +23
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sediment_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#   ifdef SUSPLOAD
#   endif
#   ifdef BEDLOAD
#   endif
#   if defined MIXED_BED || defined COHESIVE_BED
#   endif
#   ifdef AVERAGES
#   ifdef SUSPLOAD
#   endif
#   ifdef BEDLOAD
#   endif
#   endif
#   ifdef SUSPLOAD
#   endif
#   ifdef BEDLOAD
#   endif
#   if defined MIXED_BED || defined COHESIVE_BED
#   endif
# elif defined MUSTANG
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='sediments_mustang'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif /* SEDIMENT */
#endif
#if defined SUBSTANCE && !defined MUSTANG
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='substance'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BBL
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bbl_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif /* BBL */
#ifdef MRL_WCI
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wci_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef SOLVE3D
# endif
# ifdef SOLVE3D
# endif
# ifdef AVERAGES
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wci_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef SOLVE3D
#  endif
#  ifdef SOLVE3D
#  endif
# endif /* AVERAGES */
#endif /* MRL_WCI */
#if defined MRL_WCI || defined OW_COUPLING
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wave_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef AVERAGES
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wave_average_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif /* WKB_WWAVE || OW_COUPLING || WAVE_OFFLINE */
#if defined WKB_WWAVE
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wkb_wwave'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# ifdef WAVE_BREAK_TG86
# elif defined WAVE_BREAK_TG86A
# elif defined WAVE_BREAK_CT93
# else
# endif
# ifdef WAVE_ROLLER
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='wkb_roller'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# endif
#endif /* WKB_WWAVE */
#ifdef ANA_PSOURCE
# ifdef PSOURCE_NCFILE
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='psource_ncfile'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  ifdef PSOURCE_NCFILE_TS
#  else
#   ifdef TRACERS
#   endif /* TRACERS */
#  endif /* PSOURCE_NCFILE_TS  */
#  ifdef MPI
#  endif
#  ifdef TRACERS
#  endif /* TRACERS */
#  ifndef PSOURCE_NCFILE_TS
#   ifdef TRACERS
#   endif /* TRACERS */
#  endif
# else
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='psource'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#  if defined SOLVE3D && defined TRACERS
#  endif
#  ifdef MPI
#  endif
#  if defined SOLVE3D && defined TRACERS
#  endif
# endif   /*  define PSOURCE_NCFILE  */
#endif /* ANA_PSOURCE */
#if defined BHFLUX || ( defined BWFLUX && defined SALINITY)
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bottom_forcing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
# if defined MPI && defined PARALLEL_FILES && defined MPI_NOLAND
# endif
# if defined MPI && defined PARALLEL_FILES && !defined MPI_NOLAND
# endif
#endif  /* BHFLUX || (SALINITY && BWFLUX )  */
#ifdef ONLINE   /* JG Online */
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='online'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif /* ONLINE */
#ifdef LOGFILE
#endif
#ifdef MPI
#endif
      return
  99  MPI_master_only write(stdout,'(/1x,A,A/14x,A)')
     &  'SETUP_KWDS ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
