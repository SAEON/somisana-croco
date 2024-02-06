      program ncrename
      implicit none
      logical found
      character*80 ncname, varname, newname,name
      integer*4 stdout, n, ierr, ncid, varid,  nvars
     &           , lstr, lvar, lnew, lenstr, iargc
      parameter (stdout=6)
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_64bit_data
      integer*4 nf_cdf5
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_format_64bit_offset
      integer*4 nf_format_64bit_data
      integer*4 nf_format_cdf5
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_int64
      external        nf_put_att_int64
      integer*4         nf_get_att_int64
      external        nf_get_att_int64
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_compact
      parameter (nf_compact = 2)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_szip
      external nf_def_var_szip
      integer*4 nf_inq_var_szip
      external nf_inq_var_szip
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_def_var_filter
      external nf_def_var_filter
      integer*4 nf_inq_var_filter
      external nf_inq_var_filter
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      integer*4 nf_set_log_level
      external nf_set_log_level
      n=iargc()
      if (n.ne.3) then
        write(stdout,'(/1x,A,1x,A/32x,A/)') 'Usage of ncrename',
     &             'should be:', 'ncrename file.nc name newname'
        stop
      endif
      call getarg(1,ncname)
      call getarg(2,varname)
      call getarg(3,newname)
      lstr=lenstr(ncname)
      lvar=lenstr(varname)
      lnew=lenstr(newname)
      found=.false.
      lstr=lenstr(ncname)
      ierr=nf_open (ncname(1:lstr), nf_write, ncid)
      if (ierr .ne. nf_noerr) then
        write(stdout,'(/8x,A,1x,A,A/)') 'Cannot open netCDF file',
     &                                          ncname(1:lstr),'.'
        goto 100
      endif
      ierr=nf_redef (ncid)
      if (ierr .ne. nf_noerr) then
        write(stdout,'(/8x,A,1x,A,1x,A,A/)') 'Cannot switch to',
     &    'redefinition mode for netCDF file', ncname(1:lstr),'.'
        goto 99
      endif
      ierr=nf_inq_varid (ncid, varname(1:lvar), varid)
      if (ierr .ne. nf_noerr) goto 1
      ierr=nf_rename_var (ncid, varid, newname(1:lnew))
      if (ierr .eq. nf_noerr) then
        write(stdout,'(/8x,5A/)') 'Renamed variable ''',
     &  varname(1:lvar), ''' into ''', newname(1:lnew),'''.'
        found=.true.
      else
        write(stdout,'(/8x,6A/8x,A,I3/)') 'Cannot rename ',
     &     'variable ''',  varname(1:lvar), ''' into ''',
     &     newname(1:lnew),'''.', 'netCDF error status =',i err
        if (ierr. eq. nf_enameinuse) then
          write(stdout,'(8x,A)') 'This name is already in use.'
        endif
      endif
      goto 98
  1   ierr=nf_inq_dimid (ncid, varname(1:lvar), varid)
      if (ierr .ne. nf_noerr) goto 2
      ierr=nf_rename_dim (ncid, varid, newname(1:lnew))
      if (ierr .eq. nf_noerr) then
        write(stdout,'(/8x,5A/)') 'Renamed dimension ''',
     &    varname(1:lvar), ''' into ''', newname(1:lnew), '''.'
        found=.true.
      else
        write(stdout,'(/8x,6A/8x,A,I3/)') 'Cannot rename ',
     &   'dimension ''', varname(1:lvar), ''' into ''',
     &    newname(1:lnew), '''.', 'netCDF error status =', ierr
        if (ierr. eq. nf_enameinuse) then
          write(stdout,'(8x,A/)') 'This name is already in use.'
        endif
      endif
      goto 98
  2   ierr=nf_rename_att (ncid, nf_global, varname(1:lvar),
     &                                     newname(1:lnew))
      if (ierr .eq. nf_noerr) then
        write(stdout,'(/8x,5A/)') 'Renamed global attribute ''',
     &     varname(1:lvar), ''' into ''', newname(1:lnew), '''.'
        found=.true.
        goto 98
      else
        ierr=nf_inq_nvars (ncid, nvars)
        if (ierr .eq. nf_noerr) then
          do varid=1,nvars
            ierr=nf_rename_att (ncid, varid, varname(1:lvar),
     &                                       newname(1:lnew))
            if (ierr .eq. nf_noerr) then
              ierr=nf_inq_varname (ncid, varid, name)
              n=lenstr (name)
              write(stdout,'(/8x,7A/)') 'Renamed attribute ''',
     &          varname(1:lvar), ''' into ''', newname(1:lnew),
     &                   ''' for variable ''', name(1:n), '''.'
              found=.true.
            endif
          enddo
        else
          write(stdout,'(/8x,5A/)') 'Cannot determine number of ',
     &            'variables in netCDF file ', ncname(1:lstr), '.'
        endif
      endif
      if (.not.found) write(stdout,'(/8x,6A/)') 'Cannot find ',
     &      'object ''', varname(1:lvar), ''' in netCDF file ',
     &                                      ncname(1:lstr), '.'
  98  ierr=nf_enddef (ncid)
      if (ierr. ne. nf_noerr) then
        write(stdout,'(/8x,A,1x,A,1x,A,A/)') 'Cannot switch',
     &    'to data mode for netCDF file', ncname(1:lstr),'.'
      endif
  99  ierr=nf_close (ncid)
 100  stop
      end
