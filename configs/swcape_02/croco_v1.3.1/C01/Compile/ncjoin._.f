      program ncjoin
      implicit none
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
      integer*4, parameter :: maxdims=32, maxvars=96
      integer*4, parameter :: part_XI=1,  part_ETA=2, part_2D=3
      logical complete, clean_set,  digit, var_mask,    lnewvar
      integer*4 nargs, nnodes,  size_XI,  XI_rho, id_xi_rho,  id_xi_u,
     &        arg,   node,    size_ETA, ETA_rho,id_eta_rho, id_eta_v,
     &        ierr,  maxnodes,size_S,   tsize,  unlimdimid, rec,
     &        ntest, nctarg,  ndims,    size,   code_size,  lvar,
     &        nvars, ngatts,  varatts,  size1,  code_size_bak,
     &        i,j,k, is,ie,   lncn,     ltrg,   lstr, lbak, lenstr,
     &        LLm, MMm
      character(len=8) sffx, sffx_bak
      character(len=32) vname(maxvars), dimname(maxdims)
      character(len=64) nctestname, nctargname, root, root_bak, string
      character(len=64), dimension(:), allocatable :: ncname
      integer*4, dimension(:), allocatable:: ncid, xi_start, eta_start
      integer*4, dimension(:,:), allocatable :: vid, dimsize
      logical, dimension(:), allocatable :: western_edge,
     &         eastern_edge, southern_edge, northern_edge
      logical series(maxvars)
      integer*4, dimension(maxvars) :: varid,   vnode,  vdims,
     &                               vartype, part_type
      integer*4, dimension(maxdims) :: dimid,   ldim,   ibuff,
     &                               start,   count,  start1
      integer*4, dimension(maxdims,maxvars) :: dimids
      integer*4 max_buff_size, alloc_buff_size
      real*8, allocatable, dimension(:) :: buff
      character*2000 :: buff_str
      integer*4 count1(maxdims)
      integer*4 max_bfr_out, alloc_bfr_out
      real*8, allocatable, dimension(:) :: bfr_out
      logical del_part_files
      character(len=128) rmcmd
      real*4 tstart, RUN_time, CPU_time(2)
      integer*4 iclk(2), nclk, clk_rate, clk_max, iclk_init
      integer*8 net_read_size, net_wrt_size, net_fcrt_clk,
     &          net_read_clk,  net_wrt_clk,  net_assm_clk,
     &          net_sync_clk,  net_gray_clk, inc_clk
      real*8 ReadSize, ReadTime, WrtSize,  WrtTime,
     &       FcrtTime, AssmTime, SyncTime, GrayTime
      integer*8 net_rmcmd_clk
      net_fcrt_clk = 0
      call etime(CPU_time, tstart)
      nclk=1
      call system_clock (iclk(nclk), clk_rate, clk_max)
      iclk_init=iclk(nclk)
      net_read_clk=0
      net_read_size=0
      net_wrt_size=0
      net_wrt_clk=0
      net_sync_clk=0
      net_assm_clk=0
      net_gray_clk=0
      del_part_files=.false.
      net_rmcmd_clk=0
      ntest=-1
      maxnodes=-1
      max_buff_size=0
      alloc_buff_size=0
      max_bfr_out=0
      alloc_bfr_out=0
      nargs=iargc()
      arg=0
  1   nnodes=-1
      root_bak(1:1)=' '
      sffx_bak(1:1)=' '
      code_size_bak=-1
  2   arg=arg+1
        call getarg(arg,nctestname)
        lncn=lenstr(nctestname)
        if (arg.eq.1 .and. (lncn.eq.2 .and. nctestname(1:2).eq.'-d'
     &     .or. lncn.eq.8 .and. nctestname(1:8).eq.'--delete') ) then
           write(*,'(/1x,2A/)') '>>>> Flag to delete partial files ',
     &                                                  'is raised.'
           del_part_files=.true.
          goto 2
        endif
        if (ntest.ne.-1) then
          ierr=nf_close(ntest)
          ntest=-1
        endif
        ierr=nf_open (nctestname, nf_nowrite, ntest)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_att (ntest, nf_global, 'partition_ucla', i, lvar)
          if (ierr .eq. nf_noerr) then
            if (i.eq.nf_int .and. lvar.eq.6) then
              ierr=nf_get_att_int (ntest,nf_global, 'partition_ucla',
     &           ibuff)
              if (ierr .eq. nf_noerr) then
                if (nnodes.eq.-1) then
                  nnodes=ibuff(2)
                  if (nnodes.gt.maxnodes) then
                    maxnodes=nnodes
                    if (allocated(ncid)) then
                      deallocate(dimsize)
                      deallocate(vid)
                      deallocate(northern_edge)
                      deallocate(southern_edge)
                      deallocate(eastern_edge)
                      deallocate(western_edge)
                      deallocate(eta_start)
                      deallocate(xi_start)
                      deallocate(ncname)
                      deallocate(ncid)
                    endif
                    allocate (ncid(0:nnodes-1))
                    allocate (ncname(0:nnodes-1))
                    allocate (xi_start(0:nnodes-1))
                    allocate (eta_start(0:nnodes-1))
                    allocate (western_edge(0:nnodes-1))
                    allocate (eastern_edge(0:nnodes-1))
                    allocate (southern_edge(0:nnodes-1))
                    allocate (northern_edge(0:nnodes-1))
                    allocate (vid(maxvars,0:nnodes-1))
                    allocate (dimsize(maxdims,0:nnodes))
                  endif
                  complete=.false.
                  do node=0,nnodes-1
                    ncid(node)=-1
                    xi_start(node)=-1
                    eta_start(node)=-1
                  enddo
                elseif (nnodes.ne.ibuff(2)) then
                  write(*,'(/1x,2A,I4/14x,3A/14x,A,I4,4x,A/)')
     &                 '### WARNING: Number of MPI nodes in global ',
     &                 'attribute ''partition'', nnodes =', ibuff(2),
     &                 'in netCDF file ''',       nctestname(1:lncn),
     &                 ''' contradicts that from the initial',
     &                 'file in the sequence, nnodes =',      nnodes,
     &                                   ' ==> The file is ignored.'
                  arg=arg-1
                  goto 5
                endif
                node=ibuff(1)
                if (ncid(node).ne.-1) then
                  write(*,'(/1x,2A,I4,1x,A)') '### ERROR: netCDF ID ',
     &                         'for file corresponding to MPI-node =',
     &                                   node,  'is already in use.'
                  stop
                endif
                if (ncid(node).eq.-1 .and. xi_start(node).eq.-1
     &                         .and.  eta_start(node).eq.-1) then
                  ncid(node)=ntest
                  ncname(node)=nctestname
                  xi_start(node)=ibuff(3)
                  eta_start(node)=ibuff(4)
                  digit=.false.
                  is=0
                  ie=0
                  i=lncn+1
                  do while (is.eq.0 .and. i.gt.1)
                    i=i-1
                    if (nctestname(i:i).ge.'0' .and.
     &                  nctestname(i:i).le.'9') then
                      if (.not.digit) then
                        if (i.lt.lncn) then
                          if (nctestname(i+1:i+1).eq.'.') then
                            ie=i
                            digit=.true.
                          endif
                        else
                          ie=i
                          digit=.true.
                        endif
                      endif
                    elseif (digit .and. nctestname(i:i).eq.'.') then
                      digit=.false.
                      is=i+1
                    endif
                  enddo
                  if (is.gt.0 .and. ie.ge.is) then
                    root=nctestname(1:is-1)
                    if (ie.lt.lncn) then
                      sffx=nctestname(ie+1:lncn)
                    else
                      sffx(1:1)=' '
                    endif
                    k=0
                    do i=is,ie
                      k=10*k + ichar(nctestname(i:i))-48
                    enddo
                    code_size=ie-is+1
                  else
                    write(*,'(/1x,3A/)')        '### ERROR: Cannot ',
     &                  'determine MPI node number from filename ''',
     &                                     nctestname(1:lncn), '''.'
                  endif
                  ierr=nf_noerr
                  if (root_bak(1:1).eq.' ') then
                    root_bak=root
                  else
                    lvar=lenstr(root)
                    lbak=lenstr(root_bak)
                    if (lvar.ne.lbak .or. root.ne.root_bak) then
                      ierr=ierr+1
                      write(*,'(/8x,6A/17x,3A/)') 'WARNING: file ''',
     &                     nctestname(1:lncn),   ''' has different ',
     &                    'root name ''',  root(1:lvar),   ''' than',
     &                    'previously found root name ''',
     &                     root_bak(1:lbak), ''' from the same set.'
                    endif
                  endif
                  if (sffx_bak(1:1).eq.' ') then
                    sffx_bak=sffx
                  else
                    lvar=lenstr(sffx)
                    lbak=lenstr(sffx_bak)
                    if (lvar.ne.lbak .or. sffx.ne.sffx_bak) then
                      ierr=ierr+1
                      write(*,'(/8x,7A/17x,3A/)')       'WARNING: ',
     &                  'file ''',  nctestname(1:lncn),   ''' has ',
     &                  'different suffix name ''',    sffx(1:lvar),
     &                  ''' than','previously found suffix name ''',
     &                  sffx_bak(1:lbak),   ''' from the same set.'
                    endif
                  endif
                  if (code_size_bak.eq.-1) then
                    code_size_bak=code_size
                  elseif (code_size .ne. code_size_bak) then
                    ierr=ierr+1
                    write(*,'(/8x,A,I2,1x,A/17x,3A,I2,A/)')
     &              'WARNING: number of digits in MPI node segment',
     &               code_size, 'in filename', '''',
     &               nctestname(1:lncn),
     &                ''' is different than previously determined',
     &                                     code_size_bak, '.'
                  endif
                  if (k.ne.node) then
                    ierr=ierr+1
                    write(*,'(/8x,3A,I3/17x,2A/17x,A,I3,A/)')
     &                   'WARNING: file ''', nctestname(1:lncn),
     &                   ''' belongs to different MPI node',   node,
     &                   '(as determined from its global attribute',
     &                   '''partition'')', 'than node', k,
     &                   ' determined from to the file name.'
                  endif
                  if (ierr.ne.nf_noerr) goto 97
                else
                  arg=arg-1
                  goto 5
                endif
              else
                write(*,'(/1x,2A/14x,3A/)')   '### WARNING: Cannot ',
     &           'aquire global attribute ''partition'' from netCDF',
     &                                 'file ''', nctestname(1:lncn),
     &                               '''. ==> This file is ignored.'
              endif
            else
              write(*,'(/1x,2A/14x,3A/)')  '### WARNING: Wrong type ',
     &                'or size of global attribute ''partition'' in ',
     &                           'netCDF file ''', nctestname(1:lncn),
     &                                '''. ==> This file is ignored.'
            endif
          else
            write(*,'(/1x,3A/)') '### WARNING: ''', nctestname(1:lncn),
     &      ''' is not a partial netCDF file: ==> The file is ignored.'
          endif
        else
          write(*,'(/1x,4A/14x,A/)')   '### WARNING: Cannot open ''',
     &                   nctestname(1:lncn), ''' as a netCDF file: ',
     &                nf_strerror(ierr), ' ==> The file is ignored.'
        endif
   5    continue
        if (nnodes.gt.0) then
          complete=.true.
          do node=0,nnodes-1
            if (ncid(node).lt.0) complete=.false.
          enddo
        endif
      if (.not.complete .and. arg.lt.nargs) goto 2
      if (complete) then
        lncn=lenstr(ncname(0))
        write(*,'(2(1x,A,I4),1x,A,2x,A,2I5)')   'Processing set of ',
     &                         nnodes, 'files', 0, ncname(0)(1:lncn),
     &                          'i,jSW =', xi_start(0), eta_start(0)
        do node=1,nnodes-1
          if (node.lt.16 .or. (nnodes.gt.16 .and.
     &                         node.eq.nnodes-1 )) then
            write(*,'(29x,I4,1x,A,2x,A,2I5)') node,
     &                  ncname(node)(1:lncn), 'i,jSW =',
     &                  xi_start(node), eta_start(node)
          elseif (nnodes.gt.16 .and. node.lt.18) then
            write(*,'(24x,A)') '.................................'
          endif
        enddo
        if (ntest.ne.-1) then
          ierr=nf_close(ntest)
          ntest=-1
        endif
        do node=0,nnodes-1
          ncid(node)=-1
        enddo
      elseif (arg.lt.nargs) then
        goto 1
      else
        write(*,*) 'stop at 466'
        stop
      endif
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
        do node=0,nnodes-1
          lncn=lenstr(ncname(node))
          if (ncid(node).eq.-1) ierr=nf_open (ncname(node),
     &                               nf_nowrite, ncid(node))
          if (ierr .eq. nf_noerr) then
            ierr=nf_inq (ncid(node), ibuff(1), ibuff(2),
     &                               ibuff(3), ibuff(4))
            if (ierr .ne. nf_noerr) then
              write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot make ',
     &                        'general inquiry into netCDF file ''',
     &               ncname(node)(1:lncn), '''.', nf_strerror(ierr)
              goto 97
            elseif (ibuff(1) .gt. maxdims) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of dimensions', ibuff(1), 'in netCDF file ''',
     &          ncname(node)(1:lncn),    '''',   'exceeds limit.  ',
     &             'Increase parameter maxdims in file "ncjoin.F".'
              goto 97
            elseif (ibuff(2) .gt. maxvars) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of variables',  ibuff(2), 'in netCDF file ''',
     &          ncname(node)(1:lncn),  '''',      'exceeds limit. ',
     &             'Increase parameter maxvars in file "ncjoin.F".'
              goto 97
            elseif (node.eq.0) then
              ndims=ibuff(1)
              ngatts=ibuff(3)
              unlimdimid=ibuff(4)
            else
              if (ibuff(1) .ne. ndims) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                    'file ''', ncname(node)(1:lncn), ''' has ',
     &                    'different number of dimensions than ''',
     &                                       ncname(0)(1:lstr), '''.'
                ierr=ierr+1
              endif
              if (ibuff(3) .ne. ngatts) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                   'file ''',  ncname(node)(1:lncn), ''' has ',
     &               'different number of global attributes than ''',
     &                                       ncname(0)(1:lstr),'''.'
                ierr=ierr+1
              endif
              if (ibuff(4) .ne. unlimdimid) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                    'file ''', ncname(node)(1:lncn), ''' has ',
     &                'different ID for unlimited dimension than ''',
     &                                      ncname(0)(1:lstr), '''.'
                ierr=ierr+1
              endif
              if (ierr .ne. nf_noerr) goto 97
            endif
            do j=1,ibuff(1)
              ierr=nf_inq_dim (ncid(node),j,string,dimsize(j,node))
              if (ierr .eq. nf_noerr) then
                lstr=lenstr(string)
                if (node.eq.0) then
                  ldim(j)=lstr
                  dimname(j)=string(1:lstr)
                elseif (lstr.ne.ldim(j) .or. string(1:lstr).ne.
     &                                 dimname(j)(1:ldim(j)) )then
                  write(*,'(/1x,2A,I3,3A/12x,6A/12x,3A/)')    '### ',
     &                 'ERROR: Name of dimension #', j, ', named ''',
     &                  string(1:lstr),   ''' in netCDF file',  '''',
     &                  ncname(node)(1:lncn),   ''' does not match ',
     &                 'name ''',  dimname(j)(1:ldim(j)), ''' with ',
     &                 'the corresponding name from netCDF file ''',
     &                                  ncname(0)(1:lncn),    '''.'
                  goto 97
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &            'Cannot determine name and size of dimension #', j,
     &            'in netCDF file ''',   ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                goto 97
              endif
            enddo
            if (node.eq.0) nvars=0
            do i=1,ibuff(2)
              ierr=nf_inq_varname (ncid(node), i, string)
              if (ierr .eq. nf_noerr) then
                lstr=lenstr(string)
                ierr=nf_inq_varndims (ncid(node), i, k)
                if (ierr .eq. nf_noerr) then
                  lnewvar=.true.
                  do j=1,nvars
                    lvar=lenstr(vname(j))
                    if (lstr.eq.lvar .and. string(1:lstr)
     &                               .eq.vname(j)(1:lvar)) then
                      lnewvar=.false.
                      vid(j,node)=i
                      if (k.gt.vdims(j)) then
                        vdims(j)=k
                        vnode(j)=node
                      endif
                    endif
                  enddo
                  if (lnewvar) then
                    nvars=nvars+1
                    vname(nvars)=string(1:lstr)
                    vid(nvars,node)=i
                    vnode(nvars)=node
                    vdims(nvars)=k
                  endif
                else
                  write(*,'(/1x,3A,I3/12x,5A/12x,A)')  '### ERROR: ',
     &                      'Cannot determine number of dimensions ',
     &                      'for variable with id =', i,  'named ''',
     &                       string(1:lstr),  ''' in netCDF file ''',
     &                       ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                  goto 97
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &                 'Cannot determine name of variable with id =',
     &                  i, 'in netCDF file ''', ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
                goto 97
              endif
            enddo
            if (node.gt.0) then
              ierr=nf_close(ncid(node))
              ncid(node)=-1
            endif
          else
            write(*,'(/1x,A,1x,3A/14x,A)')    '### ERROR: Cannot ',
     &                 'open netCDF file ''', ncname(node)(1:lncn),
     &                                    '''.', nf_strerror(ierr)
            goto 97
          endif
        enddo
        id_xi_rho=0
        id_eta_rho=0
        id_xi_rho=0
        id_eta_v=0
        XI_rho=0
        ETA_rho=0
        size_XI=1
        size_ETA=1
        size_S=1
        tsize=1
        LLm=ibuff(5)
        MMm=ibuff(6)
        id_xi_u = 0
        do i=1,ndims
          dimsize(i,nnodes)=0
          lvar=lenstr(dimname(i))
          if (lvar.eq.6 .and. dimname(i)(1:lvar).eq.'xi_rho') then
            id_xi_rho=i
            do node=0,nnodes-1
                dimsize(i,nnodes)=LLm+2
                size_XI=max(size_XI,dimsize(i,node))
                XI_rho=max(XI_rho, dimsize(i,nnodes))
            enddo
          elseif (lvar.eq.4 .and.dimname(i)(1:lvar).eq.'xi_u') then
            id_xi_u=i
            do node=0,nnodes-1
              if (xi_start(node).gt.1) then
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                 dimsize(i,node) +xi_start(node)-2 )
              else
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                                  dimsize(i,node) )
              endif
                dimsize(i,nnodes)= LLm +1
                size_XI=max(size_XI,dimsize(i,node))
                XI_rho=max(XI_rho, dimsize(i,nnodes)+1)
            enddo
          elseif (lvar.eq.7.and. dimname(i)(1:lvar).eq.'eta_rho') then
            id_eta_rho=i
            do node=0,nnodes-1
              ETA_rho=max(ETA_rho, dimsize(i,nnodes))
                 dimsize(i,nnodes)=Mmm+2
                size_ETA=max(size_ETA,dimsize(i,node))
                ETA_rho=max(ETA_rho, dimsize(i,nnodes))
           enddo
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_v') then
            id_eta_v=i
            do node=0,nnodes-1
              if (eta_start(node).gt.1) then
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                 dimsize(i,node) +eta_start(node)-2 )
              else
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                                   dimsize(i,node))
              endif
                dimsize(i,nnodes)=Mmm +1
                size_ETA=max(size_ETA,dimsize(i,node))
                ETA_rho=max(ETA_rho, dimsize(i,nnodes)+1)
            enddo
          else
            dimsize(i,nnodes)=dimsize(i,0)
            do node=1,nnodes-1
              if (dimsize(i,0).ne.dimsize(i,node)) then
                lncn=lenstr(ncname(node))
                write(*,'(/1x,A,I3,3A,I4,1x,A/12x,4A/12x,3A,I4,A/)')
     &                 '### ERROR: Nonpartitionable dimension #',  i,
     &                 ' named ''', dimname(i)(1:lvar), ''', size =',
     &                  dimsize(i,node),   'in netCDF',    'file ''',
     &                  ncname(node)(1:lncn),    ''' has different ',
     &                 'size than the corresponding',
     &                 'dimension from file ''',   ncname(0)(1:lncn),
     &                     ''', which has size =', dimsize(i,0), '.'
                goto 97
              endif
            enddo
            if (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'s_rho') then
              size_S=max(size_S, dimsize(i,0))
            elseif (lvar.eq.3.and.dimname(i)(1:lvar).eq.'s_w') then
              size_S=max(size_S, dimsize(i,0))
            endif
          endif
          if (i.eq. unlimdimid) then
            tsize=dimsize(i,nnodes)
            dimsize(i,nnodes)=nf_unlimited
          endif
        enddo
        do node=0,nnodes-1
          western_edge(node)=.true.
          eastern_edge(node)=.true.
          southern_edge(node)=.true.
          northern_edge(node)=.true.
          if (xi_start(node).gt.1) then
            western_edge(node)=.false.
          endif
          if (id_xi_rho.gt.0) then
            if ( xi_start(node)+dimsize(id_xi_rho,node)
     &          .lt.XI_rho ) eastern_edge(node)=.false.
          endif
          if (id_xi_u.gt.0) then
            if ( xi_start(node)+dimsize(id_xi_u,node)
     &          .lt.XI_rho ) eastern_edge(node)=.false.
          endif
          if (eta_start(node).gt.1) then
            southern_edge(node)=.false.
          endif
          if (id_eta_rho.gt.0) then
            if ( eta_start(node)+dimsize(id_eta_rho,node)
     &           .lt.ETA_rho ) northern_edge(node)=.false.
          endif
          if (id_eta_v.gt.0) then
            if ( eta_start(node)+dimsize(id_eta_v,node)
     &          .lt.ETA_rho ) northern_edge(node)=.false.
          endif
        enddo
        max_buff_size=size_XI*size_ETA*size_S
        i=lenstr(root_bak)
        if (sffx_bak(1:1).ne.' ') then
          j=lenstr(sffx_bak)
          if (root_bak(i:i).eq.'.' .and. sffx_bak(1:1).eq.'.') then
            nctargname=root_bak(1:i)/ /sffx_bak(2:j)
          else
            nctargname=root_bak(1:i)/ /sffx_bak(1:j)
          endif
        else
          nctargname=root_bak(1:i)
        endif
        ltrg=lenstr(nctargname)
        ierr=nf_create (nctargname(1:ltrg), nf_clobber+nf_64bit_offset,
     &                                                          nctarg)
        if (ierr .eq. nf_noerr) then
          write(*,'(/1x,3A)')  'Created netCDF file ''',
     &                        nctargname(1:ltrg), '''.'
        else
          write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot create ',
     &                          'netCDF file ''', nctargname(1:ltrg),
     &                                      '''.', nf_strerror(ierr)
          goto 97
        endif
        size_XI=1
        size_ETA=1
        size_S=1
        do i=1,ndims
          lvar=lenstr(dimname(i))
          ierr=nf_def_dim (nctarg, dimname(i)(1:lvar),
     &                     dimsize(i,nnodes), dimid(i))
          if (ierr .eq. nf_noerr) then
            if (dimid(i) .eq. i) then
              if (dimname(i)(1:3) .eq. 'xi_') then
                size_XI=max(size_XI, dimsize(i,nnodes))
              elseif (dimname(i)(1:4) .eq. 'eta_') then
                size_ETA=max(size_ETA, dimsize(i,nnodes))
              elseif (dimname(i)(1:5) .eq. 's_rho' .or.
     &                dimname(i)(1:3) .eq. 's_w') then
                size_S=max(size_S, dimsize(i,nnodes))
              endif
            else
              write(*,'(/1x,2A,I3,1x,5A/12x,2A,I3,A/)')  '### ERROR: ',
     &        'id =', dimid(i), 'for dimension ''', dimname(i)(1:lvar),
     &        ''' from netCDF file ''',    nctargname(1:ltrg),    '''',
     &                    'differs from ', 'the original id =', i, '.'
              goto 97
            endif
          else
           write(*,'(/1x,4A/12x,A/)')      '### ERROR: Cannot define ',
     &    'dimension ''', dimname(i)(1:lvar), '''.', nf_strerror(ierr)
              goto 97
          endif
        enddo
        max_bfr_out=size_XI*size_ETA*size_S
        lncn=lenstr(ncname(0))
        if (ncid(0).eq.-1) ierr=nf_open (ncname(0), nf_nowrite,
     &                                                 ncid(0))
        if (ierr .eq. nf_noerr) then
          do i=1,ngatts
            ierr=nf_inq_attname (ncid(0), nf_global, i, string)
            if (ierr. eq. nf_noerr) then
              lvar=lenstr(string)
              if (string(1:lvar) .ne. 'partition') then
                ierr=nf_copy_att (ncid(0), nf_global, string(1:lvar),
     &                                             nctarg, nf_global)
                if (ierr .ne. nf_noerr) then
                  write(*,'(/1x,4A/12x,3A/12x,A)')     '### ERROR: ',
     &             'Cannot copy global attribute ''', string(1:lvar),
     &             ''' into netCDF',  'file ''',  nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
                  goto 97
                endif
              endif
            else
            write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: Cannot',
     &                    ' determine name of global attribute #', i,
     &                     'from netCDF file ''',  ncname(0)(1:lncn),
     &                                     '''.',  nf_strerror(ierr)
              goto 97
            endif
          enddo
        else
          write(*,'(/1x,A,1x,3A/14x,A)')     '### ERROR: Cannot open ',
     &   'netCDF file ''', ncname(0)(1:lncn), '''.', nf_strerror(ierr)
          goto 97
        endif
        do i=1,nvars
          node=vnode(i)
          lncn=lenstr(ncname(node))
          if (ncid(node).eq.-1) ierr=nf_open  (ncname(node),
     &                               nf_nowrite, ncid(node))
          if (ierr .eq. nf_noerr) then
            ierr=nf_inq_var (ncid(node), vid(i,node), vname(i),
     &                vartype(i), vdims(i), dimids(1,i),  varatts)
            if (ierr .eq. nf_noerr) then
              lvar=lenstr(vname(i))
              ierr=nf_def_var (nctarg, vname(i)(1:lvar),vartype(i),
     &                            vdims(i), dimids(1,i), varid(i))
              if (ierr .eq. nf_noerr) then
                do j=1,varatts
                  ierr=nf_inq_attname (ncid(node), vid(i,node),
     &                                              j, string)
                  if (ierr .eq. nf_noerr) then
                    lstr=lenstr(string)
                    ierr=nf_copy_att (ncid(node), vid(i,node),
     &                       string(1:lstr), nctarg, varid(i))
                    if (ierr. ne. nf_noerr) then
                      write(*,'(/1x,2A,I3,3A/12x,4A)')   '### ERROR: ',
     &                 'Cannot copy attribute #', j,' for variable ''',
     &                  vname(i)(1:lvar),  ''' into netCDF', 'file ''',
     &                  nctargname(1:ltrg), '''.  ', nf_strerror(ierr)
                      goto 97
                    endif
                  else
                    write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: ',
     &                             'Cannot get name of attribute #', j,
     &                      'for variable ''', vname(i)(1:lvar), '''.',
     &                                               nf_strerror(ierr)
                    goto 97
                  endif
                enddo
              else
                write(*,'(/8x,4A/)') 'ERROR: Cannot define ',
     &                  'variable ''', vname(i)(1:lvar), '''.'
                goto 97
              endif
            else
            write(*,'(/8x,2A/15x,A,I3,1x,3A/)')  '### ERROR: Cannot ',
     &        'determine name, type and attributes for variable #', i,
     &            'from netCDF file ''', ncname(node)(1:lncn),  '''.'
              goto 97
            endif
            series(i)=.false.
            part_type(i)=0
            do j=1,vdims(i)
              if (dimids(j,i).eq.id_xi_rho .or.
     &            dimids(j,i).eq.id_xi_u) then
                part_type(i)=part_type(i)+1
              elseif (dimids(j,i).eq.id_eta_rho .or.
     &                dimids(j,i).eq.id_eta_v) then
                part_type(i)=part_type(i)+2
              elseif (dimids(j,i).eq.unlimdimid) then
                series(i)=.true.
              endif
            enddo
            if (node.gt.0) then
              ierr=nf_close(ncid(node))
              ncid(node)=-1
            endif
          else
            write(*,'(/1x,A,1x,3A/12x,A)')  '### ERROR: Cannot open ',
     &                  'netCDF file ''', ncname(node)(1:lncn), '''.',
     &                                              nf_strerror(ierr)
            goto 97
          endif
        enddo
        ierr=nf_enddef (nctarg)
        if (max_bfr_out .gt. alloc_bfr_out) then
          if (allocated(bfr_out)) deallocate(bfr_out)
        endif
        if (max_buff_size .gt. alloc_buff_size) then
          if (allocated(buff)) deallocate(buff)
          allocate(buff(max_buff_size))
          alloc_buff_size=max_buff_size
          write(*,*) 'allocated "buff" with  max_buff_size =',
     &                                         max_buff_size
        endif
        if (max_bfr_out .gt. alloc_bfr_out) then
          allocate(bfr_out(max_bfr_out))
          alloc_bfr_out=max_bfr_out
          write(*,*) 'allocated "bfr_out" with ',
     &               'max_bfr_out =', max_bfr_out
        endif
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_fcrt_clk=net_fcrt_clk+inc_clk
        do rec=1,tsize
          if (tsize.gt.1) then
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate,clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_gray_clk=net_gray_clk+inc_clk
            write(*,'(F8.1,1x,A,I8,1x,A,I8,1x,A)')
     &         dble(iclk(nclk)-iclk_init)/dble(clk_rate),
     &     'Processing record', rec, 'out of', tsize,  '...'
          endif
          do i=1,nvars
            if (rec.eq.1 .or. series(i)) then
              if (part_type(i).eq.0 .and. .not.series(i)) then
                lvar=lenstr(vname(i))
                write(*,'(16x,3A)') 'Copy scalar variable: ''',
     &                               vname(i)(1:lvar), '''...'
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_var_text  (ncid(0), vid(i,0), buff_str)
                elseif (vartype(i) .eq. nf_byte) then
                  ierr=nf_get_var_int1   (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_short) then
                  ierr=nf_get_var_int2   (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_var_int    (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_var_real   (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_var_double (ncid(0), vid(i,0), buff)
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A/)') '### ERROR: scalar variable ',
     &              '''', vname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (ierr .eq. nf_noerr) then
                  if (vartype(i) .eq. nf_char) then
                    ierr=nf_put_var_text  (nctarg, varid(i), buff_str)
                  elseif (vartype(i) .eq. nf_byte) then
                    ierr=nf_put_var_int1   (nctarg, varid(i), buff)
                  elseif (vartype(i) .eq. nf_short) then
                    ierr=nf_put_var_int2   (nctarg, varid(i), buff)
                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_put_var_int    (nctarg, varid(i), buff)
                  elseif (vartype(i) .eq. nf_float) then
                    ierr=nf_put_var_real   (nctarg ,varid(i), buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_put_var_double (nctarg, varid(i), buff)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/1x,4A/12x,3A/12x,A)')  '### ERROR: ',
     &                            'Cannot write scalar variable ''',
     &                           vname(i)(1:lvar), ''' into netCDF',
     &                                'file ''', nctargname(1:ltrg),
     &                                    '''.',  nf_strerror(ierr)
                    goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot ',
     &                 'read scalar variable ''', vname(i)(1:lvar),
     &                                    '''.', nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_type(i).eq.0) then
                lvar=lenstr(vname(i))
                write(*,'(16x,3A)') 'Copy non-partitioned array: ''',
     &                                    vname(i)(1:lvar), '''...'
                size=1
                do j=1,vdims(i)
                  if (dimids(j,i).eq.unlimdimid) then
                    start(j)=rec
                    count(j)=1
                  else
                    start(j)=1
                    count(j)=dimsize(dimids(j,i),0)
                  endif
                  size=size*count(j)
                enddo
                if (vartype(i) .eq. nf_char .or.
     &              vartype(i) .eq. nf_byte) then
                  size=size*1
                elseif (vartype(i) .eq. nf_short) then
                  size=size*2
                elseif (vartype(i) .eq. nf_int .or.
     &                  vartype(i) .eq. nf_float) then
                  size=size*4
                elseif (vartype(i) .eq. nf_double) then
                  size=size*8
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,3A/)')  '### ERROR: variable ''',
     &                  vname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (size .gt. 8*max_buff_size) then
                  if (allocated(buff)) deallocate(buff)
                  max_buff_size=(size+7)/8
                  allocate(buff(max_buff_size))
                  write(*,*) 'allocated "buff" with max_buff_size =',
     &                                              max_buff_size
                endif
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_vara_text  (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_byte) then
                  ierr=nf_get_vara_int1   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_short) then
                  ierr=nf_get_vara_int2   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_vara_int    (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_vara_real   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_vara_double (ncid(0), vid(i,0),
     &                                    start,count, buff)
                endif
                if (ierr .eq. nf_noerr) then
                  if (vartype(i) .eq. nf_char) then
                    ierr=nf_put_vara_text  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_byte) then
                    ierr=nf_put_vara_int1  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_short) then
                    ierr=nf_put_vara_int2  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_put_vara_int   (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_float) then
                    ierr=nf_put_vara_real  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_put_vara_double(nctarg, varid(i),
     &                                    start,count, buff)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/8x,4A,I3/15x,3A/)')    '### ERROR: ',
     &                 'Cannot write variable ''', vname(i)(1:lvar),
     &              ''' for time record',rec, 'into netCDF file ''',
     &                  nctargname(1:ltrg),'''.', nf_strerror(ierr)
                    goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A,I3,A/15x,A/)')     '### ERROR: ',
     &              'Cannot read variable ''',     vname(i)(1:lvar),
     &              ''' for time record',rec,'.', nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_type(i).gt.0) then
                lvar=lenstr(vname(i))
                write(*,'(16x,2A,I3,1x,3A)')  'Assembly partitioned ',
     &                       'array type', part_type(i),  'name = ''',
     &                                         vname(i)(1:lvar), ''''
                bfr_out=0.D0
                do node=0,nnodes-1
                  var_mask=.true.
                  if (part_type(i).eq.1 .and. lvar.gt.6) then
                    if (vname(i)(lvar-5:lvar).eq.'_south'
     &                        .and. southern_edge(node)) then
                      var_mask=.true.
                    elseif (vname(i)(lvar-5:lvar).eq.'_north'
     &                           .and. northern_edge(node)) then
                      var_mask=.true.
                    endif
                  elseif (part_type(i).eq.2 .and. lvar.gt.5) then
                    if (vname(i)(lvar-4:lvar).eq.'_west'
     &                                .and. western_edge(node)) then
                      var_mask=.true.
                    elseif (vname(i)(lvar-4:lvar).eq.'_east'
     &                          .and. eastern_edge(node)) then
                      var_mask=.true.
                    endif
                  elseif (part_type(i).eq.3) then
                    var_mask=.true.
                  endif
                  if (var_mask) then
                    size=1
                    size1=1
                    do j=1,vdims(i)
                      k=dimids(j,i)
                      if (k.eq.id_xi_rho  .or. k.eq.id_xi_u  .or.
     &                    k.eq.id_eta_rho .or. k.eq.id_eta_v) then
                        start(j)=1
                        count(j)=dimsize(k,node)
                        if (k.eq.id_xi_rho) then
                          start1(j)=xi_start(node)
                          count1(j)=XI_rho
                        elseif (k.eq.id_xi_u) then
                          start1(j)=max(xi_start(node)-1,1)
                          count1(j)=XI_rho-1
                        elseif (k.eq.id_eta_rho) then
                          start1(j)=eta_start(node)
                          count1(j)=ETA_rho
                        elseif (k.eq.id_eta_v) then
                          start1(j)=max(eta_start(node)-1,1)
                          count1(j)=ETA_rho-1
                        endif
                      elseif (k.eq.unlimdimid) then
                        start(j)=rec
                        count(j)=1
                        start1(j)=rec
                        count1(j)=1
                      else
                        start(j)=1
                        count(j)=dimsize(k,nnodes)
                        start1(j)=1
                        count1(j)=count(j)
                      endif
                      size=size*count(j)
                      size1=size1*count1(j)
                    enddo
                    if (vartype(i) .eq. nf_char .or.
     &                  vartype(i) .eq. nf_byte) then
                      size=size*1
                      size1=size1*1
                    elseif (vartype(i) .eq. nf_short) then
                      size=size*2
                      size1=size1*2
                    elseif (vartype(i) .eq. nf_int .or.
     &                      vartype(i) .eq. nf_float) then
                      size=size*4
                      size1=size1*4
                    elseif (vartype(i) .eq. nf_double) then
                      size=size*8
                      size1=size1*8
                    else
                      lvar=lenstr(vname(i))
                      write(*,'(/8x,4A/)') '### ERROR: variable ''',
     &                     vname(i)(1:lvar), ''' has unknown type.'
                      goto 97
                    endif
                    if (size .gt. 8*max_buff_size) then
                      if (allocated(buff)) deallocate(buff)
                      max_buff_size=(size+7)/8
                      allocate(buff(max_buff_size))
                      write(*,*) 'allocated "buff" with ',
     &                   'max_buff_size =', max_buff_size
                    endif
                    if (ncid(node).eq.-1) ierr=nf_open (ncname(node),
     &                                       nf_nowrite, ncid(node))
                    if (ierr.eq.nf_noerr) then
                      nclk=3-nclk
                      call system_clock (iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      net_gray_clk=net_gray_clk+inc_clk
                      if (vartype(i) .eq. nf_char) then
                        ierr=nf_get_vara_text (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_byte) then
                        ierr=nf_get_vara_int1 (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_short) then
                        ierr=nf_get_vara_int2 (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_int) then
                        ierr=nf_get_vara_int  (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_float) then
                        ierr=nf_get_vara_real (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_double) then
                        ierr=nf_get_vara_double(ncid(node),vid(i,node),
     &                                             start, count, buff)
                      endif
                      if (ierr.eq.nf_noerr) then
                        net_read_size=net_read_size+size
                        nclk=3-nclk
                        call system_clock (iclk(nclk),clk_rate,clk_max)
                        inc_clk=iclk(nclk)-iclk(3-nclk)
                        net_read_clk=net_read_clk+inc_clk
                      else
                        lvar=lenstr(vname(i))
                        lncn=lenstr(ncname(node))
                        write(*,'(/1x,4A,I3/15x,3A/15x,A/)')  '### ',
     &                              'ERROR: Cannot read variable ''',
     &                   vname(i)(1:lvar), ''' for time record', rec,
     &                  'from netCDF file ''',  ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
                        goto 97
                      endif
                      if (node.gt.0) then
                        ierr=nf_close(ncid(node))
                        ncid(node)=-1
                      endif
                    else
                      lncn=lenstr(ncname(node))
                      write(*,'(/1x,A,1x,3A/14x,A)')   '### ERROR: ',
     &                                  'Cannot open netCDF file ''',
     &                ncname(node)(1:lncn), '''.', nf_strerror(ierr)
                      goto 97
                    endif
                    if (size1 .gt. 8*max_bfr_out) then
                      if (allocated(bfr_out)) deallocate(bfr_out)
                      max_bfr_out=(size1+7)/8
                      allocate(bfr_out(max_bfr_out))
                      write(*,*) 'allocated "bfr_out" with ',
     &                           'max_bfr_out =', max_bfr_out
                    endif
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_gray_clk=net_gray_clk+inc_clk
                    if (vartype(i) .eq. nf_char) then
                      call assembly_text (buff, count,
     &                                 bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_byte) then
                      call assembly_byte (buff, count,
     &                                 bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_short) then
                      call assembly_int2 (buff, count,
     &                                    bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_int) then
                      call assembly_int  (buff, count,
     &                                    bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_float) then
                      call assembly_real (buff, count,
     &                                    bfr_out,
     &                                 start1,count1, vdims(i))
                     elseif (vartype(i) .eq. nf_double) then
                      call assembly_double (buff,count,bfr_out,
     &                                 start1,count1, vdims(i))
                    endif
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_assm_clk=net_assm_clk+inc_clk
                  endif
                  if (node .eq. nnodes-1) then
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_gray_clk=net_gray_clk+inc_clk
                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_byte) then
                      ierr=nf_put_vara_int1   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_short) then
                      ierr=nf_put_vara_int2   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_float) then
                      ierr=nf_put_vara_real   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_put_vara_double (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    endif
                    if (ierr.eq.nf_noerr) then
                      net_wrt_size=net_wrt_size+size1
                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      net_wrt_clk=net_wrt_clk+inc_clk
                    else
                      lvar=lenstr(vname(i))
                      lncn=lenstr(vname(i))
                      write(*,'(/1x,3A,I3/12x,3A/12x,A/)')
     &                         '### ERROR: Cannot write variable ''',
     &                     vname(i)(1:lvar),''' for time record',rec,
     &                    'into netCDF file ''',  nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
                      goto 97
                    endif
                  endif
                enddo
              endif
            endif
          enddo
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate, clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          net_gray_clk=net_gray_clk+inc_clk
          ierr=nf_sync (nctarg)
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate, clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          net_sync_clk=net_sync_clk+inc_clk
        enddo
        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false.
  98    write(*,*) 'closing files...'
        do node=0,nnodes-1
          if (ncid(node).ne.-1) then
            ierr=nf_close(ncid(node))
            ncid(node)=-1
          endif
        enddo
        write(*,*) '...........input'
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
        ierr=nf_close (nctarg)
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_sync_clk=net_sync_clk+inc_clk
        write(*,*) '...........output'
        if (del_part_files) then
          if (clean_set) then
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate, clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_gray_clk=net_gray_clk+inc_clk
            write(*,'(/1x,A)') 'Deleting partial files...'
            do node=0,nnodes-1
              rmcmd='/bin/rm -f '/ /ncname(node)
              lstr=lenstr(rmcmd)
              if (node.lt.16 .or. (nnodes.gt.16 .and.
     &                         node.eq.nnodes-1 )) then
                write(*,'(27x,3A)')  '''', rmcmd(1:lstr), ''''
              elseif (nnodes.gt.16 .and. node.lt.18) then
                write(*,'(24x,A)') '.................................'
              endif
              call system (rmcmd(1:lstr))
            enddo
            write(*,*)
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate, clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_rmcmd_clk=net_rmcmd_clk+inc_clk
          else
            write(*,'(/1x,2A/)')  '### ERROR: Not removing ',
     &                     'partial files because of errors.'
          endif
        endif
      if (arg .lt. nargs)  goto 1
      call etime(CPU_time, RUN_time)
      RUN_time=RUN_time-tstart
      write(*,'(/3(1x,A,F11.2,1x))') 'CPU_time:  run =', RUN_time,
     &                   'usr =', CPU_time(1),  'sys =', CPU_time(2)
      if (clk_rate.gt.0) then
        ReadSize=1.0D-6*net_read_size
        WrtSize=1.0D-6*net_wrt_size
        ReadTime=net_read_clk/dble(clk_rate)
        AssmTime=net_assm_clk/dble(clk_rate)
        WrtTime = net_wrt_clk/dble(clk_rate)
        SyncTime=net_sync_clk/dble(clk_rate)
        FcrtTime=net_fcrt_clk/dble(clk_rate)
        write(*,'(/1x,A,22x,F12.2,1x,A)') 'Analysis/file creation :',
     &                                               FcrtTime, 'sec'
        write(*,'(8x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &         'Total data read :', ReadSize, 'MBytes in',  ReadTime,
     &                          'sec (', ReadSize/ReadTime, 'MB/sec)'
        write(*,'(5x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &      'Total data written :', WrtSize,  'MBytes in',   WrtTime,
     &                          'sec (',  WrtSize/WrtTime,  'MB/sec)'
        write(*,'(5x,A,22x,F12.2,1x,A)')      'Data assembly time :',
     &                                               AssmTime, 'sec'
        write(*,'(2x,A,22x,F12.2,1x,A)')   'Output file sync time :',
     &                                               SyncTime, 'sec'
        if (del_part_files) then
          write(*,'(1x,A,22x,F12.2,1x,A)') 'Removing partial files :',
     &                          net_rmcmd_clk/dble(clk_rate), 'sec'
        endif
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
        GrayTime=net_gray_clk/dble(clk_rate)
        write(*,'(14x,A,22x,F12.2,1x,A)') 'Gray time :',GrayTime,'sec'
        inc_clk=iclk(nclk)-iclk_init
        write(*,'(47x,A/12x,A,11x,F12.2,1x,A/)') '------------------',
     &   'Elapsed wall-clock time:', inc_clk/dble(clk_rate), 'sec'
      endif
      stop
      end
      subroutine assembly_text (buff, count,  BFr_out,
     &                          start1, count1, vdims)
      implicit none
      character(len=1) buff(*), bfr_out(*)
      integer*4 vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)
      istr=start1(1)
      imax1=count1(1)
      if (ndims.gt.1) then
        jmax=count(2)
        jstr=start1(2)
        jmax1=count1(2)
      else
        jstr=1
        jmax=1
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine assembly_byte (buff,  count, bfr_out,
     &                          start1, count1, vdims)
      implicit none
      integer(kind=1) ::  buff(*), bfr_out(*)
      integer*4 vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)
      istr=start1(1)
      imax1=count1(1)
      if (ndims.gt.1) then
        jmax=count(2)
        jstr=start1(2)
        jmax1=count1(2)
      else
        jstr=1
        jmax=1
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine assembly_int2 (buff, count,  bfr_out,
     &                          start1, count1, vdims)
      implicit none
      integer(kind=2) :: buff(*), bfr_out(*)
      integer*4 vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)
      istr=start1(1)
      imax1=count1(1)
      if (ndims.gt.1) then
        jmax=count(2)
        jstr=start1(2)
        jmax1=count1(2)
      else
        jstr=1
        jmax=1
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine assembly_int  (buff, count,  bfr_out,
     &                          start1, count1, vdims)
      implicit none
      integer(kind=4) :: buff(*), bfr_out(*)
      integer*4 vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)
      istr=start1(1)
      imax1=count1(1)
      if (ndims.gt.1) then
        jmax=count(2)
        jstr=start1(2)
        jmax1=count1(2)
      else
        jstr=1
        jmax=1
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine assembly_real (buff, count,  bfr_out,
     &                          start1, count1, vdims)
      implicit none
      real(kind=4) :: buff(*), bfr_out(*)
      integer*4 vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)
      istr=start1(1)
      imax1=count1(1)
      if (ndims.gt.1) then
        jmax=count(2)
        jstr=start1(2)
        jmax1=count1(2)
      else
        jstr=1
        jmax=1
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine assembly_double (buff,count, bfr_out,
     &                          start1, count1, vdims)
      implicit none
      real(kind=8) :: buff(*), bfr_out(*)
      integer*4 vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1
      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)
      istr=start1(1)
      imax1=count1(1)
      if (ndims.gt.1) then
        jmax=count(2)
        jstr=start1(2)
        jmax1=count1(2)
      else
        jstr=1
        jmax=1
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif
      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end
