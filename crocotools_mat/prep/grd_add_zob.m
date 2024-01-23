% Add the roughness length scale to the grid file
% This can be spatially varing, but for now we just assign it to a 
% constant value
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts

clear all
close all
crocotools_param
zob_const = 0.001;

nc=netcdf(grdname,'write');

h  = nc{'h'}(:);
zob_var = h.*0 + zob_const;

nc{'z0b'} = ncdouble('eta_rho', 'xi_rho');
nc{'z0b'}.long_name = ncchar('roughness length scale');
nc{'z0b'}.long_name = 'roughness length scale';
nc{'z0b'}.units = ncchar('meter');
nc{'z0b'}.units = 'meter';

nc{'z0b'}(:)=zob_var;

close(nc)

