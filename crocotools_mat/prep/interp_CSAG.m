function interp_CSAG(CSAG_dir,Roa,interp_method,...
                     lon1,lat1,tin,indx_var,indx_rain,...
		     nc_frc,nc_blk,lon,lat,angle)

% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%
% Based on interp_CFSR.m
%
%
% 1: Air temperature: Convert from Kelvin to Celsius
%
vname='tas';
nc=netcdf([CSAG_dir,'wrfout_d03_',vname,'_regridded_sub.nc']);
tair=squeeze(nc{vname}(indx_var(tin),:,:));
close(nc);
tair=get_missing_val(lon1,lat1,tair,nan,Roa,nan);
tair=tair-273.15;
tair=interp2(lon1,lat1,tair,lon,lat,interp_method);
%
% 2: Relative humidity: Convert from % to fraction
%
% Get Specific Humidity [Kg/Kg]
%
vname='huss';
nc=netcdf([CSAG_dir,'wrfout_d03_',vname,'_regridded_sub.nc']);
shum=squeeze(nc{vname}(indx_var(tin),:,:));
close(nc);
shum=get_missing_val(lon1,lat1,shum,nan,Roa,nan);
shum=interp2(lon1,lat1,shum,lon,lat,interp_method);
%
% computes specific humidity at saturation (Tetens  formula)
% (see air_sea tools, fonction qsat)
%
rhum=shum./qsat(tair);
%
% 3: Precipitation rate: Convert from mm/hr to cm/day
%
vname='rain';
nc=netcdf([CSAG_dir,'wrfout_d03_',vname,'_regridded_sub.nc']);
prate=squeeze(nc{vname}(indx_rain(tin),:,:));
close(nc);
prate=get_missing_val(lon1,lat1,prate,nan,Roa,nan);
prate=prate*0.1*24.0;
prate=interp2(lon1,lat1,prate,lon,lat,interp_method);
prate(abs(prate)<1.e-4)=0;
%
% 4: Net shortwave flux: [W/m^2]
%      CROCO convention: downward = positive
%
% Downward solar shortwave
%
vname='rsds';
nc=netcdf([CSAG_dir,'wrfout_d03_',vname,'_regridded_sub.nc']);
dswrf=squeeze(nc{vname}(indx_var(tin),:,:));
close(nc);
dswrf=get_missing_val(lon1,lat1,dswrf,nan,Roa,nan);
%  
% Upward solar shortwave - not provided by CSAG so we have to estimate an
% albedo to calculate the net radiation. Comparing with CFSR downward and
% upward short wave radiation fluxes, an albedo of 6% is considered
% reasonable and in good agreement with CFSR values (when downward
% radiation is highest i.e. excluding early morning and late evening)
%
albedo=6./100;
%  
%  Net solar shortwave radiation  
%
radsw=(1-albedo)*dswrf;
%----------------------------------------------------  
% GC le 31 03 2009
%  radsw is NET solar shortwave radiation
%  no more downward only solar radiation
% GC  bug fix by F. Marin IRD/LEGOS
%-----------------------------------------------------
radsw=interp2(lon1,lat1,radsw,lon,lat,interp_method);
radsw(radsw<1.e-10)=0;
%
% 5: Net outgoing Longwave flux:  [W/m^2]
%      CROCO convention: positive upward (opposite to nswrf !!!!)
%
% Get the net longwave flux [W/m^2]
%
%  5.1 get the downward longwave flux [W/m^2]
%
vname='rlds';
nc=netcdf([CSAG_dir,'wrfout_d03_',vname,'_regridded_sub.nc']);
dlwrf=squeeze(nc{vname}(indx_var(tin),:,:));
close(nc);
dlwrf=get_missing_val(lon1,lat1,dlwrf,nan,Roa,nan);
radlw_in=interp2(lon1,lat1,dlwrf,lon,lat,interp_method);
%
%  5.2 get the upward longwave flux [W/m^2]
%
%  This parameter is not provided by CSAG, but we don't need it as we
%  define BULK_LW for all runs, which means that the outward longwave flux
%  is calculated using the model SST and Stefan's law
%
% 6: Wind & Wind stress [m/s]
%
vname='u10';
nc=netcdf([CSAG_dir,'wrfout_d03_wind_10m_regridded_sub.nc']);
uwnd=squeeze(nc{vname}(indx_var(tin),:,:));
uwnd=get_missing_val(lon1,lat1,uwnd,nan,Roa,nan);
uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
%
vname='v10';
vwnd=squeeze(nc{vname}(indx_var(tin),:,:));
close(nc)
vwnd=get_missing_val(lon1,lat1,vwnd,nan,Roa,nan);
vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
%
% Compute the stress
%
wspd=sqrt(uwnd.^2+vwnd.^2);
[Cd,uu]=cdnlp(wspd,10.);
rhoa=air_dens(tair,rhum*100);
tx=Cd.*rhoa.*uwnd.*wspd;
ty=Cd.*rhoa.*vwnd.*wspd;
%
% Rotations on the CROCO grid
%
cosa=cos(angle);
sina=sin(angle);
%
sustr=rho2u_2d(tx.*cosa+ty.*sina);
svstr=rho2v_2d(ty.*cosa-tx.*sina);
%
% uwnd et vwnd sont aux points 'rho'
%
u10=rho2u_2d(uwnd.*cosa+vwnd.*sina);
v10=rho2v_2d(vwnd.*cosa-uwnd.*sina);
%
% Fill the CROCO files
%
if ~isempty(nc_frc)
  nc_frc{'sustr'}(tin,:,:)=sustr;
  nc_frc{'svstr'}(tin,:,:)=svstr;
end
if ~isempty(nc_blk)
  nc_blk{'tair'}(tin,:,:)=tair;
  nc_blk{'rhum'}(tin,:,:)=rhum;
  nc_blk{'prate'}(tin,:,:)=prate;
  nc_blk{'wspd'}(tin,:,:)=wspd;
%  nc_blk{'radlw'}(tin,:,:)=radlw;
  nc_blk{'radlw_in'}(tin,:,:)=radlw_in;
  nc_blk{'radsw'}(tin,:,:)=radsw;
  nc_blk{'uwnd'}(tin,:,:)=u10;
  nc_blk{'vwnd'}(tin,:,:)=v10;
%  nc_blk{'sustr'}(tout,:,:)=sustr;
%  nc_blk{'svstr'}(tout,:,:)=svstr;
end


