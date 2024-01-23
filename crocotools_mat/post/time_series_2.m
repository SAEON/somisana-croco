function ts=time_series_2(hisfile,gridfile,lon0,lat0,vname,vlevel,coef,ref_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  As per time_series.m in croco_tools, only we write return a timeseries
%  object
%
%  Get a vertical Profile
%
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Defaults values
%
if nargin < 1
  error('You must specify a file name')
end
if nargin < 2
  gridfile=hisfile;
  disp(['Default grid name: ',gridfile])
end
if nargin < 3
  lon0=[];
end
if nargin < 4
  lat0=[];
end
if nargin < 5
  vname='temp';
  disp(['Default variable to plot: ',vname])
end
if nargin < 6
  vlevel=-10;
  disp(['Default vertical level: ',num2str(vlevel)])
end
if nargin < 7
  coef=1;
  disp(['Default coef: ',num2str(coef)])
end
%
% Get default values
%
if isempty(gridfile)
  gridfile=hisfile;
end
if vname(1)=='u'
  [lat,lon,mask]=read_latlonmask(gridfile,'u');
elseif vname(1)=='v'
  [lat,lon,mask]=read_latlonmask(gridfile,'v');
else
  [lat,lon,mask]=read_latlonmask(gridfile,'r');
end
if isempty(lon0) | isempty(lat0)
  lat0=mean(mean(lat));
  lon0=mean(mean(lon));
end
%
% Find j,i indices for the profile
%
disp(['lon0 = ',num2str(lon0),' - lat0 = ',num2str(lat0)])
[J,I]=find((lat(1:end-1,1:end-1)<=lat0 & lat(2:end,2:end)>lat0 &...
            lon(2:end,1:end-1)<=lon0 & lon(1:end-1,2:end)>lon0)==1);
if isempty(I) |  isempty(J)
  disp('Warning: profile place not found')
  [M,L]=size(lon);
  I=round(L/2);
  J=round(M/2);
end
% GGF hack edit needed in the case where two points are found
J=J(1);
I=I(1);
disp(['I = ',int2str(I),' J = ',int2str(J)])
lon1=lon(J,I);
lat1=lat(J,I);
disp(['lon1 = ',num2str(lon1),' - lat1 = ',num2str(lat1)])
%
% get the vertical levels
%
nc=netcdf(hisfile);
if strcmp(vname,'zeta') | strcmp(vname,'ubar') | strcmp(vname,'vbar') | ...
   strcmp(vname,'sustr') | strcmp(vname,'svstr') | strcmp(vname,'shflux') | ...
   strcmp(vname,'swflux') | strcmp(vname,'shflx_rsw') | strcmp(vname,'swrad') | ...
   strcmp(vname,'shflx_rlw') | strcmp(vname,'shflx_sen') | strcmp(vname,'shflx_lat') | ...
   strcmp(vname,'sst_skin')
  var=coef*squeeze(nc{vname}(:,J,I));
elseif vlevel>0
%
% Read the variable
%
  if vname(1)=='*'
    if strcmp(vname,'*Ke')
      u=mean(squeeze(nc{'u'}(:,vlevel,J,I-1:I)),2);
      v=mean(squeeze(nc{'v'}(:,vlevel,J-1:J,I)),2);
      var=coef.*0.5.*(u.^2+v.^2);
    elseif strcmp(vname,'*Speed')
      u=mean(squeeze(nc{'u'}(:,vlevel,J,I-1:I)),2);
      v=mean(squeeze(nc{'v'}(:,vlevel,J-1:J,I)),2);
      var=coef.*sqrt(u.^2+v.^2);
    elseif strcmp(vname,'*Rho_pot')
      temp=squeeze(nc{'temp'}(:,vlevel,J,I));
      salt=squeeze(nc{'salt'}(:,vlevel,J,I));
      var=coef*rho_pot(temp,salt);
    elseif strcmp(vname,'*Chla')
      sphyto=squeeze(nc{'SPHYTO'}(:,vlevel,J,I));
      lphyto=squeeze(nc{'LPHYTO'}(:,vlevel,J,I));
      theta_m  =0.020;
      CN_Phyt  = 6.625;
      var=coef*theta_m*(sphyto+lphyto)*CN_Phyt*12.;
      var(var<=0)=NaN;
    else
      disp('Sorry not implemented yet')
      return
    end 
  else
    var=coef*squeeze(nc{vname}(:,vlevel,J,I));
  end
else
  ng=netcdf(gridfile);
  if vname(1)=='u'
    zeta=mean(squeeze(nc{'zeta'}(:,J,I:I+1)),2);
    h=mean(squeeze(ng{'h'}(J,I:I+1)),2);
  elseif vname(1)=='v'
    zeta=mean(squeeze(nc{'zeta'}(:,J:J+1,I)),2);
    h=mean(squeeze(ng{'h'}(J:J+1,I)),1);
  else
    zeta=squeeze(nc{'zeta'}(:,J,I));
    h=squeeze(ng{'h'}(J,I));
  end
  close(ng)
  theta_s=nc.theta_s(:);
  if (isempty(theta_s))
%    disp('Rutgers version')
    theta_s=nc{'theta_s'}(:);
    theta_b=nc{'theta_b'}(:);
    Tcline=nc{'Tcline'}(:);
  else 
%    disp('UCLA version');
    theta_b=nc.theta_b(:);
    Tcline=nc.Tcline(:);
  end
  if (isempty(Tcline))
%    disp('UCLA version 2');
    hc=nc.hc(:);
  else
    hmin=min(min(h));
    hc=min(hmin,Tcline);
  end
  Nr=length(nc('s_rho'));
  s_coord=1;
  VertCoordType = nc.VertCoordType(:);
  if isempty(VertCoordType),
    vtrans=nc{'Vtransform'}(:);
    if ~isempty(vtrans),
      s_coord=vtrans;
    end
  elseif VertCoordType=='NEW', 
   s_coord=2;
  end;
  if s_coord==2,
   hc=Tcline;
  end
%
%   Read the variable
%
  if vname(1)=='*'
    if strcmp(vname,'*Ke')
      u=mean(squeeze(nc{'u'}(:,:,J,I-1:I)),3);
      v=mean(squeeze(nc{'v'}(:,:,J-1:J,I)),3);
      var2=coef.*0.5.*(u.^2+v.^2);
    elseif strcmp(vname,'*Speed')
      u=mean(squeeze(nc{'u'}(:,:,J,I-1:I)),3);
      v=mean(squeeze(nc{'v'}(:,:,J-1:J,I)),3);
      var2=coef.*sqrt(u.^2+v.^2);
    elseif strcmp(vname,'*Rho_pot')
      temp=squeeze(nc{'temp'}(:,:,J,I));
      salt=squeeze(nc{'salt'}(:,:,J,I));
      var2=coef*rho_pot(temp,salt);
    elseif strcmp(vname,'*Chla')
      sphyto=squeeze(nc{'SPHYTO'}(:,:,J,I));
      lphyto=squeeze(nc{'LPHYTO'}(:,:,J,I));
      theta_m  =0.020;
      CN_Phyt  = 6.625;
      var2=coef*theta_m*(sphyto+lphyto)*CN_Phyt*12.;
      var2(var2<=0)=NaN;
    else
      disp('Sorry not implemented yet')
      return
    end 
  else
    var2=coef*squeeze(nc{vname}(:,:,J,I));
  end
%
  [T,N]=size(var2);
  if N==Nr+1
    type='w';
  else
    type='r';
  end
%
  var=0*(1:T);
%
  for l=1:T
    Z=squeeze(zlevs(h,squeeze(zeta(l,:,:)),theta_s,theta_b,hc,Nr,type,s_coord));
    var(l)=interp1(Z,var2(l,:),vlevel);
  end
end
%
% Get the time
%
time=nc{'scrum_time'}(:);
time=time/(24*3600)+ref_time;
close(nc)
%
ts=timeseries(var,datestr(time));

return

