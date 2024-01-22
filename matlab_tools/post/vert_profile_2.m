function [var,Z]=vert_profile_2(hisfile,gridfile,lon0,lat0,vname,tindex,coef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get a vertical Profile
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  this is as per vert_profile() in croco_tools, but edited so that
%  the profile and z levels are output by the function rather than plotted
%
%  Also modified to write out u and v relative to true north, rather than
%  grid north
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
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
%  Updated 02-Nov-2006 by Pierrick Penven (Yorig)
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
  tindex=1;
  disp(['Default time index: ',num2str(tindex)])
end
if nargin < 7
  coef=1;
  disp(['Default coef: ',num2str(coef)])
end
if nargin < 8
  Yorig=NaN;
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
%disp(['lon0 = ',num2str(lon0),' - lat0 = ',num2str(lat0)])
[J,I]=find((lat(1:end-1,1:end-1)<=lat0 & lat(2:end,2:end)>lat0 &...
            lon(2:end,1:end-1)<=lon0 & lon(1:end-1,2:end)>lon0)==1);
if isempty(I) |  isempty(J)
  %disp('Warning: profile place not found')
  %[M,L]=size(lon);
  %I=round(L/2);
  %J=round(M/2);
  
  % rather just return an array of nans
  nc=netcdf(hisfile);
  Nr=length(nc('s_rho'));
  close(nc);
  var=nan(Nr,1);
  Z=nan(Nr,1);
  return
  
end
% GGF hack edit needed in the case where two points are found
J=J(1);
I=I(1);
%disp(['I = ',int2str(I),' J = ',int2str(J)])
%lon1=lon(J,I);
%lat1=lat(J,I);
%disp(['lon1 = ',num2str(lon1),' - lat1 = ',num2str(lat1)])
%
% get the vertical levels
%
ng=netcdf(gridfile);
nc=netcdf(hisfile);
if vname(1)=='u'
  zeta=mean(squeeze(nc{'zeta'}(tindex,J,I:I+1)),2);
  h=mean(squeeze(ng{'h'}(J,I:I+1)),2);
  angle=mean(squeeze(nc{'angle'}(J,I:I+1)),2);
elseif vname(1)=='v'
  zeta=mean(squeeze(nc{'zeta'}(tindex,J:J+1,I)),1);
  h=mean(squeeze(ng{'h'}(J:J+1,I)),1);
  angle=mean(squeeze(nc{'angle'}(J:J+1,I)),1);
else
  zeta=squeeze(nc{'zeta'}(tindex,J,I));
  h=squeeze(ng{'h'}(J,I));
  angle=squeeze(nc{'angle'}(J,I));
end
close(ng)
theta_s=nc.theta_s(:);
if (isempty(theta_s))
%  disp('Rutgers version')
  theta_s=nc{'theta_s'}(:);
  theta_b=nc{'theta_b'}(:);
  Tcline=nc{'Tcline'}(:);
else 
%  disp('UCLA version');
  theta_b=nc.theta_b(:);
  Tcline=nc.Tcline(:);
end
if (isempty(Tcline))
%  disp('UCLA version 2');
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
% Read the variable
%
if vname(1)=='*'
  if strcmp(vname,'*Ke')
    u=mean(squeeze(nc{'u'}(tindex,:,J,I-1:I)),2);
    v=mean(squeeze(nc{'v'}(tindex,:,J-1:J,I)),2);
    var=coef.*0.5.*(u.^2+v.^2);
  elseif strcmp(vname,'*Speed')
    u=mean(squeeze(nc{'u'}(tindex,:,J,I-1:I)),2);
    v=mean(squeeze(nc{'v'}(tindex,:,J-1:J,I)),2);
    var=coef.*sqrt(u.^2+v.^2);
  elseif strcmp(vname,'*Rho')
    temp=squeeze(nc{'temp'}(tindex,:,J,I));
    salt=squeeze(nc{'salt'}(tindex,:,J,I));
    z=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,Nr,'r',s_coord));
    var=coef*rho_eos(temp,salt,z);
  elseif strcmp(vname,'*Rho_pot')
    temp=squeeze(nc{'temp'}(tindex,:,J,I));
    salt=squeeze(nc{'salt'}(tindex,:,J,I));
    var=coef*rho_pot(temp,salt);
  elseif strcmp(vname,'*Chla')
    sphyto=squeeze(nc{'SPHYTO'}(tindex,:,J,I));
    lphyto=squeeze(nc{'LPHYTO'}(tindex,:,J,I));
    theta_m  =0.020;
    CN_Phyt  = 6.625;
    var=coef*theta_m*(sphyto+lphyto)*CN_Phyt*12.;
    var(var<=0)=NaN;
  else
    disp('Sorry not implemented yet')
    return
  end 
elseif vname(1)=='u' || vname(1)=='v'
  var_u=coef*squeeze(nc{'u'}(tindex,:,J,I));  
  var_v=coef*squeeze(nc{'v'}(tindex,:,J,I)); 
  %
  %  Rotation to get back angles relative to true north
  %
  cosa = cos(angle);
  sina = sin(angle);
  if vname(1)=='u'
  var = var_u.*cosa - var_v.*sina;
  else
  var = var_v.*cosa + var_u.*sina;
  end
else
  var=coef*squeeze(nc{vname}(tindex,:,J,I));
end
%
N=length(var);
if N==Nr+1
  type='w';
else
  type='r';
end
Z=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,Nr,type,s_coord));
close(nc)
%
% Get the date
%
% [day,month,year,imonth,thedate]=...
% get_date(hisfile,tindex,Yorig);
% %
% plot(var,Z,'k')
% hold on
% plot(var,Z,'r.')
% hold off
% ylabel('Depth [m]')
% xlabel([vname,' - ',thedate])
return

