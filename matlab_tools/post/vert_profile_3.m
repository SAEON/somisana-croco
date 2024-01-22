function [var,Z]=vert_profile_3(diafile,hisfile,gridfile,lon0,lat0,vname,tindex,coef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get a vertical Profile
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  this is as per vert_profile() in croco_tools, but edited so that
%  the profile and z levels are output by the function rather than plotted
%  rev *_3 edited to allow for extracting from a diagnostics file- edit
%  needed as zeta variable is in the hisfile, but not the diafile, so we
%  need to read in both
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
  disp('Warning: profile place not found')
  [M,L]=size(lon);
  I=round(L/2);
  J=round(M/2);
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
nd=netcdf(diafile);
var=coef*squeeze(nd{vname}(tindex,:,J,I));

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

