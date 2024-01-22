function extract_GLORYS_Mydata(HYCOM_dir,HYCOM_prefix,path,Y,M,...
                      lonT,latT,lonU,latU,lonV,latV,depth,time,...
                      trange,krange,jrange,...
                      i1min,i1max,i2min,i2max,i3min,i3max,...
                      Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract a subset from a local GLORYS file
% Write it in a local file (keeping the classic
% SODA netcdf format) 
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
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated   6-Sep-2006 by Pierrick Penven
%  Updated : 23-Oct-2009 by Gildas Cambon Adapatation for the special tratment 
%            for the bry file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(['    Download GLORYS for ',num2str(Y),...
      ' - ',num2str(M)])

% Initialise output variables here
% we need to first get dimensions of the raw (i.e. not subset) input grid 
% for initialising of output arrays before subsetting
glorys_day=datenum(Y,M,1)-datenum(1990,1,1);
glorysfile=[path,'algoa_',num2str(glorys_day,'%06d'),'.nc'];
nc=netcdf(glorysfile);
NxT=length(nc{'lonT'}(:));
NyT=length(nc{'latT'}(:));
NxU=length(nc{'lonU'}(:));
NyU=length(nc{'latU'}(:));
NxV=length(nc{'lonV'}(:));
NyV=length(nc{'latV'}(:));
close(nc)
Nz=length(krange);
Nt=length(trange);
ssh=NaN(Nt,NyT,NxT);
u=NaN(Nt,Nz,NyU,NxU);
v=NaN(Nt,Nz,NyV,NxV);
temp=NaN(Nt,Nz,NyT,NxT);
salt=NaN(Nt,Nz,NyT,NxT);
  
for t=1:Nt
  
    %
    % Get the dataset attributes
    %
    %disp('Get the dataset attributes')
    glorys_day=datenum(Y,M,t)-datenum(1990,1,1);
    glorysfile=[path,'algoa_',num2str(glorys_day,'%06d'),'.nc'];
    nc=netcdf(glorysfile);
    
    %
    % Get SSH
    %
    %disp('    ...SSH')
    %
    ssh_t=nc{'ssh'}(:);
    ssh_t(ssh_t>999)=NaN;
    ssh(t,:,:)=ssh_t;

    % Get U
    %
    %disp('    ...U')
    %
    u_t=nc{'u'}(:);
    u_t(u_t>999)=NaN;
    u(t,:,:,:)=u_t;
    %
    % Get V
    %
    %disp('    ...V')
    %
    v_t=nc{'v'}(:);
    v_t(v_t>999)=NaN;
    v(t,:,:,:)=v_t;
    
    %
    % Get TEMP
    %
    %disp('    ...TEMP')
    temp_t=nc{'temp'}(:);
    temp_t(temp_t>999)=NaN;
    temp(t,:,:,:)=temp_t;
    %
    % Get SALT
    %
    %disp('    ...SALT')
    salt_t=nc{'salt'}(:);
    salt_t(salt_t>999)=NaN;
    salt(t,:,:,:)=salt_t;

end

%
% Subset the data
%
if ~isempty(i2min)
  irange=i2min:i2max;
  ssh_sub=ssh(trange,jrange,irange);
  u_sub=u(trange,krange,jrange,irange);
  v_sub=v(trange,krange,jrange,irange);
  temp_sub=temp(trange,krange,jrange,irange);
  salt_sub=salt(trange,krange,jrange,irange);
end
%
if ~isempty(i1min)
  irange=i1min:i1max;
  
  ssh_sub0=ssh(trange,jrange,irange);
  ssh_sub=cat(ndims(ssh_sub),ssh_sub0,ssh_sub);
  
  u_sub0=u(trange,krange,jrange,irange);
  u_sub=cat(ndims(u_sub),u_sub0,u_sub);
  
  v_sub0=v(trange,krange,jrange,irange);
  v_sub=cat(ndims(v_sub),v_sub0,v_sub);
  
  temp_sub0=temp(trange,krange,jrange,irange);
  temp_sub=cat(ndims(temp_sub),temp_sub0,temp_sub);
  
  salt_sub0=salt(trange,krange,jrange,irange);
  salt_sub=cat(ndims(salt_sub),salt_sub0,salt_sub);
  
end
%
if ~isempty(i3min)
  irange=i3min:i3max;
  
  ssh_sub0=ssh(trange,jrange,irange);
  ssh_sub=cat(ndims(ssh_sub),ssh_sub,ssh_sub0);
  
  u_sub0=u(trange,krange,jrange,irange);
  u_sub=cat(ndims(u_sub),u_sub,u_sub0);
  
  v_sub0=v(trange,krange,jrange,irange);
  v_sub=cat(ndims(v_sub),v_sub,v_sub0);
  
  temp_sub0=temp(trange,krange,jrange,irange);
  temp_sub=cat(ndims(temp_sub),temp_sub,temp_sub0);
  
  salt_sub0=salt(trange,krange,jrange,irange);
  salt_sub=cat(ndims(salt_sub),salt_sub,salt_sub0);
  
end

%
% Create the GLORYS file
%
create_GLORYS_Mydata([HYCOM_dir,HYCOM_prefix,'Y',num2str(Y),'M',num2str(M),'.cdf'],...
            lonT,latT,lonU,latU,lonV,latV,depth,time,...
            temp_sub,salt_sub,u_sub,v_sub,ssh_sub,[],[],Yorig)
%
return
