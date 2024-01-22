function extract_BLUELINK_Mydata(BLUELINK_dir,BLUELINK_prefix,path,year,month,...
                      lonT,latT,lonUV,latUV,depth,time,...
                      trange,krange,jrange,...
                      i1min,i1max,i2min,i2max,i3min,i3max,...
                      Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract a subset from local BLUELINK files
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
disp(['    Download BLUELINK for ',num2str(year),...
      ' - ',num2str(month)])
%
% Get the dataset attributes
%
disp('Get the dataset attributes')
%nc=netcdf(bluelinkfile);

%
% Get SSH
%
disp('    ...SSH')
var='eta_t';
bluelinkfile=[path,var,'_',num2str(year),'_',num2str(month),'.nc'];
nc=netcdf(bluelinkfile);
ssh=nc{var}(:);
close(nc);
scale = ncreadatt(bluelinkfile,var,'scale_factor');
offset = ncreadatt(bluelinkfile,var,'add_offset');
missingval = ncreadatt(bluelinkfile,var,'missing_value');
if missingval<0
  ssh(ssh<=(0.999*missingval))=NaN;
else
  ssh(ssh>=(0.999*missingval))=NaN;
end
ssh=ssh.*scale+offset;
ssh=double(ssh);

%
% Get TAUX
%
% not in the roms file
%
% Get TAUY
%
% not in the roms file
%
% Get U
%
disp('    ...U')
var='u';
bluelinkfile=[path,var,'_',num2str(year),'_',num2str(month),'.nc'];
nc=netcdf(bluelinkfile);
u=nc{var}(:);
close(nc);
scale = ncreadatt(bluelinkfile,var,'scale_factor');
offset = ncreadatt(bluelinkfile,var,'add_offset');
missingval = ncreadatt(bluelinkfile,var,'missing_value');
if missingval<0
  u(u<=(0.999*missingval))=NaN;
else
  u(u>=(0.999*missingval))=NaN;
end
u=u.*scale+offset;
u=double(u);

%
% Get V
%
disp('    ...V')
var='v';
bluelinkfile=[path,var,'_',num2str(year),'_',num2str(month),'.nc'];
nc=netcdf(bluelinkfile);
v=nc{var}(:);
close(nc);
scale = ncreadatt(bluelinkfile,var,'scale_factor');
offset = ncreadatt(bluelinkfile,var,'add_offset');
missingval = ncreadatt(bluelinkfile,var,'missing_value');
if missingval<0
  v(v<=(0.999*missingval))=NaN;
else
  v(v>=(0.999*missingval))=NaN;
end
v=v.*scale+offset;
v=double(v);

%
% Get TEMP
%
disp('    ...TEMP')
var='temp';
bluelinkfile=[path,var,'_',num2str(year),'_',num2str(month),'.nc'];
nc=netcdf(bluelinkfile);
temp=nc{var}(:);
close(nc);
scale = ncreadatt(bluelinkfile,var,'scale_factor');
offset = ncreadatt(bluelinkfile,var,'add_offset');
missingval = ncreadatt(bluelinkfile,var,'missing_value');
if missingval<0
  temp(temp<=(0.999*missingval))=NaN;
else
  temp(temp>=(0.999*missingval))=NaN;
end
temp=temp.*scale+offset;
temp=double(temp);

%
% Get SALT
%
disp('    ...SALT')
var='salt';
bluelinkfile=[path,var,'_',num2str(year),'_',num2str(month),'.nc'];
nc=netcdf(bluelinkfile);
salt=nc{var}(:);
close(nc);
scale = ncreadatt(bluelinkfile,var,'scale_factor');
offset = ncreadatt(bluelinkfile,var,'add_offset');
missingval = ncreadatt(bluelinkfile,var,'missing_value');
if missingval<0
  salt(salt<=(0.999*missingval))=NaN;
else
  salt(salt>=(0.999*missingval))=NaN;
end
salt=salt.*scale+offset;
salt=double(salt);

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
% Create the BLUELINK file
%
create_BLUELINK_Mydata([BLUELINK_dir,BLUELINK_prefix,'Y',num2str(year),'M',num2str(month),'.cdf'],...
            lonT,latT,lonUV,latUV,lonUV,latUV,depth,time,...
            temp_sub,salt_sub,u_sub,v_sub,ssh_sub,[],[],Yorig)
%
return
