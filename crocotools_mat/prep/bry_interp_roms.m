function bry_interp_roms(bryname,lon,lat,hisname,...
                    dataname,vname,obcndx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  largely based on bry_interp() function in roms_tools
%
%  edited for use in 1-way offline nesting from existing roms run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[M,L]=size(lon);
%
% set the default value if no data
%
% default=NaN;
%
% Read in the datafile
%
nc=netcdf(hisname,'r');
lon_rho_his=nc{'lon_rho'}(:);
lat_rho_his=nc{'lat_rho'}(:);
lon_u_his=nc{'lon_u'}(:);
lat_u_his=nc{'lat_u'}(:);
lon_v_his=nc{'lon_v'}(:);
lat_v_his=nc{'lat_v'}(:);
if (dataname == 'u' | dataname == 'ubar')
    lon_his=lon_u_his;
    lat_his=lat_u_his;
elseif (dataname == 'v' | dataname == 'vbar')
    lon_his=lon_v_his;
    lat_his=lat_v_his;
else
    lon_his=lon_rho_his;
    lat_his=lat_rho_his;
end

Nz=length(nc('s_rho'));
data_his=nc{dataname}(:);
time=nc{'time'}(:);
tlen=length(time);

%
% get the boundary position
%
if obcndx==1
%
% Southern boundary
% 
  iroms=(1:L);
  jroms=1;
elseif obcndx==2
%
% Eastern boundary
% 
  iroms=L;
  jroms=(1:M);
elseif obcndx==3
%
% Northern boundary
% 
  iroms=(1:L);
  jroms=M;
elseif obcndx==4
%
% Western boundary
% 
  iroms=1;
  jroms=(1:M);
end
%
lon=lon(jroms,iroms);
lat=lat(jroms,iroms);
%
% get a data subgrid (dependant of the OBC used)
%
%dl=1.6;  %-> good for WOA2009
%dl=2.0;   %-> good for CARS2009  

% lonmin=min(min(lon))-dl;
% lonmax=max(max(lon))+dl;
% latmin=min(min(lat))-dl;
% latmax=max(max(lat))+dl;
% %
% j=find(Y>=latmin & Y<=latmax);
% i1=find(X-360>=lonmin & X-360<=lonmax);
% i2=find(X>=lonmin & X<=lonmax);
% i3=find(X+360>=lonmin & X+360<=lonmax);
% x=cat(1,X(i1)-360,X(i2),X(i3)+360);
% y=Y(j);
%
% Open the Z-boundary file
%
nc_bry=netcdf(bryname,'write');
%
% Check the time
%
%tbry=nc_bry{'bry_time'}(:); 

% loop on time
%
%missval=nc{dataname}.missing_value(:);
for l=1:tlen
%for l=1:1
  %disp(['time index: ',num2str(l),' of total: ',num2str(tlen)])
  dims=size(lon);
  if (ndims(data_his)==4) 
      for k=1:Nz
          data=squeeze(data_his(l,k,:,:)); % 2D data for this level and timestep
          if dims(1)==1
            datasgrid(k,:)=interp2(lon_his,lat_his,data,lon,lat,'spline');
          else
            datasgrid(k,:)=(interp2(lon_his,lat_his,data,lon,lat,'spline'))';  
          end 
      end
  end
  if (ndims(data_his)==3)
      data=squeeze(data_his(l,:,:)); % 2D data for this timestep
      if dims(1)==1
        datasgrid(:)=interp2(lon_his,lat_his,data,lon,lat,'spline');
      else
        datasgrid(:)=(interp2(lon_his,lat_his,data,lon,lat,'spline'))';  
      end 
  end
  nc_bry{vname}(l,:,:,:)=datasgrid;
end
close(nc_bry);
close(nc);
return

