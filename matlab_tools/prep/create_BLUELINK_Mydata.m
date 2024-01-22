function create_BLUELINK_Mydata(fname,lonT,latT,lonU,latU,lonV,latV,depth,time,...
                     temp,salt,u,v,ssh,taux,tauy,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Create the OGCM file
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
%  Updated    6-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
missval=NaN;
disp('    Create the OGCM file')
nc=netcdf(fname,'clobber');
redef(nc);
nc('lonT')=length(lonT);
nc('latT')=length(latT);
nc('lonU')=length(lonU);
nc('latU')=length(latU);
nc('lonV')=length(lonV);
nc('latV')=length(latV);
nc('depth')=length(depth);
nc('time')=length(time);
nc{'temp'}=ncfloat('time','depth','latT','lonT') ;
nc{'temp'}.long_name=ncchar('TEMPERATURE');
nc{'temp'}.long_name='TEMPERATURE';
nc{'temp'}.units=ncchar('deg. C');
nc{'temp'}.units='deg. C';
nc{'temp'}.missing_value=missval;
nc{'salt'}=ncfloat('time','depth','latT','lonT') ;
nc{'salt'}.long_name=ncchar('SALINITY');
nc{'salt'}.long_name='SALINITY';
nc{'salt'}.units=ncchar('ppt');
nc{'salt'}.units='ppt';
nc{'salt'}.missing_value=missval;
nc{'u'}=ncfloat('time','depth','latU','lonU') ;
nc{'u'}.long_name=ncchar('ZONAL VELOCITY');
nc{'u'}.long_name='ZONAL VELOCITY';
nc{'u'}.units=ncchar('m/sec');
nc{'u'}.units='m/sec';
nc{'u'}.missing_value=missval;
nc{'v'}=ncfloat('time','depth','latV','lonV') ;
nc{'v'}.long_name=ncchar('MERIDIONAL VELOCITY');
nc{'v'}.long_name='MERIDIONAL VELOCITY';
nc{'v'}.units=ncchar('m/sec');
nc{'v'}.units='m/sec';
nc{'v'}.missing_value=missval;
nc{'ubar'}=ncfloat('time','latU','lonU') ;
nc{'ubar'}.long_name=ncchar('ZONAL BAROTROPIC VELOCITY');
nc{'ubar'}.long_name='ZONAL BAROTROPIC VELOCITY';
nc{'ubar'}.units=ncchar('m/sec');
nc{'ubar'}.units='m/sec';
nc{'ubar'}.missing_value=missval;
nc{'vbar'}=ncfloat('time','latV','lonV') ;
nc{'vbar'}.long_name=ncchar('MERIDIONAL BAROTROPIC VELOCITY');
nc{'vbar'}.long_name='MERIDIONAL BAROTROPIC VELOCITY';
nc{'vbar'}.units=ncchar('m/sec');
nc{'vbar'}.units='m/sec';
nc{'vbar'}.missing_value=missval;
nc{'taux'}=ncfloat('time','latU','lonU') ;
nc{'taux'}.long_name=ncchar('TAU_X');
nc{'taux'}.long_name='TAU_X';
nc{'taux'}.units=ncchar('N.m-2');
nc{'taux'}.units='N.m-2';
nc{'taux'}.missing_value=missval;
nc{'tauy'}=ncfloat('time','latV','lonV') ;
nc{'tauy'}.long_name=ncchar('TAU_Y');
nc{'tauy'}.long_name='TAU_Y';
nc{'tauy'}.units=ncchar('N.m-2');
nc{'tauy'}.units='N.m-2';
nc{'tauy'}.missing_value=missval;
nc{'ssh'}=ncfloat('time','latT','lonT') ;
nc{'ssh'}.long_name=ncchar('SEA LEVEL HEIGHT');
nc{'ssh'}.long_name='SEA LEVEL HEIGHT';
nc{'ssh'}.units=ncchar('m');
nc{'ssh'}.units='m';
nc{'ssh'}.missing_value=missval;
nc{'lonT'}=ncdouble('lonT') ;
nc{'lonT'}.units=ncchar('degrees_east');
nc{'lonT'}.units='degrees_east';
nc{'latT'}=ncdouble('latT') ;
nc{'latT'}.units=ncchar('degrees_north');
nc{'latT'}.units='degrees_north';
nc{'lonU'}=ncdouble('lonU') ;
nc{'lonU'}.units=ncchar('degrees_east');
nc{'lonU'}.units='degrees_east';
nc{'latU'}=ncdouble('latU') ;
nc{'latU'}.units=ncchar('degrees_north');
nc{'latU'}.units='degrees_north';
nc{'lonV'}=ncdouble('lonV') ;
nc{'lonV'}.units=ncchar('degrees_east');
nc{'lonV'}.units='degrees_east';
nc{'latV'}=ncdouble('latV') ;
nc{'latV'}.units=ncchar('degrees_north');
nc{'latV'}.units='degrees_north';
nc{'depth'}=ncdouble('depth') ;
nc{'depth'}.units=ncchar('meters');
nc{'depth'}.units='meters';
nc{'time'}=ncdouble('time') ;
eval(['nc{''time''}.units = ncchar(''days since 1-Jan-',num2str(Yorig),' 00:00:0.0'');'])
eval(['nc{''time''}.units = ''days since 1-Jan-',num2str(Yorig),' 00:00:0.0'';'])
endef(nc);
%
% File the file
%
disp('    Fill the OGCM file')
nc{'depth'}(:)=depth;
nc{'latT'}(:)=latT;
nc{'lonT'}(:)=lonT;
nc{'latU'}(:)=latU;
nc{'lonU'}(:)=lonU;
nc{'latV'}(:)=latV;
nc{'lonV'}(:)=lonV;
%
for tndx=1:length(time)
%
nc{'time'}(tndx)=time(tndx);
%
if length(time)==1
  nc{'ssh'}(tndx,:,:)=ssh;
  %nc{'taux'}(tndx,:,:)=taux;
  %nc{'tauy'}(tndx,:,:)=tauy;
  u1=u;
  v1=v;
  nc{'u'}(tndx,:,:,:)=u1;
  nc{'v'}(tndx,:,:,:)=v1;
  nc{'temp'}(tndx,:,:,:)=temp;
  nc{'salt'}(tndx,:,:,:)=salt;
else
  nc{'ssh'}(tndx,:,:)=squeeze(ssh(tndx,:,:));
  %nc{'taux'}(tndx,:,:)=squeeze(taux(tndx,:,:));
  %nc{'tauy'}(tndx,:,:)=squeeze(tauy(tndx,:,:));
  u1=squeeze(u(tndx,:,:,:));
  v1=squeeze(v(tndx,:,:,:));
  nc{'u'}(tndx,:,:,:)=u1;
  nc{'v'}(tndx,:,:,:)=v1;
  nc{'temp'}(tndx,:,:,:)=squeeze(temp(tndx,:,:,:));
  nc{'salt'}(tndx,:,:,:)=squeeze(salt(tndx,:,:,:));
end
%
% Compute the barotropic velocities
%
masku=isfinite(u1);
maskv=isfinite(v1);
u1(isnan(u1))=0;
v1(isnan(v1))=0;
dz=gradient(depth);
NZ=length(depth);
du=0*squeeze(u1(1,:,:));
zu=du;
dv=0*squeeze(v1(1,:,:));
zv=dv;
for k=1:NZ
  du=du+dz(k)*squeeze(u1(k,:,:));
  zu=zu+dz(k)*squeeze(masku(k,:,:));
  dv=dv+dz(k)*squeeze(v1(k,:,:));
  zv=zv+dz(k)*squeeze(maskv(k,:,:));
end
du(zu==0)=NaN;
dv(zv==0)=NaN;
zu(zu==0)=NaN;
zv(zv==0)=NaN;
ubar=du./zu;
vbar=dv./zv;
%
nc{'ubar'}(tndx,:,:)=ubar;
nc{'vbar'}(tndx,:,:)=vbar;
%
end
%
close(nc)
%
return
