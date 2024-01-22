%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  Build a CROCO boundary file
%
%  Interpolate surface elevation, temperature, salinity and velocity from a
%  HYCOM simulation to get boundary conditions for
%  CROCO (boundary netcdf file) .
%
%  Data source : HYCOM
%    https://hycom.org/
%    http://tds.hycom.org/thredds/catalog.html
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
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
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param
%
bry_prefix=[bry_prefix,'_HYCOM_'];
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%

%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname,'r');
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
hmax=max(max(nc{'h'}(:)));
close(nc);

%
% Loop through all months
%
for Y=Ymin:Ymax
    if Y==Ymin
        mo_min=Mmin;
    else
        mo_min=1;
    end
    if Y==Ymax
        mo_max=Mmax;
    else
        mo_max=12;
    end
    for M=mo_min:mo_max
        disp(' ')
        disp(['Processing  year ',num2str(Y),...
            ' - month ',num2str(M)])
        disp(' ')
        
        %
        % open hycom data file for this month
        hycom_data=[hycom_dir,num2str(year),'_',num2str(month),'.nc'];
        nc_hycom=netcdf(hycom_data,'r');
        hycom_time=nc_hycom{'time'}(:); % CHECK WHAT TIME IS ACTUALLY CALLED!!
        % NEED TO CONVERT hycom_time TO SECONDS SINCE CROCO TIME ORIGIN
        
        % create the boundary file for this month
        bryname=[blk_prefix,'Y',num2str(Y),...
                'M',num2str(M),nc_suffix];
        create_bryfile(bryname,grdname,CROCO_title,obc,...
                 theta_s,theta_b,hc,N,...
                 hycom_time,0,'clobber',vtransform);
        
        
    end
end








%
% Create the boundary file in Z-coordinates
%
if (makeZbry)
  disp(' ')
  disp(' Create the boundary Z-file...')
%
% get Z
%
    nc=netcdf(temp_ann_data,'r');
  Z=nc{'Z'}(:);
  kmax=max(find(Z<hmax))-1;
  Z=Z(1:kmax);
  close(nc)
  create_bry_Z(Zbryname,grdname,CROCO_title,obc,...
                Z,woa_time,woa_cycle,'clobber');
  disp(' ')
  disp(' Horizontal extrapolations')
%
% Loop on the lateral boundaries 
%
  for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
	suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
	suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
	suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
	suffix='_west';
      end
      disp('  Temperature...')
      bry_interp(Zbryname,lon,lat,temp_month_data,temp_ann_data,...
               'temperature',['temp',suffix],obcndx,Roa);
      disp('  Salinity...')
      bry_interp(Zbryname,lon,lat,salt_month_data,salt_ann_data,...
               'salinity',['salt',suffix],obcndx,Roa);        
    end
  end
end

%
% Vertical interpolations 
%
if (makebry)
  disp(' ')
  disp(' Vertical interpolations')

%
% Loop on the lateral boundaries 
%
  for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
	suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
	suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
	suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
	suffix='_west';
      end
      disp(' ')
      disp('  Temperature...')
      vinterp_bry(bryname,grdname,Zbryname,['temp',suffix],obcndx);
      disp(' ')
      disp('  Salinity...')
      vinterp_bry(bryname,grdname,Zbryname,['salt',suffix],obcndx);
      if (insitu2pot)
        disp(' ')
        disp('  Compute potential temperature from in-situ...')
        getpot_bry(bryname,grdname,obcndx)
      end

%
% Geostrophy
%
      disp(' ')
      disp('  Compute geostrophic currents')
      geost_currents_bry(bryname,grdname,Zbryname,frcname,zref,obcndx)
    end
  end

%
% Remove avg SSH
%
  rmavgssh(bryname,grdname,obc)

end



%----------------------------------------------------------------------------
% Make a few plots
%----------------------------------------------------------------------------
if makeplot==1
  disp(' ')
  disp(' Make a few plots...')
  test_bry(bryname,grdname,'temp',1,obc)
  figure
  test_bry(bryname,grdname,'salt',1,obc)
  figure
  test_bry(bryname,grdname,'u',1,obc)
  figure
  test_bry(bryname,grdname,'v',1,obc)
  figure
  test_bry(bryname,grdname,'temp',6,obc)
  figure
  test_bry(bryname,grdname,'salt',6,obc)
  figure
  test_bry(bryname,grdname,'u',6,obc)
  figure
  test_bry(bryname,grdname,'v',6,obc)
end
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
