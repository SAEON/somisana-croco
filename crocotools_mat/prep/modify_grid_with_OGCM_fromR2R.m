%--------------------------------------------------------------
%
%  Modify child grid topography such that it matches the interpolated
%  parent topography at the boundaries.
%
%  This script is for use with make_croco2croco.
%
%   (c) 2007 Jeroen Molemaker, UCLA
% 	modified by Serena Illig and Florian MONETTI, IRD (2010)
%--------------------------------------------------------------
%
  clear all
  crocotools_param

%
% CROCO parent and child grid files
  cgrid = grdname;
  cgrid2=[grdname(1:end-3),'_',OGCM_prefix,'modified.nc'];
              

% Copy cgrid into cgrid2 file 
  copyfile(cgrid,cgrid2,'f');

% Only match to parent topography on open boundaries
  obcflag              = [1 0 1 1];      % open boundaries flag (1=open , [S E N W])
%
% End user-defined----------------------------------------------
% Get minimal parent subgrid bounds
   %limits=new_bry_subgrid(pgrid,cgrid);
   %imin=limits(1); imax=limits(2);
   %jmin=limits(3); jmax=limits(4);

% Get topography data from childgrid
  nc = netcdf(cgrid,'nowrite');
  hc   = nc{'hraw'}(:);
  lonc = nc{'lon_rho'}(:);  
  latc = nc{'lat_rho'}(:);
  maskr = nc{'mask_rho'}(:);
  h__   = nc{'h'}(:);
  close(nc);
  [Mc,Lc]=size(hc);
  min_hc = min(min(hc));
  
  contourf(hc)

% Estimate parent grid topo with TEMP data 
  load('/hdd/DATA/SODA_3.4.1/Grid/e3t.50serena');
  
  % 1. Interpolate SODA to CROCO grid @ open OB
  [OGCM_dir,OGCM_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),'.nc']
  nc=netcdf([OGCM_dir,OGCM_prefix,'Y',num2str(Ymin),'M',num2str(Mmin),'.nc']);
 % lonT=nc{'lonT'}(:);
 % latT=nc{'latT'}(:);
 % Z=nc{'depth'}(:);
 % NZ=length(Z);
 % tin=1;
  lonT=nc{'lonT'}(:);
  latT=nc{'latT'}(:);
  Z=nc{'depth'}(:);
  NZ=length(Z);
  tin=1;  
  
  for k=1:NZ
    temp(k,:,:)=ext_data_OGCM_nofill(nc,lonT,latT,'temp',tin,lonc,latc,k, Roa,interp_method);
  end
  
  for i=1:length(lonc(1,:))
  for j=1:length(latc(:,1))    
        hp(j,i)=NaN;
        for k=1:NZ
            if (isfinite(temp(k,j,i)))
                hp(j,i)=Z(k)+e3t(k)/2.;
            end
        end 
  end
  end
  
  lonp=lonc;
  latp=latc;
  maskp=hp-hp+1;
  maskp(isnan(maskp))=0.;

 
  %nc    = netcdf(pgrid,'nowrite');
 % hp    = squeeze(nc{'h'}(jmin:jmax,imin:imax));
  %lonp  = squeeze(nc{'lon_rho'}(jmin:jmax,imin:imax));
  %latp  = squeeze(nc{'lat_rho'}(jmin:jmax,imin:imax));
 % maskp = squeeze(nc{'mask_rho'}(jmin:jmax,imin:imax));
  %close(nc);
 [Mp,Lp]=size(hp);
%  In order to most closely match the surface area of the boundary, we'll
%  set the depth of the parent grid to the minimum depth of the child where
%  ever it is masked
   hp(~maskp) = min_hc;

% Get interpolation coefficient to go to (lonc,latc).
  [elem,coef] = get_tri_coef(lonp,latp,lonc,latc,maskp);

%% parent grid topo at child locations
  hpi = sum(coef.*hp(elem),3);


  dist = zeros(Mc,Lc,4);
  for i = 1:Mc   %% north south
   for j = 1:Lc    %% east west
     dist(i,j,1) =      i/Mc + (1-obcflag(1))*1e6; % South
     dist(i,j,2) = (Lc-j)/Lc + (1-obcflag(2))*1e6; % East
     dist(i,j,3) = (Mc-i)/Mc + (1-obcflag(3))*1e6; % North
     dist(i,j,4) =      j/Lc + (1-obcflag(4))*1e6; % West
   end
  end
  dist = min(dist,[],3);

  alpha = 0.5*tanh(50*(dist-0.06))+0.5; %% Feel free to play with this function.
% alpha = 1.0 - 0.5*(cos(pi*dist.^2)+1);



  hcn = alpha.*hc + (1-alpha).*hpi;
  
  
  %%%
  hcns=smoothgrid(hcn,maskr,hmin,hmax_coast,hmax,...
             rtarget,n_filter_deep_topo,n_filter_final);
  %%%
  
  
  nc = netcdf(cgrid2,'write');
  nc{'h'}(:) = hcns;
  nc{'hraw'}(:) = hcn;
  close(nc);

%% Visualize the modification
warning off
  sc0 = min(min(hcn));
  sc1 = max(max(hcn));
  figure(1)
  subplot(2,2,1)
  pcolor(lonc,latc,hpi);caxis([sc0 sc1]);colorbar;shading flat
  title('Interpolated Parent Topo')
  subplot(2,2,2)
  pcolor(lonc,latc,hcns./maskr);caxis([sc0 sc1]);colorbar;shading flat

  title('Boundary Smoothed Child Topo')
  subplot(2,2,3)
  pcolor(lonc,latc,hcns-hpi);colorbar;shading flat
  title('Difference between Parent and smoothed child Topo');
  subplot(2,2,4)
  pcolor(lonc,latc,alpha);colorbar;shading flat
  title('Parent/Child transition function');
warning on
%% R-factors

  r1 = 0.*hc;
  r2 = 0.*hc;
  r1(1:end-1,:) = 0.5*(hc(2:end,:) - hc(1:end-1,:))./( hc(2:end,:) + hc(1:end-1,:) );
  r2(:,1:end-1) = 0.5*(hc(:,2:end) - hc(:,1:end-1))./( hc(:,2:end) + hc(:,1:end-1) );
  rc= 2*max(abs(r1),abs(r2));

  r1n = 0.*hcns;
  r2n = 0.*hcns;
  r1n(1:end-1,:) = 0.5*(hcns(2:end,:) - hcns(1:end-1,:))./( hcns(2:end,:) + hcns(1:end-1,:) );
  r2n(:,1:end-1) = 0.5*(hcns(:,2:end) - hcns(:,1:end-1))./( hcns(:,2:end) + hcns(:,1:end-1) );
  rn = 2*max(abs(r1n),abs(r2n));

  r1p = 0.*hp;
  r2p = 0.*hp;
  r1p(1:end-1,:) = 0.5*(hp(2:end,:) - hp(1:end-1,:))./( hp(2:end,:) + hp(1:end-1,:) );
  r2p(:,1:end-1) = 0.5*(hp(:,2:end) - hp(:,1:end-1))./( hp(:,2:end) + hp(:,1:end-1) );
  rp = 2*max(abs(r1p),abs(r2p));

  r1pi = 0.*hpi;
  r2pi = 0.*hpi;
  r1pi(1:end-1,:) = 0.5*(hpi(2:end,:) - hpi(1:end-1,:))./( hpi(2:end,:) + hpi(1:end-1,:) );
  r2pi(:,1:end-1) = 0.5*(hpi(:,2:end) - hpi(:,1:end-1))./( hpi(:,2:end) + hpi(:,1:end-1) );
  rpi = 2*max(abs(r1pi),abs(r2pi));

  figure(2)
  sc0 = 0;
  sc1 = 0.30;
  subplot(2,2,1)
  pcolor(lonp,latp,rp);caxis([sc0 sc1]);colorbar;shading flat
  title('r-factor parent Topo')

  subplot(2,2,2)
  pcolor(lonc,latc,rpi);caxis([sc0 sc1]);colorbar;shading flat
  title('r-factor interpolated parent Topo')

  subplot(2,2,3)
  pcolor(lonc,latc,rc);caxis([sc0 sc1]);colorbar;shading flat
  title('r-factor original Child Topo')

  subplot(2,2,4)
  pcolor(lonc,latc,rn);caxis([sc0 sc1]);colorbar;shading flat
  title('r-factor blended  Child Topo')
