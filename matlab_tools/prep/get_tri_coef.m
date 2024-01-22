function [elem,coef,nnel] = get_tri_coef(lonp,latp,lonc,latc,maskp)

%
% Inputs:
%          parent lon and lat 2d arrays (lonp, latp)
%          child lon and lat 2d arrays (lonc, latc)
%
% Ouputs:
%          elem - pointers to 2d gridded data (at lonp,latp locations) from
%                 which the interpolation is computed (3 for each child point)
%          coef - linear interpolation coefficients
%          nnel - 2d pointer to the nearest neighbor in the horizontal
%
%  Use:
%          To subsequently interpolate data from Fp(lonp,latp) to Fc(lonc,latc), the following
%          will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
%          griddata calls. Since it avoids repeated triangulations and tsearches (that are done
%          with every call to griddata) it should be much faster.
%
%   (c) 2007, Jeroen Molemaker, UCLA 
%-------------------------------------------------------------------------


 [Mp,Lp] = size(lonp);
 [Mc,Lc] = size(lonc);

%
%%  Project lon, lat with a gnomonic projection for accurate distances.
%
   width1 = max(max(lonp))-min(min(lonp));
   width2 = max(max(latp))-min(min(latp));
   width3 = max(max(lonc))-min(min(lonc));
   width4 = max(max(latc))-min(min(latc));
   width = max([width1 width2 width3 width4]);
   if width<100
    lon0 = mean(mean(mean(lonc)));
    lat0 = mean(mean(mean(latc)));
    [xp,yp] = gnomonic(lonp,latp,lon0,lat0);
    [xc,yc] = gnomonic(lonc,latc,lon0,lat0);
   else
    disp('too big for gnomonic projection')
    xp = lonp; yp = latp;
    xc = lonc; yc = latc;
   end

  Xp    = [reshape(xp,Mp*Lp,1) reshape(yp,Mp*Lp,1) ];
  Xc    = [reshape(xc,Mc*Lc,1) reshape(yc,Mc*Lc,1) ];

  tri     = delaunay(xp,yp);
  [tn,pn] = tsearchn(Xp,tri,Xc);

% Fix to deal with child points that are outside parent grid (those points should be masked!)
  if (length(tn(~isfinite(tn)))>0);
    disp('Warning in get_tri_coef: outside point(s) detected.');
    [xc,yc] = fix_outside_child(xc,yc,tn);
    Xc      = [reshape(xc,Mc*Lc,1) reshape(yc,Mc*Lc,1) ];
    [tn,pn] = tsearchn(Xp,tri,Xc);
  end;

  elem = reshape(tri(tn,:),Mc,Lc,3);
  coef = reshape(       pn,Mc,Lc,3);

  xpm = xp;
  ypm = yp;
  xpm(find(maskp==0)) = xpm(find(maskp==0)) + 1e4;  %% move the masked points far away
  ypm(find(maskp==0)) = ypm(find(maskp==0)) + 1e4;  %% move the masked points far away

  if 0
   figure(1)
   plot(xp,yp,'.b')
   figure(2)
   plot(xpm,ypm,'.r')
   error('debugging in get_tri_coef')
  end

  %%  We need a re-triangulation here.
  trim  = delaunay(xpm,ypm);
  Xpm   = [reshape(xpm,Mp*Lp,1) reshape(ypm,Mp*Lp,1) ];
  nnel  = dsearchn(Xpm,trim,Xp);                    %% find the nearest non-masked neighbor
  nnel = reshape(     nnel,Mp,Lp  );

 if 0
  figure(1);triplot(tri,xp,yp)
  error('debugging in get_tri_coef')
 end

 if 0
  mask = zeros(1,Mp,Lp);
  mask(1,:,:) =maskp;
  figure(1);pcolor(lonp,latp,squeeze(mask));colorbar;shading flat
  fmask = fillmask(mask,1,maskp,nnel);
  figure(2);pcolor(lonp,latp,squeeze(fmask));colorbar;shading flat
  error('debugging in get_tri_coef')
 end
