function [lon,lat]=read_d3d_grd(grd_file)
% This script reads a deld3d curvilinear grid and makes sure the south-west corner has 
% index (1,1), and also ordered (eta,xi) as needed by croco 

% read the grid using the openearthtools wlgrid.m function (copied into this
% repo)
grd = wlgrid('read', grd_file);
lon=grd.X;
lat=grd.Y;

[n1,n2]=size(lon);

disp(['delft3d grid has dimensions (',num2str(n1),',', num2str(n2),')']);

% Delft3d convention is to have the grid indices ordered (xi,eta)
% so we transpose them to be consistent with CROCO convention of (eta,xi)
lon=lon';
lat=lat';

disp('check if we need to flip dimensions');
% check if we need to flip the dimensions to ensure convension for specifying
% which sides represent N, S, E and W.
% south-west corner should have index (1,1)
%
% check the north south dimension. First index of the eta dimension
% should further south than the top index (even with a tilted grid)
if nanmean(lat(1,:))>nanmean(lat(end,:))
    lon=flip(lon,1);
    lat=flip(lat,1);
    disp('eta dimension was flipped');
else
    disp('eta dimension was not flipped');
end
%
% check the west east dimension. First index of the xi dimension
% should be west of the last index (even with a tilted grid)
if nanmean(lon(:,1))>nanmean(lon(:,end))
    lon=flip(lon,2);
    lat=flip(lat,2);
    disp('xi dimension was flipped');
else
    disp('xi dimension was not flipped');
end

end

