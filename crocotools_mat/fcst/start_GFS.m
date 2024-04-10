%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Add the paths of the different toolboxes
% 
%  this file is here purely for the case of when we need to run reformat_GFS.m 
%  via our matlab docker image. So paths are relative to the docker image 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Add the paths of the different toolboxes'])
% paths are relative to the docker image used as part of somisana ops
% adjust for local development
addpath(genpath('/somisana-croco/crocotools_mat/misc')); % path inside our docker image
addpath('/somisana-croco/crocotools_mat/fcst/'); % path inside our docker image

