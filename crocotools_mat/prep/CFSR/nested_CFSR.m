% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% this script loops through all months as defined in crocotools_param and 
% generates a nested bulk file for each month, bypassing the nestgui interface

clear all
close all
crocotools_param

nc_suffix=['.nc.',num2str(level)];
grdname=[grdname,'.',num2str(level)];
%[folder,grd,ext] = fileparts(grdname);

% will use values in crocotools_param if below is commented
% Ymin=2009;
% Ymax=2009;
% Mmin=9;
% Mmax=9;

for i=Ymin:Ymax
    
    if i==Ymin
        M_start=Mmin;
    else
        M_start=1;
    end
    if i==Ymax
        M_end=Mmax;
    else
        M_end=12;
    end
    
    for j=M_start:M_end
        
        blkname_p=[folder_p,'/blk_CFSR_Y',num2str(i),'M',num2str(j),'.nc'];
        [folder_p,blk_p,ext_p] = fileparts(blkname_p);
        
        blkname=[blk_p,nc_suffix];
        
        nested_bulk(grdname,blkname_p,blkname)
                
    end
end

