function [u_bry,v_bry,ubar_bry,vbar_bry,...
          temp_bry,salt_bry]=vinterp_OGCM_bry(zr_bry,zu_bry,zv_bry,...
                                              dzr_bry,dzu_bry,dzv_bry,...
                                              u_bry,v_bry,ubar_bry,vbar_bry,...
                                              temp_bry,salt_bry,...
                                              N,Z,conserv);
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%
%
% Perform the vertical interpolations for the bry_file with OGCM for each boundary
%
%
%
%
zr_bry=squeeze(zr_bry);
zu_bry=squeeze(zu_bry);
zv_bry=squeeze(zv_bry);
dzr_bry=squeeze(dzr_bry);
dzu_bry=squeeze(dzu_bry);
dzv_bry=squeeze(dzv_bry);
%
% Add a level on top and bottom with no-gradient
%
u_bry=cat(1,u_bry(1,:),u_bry);
u_bry=cat(1,u_bry,u_bry(end,:));
v_bry=cat(1,v_bry(1,:),v_bry);
v_bry=cat(1,v_bry,v_bry(end,:));
temp_bry=cat(1,temp_bry(1,:),temp_bry);
temp_bry=cat(1,temp_bry,temp_bry(end,:));
salt_bry=cat(1,salt_bry,salt_bry(end,:));
salt_bry=cat(1,salt_bry(1,:),salt_bry);
% 
% Perform the vertical interpolations 
%
u_bry=squeeze(ztosigma_1d(flipdim(u_bry,1),zu_bry,flipud(Z)));
v_bry=squeeze(ztosigma_1d(flipdim(v_bry,1),zv_bry,flipud(Z)));
temp_bry=squeeze(ztosigma_1d(flipdim(temp_bry,1),zr_bry,flipud(Z)));
salt_bry=squeeze(ztosigma_1d(flipdim(salt_bry,1),zr_bry,flipud(Z)));
%
% Correct the horizontal transport 
% i.e. remove the interpolated tranport and add 
%      the OGCM transport
%
if conserv==1
  u_bry=u_bry-squeeze(tridim(squeeze(sum(u_bry.*dzu_bry)./sum(dzu_bry)),N));
  v_bry=v_bry-squeeze(tridim(squeeze(sum(v_bry.*dzv_bry)./sum(dzv_bry)),N));
  u_bry =u_bry + squeeze(tridim(ubar_bry,N));
  v_bry = v_bry + squeeze(tridim(vbar_bry,N));
end
%
% Barotropic velocities
%
ubar_bry=squeeze(sum(u_bry.*dzu_bry)./sum(dzu_bry));
vbar_bry=squeeze(sum(v_bry.*dzv_bry)./sum(dzv_bry));
%
return
