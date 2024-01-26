function var = read_CFSR_hourly(ifile,vname,...
                      tndx,...
                      i1min,i1max,i2min,i2max,i3min,i3max,...
                      jmin,jmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract a subset from downloaded CFSR data
% 
% From extract_NCEP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get the variable 2D subset (take care of greenwitch)
%
nc=netcdf(ifile);
numdimvar0=ncsize(nc{vname});
numdimvar=length(numdimvar0);
if ~isempty(i1min)
    if ( numdimvar > 3 )
      var1=squeeze(nc{vname}(tndx,1,jmin:jmax,i1min:i1max));
    else
      var1=squeeze(nc{vname}(tndx,jmin:jmax,i1min:i1max));
    end
else
    var1=[];
end
%
if ~isempty(i2min)
    if ( numdimvar > 3 )
      var2=squeeze(nc{vname}(tndx,1,jmin:jmax,i2min:i2max));
    else
      var2=squeeze(nc{vname}(tndx,jmin:jmax,i2min:i2max));
    end
else
    var2=[];
end
%  
if ~isempty(i3min)
    if ( numdimvar > 3 )
      var3=squeeze(nc{vname}(tndx,1,jmin:jmax,i3min:i3max));  
    else
      var3=squeeze(nc{vname}(tndx,jmin:jmax,i3min:i3max));    
    end
else
    var3=[];
end
%
var=cat(3,var1,var2,var3);
%
% North-South inversion
%
if (length(size(var))==2)
    var=flipdim(var,1);
elseif (length(size(var))==3)
    var=flipdim(var,2);
end
%
missing_value = ncreadatt(ifile,vname,'_FillValue');
if isempty(missing_value)
    missing_value=3.399999952144364e+38;
end
close(nc)
%
% Correct the variable
%
var(var==missing_value)=NaN;

return

