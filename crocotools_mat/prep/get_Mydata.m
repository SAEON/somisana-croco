function var=get_Mydata(path,fname,vname,trange,krange,jrange,...
                    i1min,i1max,i2min,i2max,i3min,i3max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  var=get_Mydata(path,fname,vname,trange,krange,jrange,...
%             i1min,i1max,i2min,i2max,i3min,i3max)
%  
%  script analagous to getdap.m, but when local file is used
%  Download a data subsets from a OPENDAP server.
%
%  Take care of the greenwitch meridian
%  (i.e. get 3 subgrids defined by i1min,i1max,i2min,i2max,i3min,i3max
%  and concatenate them).
%
%  In case of network failure, the program resend the OPENDAP 
%  query until it works (no more than 100 times though...)
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
%  Copyright (c) 2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
url=[path,fname];
%
var=[];
%
if ~isempty(i2min)
  irange=['[',num2str(i2min),':',num2str(i2max),']'];
  var=readdap(url,vname,[trange,krange,jrange,irange]);
end
%
if ~isempty(i1min)
  irange=['[',num2str(i1min),':',num2str(i1max),']'];
  var0=readdap(url,vname,[trange,krange,jrange,irange]);
  var=cat(ndims(var),var0,var);
end
%
if ~isempty(i3min)
  irange=['[',num2str(i3min),':',num2str(i3max),']'];
  var0=readdap(url,vname,[trange,krange,jrange,irange]);
  var=cat(ndims(var),var,var0);
end
%
return
