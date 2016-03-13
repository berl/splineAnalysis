%     Copyright 2016 Brian R. Long b.r.long.08@gmail.com
%
%
%
%     If you use this software, please consider citing our related paper:
%  B.R. Long and T.Q. Vu 'Spatial structure and diffusive dynamics from single-particle
% trajectories using spline analysis'
% Biophys J. 2010 Apr 21;98(8):1712-21.
%
%     This file is part of splineAnalysis.
%
%     splineAnalysis is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     splineAnalysis is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with splineAnalysis.  If not, see <http://www.gnu.org/licenses/>.


function out = combinatron2016(trajectory, splineParam)

%  input:   trajectory          single trajectory in crocker&grier format
%                               (see splineanalysis2016.m for details)
%           splineParam         struct with various parameters 
%
%
%  out:     struct with 3 fields:
%                       a, b  from shapefinder3
%                       rottraj  the rotated trajectory from preshape10
preout = preshape2016(trajectory, splineParam);


if ~isstruct(preout)==1
    out = preout;
    return
end



r1= splineParam.radiusFactor*mean(preout.ystd);
r2 = splineParam.radiusRatio*r1;
[out.a out.b ]= shapefinder3(preout.rottraj, r1,r2, preout.leftmost, splineParam);


%%%%%%  IMPORTANT: add a point to each end of the array that will
%%%%%%  be just an extension of the last  (or first) displacement

[checkarows, checkacol] = size(out.a);
out.a = [ 0, 0,1 ; out.a(out.a(:,3)~=0,1:3); 0, 0,1; 0, 0 ,1];

out.a(1,1) = out.a(2,1)-(out.a(3,1)-out.a(2,1));
out.a(1,2) = out.a(2,2)-(out.a(3,2)-out.a(2,2));

out.a(end-1,1) = out.a(end-2,1)+ (out.a(end-2,1)-out.a(end-3,1));
out.a(end-1,2) = out.a(end-2,2)+ (out.a(end-2,2)-out.a(end-3,2));
% ADD a second extra point to the 'end'  to really make sure the spline
% keeps going past the end of the trajectory

out.a(end,1)= out.a(end-1,1)+(out.a(end-2,1)-out.a(end-3,1));
out.a(end,2)= out.a(end-1,2)+(out.a(end-2,2)-out.a(end-3,2));
out.rottraj = preout.rottraj;




