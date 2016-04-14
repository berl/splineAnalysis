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




function out = preshape2016(originaltraj, splineParam)

% function to establish the two radius parameters needed by
% shapefinder3.

% input                 traj            a trajectory. columns 1 and 2 are x
% and y, respectively
% output                out             struct with several fields for
% characterizing the input trajectory, and the trajectory rotated along the
% direction of maximum extent



%so I'm going to run through 180 degrees in 1 degree increments to determine the
%  largest maximum extent of the trajectory.  this direction i'll call the
%  semimajor direction


% then, i'm going to bin the semimajor direction into 10 (or so) segments
n =splineParam.nTrajShapePoints;



%%%%%%%%%%%%  Some important filtering:
%  0. see below

%   1. eliminate trajectories with clear pixel biasing:  more than 8% of
%   positions have a .5 or 0 decimal place:
%  2. throw out trajectories  that don't move more than 2
%   pixels in either direction.
%    3. ignore trajectories that don't have a 'linear shape' as defined
%    below


%%%%%%%%%%%%%%%%% JUNE 2010 MODIFICATION:

%%%%%  ADDING FILTER HERE to throw out trajectories with less than 100 points.
%%%%%   100 is a reasonable cutoff- if every 4th point is used, there are
%%%%%   25 independent values for MSD (or whatever) and 1/sqrt(25) = 20%, so all MSD statistical errors will be less than 20%
% this means that MSD points 1:4 should provide useful information, not
% statistically insignificant data

%0.
if numel(originaltraj(:,1)) < splineParam.minTrajLength
    out = 15;
    return
end



% 1. check for pixel biasing: ignores trajectories with more than 8% of
% positions with .5 or .0 decimal place
%|||||||||||||||||  CRITICAL NOTE:  ||||||||||||
%  this _MUST_ be done at a
%  point in the analysis when trajectory is IN PIXELS!!!
pixelbx = mod(originaltraj(:,1),1);
pixelby = mod(originaltraj(:,2),1);
num0x=(sum(pixelbx==0));
num0y=sum(pixelby==0);
num5x = sum(pixelbx==0.5);
num5y = sum(pixelby==0.5);
totalpixelb = num0x+num0y+num5x+num5y;
totalpointspix = (2*numel(pixelbx));
if (sum(pixelbx==0)+sum(pixelbx==0.5)+sum(pixelby==0)+sum(pixelby==0.5))/(2*numel(pixelbx)) >.08
    out = 1;
    %pause(5)
    return
end

%%%%%%%%%%%%%%%%%% 2.  ANOTHER FILTER:
%  now throwing out trajectories of particles that don't move more than .5
%  UNIT in both x and y directions.
% for june 2010 analysis, this isnt much of a filter: i'm actually looking
% for immobile particles, so this is really a strict case of zero movement-
% stdev of around +/-50nm covers 2/3 of a pixel, and for N >100, there should be
% substantial population of points beyond +/-1 std dev


%|||||||||||||||||  CRITICAL NOTE:  ||||||||||||
% this should be done at a
%  point in the analysis when trajectory is IN PIXELS!!!, otherwise it's a
%  restriction on overall size of the trajectory
if abs(max(originaltraj(:,1))-min(originaltraj(:,1))) < .5 && abs( max(originaltraj(:,2))-min(originaltraj(:,2))) <.5
    
    out = 2;
    return
end



% now move the center of mass of the trajectory to the origin. not
% strictly necessary, but it makes debugging plots a bit easier.
traj(:,1) = originaltraj(:,1)-mean(originaltraj(:,1));
traj(:,2) = originaltraj(:,2)-mean(originaltraj(:,2));




for k =1:180
    
    % rotangle function takes angles in degrees
    rotmatrixk = rotangle(k);
    rottraj=traj*rotmatrixk;
    widthx(k) = max(rottraj(:,1))-min(rottraj(:,1));
    widthy(k) = max(rottraj(:,2))-min(rottraj(:,2));
end
%hold off

maxangle = find(widthx == max(widthx));

maxrotmatrix = rotangle(maxangle);
maxtraj = traj*maxrotmatrix;


%%%% maxtraj is now the original trajectory, centered on the origin
%%%% and rotated so its semimajor axis is along the X axis.

%%% divide the semimajor axis into n bins
if ~splineParam.noPlots
    figure
end
binsize = (max(maxtraj(:,1))-min(maxtraj(:,1)))/n;
for j =1:n
    
    jthbinmin = min(maxtraj(:,1))+(j-1)*binsize;
    jthbinmax = jthbinmin+binsize;
    trajsectionj = maxtraj(maxtraj(:,1)>=jthbinmin & maxtraj(:,1) <= jthbinmax,:);
    if numel(trajsectionj)==0
        continue
    end
    
    
    out.ymean(j) = mean(trajsectionj(:,2));
    out.xmean(j) = mean(trajsectionj(:,1));
    out.ystd(j) = std(trajsectionj(:,2));
    max(trajsectionj(:,2))
    min(trajsectionj(:,2))
    out.width(j) = max(trajsectionj(:,2))-min(trajsectionj(:,2));
    
    if ~splineParam.noPlots
        hold all, plot(trajsectionj(:,1), trajsectionj(:,2)), axis equal
    end
    
end
    if ~splineParam.noPlots
hold off
    end
out.badstdy = std(maxtraj(:,2));
out.badwidthy = max(maxtraj(:,2))-min(maxtraj(:,2));



%%%% latent 3rd trajectory screen/filter here.  
%  
%   if the ratio of xextent to the  yextent of the trajectory  is
%   less than 1.75, the trajectory isn't really curvilinear, so
%   ignore it:


%   checkaspectratio = aspectratio(maxtraj);
%
%     if  checkaspectratio.aspectrm <1.75
%      %if mean(out.ystd)/(max(maxtraj(:,1))-min(maxtraj(:,1))) > 0.10
%          out = checkaspectratio.aspectrm ;
%          return
%      end
%



%%%%%%%%%VERY IMPORTANT STEP BELOW:  replacing x and y data with
%%%%%%%%% CENTERED (on origin) and ROTATED x and y data.
out.rottraj = [maxtraj, originaltraj(:,3:end)];

out.leftmost = maxtraj(maxtraj(:,1)==min(maxtraj(:,1)),1:2);
if numel(out.leftmost) > 2
    out.leftmost = out.leftmost(1,:);  %in trajectories with pixel binning problems, there can be more than
    %one point as 'leftmost'
end

%%%%%% rotating trajectories so they're vertically oriented for subsequent
%%%%%% plots 
out.rottraj(:,1:2) = fliplr(out.rottraj(:,1:2));

out.leftmost = fliplr(out.leftmost);