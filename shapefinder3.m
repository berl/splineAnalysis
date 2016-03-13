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




function [comout, locout] = shapefinder3(traj, rshape, rshape2, firstpoint, splineParam)

% %   inputs:     traj            a trajectory
%                 rshape          the radius over which averages will be taken
%                 rdisp           the distance between locout points
%
%
%   output:      comout           points approximating the shape of the
%                                 trajectory, including center-of-mass
%                                 adjustments.
%                   locout        same idea without the center-of-mass
%                                 corrections


%%%%%%%%%%%%%%  the first location is now an input,firstpoint, inferred from the output of the
%%%%%%%%%%%%%%  preshape function.

% initializing a couple boolean values here:
endit = 0;
expando =0;




if nargin==4
    locout(1,1:2)=firstpoint;
else
    
    %  point 1  is a little unusual.  the first locout point is just the
    %  first point of the trajectory
    locout(1,1:2) =[traj(1,1) traj(1,2)];
end


loc1 = locout(1,1:2);



%  find the distance of each point from the first point in the trajectory:
distance1 = (traj(:,1)-loc1(1)).*(traj(:,1)-loc1(1))+(traj(:,2)-loc1(2)).*(traj(:,2)-loc1(2));
%now find the center of mass within rshape of locout(1,1:2):
comout(1,1:2) =  mean(traj( distance1 < rshape^2 ,1:2),1);
comout(1,3) =1;


%%%%%%% now run through a loop stepping along the trajectory. 
%%%%%%% This is arranged as a for loop, but functions a bit more like a while loop...subsequent
%%%%%%% versions could be an open-ended structure
if ~splineParam.noPlots
    fig1 = figure;
end


for i = 2:(4*splineParam.nTrajShapePoints)
    
    %%%%  first, calculate the vectors between the com location (i-1) and each traj point
    vectorsi(:,1) = traj(:,1)-comout(i-1,1);
    vectorsi(:,2) = traj(:,2)-comout(i-1,2);
    if expando
        rshape2 = 2*rshape2;
        endit = 1;
    end
    
    %%% now look at the particles in the annulus rshape < r < rshape2
    avectorsi =  vectorsi( (sum(vectorsi.^2,2)>rshape^2)&(sum(vectorsi.^2,2)<=rshape2^2),:);
    
    apoints =        traj( (sum(vectorsi.^2,2)>rshape^2)&(sum(vectorsi.^2,2)<=rshape2^2),:);
    if ~splineParam.noPlots
        figure(fig1)
        plot(traj(:,1), traj(:,2), '.'), axis equal
        hold all, plot(apoints(:,1), apoints(:,2), 'o')
        hold all, plot(comout(i-1,1), comout(i-1,2),'* black')
        
        
        rectangle('Position',[comout(i-1,1)-rshape, comout(i-1,2)-rshape, 2*rshape, 2*rshape], 'Curvature', 1)
        rectangle('Position',[comout(i-1,1)-rshape2, comout(i-1,2)-rshape2, 2*rshape2,2*rshape2], 'Curvature', 1)
    end
    
    % now I need to identify the direction of the previous point... for the
    % first time through the loop (i =2), there is no previous direction,
    % only a previous point.  so i'll have to use the direction with the
    % highest density of particles.
    
    
    %otherwise
    if i ~=2
        if ~expando
            prevector = [comout(i-1,1)-comout(i-2,1), comout(i-1,2)-comout(i-2,2)];
        else   %if we're in expando mode (rshape2->2*rshape2), use previous, as long as there is one...
            prevector = [comout(max(i-2,1),1)-comout(max(i-3,1),1), comout(max(i-2,1),2)-comout(max(i-3,1),2)] ;
        end
    else
        prevector = [1,0];  %   if there's no previous angle, just use the x axis.
    end
    
    if prevector(1)==0   % if there happened to be an expando mode point at the first point, this statement is needed
        prevector = [1, 0];
    end
    
    %%%%%% below i'm using dot and cross products to figure out the angle.
    %%%%%  there still will be a reckoning needed to return to cartesian
    %%%%%  coordinates when the actual displacement is calculated.
    
    %    now calculating the cosine of the angle between prevector and avectorsi
    %    using the dot product
    size(avectorsi(:,1)*prevector(1) + avectorsi(:,2)*prevector(2));
    size((sqrt(sum(avectorsi.^2,2)*sum(prevector.^2))));
    
    cosi = (avectorsi(:,1)*prevector(1) + avectorsi(:,2)*prevector(2))./(sqrt(sum(avectorsi.^2,2)*sum(prevector.^2)));
    %%%%% but the cosine doesn't tell you the whole story. I'll use the sign of the cross
    %%%%% product to get a sign to associate with the cosine.
    
    
    signcrossi = sign(avectorsi(:,1).*prevector(2)-avectorsi(:,2).*prevector(1));
    
    

    signcrossi(signcrossi(:) == 0)= 1; 
    %   assign a nonzero value
    %   to the sign of the crossproduct for cases when the same location is
    %   found in the annuli of sequential angle searches.
    
    
    
    
    
    angleinfoi = [cosi, signcrossi];
    nangles =50;
    subtend = splineParam.thetaRange;
    angleincrement = (1/(2*(nangles+1)))*2*pi;
    for bin = -nangles:nangles
        binanglemin = max(-pi,bin*angleincrement-subtend/2);
        binanglemax = min(bin*angleincrement + subtend/2, pi);
        nsubtend(bin+nangles+1) = sum(acos(cosi).*signcrossi > binanglemin & acos(cosi).*signcrossi <=binanglemax);
    end
    angles = (-nangles:nangles)*angleincrement;
    
    
    
    
    
    
    %%%% the first direction is just towards the maximum angle, but subsequent
    %%%% (i>2) angles should be restricted to less than condition below
    if i ==2
        anglei = mean( angles(nsubtend== max(nsubtend)));
    else
        anglei = mean( angles((abs(angles) < splineParam.thetaMax)&(nsubtend== max(nsubtend(abs(angles) < splineParam.thetaMax)))));
    end

    
    
    
    sectoranglemin = max(-pi,anglei-subtend/2);
    sectoranglemax = min(anglei + subtend/2, pi);
    checktheangles(i) = anglei;
    % the actual location set for the i-th point is the mean of the
    % particle locations within the angle 'subtend/2' of 'anglei'.
    locsubtend = apoints(acos(cosi).*signcrossi > sectoranglemin & acos(cosi).*signcrossi <=sectoranglemax,:);
    
    
    
    
    if numel(locsubtend(:,1))== 0
        if endit   %if we've already been through this once (rshape2 is already expanded), quit the loop and get out
            break
        end
        expando = 1;
        comout(i,1:2) = comout(i-1,1:2);
        comout(i,3) = 0;
        continue
    elseif numel(locsubtend(:,1))< 3
        %  if there are only 1 or 2 particles determining the location, mark
        %  the 3rd column of the comout array with a 0, in case these should be ignored
        %  in the spline fitting.
        comout(i,1) = mean(locsubtend(:,1));
        comout(i,2) = mean(locsubtend(:,2));
        comout(i,3)=0;
        if endit   %if we've already been through this once, reset rshape2 and hope that
            % things are back to normal.
            rshape2 = rshape2/2;
            endit = 0;
            expando =0;
        end
        
    else
        comout(i,3)=1;
        comout(i,1) = mean(locsubtend(:,1));
        comout(i,2) = mean(locsubtend(:,2));
        if endit
            rshape2 = rshape2/2;
            endit = 0
            expando =0;
        end
    end

    
    
    
end


