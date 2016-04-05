%     Copyright 2016 Brian R. Long b.r.long.08@gmail.com
%
%
%
%     If you use this software, please consider citing our related paper:
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


function out = splineanalysis2016(alltraj, savepath, savename, interactive, splineParam)




%  % %   %  %  INPUT
%           alltraj            2d array in the Crocker and Grier format,
%                               i.e. columns of
%                        |x | y| ... |frame number| particle number |

%                       where '...' indicates columns of other data
%                       such as brightness, eccentricity, etc.
%                       the rows are ordered by particle number and then  by frame number
%                       so particle 1's trajectory is first, followed by
%                       particle 2's  trajectory and  each particle's
%                       trajectory starts with its earliest frame and
%                       continues until the last frame where that particle is
%                       recorded.


%           savepath           a string with the full path to the location where output files
%                               will be saved.
%           savename           a string with the name prefix for files to
%                               be saved





%           interactive         boolean value where 1 uses an interactive
%                               mode in which trajectories can be
%                               classified after spline fitting.   see
%                               details below

%           splineParam         a struct with the following fields to
%                               control the free parameters for spline
%                               fitting.  If any fields are not set, (or if no variable is passed to splineanalysis2010), default values are
%                               used.

%               splineParam.thetarange          angular range used in
%                                               search for next spline guide point
%               splineParam.thetamax            the maximum angle between
%                                               sequential spline guide points
%               splineParam.radiusFactor        multaplicative factor
%                                               relating r1 to average
%                                               lateral standard deviation
%               splineParam.radiusRatio         ratio between r1 (normal
%                               mode) and r2 (expanded mode) used when looking for spline guide
%                               points when trajectory data is sparse.
%
%               splineParam.nTrajShapePoints    number of segments to
%                               divide the trajectory into when measuring the average local
%                               'width' of the trajectory
%               splineParam.nSplineCurvePoints  number of interpolated points in the
%                                               spline curve itself (interpolating between guide points)
%               splineParam.minTrajLength       minimum number of points
%                               per trajectory.  a value of 100
%



%  % %   %  %   OUTPUT

%            out            a struct with the following fields

%           out.traj        a matrix with the analyzed trajectory
%
% |x|y|...|frame number|particle number|splinex|spliney|signedDistance|deltaPerpX|deltaPerpY|deltaParallelX|deltaParallelY|distanceParallel|distancePerp|classification|

% |x|y|...|frame number|particle number|    this is the original input matrix
% |splinex|spliney|                         components of the spline curve at each point in the trajectory
% |signedDistance|                          signed distance along the spline
% |deltaPerpX|deltaPerpY|                   trajectory displacement perpendicular to the spline
% |deltaParallelX|deltaParallelY|           trajectory displacement parallel to the spline
% |distanceParallel|                        distance parallel to the spline
% |distancePerp|                            distance perpendicular to the spline
% |classification|                          user 






nofigurebuildup = 1;
if ~splineParam.noPlots
    figure
end


[rowsall, colsall] = size(alltraj);
toadd = zeros(rowsall,10);
out.traj = [alltraj, toadd];
thenums = unique(alltraj(:,colsall));

realn=0;
%%%%%%%%%%%  realn is the total number of trajectories that are actually
%%%%%%%%%%%  spline-fitted.


ignoredlist=[];
for j = 1:numel(thenums)
    thenums(j)
    whichrowsj = find(alltraj(:,colsall) == thenums(j));
    jtraj = alltraj(whichrowsj,:);
    [inrows, incols] =size(jtraj)
    
    %%%%%%%%%%%%%%%% for each trajectory, run combinatron2016, which uses preshape10 and
    %%%%%%%%%%%%%%%% shapefinder3 to determine the points through which the
    %%%%%%%%%%%%%%%% spline should be fit.  it also rotates the
    %%%%%%%%%%%%%%%% trajectories, so i'll plot the rotated version
    %%%%%%%%%%%%%%%% instead.
    
    automaticout = combinatron2016(jtraj, splineParam);
    
    if ~isstruct(automaticout)   % combinatron filters out bad trajectories and returns a 0 (not a struct) as an indicator of a bad traj.
        
        %%%%%%    need to mark this trajectory as IGNORED...
        ignoredlist(j,1) = automaticout;
        ignoredlist(j,2) = thenums(j);
        continue
    else
        
    end
    ignoredlist(j,2) = thenums(j);
    
    automaticout.a;
    jtraj = automaticout.rottraj;
    
    jinput=automaticout.a(:,1:2);
    scrsz = get(0,'ScreenSize');
    realn = realn+1;
    xypairs =jinput';
    parameterize = linspace(0,2*pi,numel(xypairs)/2);
    pp = spline(parameterize,xypairs);
    yy = ppval(pp, linspace(0,2*pi,splineParam.nSplineCurvePoints));
    splinefit = [yy(1,:)', yy(2,:)'];
    
    
%% plotting    
    if ~splineParam.noPlots
        figure('Position',[10 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]),clf, subplot(3,2,[1 5]), plot(jtraj(:,1),jtraj(:,2), '- r'), axis equal, hold on, plot(jtraj(:,1),jtraj(:,2), '. blue'), plot(jinput(:,1),jinput(:,2),'or');
        xlabel('x (pixels)');
        ylabel('y (pixels)');
        title({['File ',savepath]; ['Trajectory ', int2str(thenums(j))]})
        mainfig = gcf;
        
        
        figure(mainfig)        
        hold on, plot(yy(1,:),yy(2,:),'. black', 'MarkerSize',5), axis equal        
        figure, plot(jtraj(:,1)-min(jtraj(:,1)),jtraj(:,2)-min(jtraj(:,2)), '- r'), axis equal, %xlim([-.5 3.5]),ylim('auto'), %hold on, plot(jtraj(:,1)-min(jtraj(:,1)),jtraj(:,2)-min(jtraj(:,2)), '. blue');
        xlabel('x (pixels)','FontSize', 16);
        ylabel('y (pixels)','FontSize', 16);
        hold on, plot(yy(1,:)-min(jtraj(:,1)),yy(2,:)-min(jtraj(:,2)),'. black', 'MarkerSize',3), axis equal
        set(gca, 'FontSize', 16);
        set(findobj(findobj(gca),'type','text'), 'FontSize', 16);
        saveas(gcf,[savepath,savename,'show10.pdf'],'pdf');
        saveas(gcf,[savepath,savename,'show10.fig'],'fig');
    end
 %%   
    outjtraj = jtraj;
    [numpointsj,~] = size(outjtraj);
    jsplinesout= zeros(numpointsj-1,2);
    
    jsigndistancesout = zeros(numpointsj-1,1);
    deltaparalleljout = zeros(numpointsj-1,2);
    alongjout = zeros(numpointsj-1,1);
    transversejout = zeros(numpointsj-1,1);
    deltaperpjout = zeros(numpointsj-1,2);

    for k = 0:5  % nonoverlapping displacements sampled at logarithmic intervals 
        
        
        jtrajk = jtraj(1:(2^k):end,:);
        [numpointsk,~] = size(jtrajk);
        jsigndistancesk = zeros(numpointsk,1);
        jsplinesk = zeros(numpointsk,2);
        
        
        
        deltaperpjk = zeros(numpointsk,2);
        deltaparalleljk = zeros(numpointsk,2);
        alongjk = zeros(numpointsk,1);
        
        transversejk = zeros(numpointsk,1);
        
        for i = 1:numpointsk-1
            
            
            alldisp=    [jtrajk(i,1)- splinefit(:,1), jtrajk(i,2)- splinefit(:,2)];
            alldispnext = [jtrajk(i+1,1)- splinefit(:,1), jtrajk(i+1,2)- splinefit(:,2)];
            thepoint =  find( alldisp(:,1).*alldisp(:,1)+alldisp(:,2).*alldisp(:,2) == min(alldisp(:,1).*alldisp(:,1)+alldisp(:,2).*alldisp(:,2)));
            thepointnext = find( alldispnext(:,1).*alldispnext(:,1)+alldispnext(:,2).*alldispnext(:,2) == min(alldispnext(:,1).*alldispnext(:,1)+alldispnext(:,2).*alldispnext(:,2)));
            if numel(thepointnext)>1
                thepointnext = max(thepointnext(thepointnext~=thepoint));
            end
            dispi = alldisp(thepoint,:);
            
            
            if thepoint == 1
                tangenti = splinefit(thepoint+1,:)-splinefit(thepoint,:);
            elseif thepoint == numel(splinefit)/2
                tangenti = splinefit(thepoint,:)-splinefit(thepoint-1,:);
            else
                tangenti = splinefit(thepoint+1,:)-splinefit(thepoint-1,:);
            end
            
            if k==0
                jsigndistancesk(i) = (sign(dispi(1)*tangenti(2)-dispi(2)*tangenti(1))*sqrt(sum(dispi.*dispi)));
                jsplinesk(i,:) = splinefit(thepoint,:);
                
                %%%%%%%%%%%%%% the full (output) distances and splines are
                %%%%%%%%%%%%%% based on the \Delta t = 1 frame data (i.e. k=0)
                jsigndistancesout(i) = (sign(dispi(1)*tangenti(2)-dispi(2)*tangenti(1))*sqrt(sum(dispi.*dispi)));
                jsplinesout(i,:) = splinefit(thepoint,:) ;
            else
                jsigndistancesk(i) = (sign(dispi(1)*tangenti(2)-dispi(2)*tangenti(1))*sqrt(sum(dispi.*dispi)));
                jsplinesk(i,:) = splinefit(thepoint,:);
            end
            
            
            %%%%%%%%%%%%%  here is the displacement vector between the ith
            %%%%%%%%%%%%%  point and the i+1 point...
            deltai = jtrajk(i+1,1:2) - jtrajk(i,1:2);
            %%%%%%%%%%%%%  and now the spline curve vector between the same two
            %%%%%%%%%%%%%  points
            deltasplinei = splinefit(thepointnext,:)- jsplinesk(i,:);
            
            
            %%%%%%%%%%%%%%%%%% 2009.04.27... there is a minor problem with NaNs
            %%%%%%%%%%%%%%%%%% here due to the rare case when the same spline fit
            %%%%%%%%%%%%%%%%%% point is the closest to 2 consecutive trajectory
            %%%%%%%%%%%%%%%%%% points.  This is why splineParam.nSplineCurvePoints 
            %%%%%%%%%%%%%%%%%% is set to 10000 for adequate resolution of
            %%%%%%%%%%%%%%%%%% the spline curve
            
            
           
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            if k==0
                deltaparalleljk(i,1:2) = (deltai*deltasplinei'/(deltasplinei*deltasplinei'))*deltasplinei;
                
                deltaperpjk(i,1) = deltai(1)-deltaparalleljk(i,1);
                deltaperpjk(i,2) = deltai(2)-deltaparalleljk(i,2);
                
                
                
                %%%%%%%%%%%%%% the full (output) parallel and perpendicular displacements are
                %%%%%%%%%%%%%% based on the \Delta t = 1 frame data (i.e. k=0)
                deltaparalleljout(i,1:2) = (deltai*deltasplinei'/(deltasplinei*deltasplinei'))*deltasplinei;
                
                deltaperpjout(i,1) = deltai(1)-deltaparalleljk(i,1);
                deltaperpjout(i,2) = deltai(2)-deltaparalleljk(i,2);
            else
                deltaparalleljk(i,1:2) = (deltai*deltasplinei'/(deltasplinei*deltasplinei'))*deltasplinei;
                
                deltaperpjk(i,1) = deltai(1)-deltaparalleljk(i,1);
                deltaperpjk(i,2) = deltai(2)-deltaparalleljk(i,2);
            end
            
            trajstep = deltai;
            splinestep = deltasplinei;
            dotcheck = deltaperpjk(i,:)*deltaparalleljk(i,:)';

            alongjk(i) = sign(thepointnext-thepoint)*(deltai*deltasplinei'/sqrt(deltasplinei*deltasplinei'));
            if k == 0
                alongjk(i) = sign(thepointnext-thepoint)*(deltai*deltasplinei'/sqrt(deltasplinei*deltasplinei'));
                aaaa = sign(deltaperpjk(i,1)*deltasplinei(2)-deltaperpjk(i,2)*deltasplinei(1))*sqrt(deltaperpjk(i,:)*deltaperpjk(i,:)');
                transversejk(i) = sign(thepointnext-thepoint)*aaaa;
                alongjout(i) = sign(thepointnext-thepoint)*(deltai*deltasplinei'/sqrt(deltasplinei*deltasplinei'));
                aaaa = sign(deltaperpjk(i,1)*deltasplinei(2)-deltaperpjk(i,2)*deltasplinei(1))*sqrt(deltaperpjk(i,:)*deltaperpjk(i,:)');
                transversejout(i) = sign(thepointnext-thepoint)*aaaa;
            else
                alongjk(i) = sign(thepointnext-thepoint)*(deltai*deltasplinei'/sqrt(deltasplinei*deltasplinei'));
                aaaa = sign(deltaperpjk(i,1)*deltasplinei(2)-deltaperpjk(i,2)*deltasplinei(1))*sqrt(deltaperpjk(i,:)*deltaperpjk(i,:)');
                transversejk(i) = sign(thepointnext-thepoint)*aaaa;
            end
            if isnan(transversejk(i))
                
            end
            
        end

        
        output=transversejk ;
        [ykpar,xkpar]  = hist(alongjk,linspace(-.6,.6,29));
        [ykperp,xkperp]= hist(transversejk,linspace(-.6,.6,29));

        
        
        
        
        if k== 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [disthisty, disthistx] =  hist(jsigndistancesk,25);
            [meanhist,widthhist] = stats(jsigndistancesk);
            [alonghist1,alonghist2] = hist(alongjk,31);
            [perphist1,perphist2]   = hist(transversejk,alonghist2);
            
            
            if ~splineParam.noPlots
                figure(mainfig), subplot(3,2, 2), plot(jtrajk(:,5), jtrajk(:,2)-jtrajk(1,2),'blue'), hold on, plot(jtrajk(:,5), jtrajk(:,1)-jtrajk(1,1), 'red'), hold off
                title('x- and y- trajectories');
                xlabel('time (frames)');
                ylabel('displacement (pixels)');
                
                
                
                
                figure(mainfig),subplot(3,2,4), bar(disthistx, disthisty)
                title({['Histogram of distances from spline   ']; [int2str(numpointsk),'  data points']})
                text('Position',[.03,.85],'String',['Width = ', num2str(2*widthhist),'pixels' ]  ,'Units','normalized','FontSize',7)
                xlabel('pixels')
                ylabel('N')
                
                figure,hold on, bar(disthistx, disthisty)
                title({['Histogram of distances from spline   ']; [int2str(numpointsk),'  data points']})
                text('Position',[.03,.85],'String',['Width = ', num2str(2*widthhist),'pixels' ]  ,'Units','normalized','FontSize',7)
                xlabel('distance from spline curve(pixels)')
                ylabel('number of positions')
                set(gca, 'FontSize', 16);
                set(findobj(findobj(gca),'type','text'), 'FontSize', 16);
                %
                saveas(gcf,[savepath,savename,'histogram10.pdf'],'pdf');
                saveas(gcf,[savepath,savename,'histogram10.fig'],'fig');
                
                
                figure(mainfig),subplot(3,2,6), plot( alonghist2,alonghist1,'*- r'), hold on, plot( perphist2,perphist1, '.- blue')
                title(['Displacements || (red) and \perp  (blue) to spline  '])
                xlabel('pixels')
                ylabel('N')

                
            end

        end
        
    end
    
    

    whichjsize = size(whichrowsj);
    rhs = size(jsplinesout(:,1));
    
    
    outjrows = whichrowsj(1:end-1);
    
    
    
    jsplinesout;
    out.traj(outjrows,incols+1) = jsplinesout(:,1);
    out.traj(outjrows,incols+2) = jsplinesout(:,2);
    out.traj(outjrows,incols+3) = jsigndistancesout;
    
    out.traj(outjrows,incols+4:incols+5) = deltaperpjout;
    out.traj(outjrows,incols+6:incols+7) = deltaparalleljout;
    out.traj(outjrows,incols+8) = alongjout;
    out.traj(outjrows,incols+9) = transversejout;
    
    
    
    % internalhist =   histanalysis(out.traj(outjrows,:));
    %
    %    gausswidthpar = num2str(internalhist(10));
    %    parstr1 = ['Gaussian width parallel = ',gausswidthpar];
    %    gfitpar = num2str(internalhist(11));
    %    parstr2 = ['Gaussian fit R^2 parallel = ', gfitpar];
    %    kurtosispar = num2str(internalhist(18));
    %    parstr3 = ['Kurtosis parallel = ',kurtosispar];
    %    widthpar = num2str(internalhist(15));
    %    parstr4 = ['Width parallel =', widthpar];
    %
    %    gausswidthperp= num2str(internalhist(4));
    %    perpstr1 = ['Gaussian width \perp = ', gausswidthperp];
    %    gfitperp = num2str(internalhist(5));
    %     perpstr2 = ['Gaussian fit R^2 \perp = ', gfitperp];
    %    kurtosisperp = num2str(internalhist(19));
    %    perpstr3 = ['Kurtosis \perp = ',kurtosisperp];
    %    widthperp = num2str(internalhist(17));
    %    perpstr4 = ['Width \perp =', widthperp];
    %
    %    figure(mainfig)
    % text('Position',[.03,.85],'String',{parstr1; parstr2; parstr3; parstr4},'Units','normalized','FontSize',7)
    % text('Position',[.97,.85],'String',{ perpstr1;perpstr2; perpstr3; perpstr4},'Units','normalized','FontSize',7, 'HorizontalAlignment','Right')
    
    if interactive
        usercheck = waitforbuttonpress;
        if usercheck
            
            usertype = str2num(get(gcf, 'CurrentCharacter'));
            if numel(usertype) ==0
                usertype = 1;
            end
        else
            usertype = 1;
        end
        if usertype == 0
            usertype =1;
        end
        
        toprint = {
            'Not used for subsequent analysis'
            'Category 1'
            'Category 2'
            'Category 3'
            };
        

        
        if ~splineParam.noPlots
        figure(mainfig),subplot(2,2,[1 3])
        xlabel('automatic')
        title(['Trajectory ', int2str(thenums(j))])
        end
        out.traj(whichrowsj,incols+10)= usertype;

        
        
        %%%%%%%%%%%%%%%%%   Usertype is set up to be zero (0) if the mouse is
        %%%%%%%%%%%%%%%%%   clicked during the histogram, and to return the NUMBER
        %%%%%%%%%%%%%%%%%   typed if a number is pressed.  .. tentatively:
        %%%%%%%%%%%%%%%%%       1   will be used for a useless
        %%%%%%%%%%%%%%%%%           trajectory (that somehow still got a spline)
        %%%%%%%%%%%%%       2   will be used for a "2-track" trajectory
        %%%%%%%%%%%%%       3   will be a circular region/geometry
        %%%%%%%%%%%%%       4   will be anything else that is worthy of a
        %%%%%%%%%%%%%           spline curve but not fitting into 2 or 3
        %savecheck = waitforbuttonpress;
        savecheck = 1;
        if savecheck
            tosave = getframe(mainfig);
            imwrite(tosave.cdata, [savepath,'Traj', int2str(thenums(j)),'out10.jpg']);
            saveas(gcf,[savepath,savename,'Traj', int2str(thenums(j)),'out10',toprint{usertype}(1:3),'.pdf'])
            saveas(gcf,[savepath,savename,'Traj', int2str(thenums(j)),'out10.fig'])
        end
        
        
        
        
        
   
        
        
    else
        
        pause(0.01)
        
        
        
        out.traj(whichrowsj,incols+10)= 10  ;
        %this assigns a value to the usertype column even when splineAnalysis2016 is run
        %without interactive mode.  The idea is to be able to identify
        %automatically-analyzed datasets.
        if ~splineParam.noPlots
            savecheck =1
            if savecheck
                tosave = getframe(mainfig);
                imwrite(tosave.cdata, [savepath,'Traj', int2str(thenums(j)),'out10.jpg']);
                saveas(gcf,[savepath,savename,'Traj', int2str(thenums(j)),'out10','AUTO','.pdf'])
                saveas(gcf,[savepath,savename,'Traj', int2str(thenums(j)),'out10.fig'])
            end
        end
        
    end
    

    clear('jsplines','jsigndistances','jtraj', 'jtrajk','whichrowsj','outjrows', 'automaticout')
    
    
    if nofigurebuildup
        close all
    end
   
    
end

out.ignoredlist = ignoredlist;
