function out = aspectratio(inputtraj)


%%%%%%%%%%%%%%this function calculates the aspect ratio ( >1 ) of a
%%%%%%%%%%%%%%trajectory inputted in [xvals, yvals] format.

out.traj = inputtraj;
[out.coms, out.sigs] =stats(out.traj);
%%%%% subtract the center of mass from each trajectory
out.traj(:,1) = out.traj(:,1) -out.coms(1);
out.traj(:,2) = out.traj(:,2) -out.coms(2);
 

%%%% now rotate the trajectory through 360 degrees and calculate the width
%%%% in x and y in each case....


for degs = 1:180
    
    inrads = degs*(2*pi/360);
    rotated = out.traj(:,1:2)*[cos(inrads) sin(inrads); -sin(inrads) cos(inrads)];
    [crap, out.widths(degs,1:2)] = stats(   rotated);
    out.minmax(degs,1:2) = [max(rotated(:,1))-min(rotated(:,1)),max(rotated(:,2))-min(rotated(:,2))];
    
end

out.aspectr = max(out.widths(:,1)./out.widths(:,2));
out.aspectrm = max(out.minmax(:,1)./out.minmax(:,2));


%%%%%%%%%%%%%%%%%  now adding a little bit of code to calculate the
%%%%%%%%%%%%%%%%%  correlation between adjacent displacements...

%%%%%%%%%%%%%%%%%  this is the mean correlation between \Delta_i and
%%%%%%%%%%%%%%%%%  \Delta_i+1...  should be zero for a true random walk,
%%%%%%%%%%%%%%%%%  although I think the expectation value probably scales
%%%%%%%%%%%%%%%%%  like 1/sqrt(N) for N steps in the trajectory.

[rs,cs] = size(out.traj);
test = sum(diff(out.traj(1:end-1,1:2)).*diff(out.traj(2:end,1:2)));
out.corrx = test(1)/(rs-1);
out.corry = test(2)/(rs-1);

