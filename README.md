# splineAnalysis

splineAnalysis is a MATLAB package used to analyze particle trajectories confined to curvilinear structures. 

# Introduction
The motion of biomolecular complexes at the subcellular level takes place in a highly hetergeneous environment. Early versions of this code were written to analyze trajectories confined to curvilinear structures, where traditional 2D mean-square-displacement analysis underestimates rates of diffusion. The current package is a cleaned up version of the original code, restructured to allow other researchers (with a little MATLAB experience) to use this tool.

# User Guide
The best place to start is by running a basic script [generateTestData.m](https://github.com/berl/splineAnalysis/blob/master/generateTestData.m) analyzing some artificial data of 5 diffusing trajectories confined to a sinusoid. This script illustrates a few things:
- How your trajectory data needs to be formulated before analysis.  The format we've chosen to use (after Crocker and Grier) is described in more detail in [splineAnalysis2016.m](https://github.com/berl/splineAnalysis/blob/master/splineanalysis2016.m). 
- Several analysis parameters are combined into the splineParam struct, and the values shown in [generateTestData.m](https://github.com/berl/splineAnalysis/blob/master/generateTestData.m) are chosen to illustrate how the code works on the example data.

## Analysis Parameters in `splineParam`
- `myParam.thetaRange = pi/4` The angle range searched forward along the trajectory curve.
- `myParam.nTrajShapePoints   = 10` The trajectory is spatially divided into this many clusters to provide starting points for the spline curve
- `myParam.minTrajLength = 100` A filter in [preshape2016.m](https://github.com/berl/splineAnalysis/blob/master/preshape2016.m) throws out trajectories shorter than this.
- `myParam.nSplineCurvePoints = 10000`  This sets the number of points in the resampled spline curve.  
- `myParam.noPlots = true` The intial version of this code included a lot of automatically-saved plots.  They're still in the code but may or may not be useful for a given data set.  See also the notes below about `interactiveMode`.
- `myParam.thetaMax = .2*pi` This is the largest angle acceptable between spline guide points.  
- `myParam.radiusFactor = 3` The initial search radius is this the average transverse trajectory width multiplied by this factor.
- `myParam.radiusRatio = 2`  The secondary search is the initial radius multiplied by this factor.

These parameter values work well for the example dataset (and for our original experimental data), but may be different for other data sets.

## InteractiveMode
During our data analysis, we utilized an interactive mode with several automatically-generated and saved figures and an option to manually classify analyzed trajectories by mouseclick or keyboard. See [splineAnalysis2016.m](https://github.com/berl/splineAnalysis/blob/master/splineanalysis2016.m) around `line 415` for details.

## Project Status 
Although this project is currently not an active area of development, I will entertain pull requests, suggestions and questions.




## Publication
 B.R. Long and T.Q. Vu 'Spatial structure and diffusive dynamics from single-particle trajectories using spline analysis'
 Biophys J. 2010 Apr 21;98(8):1712-21.
