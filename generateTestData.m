

for i = 1:5
trajTest(((i-1)*100+1):(i*100),1) = 10*rand(100,1)+i*30*ones(100,1)+[1:100]'
end
for i = 1:5
trajTest(((i-1)*100+1):(i*100),2) = 10*rand(100,1)+8*sin(.1*[1:100]'+rand(1))+mod(i,2)*3
end

myParam.thetaRange = pi/4
myParam.nTrajShapePoints   = 10
myParam.minTrajLength = 100
myParam.nSplineCurvePoints = 10000
myParam.noPlots = true
myParam.thetaMax = .2*pi
myParam.radiusFactor = 3
myParam.radiusRatio = 2

splineanalysis2016(trajTest,'exampleDirectory', 'test1',0, myParam)
