function obj=objFunc(x)

global dataForPlotting objForPlotting

centerPts = equationSet(x);  %Calculate the center points


obj=sum(std((centerPts)').^2);  


dataForPlotting=[dataForPlotting;x];
objForPlotting=[objForPlotting obj];
