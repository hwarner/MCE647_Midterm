% This is the Main Script for running the "parameterization" optimization
% algorithm.


clear; close all; clc

% Later we will need the PUMA robot's transformation matrix. Create a
% function that contains it if it does not exist.
if exist('transA5_0.m')~=2
   makeDHFunctions; 
end


if 0    %Set to 0 for our measured data, 1 for Dr. Richters Data
    theta1 = [-.379; -.227; -0.513; -0.478];
    theta2 = [-.642; -.884; .495; 0.519];
    theta3 = [.520; .914; -.759; -0.672];
    
else
    theta1 = [-0.74	-0.775	-0.844	-0.911	-0.999	-1.029	-1.02	-0.977	-0.919	-0.832	-0.754];
    theta2 = [-0.04	0.044	-1.12	-0.635	-0.636	-0.853	-0.328	-0.22	-0.796	0.381	-1.442];
    theta3 = [-0.349	-0.534	0.832	0.198	0.198	0.526	-0.076	-0.144	0.704	-0.807	1.457];
    
     %Limit which measurements are used
     %theta1=theta1(1:2:9);
     %theta2=theta2(1:2:9);
     %theta3=theta3(1:2:9);
end

numPoints=length(theta1);  %Number of measurement points

%% Setup the [initialValue lowerbound upperbound]
initAndBounds.gammaR=[3 0 3];
initAndBounds.phiA=[0 -pi pi];
initAndBounds.psiA=[-pi/2 -pi pi];
initAndBounds.betaC=[0 -pi pi];

%For gamma and betaC, we need to create entries for each point (gammaR1,
%   gammaR2....)
for k=1:numPoints
    initAndBounds.(['gammaR' num2str(k)]) = initAndBounds.gammaR;
    initAndBounds.(['betaC' num2str(k)]) = initAndBounds.betaC;
end

%% Create the symbolic variables
syms cx cy cz phiA psiA xp yp zp
x = sym('x',[numPoints,1]);
y = sym('y',[numPoints,1]);
z = sym('z',[numPoints,1]);
betaC = sym('betaC',[numPoints,1]);
gammaR = sym('gammaR',[numPoints,1]);
r = 0.33;


%% Create the equations
for i=1:numPoints
    
    c=circleParam(betaC(i),r,cx,cy,cz,phiA,psiA);  %Symbolic calculation of points on circle
    c(1)=c(1)-xp;
    c(2)=c(2)-yp;
    c(3)=c(3)-zp;
    c(1)=solve(c(1),cx);
    c(2)=solve(c(2),cy);
    c(3)=solve(c(3),cz);
    
    A5_0{i}=transA5_0(theta1(i),theta2(i),theta3(i));
    l=lineParam(A5_0{i},gammaR(i));
    
    eqSet(1,i)=subs(c(1),xp,l(1));
    eqSet(2,i)=subs(c(2),yp,l(2));
    eqSet(3,i)=subs(c(3),zp,l(3));
       
end

%% Make a function/file (equationSet.m) with the equations in it.
%It will be evalued by the objective function (objFunc.m)

%Get the list of variables sorted as the symbolic toolbox uses them
varList=symvar(eqSet);
numVar=length(varList);

%Create the equation function evaluated in the objFunc
matlabFunction(eqSet,'file','equationSet','vars',{varList});


%%  Setup and run the opmization

% Convert the structure containg the inital, upper bound, and initial bound
% into vectors
for g=1:numVar
   ib=initAndBounds.(char(varList(g)));
   initGuess(g)=ib(1);
   lb(g)=ib(2);
   ub(g)=ib(3);
end

%Set up a global to write data to for plotting
global dataForPlotting objForPlotting
dataForPlotting=[];
objForPlotting=[];

%Run optimization
opts = optimset('Display','on','MaxFunEvals',100000,'MaxIter',100000);
tic; 
[optInput objfinalValue] = fmincon(@objFunc,initGuess,[],[],[],[],lb,ub,[],opts);
display(['It took ' num2str(toc) 'seconds to perform the optimization.'])

% Get the final center of all the circles
locCenterAllCircles=equationSet(optInput);
locCenter=mean(locCenterAllCircles,2);  %Use the mean center location

% Uncomment this line to get only the final result
dataForPlotting=dataForPlotting(end,:);

%Make variables for each of the symbolics
m=unpackSymbolics(dataForPlotting,varList);

% Plot the results
M=plotResults(numPoints,dataForPlotting,m,r,A5_0);

%Output the results to the command line
locCenter, m.phiA, m.psiA,objfinalValue


