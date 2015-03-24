% Holly Warner
% MCE 647 Midterm Takehome
% plot_bot.m

clear theta1 theta2 theta3

d1 = .666;
d2 = .2435;
d3 = .0934;
a2 = .4318;
b  = .270;
k  = .0521;

L(1) = Link('d', d1, 'a', 0, 'alpha', -pi/2);
L(2) = Link('d', -d2, 'a', a2, 'alpha', 0);
L(3) = Link('d', d3, 'a', 0, 'alpha', 0);
L(4) = Link('d', k, 'a', 0, 'alpha', 0);
L(5) = Link('d', 0, 'a', b, 'alpha', 0);


theta1Array = [-0.74  -0.775 -0.844 -0.911 -0.999 -1.029 -1.02  -0.977 -0.919 -0.832 -0.754];
theta2Array = [-0.04   0.044 -1.12  -0.635 -0.636 -0.835 -0.328 -0.22  -0.796  0.381 -1.442];
theta3Array = [-0.349 -0.534  0.832  0.198  0.198  0.526 -0.076 -0.144  0.704 -0.807  1.457];

if 0 % Set to 0 for all data points, 1 for every other
    for i = 1:round(length(theta1Array)/2)
        theta1(i) = theta1Array(i*2-1);
        theta2(i) = theta2Array(i*2-1);
        theta3(i) = theta3Array(i*2-1);
    end
    
else
    
    theta1 = theta1Array;
    theta2 = theta2Array;
    theta3 = theta3Array;
    
end

% add center point
theta1(length(theta1)+1) = -0.86;
theta2(length(theta2)+1) = -1.576;
theta3(length(theta3)+1) =  1.45;

project_bot = SerialLink(L, 'name', 'project\_bot');

for i = 1:length(theta1)
    project_bot.plot([theta1(i) theta2(i) theta3(i) 0 0], 'workspace', [-2 2 -2 2 0 4])
    pause(5)
end

