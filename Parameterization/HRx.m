% calculate homogeneous transformation matrix for rotation about the x-axis

function HRx = HRx(theta)

HRx = [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];

end