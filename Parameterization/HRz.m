% calculate homogeneous transformation matrix for rotation about the z-axis

function HRz = HRz(theta)

HRz = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];

end