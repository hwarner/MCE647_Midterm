function T = hmc3_transrotZ(translation, theta)
% Generates a homogeneous transformation matrix for a 3D translation followed
% by a rotation about the Z-axis.
%
%   T = hmc3_transrotZ(translation, theta)
%
% Input:
%	translation		(3 x 1 matrix) translation vector
%	theta			(scalar) amount of rotation in radians (positive is screw towards +Z)
%
% Output:
%	T		(4x4 matrix) homogeneous transformation matrix

	T = [	cos(theta) -sin(theta)	0	translation(1); ...
			sin(theta)  cos(theta)	0	translation(2); ...
				0			0		1	translation(3); ...
				0			0		0		1		  ];

end