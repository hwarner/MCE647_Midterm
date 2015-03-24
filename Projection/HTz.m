% calculate homogeneous transformation matrix for translation along the z-axis

function HTz = HTz(d)

HTz = [1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];

end