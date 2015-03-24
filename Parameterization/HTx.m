% calculate homogeneous transformation matrix for translation along the x-axis

function HTx = HTx(d)

HTx = [1 0 0 d; 0 1 0 0; 0 0 1 0; 0 0 0 1];

end