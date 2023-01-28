function X=randomShape(center,size,numPoints)
    Xc=(rand(10,2)-0.5)*size+center; % rand control points seed
    Xc = Xc(boundary(Xc), :);
    Xc(length(Xc) - 1, :) = 2 * Xc(1, :) - Xc(2, :);
    X = bspl(Xc, numPoints);
end