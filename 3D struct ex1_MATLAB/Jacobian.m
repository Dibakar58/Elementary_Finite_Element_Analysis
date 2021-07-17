function [J,xyDeriva] = Jacobian(nodeCoord,nDeriva)
    J = nodeCoord'*nDeriva;                   
    xyDeriva = nDeriva/J;%%% xyDeriva means: d/dx; d/dy; d/dz
end
