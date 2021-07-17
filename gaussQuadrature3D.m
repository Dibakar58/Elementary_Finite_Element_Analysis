function [weights,locations] = gaussQuadrature3D(option)
    %%%% Gauss quadrature for eight­nodal hexahedron elements
    %%%% option 'complete'(2x2x2)
    %%%% option 'reduced'(1x1)
    %%%% locations: Gauss point locations
    %%%% weights: Gauss point weights
    switch option
        case 'complete'
        locations=...
          1/sqrt(3)*[ -1, -1, -1;
                       1, -1, -1;
                       1,  1, -1;
                      -1,  1, -1;
            
                      -1, -1, 1;
                       1, -1, 1;
                       1,  1, 1;
                      -1,  1, 1];
                  
        weights = [1;1;1;1;1;1;1;1]; 
        case 'reduced'
        locations =[0 0 0];
        weights = 8;
    end
end
    