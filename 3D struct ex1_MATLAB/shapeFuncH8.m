function [shape,nDeriva]= shapeFuncH8(x_L,y_L,z_L)
%%% shape function and derivatives for H8 (cubic) elements

    shape = 1/8*[(1-x_L)*(1-y_L)*(1-z_L);
                (1+x_L)*(1-y_L)*(1-z_L);
                (1+x_L)*(1+y_L)*(1-z_L);
                (1-x_L)*(1+y_L)*(1-z_L);
                
                (1-x_L)*(1-y_L)*(1+z_L);
                (1+x_L)*(1-y_L)*(1+z_L);
                (1+x_L)*(1+y_L)*(1+z_L);
                (1-x_L)*(1+y_L)*(1+z_L)];
            
    nDeriva = 1/8*[-(1-y_L)*(1-z_L),-(1-x_L)*(1-z_L), -(1-x_L)*(1-y_L);
                    (1-y_L)*(1-z_L),-(1+x_L)*(1-z_L), -(1+x_L)*(1-y_L);
                    (1+y_L)*(1-z_L), (1+x_L)*(1-z_L), -(1+x_L)*(1+y_L);
                    -(1+y_L)*(1-z_L), (1-x_L)*(1-z_L), -(1-x_L)*(1+y_L);
                    
                    -(1-y_L)*(1+z_L),-(1-x_L)*(1+z_L), (1-x_L)*(1-y_L);
                    (1-y_L)*(1+z_L),-(1+x_L)*(1+z_L), (1+x_L)*(1-y_L);
                    (1+y_L)*(1+z_L), (1+x_L)*(1+z_L), (1+x_L)*(1+y_L);
                    -(1+y_L)*(1+z_L), (1-x_L)*(1+z_L), (1-x_L)*(1+y_L)];
end
