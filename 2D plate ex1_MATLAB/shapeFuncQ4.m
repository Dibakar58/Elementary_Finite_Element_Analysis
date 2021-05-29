function [shape,nDeriva]= shapeFuncQ4(xi,eta)
    % shape function and derivatives for Q4 elements
    % shape : Shape functions
    % nderiva: derivatives with respect to xi and eta 
    % xi, eta: natural coordinates (-1 ... +1)
    shape = 1/4*[(1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];
            
    nDeriva = 1/4*[-(1-eta),    -(1-xi);
                    1-eta,      -(1+xi);
                    1+eta,      1+xi;
                    -(1+eta),   1-xi];
end
