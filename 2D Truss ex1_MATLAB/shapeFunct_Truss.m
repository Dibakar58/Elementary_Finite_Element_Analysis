function [shape,nderiv] = shapeFunct_Truss(xi)
%%%% shapeN : Shape functions N1 and N2
%%%% N_deriv: derivatives of N1 and N2 w.r.t. xi 
%%%% xi: natural coordinates (-1 ... +1)

shape =1/2*[1-xi;1+xi];

nderiv = [-1;1]/2;

end
  
