function [stiff,mass] = formStiffness3D(nDof,nE,elNode,nP,xyz,C,dens)

stiff = zeros(nDof);
mass = zeros(nDof);

%%% 2x2x2 quadrature
[gaussWt,gaussLoc] = gaussQuadrature3D('complete');

for e = 1:nE                           
  id = elNode(e,:); 
  elDof = [id, id+nP, id+2*nP];
  ndof = length(id);%%% ndof = 8 for H8
  
  %%% Loop for Gauss point
  for q = 1:size(gaussWt,1)                      
    GaussPoint = gaussLoc(q,:);                                                     
    x_L = GaussPoint(1);
    y_L = GaussPoint(2);
    z_L = GaussPoint(3);
    
    %%% shape functions and derivatives
    [shape,nDeriva]= shapeFuncH8(x_L,y_L,z_L);

    %%% Jacobian matrix, inverse of Jacobian    
    [J,xyzDeriv] = Jacobian(xyz(id,:),nDeriva);%%% 1
    
    %%% B matrix (Linear strain - displacement matrix)%%% ???
    B = zeros(6,3*ndof);%%% B(6x24)
    B(1,1:ndof)             = xyzDeriv(:,1)';
    B(2,ndof+1:(2*ndof))    = xyzDeriv(:,2)';
    B(3,2*ndof+1:(3*ndof))  = xyzDeriv(:,3)';
    
    B(4,1:ndof)             = xyzDeriv(:,2)';
    B(4,ndof+1:(2*ndof))    = xyzDeriv(:,1)';
    
    B(5,ndof+1:(2*ndof))    = xyzDeriv(:,3)';
    B(5,2*ndof+1:(3*ndof))  = xyzDeriv(:,2)';
    
    B(6,2*ndof+1:(3*ndof))  = xyzDeriv(:,1)';
    B(6,1:ndof)             = xyzDeriv(:,3)';
    
    %%% stiffness matrix
    stiff(elDof,elDof)= stiff(elDof,elDof)+ B'*C*B*gaussWt(q)*det(J);    
    
    %%%% mass matrix
    mass(id,id) = mass(id,id) + shape*shape'*dens*gaussWt(q)*det(J);
    mass(id+nP,id+nP)= mass(id+nP,id+nP)+ shape*shape'*dens*gaussWt(q)*det(J);
    mass(id+2*nP,id+2*nP)= mass(id+2*nP,id+2*nP)+ shape*shape'*dens*gaussWt(q)*det(J);
  end
end
