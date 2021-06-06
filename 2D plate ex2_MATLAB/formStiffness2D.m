function [stiff,mass] = formStiffness2D(nDof,nel,elNode,nn,nodeCoord,C,dens,thick)
% compute stiffness matrix (and mass matrix)
% for plane stress Q4 elements
stiff = zeros(nDof);
mass = zeros(nDof);
%%% 2 by 2 quadrature
[gaussWt,gaussLoc] = gaussQuadrature('complete');
for e = 1:nel                           
  id = elNode(e,:); 
  elDof = [id id+nn];
  ndof = length(id);%%% ndof = 4 for Q4
  %%% loop for Gauss point
  for q = 1:size(gaussWt,1)                      
    GaussPoint = gaussLoc(q,:);                                                     
    xi = GaussPoint(1);
    eta = GaussPoint(2);
    %%% shape functions and derivatives
    [shape,nDeriva] = shapeFuncQ4(xi,eta);
    %%% Jacobian matrix, inverse of Jacobian    
    [J,xyDeriva] = Jacobian(nodeCoord(id,:),nDeriva);
    %%%  BL matrix (Linear strain - displacement matrix)
    B = zeros(3,2*ndof);%%% B(3x8)
    B(1,1:ndof)             = xyDeriva(:,1)';
    B(2,ndof+1:(2*ndof))    = xyDeriva(:,2)';
    B(3,1:ndof)             = xyDeriva(:,2)';
    B(3,ndof+1:(2*ndof))    = xyDeriva(:,1)';
    %%% stiffness matrix
    stiff(elDof,elDof)= stiff(elDof,elDof)+ B'*C*B*thick*gaussWt(q)*det(J);    
    %%%% mass matrix
    mass(id,id) = mass(id,id) + shape*shape'*dens*thick*gaussWt(q)*det(J);
    mass(id+nn,id+nn)= mass(id+nn,id+nn)+ shape*shape'*dens*thick*gaussWt(q)*det(J);
  end
end
