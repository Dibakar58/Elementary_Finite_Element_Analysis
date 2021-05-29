function PlotMesh(coordinates,nodes)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Finite Element Method Mesh
% Synopsis :
%           PlotMesh(coordinates,nodes) 
%--------------------------------------------------------------------------

nel = length(nodes) ;                  % number of elements
nnel = size(nodes,2);                % number of nodes per element

% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;

for iel=1:nel   
     for i=1:nnel
         ndi= nodes(iel,i);         % extract connected node for (iel)-th element
         X(i,iel)=coordinates(ndi,1);    % extract x value of the node
         Y(i,iel)=coordinates(ndi,2);    % extract y value of the node
     end
end
    
f1 = figure ;
set(f1,'name','Mesh','numbertitle','off','Color','w') ;
fill(X,Y,'w')   
axis off ;