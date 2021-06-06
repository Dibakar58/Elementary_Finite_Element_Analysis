function [coordinates,nodes,nel,nnode] = MeshRectangular(Lx,Ly,nx,ny)
% To Mesh a membrane using 4 noded Elements        |    

nel = nx*ny ;        % Total Number of Elements in the Mesh
nnel = 4 ;           % Number of nodes per Element
% Number of nodes on the Length and Breadth
npx = nx+1 ;
npy = ny+1 ;
nnode = npx*npy ;      % Total node number
% Discretizing the Length and Breadth of the membrane
nx = linspace(0,Lx,npx) ;
ny = linspace(0,Ly,npy) ;
[xx,yy] = meshgrid(nx,ny);
coordinates = [xx(:), yy(:)];

% To get the Nodal Connectivity Matrix
NodeNo = 1:nnode ;
nodes = zeros(nel,nnel) ;

% small code for selecting nodal point number for specific element
    NodeNo = reshape(NodeNo,npy,npx);
    nodes(:,1) = reshape(NodeNo(1:npy-1,1:npx-1),nel,1);
    nodes(:,2) = reshape(NodeNo(1:npy-1,2:npx),nel,1);
    nodes(:,3) = reshape(NodeNo(2:npy,2:npx),nel,1);
    nodes(:,4) = reshape(NodeNo(2:npy,1:npx-1),nel,1);

% Plotting the Finite Element Mesh
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
% Extract X,Y coordinates for the (iel)-th element
  for iel = 1:nel
      X(:,iel) = coordinates(nodes(iel,:),1) ;
      Y(:,iel) = coordinates(nodes(iel,:),2) ;
  end
nnode = size(coordinates,1) ;        
nel = size(nodes,1) ;                
