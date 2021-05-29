function PlotFieldonMesh(coordinates,nodes,vari,varimin,varimax)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose:
%         To plot the profile of a component on mesh
% Synopsis :
%           ProfileonMesh(coordinates,nodes,component)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [X Y ] 
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]      
%           component - The components whose profile to be plotted
%           -----> components = a column vector in the order of node
%                               numbers
%--------------------------------------------------------------------------
%%%% Setting color map
% % c1 = zeros(6,1);c2 = ones(6,1);c33 = 0:0.2:1;c3 = c33';c44 = 1:-0.2:0;
% % c55 = 0.8:-0.2:0;c4 = c44';cc1 = [c1,c3,c4];cc2 = [c3,c2,c1];
% % cc3 = [ones(5,1),c55',zeros(5,1)];cmap = [cc1;cc2;cc3];

cmap = [0 0 1;0 0.25 0.75; 0 0.5 0.5;0 0.75 0.25; 0 1 0;0.25 1 0;0.5 1 0;0.75 1 0;
    1 1 0;1 0.75 0;1 0.5 0;1 0.25 0;1 0 0];

nel = length(nodes) ;                  % number of elements
nnel = size(nodes,2);                % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
profile = zeros(nnel,nel) ;
%
for iel=1:nel   
     for i=1:nnel
         nd(i)= nodes(iel,i);         % extract connected node for (iel)-th element
         X(i,iel)=coordinates(nd(i),1);    % extract x value of the node
         Y(i,iel)=coordinates(nd(i),2);    % extract y value of the node
     end   
     profile(:,iel) = vari(nd') ;         % extract component value of the node 
end
% % varimax = max(vari);varimin = min(vari);
% Plotting the FEM mesh and profile of the given component
% % fh = figure ;
% % set(fh,'name','Postprocessing','numbertitle','off') ;
fill(X,Y,profile,'EdgeColor','none');
% % axis off ;
% Colorbar Setting
% % %      SetColorbar;%% original idea of the authors
ncolor = size(cmap,1);
varistep = (varimax-varimin)/(ncolor);
colormap(cmap);
view(2);
caxis([varimin,varimax]);
axis equal;
h = colorbar;
set(h,'ytick',varimin:varistep:varimax);
set(gca,'FontSize',14);
end

              
         
 
   
     
       
       

