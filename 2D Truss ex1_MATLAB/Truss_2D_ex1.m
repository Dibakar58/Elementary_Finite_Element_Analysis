clear
%% %%% Step 1. Material properties and geometricall properties
E1 = 2e11;E2 = 70e9;%%% Elastic modulus

A1 = 0.1; A2 = 0.2;%%%% Cross-sectional areas
A = [A1;A2];

L1 = 1; L2 = 1;%%%% bars' lengths

%%%%%%%% Stress-strain relation: S = C.E
C1 = E1; C2 = E2;
C = [C1;C2];

%% %%%% Step 2. Create mesh
%%%%%%% 2.1. generate nodal coordinates
xy1 = [0, 0];
xy2 = [L1*cos(pi/4), L1*sin(pi/4)];
xy3 = [0, (L1+L2)*sin(pi/4)];

xy = [xy1; xy2; xy3];%%% combine all nodes

nP = size(xy,1);%%%% number of nodes: nP = 3

%%%%%%% 2.2. nodal conectivity for elements
nE = nP-1;%%% number of elements
eNodes = zeros(nE,2);
eNodes(1,:) = [1,2];%%%% element 1 has nodes 1 and 2
eNodes(2,:) = [2,3];%%%% element 2 has nodes 2 and 3


% % %%%%%% check geometry:
% % figure
% % plot(xy(:,1),xy(:,2),'k:','LineWidth',2);
% % axis equal

%% %%%%%% Step 3. Calculate stiffness matrix: K
nDof = 2*nP;
K = zeros(nDof,nDof);%%% Gloabl stiffness K has size of nPxnP

for i = 1:nE
    %%%% nodes of element
    eNodei = eNodes(i,:);
    
    eDof = [eNodei;eNodei+nP];
    nnDof = length(eNodei);
    
    %%% Length of element
    Le = sqrt((xy(eNodei(2),1)-xy(eNodei(1),1))^2 + (xy(eNodei(2),2)-xy(eNodei(1),2))^2);
    cos = (xy(eNodei(2),1)-xy(eNodei(1),1))/Le;
    sin = (xy(eNodei(2),2)-xy(eNodei(1),2))/Le;
    
    detJ = Le/2;
    invJ = 1/detJ;
    
    %%%%% Central gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% B matrix
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %%%%% Element stiffness matrix
    Ke = A(i,1)*(B'*C(i,1)*B)*detJ*gaussWt;
    
    %%%%% Transformation matrix
    Te = [cos, sin, 0,   0;
           0,   0, cos, sin];
    KeG = Te'*Ke*Te;
    
    K(eDof,eDof) = K(eDof,eDof) + KeG;
end

%% %%% Boundary conditions
fixedP = [1;3];
fixedDof = [fixedP;fixedP+nP];%%%% fixedDof = [u1;u3;v1;v3]

%% %%% Loading condition
F1 = 5e9; F2 = 3e9;

force = zeros(nDof,1);

force(2,1) = F1;%%%%% at node 2: Fx = F1
force(2+nP,1) = F2;%%%% at node 2: Fy = F2

%% %%%% Solution
disp = solution(nDof,fixedDof,K,force);

save('Truss_2D_Ex2.mat');


%% %%%% Post processing
u_dof = (1:nP)';
v_dof = u_dof + nP;
u = disp(u_dof,1);
v = disp(v_dof,1);

scale = 1;%%% true scale
xy_newFEM = xy + scale*[u,v];%%% deformed configuration


%%%%%%%% Now, we need to impart ANSYS results
error_u = (u(2,:)-dispANSYS(2,1))/dispANSYS(2,1)*100;
error_v = (v(2,:)-dispANSYS(2,2))/dispANSYS(2,2)*100;

xy_ansys = NLIST(:,1:2);
xy_newANSYS = xy_ansys + scale*dispANSYS;


%% %%%%%% Plot deformed configurations in ANSYS and FEM code
figure
plot(xy(:,1),xy(:,2),'k:','LineWidth',2);
hold on
plot(xy_newFEM(:,1),xy_newFEM(:,2),'b--','LineWidth',2);
hold on
plot(xy_newANSYS(:,1),xy_newANSYS(:,2),'r:','LineWidth',2);

xlabel('x (m)');
ylabel('displacement: u (m)');
grid on
view(2)
legend('Undeformed shape', 'Deformed shape - FEA code','Deformed shape - ANSYS')
legend('boxoff')
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')

%% %%%%%% Calculate strain and stress
strain_ele = zeros(nE,1);
stress_ele = zeros(nE,1);

strain_node = zeros(nE,2);
stress_node = zeros(nE,2);

for i = 1:nE
    %%%% nodes of element
    eNodei = eNodes(i,:);
    
    eDof = [eNodei;eNodei+nP];
    nnDof = length(eNodei);
    
    %%% Length of element
    Le = sqrt((xy(eNodei(2),1)-xy(eNodei(1),1))^2 + (xy(eNodei(2),2)-xy(eNodei(1),2))^2);
    cos = (xy(eNodei(2),1)-xy(eNodei(1),1))/Le;
    sin = (xy(eNodei(2),2)-xy(eNodei(1),2))/Le;
    
    detJ = Le/2;
    invJ = 1/detJ;
    
    %%%%% Central gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% Strain-displacement matrix: e = B*u 
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %%%%% Transformation matrix
    Te = [cos, sin, 0,   0;
           0,   0, cos, sin];
       
    %%%%% Transform global displacements back to local displacements
    eDof_v = [eDof(:,1);eDof(:,2)];%%%% see slide 4 for the order of displacement vector
    dispL = Te*disp(eDof_v,1);
    %%%%% strain for element and nodes
    ei = B*dispL;
    
    strain_ele(i,:) = ei;
    strain_node(i,:) = ei;
    
    %%%%% stress for elements and nodes
    Ci = C(i,1);
    
    stress_ele(i,:) = Ci*ei;
    stress_node(i,:) = Ci*ei;
end

%%%%%%%%% FINISH!! %%%%%%%%%%%%%%
%%%%%%%%% Please save data if you want to use it later %%%%%%%%%%
