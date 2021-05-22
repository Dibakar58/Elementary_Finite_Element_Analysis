clear

%% %%% Step 1. Material properties and geometrical properties
E0 = 2e11;%%% Elastic modulus

A0 = 0.1;%%%% Cross-sectional area
L0 = 0.5;%%%% Length of each truss

nE = 16; %%% 16 elements
%%%% Create a vector store cross-sectional areas for all trusses
A = A0*ones(nE,1);
%%%% Create a vector store lengths for all trusses
L = L0*ones(nE,1);

%%%%%%%% Stress-strain relation: S = C.E
C = E0*ones(nE,1);

%% %%%% Step 2. Create mesh
%%%%%%% 2.1. generate nodal coordinates
P1 = [0, 0];
P2 = [L0, 0];
P3 = [2*L0, 0];
P4 = [3*L0, 0];

P5 = [2*L0 + L0*cos(pi/3), L0*sin(pi/3)];
P6 = [L0 + L0*cos(pi/3), L0*sin(pi/3)];
P7 = [L0*cos(pi/3), L0*sin(pi/3)];

P8 = [L0, 2*L0*sin(pi/3)];
P9 = [2*L0, 2*L0*sin(pi/3)];

%%%%%% nodal coordinates: xy
xy = [P1; P2; P3; P4; P5; P6; P7; P8; P9];
nP = size(xy,1);%%%% number of nodes: nP = 9

%%%%%%% 2.2. Nodal conectivity for elements
eNodes = zeros(nE,2);
eNodes(1,:) = [1,2];%%%% element 1: 1-2
eNodes(2,:) = [2,3];%%%% element 2: 2-3
eNodes(3,:) = [3,4];%%%% element 3: 3-4
eNodes(4,:) = [7,6];%%%% element 4: 7-6
eNodes(5,:) = [6,5];%%%% element 5: 6-5
eNodes(6,:) = [8,9];%%%% element 6: 8-9
eNodes(7,:) = [1,7];%%%% element 7: 1-7
eNodes(8,:) = [7,2];%%%% element 8: 7-2
eNodes(9,:) = [2,6];%%%% element 9: 2-6
eNodes(10,:) = [6,3];%%%% element 10: 6-3
eNodes(11,:) = [3,5];%%%% element 11: 3-5
eNodes(12,:) = [5,4];%%%% element 12: 5-4
eNodes(13,:) = [7,8];%%%% element 13: 7-8
eNodes(14,:) = [8,6];%%%% element 14: 8-6
eNodes(15,:) = [6,9];%%%% element 15: 6-9
eNodes(16,:) = [9,5];%%%% element 16: 9-5


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
    
    %%%%% Strain-displacment matrix: B
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %%%%% Element stiffness matrix in local coordinate system
    Ke = A(i,1)*(B'*C(i,1)*B)*detJ*gaussWt;
    
    %%%%% Transformation matrix
    Te = [cos, sin, 0,   0;
           0,   0, cos, sin];
       
    %%%%% Element stiffness matrix in the global coordinate system
    KeG = Te'*Ke*Te;
    
    %%%%% Global stiffness matrix
    K(eDof,eDof) = K(eDof,eDof) + KeG;
end

%% %%% Boundary conditions
fixedDof = [1;1+nP;4+nP];%%%%% u1 = v1 = 0; v4 = 0

%% %%% Loading conditions
F1 = -1e9; F2 = -1e9;%%% forces are in negatve y-direction
force = zeros(nDof,1);

force(8+nP,1) = F1;%%%% force(8+nP,1) means Fy at node 8
force(9+nP,1) = F2;%%%% force(9+nP,1) means Fy at node 9

%% %%%% Solution
disp = solution(nDof,fixedDof,K,force);



%% %%%% Post processing
u_dof = (1:nP)';
v_dof = u_dof + nP;

scale = 1;
xy_newFEM = xy + scale*[disp(u_dof,1),disp(v_dof,1)];

xy_ansys = NLIST(:,1:2);
xy_newANSYS = xy_ansys + scale*dispANSYS;

%%%%%%%% Plot deformed shape
figure
for i = 1:nE
    eNodei = eNodes(i,:);
    h1 = plot(xy(eNodei,1),xy(eNodei,2),'k:','LineWidth',2);
    hold on
    h2 = plot(xy_newFEM(eNodei,1),xy_newFEM(eNodei,2),'b--','LineWidth',2);
    hold on
    eNodeiANSYS = ELIST(i,7:8);
    h3 = plot(xy_newANSYS(eNodeiANSYS,1),xy_newANSYS(eNodeiANSYS,2),'r:','LineWidth',2);
    hold on
end
axis([0 4*L0 -0.2 1]);
xlabel('x (m)');
ylabel('y (m)');
grid on
view(2)
legend('Undeformed', 'Deformed-FEM code')%,'Deformed - ANSYS')
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
    
    %%%%% B matrix
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %%%%% Transformation matrix
    Te = [cos, sin, 0,   0;
           0,   0, cos, sin];
       
    %%%%% Transform glocal displacements back to local displacements
    eDof_v = [eDof(:,1);eDof(:,2)];% order of displacement vector
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

%%%%%% Find element that has maximum strain, stress
id_max = find(strain_ele(:,1)==max(strain_ele(:,1)));
% % % id_max = find(stress_ele(:,1)==max(stress_ele(:,1)));%%% give the same result
e_max = strain_ele(id_max,1);
s_max = stress_ele(id_max,1);

%%%%%%%%% FINISH!! %%%%%%%%%%%%%%
%%%%%%%%% Please save data if you want to use it later %%%%%%%%%%
