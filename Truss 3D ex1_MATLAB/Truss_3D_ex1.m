clear

%% %%% Step 1. Material properties and geometrical properties
E0 = 2e11;%%% Elastic modulus

%%%%%%% Cross-sectional areas
A1 = 0.1;
A2 = A1;
A3 = 2*A1;

%%%%%%% Lengths of trusses
L1 = 1;
L2 = 1;
L3 = 1;

nE = 3; %%% nE = 3 : number of elements
A = zeros(nE,1);
A(1,1) = A1;A(2,1) = A2;A(3,1) = A3;

%%%%%%%% Stress-strain relation: S = C.E
C = E0*ones(nE,1);


%% %%%% Step 2. Create mesh
%%%%%%% 2.1. generate nodal coordinates
P1 = [L1, 0, 0];
P2 = [L1, L2, 0];
P3 = [0, L2, 0];
P4 = [0, 0, L3];


%%%%%% node coordinates: xyz
xyz = [P1; P2; P3; P4];
nP = size(xyz,1);%%%% number of nodes: nP = 4

%%%%%%% 2.2. Nodal conectivity for elements
eNodes = zeros(nE,2);
eNodes(1,:) = [1,2];%%%% element 1: 1-2
eNodes(2,:) = [2,3];%%%% element 2: 2-3
eNodes(3,:) = [2,4];%%%% element 3: 2-4


%% %%%%%% Step 3. Calculate stiffness matrix: K
nDof = 3*nP;%%%% each node has 3 dofs => total Dof number = 3*nP
K = zeros(nDof,nDof);%%% Gloabl stiffness K has size of nPxnP

for i = 1:nE
    %%%% nodes of element
    eNodei = eNodes(i,:);
    
    eDof = [eNodei;eNodei+nP;eNodei+2*nP];
    nnDof = length(eNodei);
    
    %%% Length of element
    Le = sqrt((xyz(eNodei(2),1)-xyz(eNodei(1),1))^2 + ...
        (xyz(eNodei(2),2)-xyz(eNodei(1),2))^2 + (xyz(eNodei(2),3)-xyz(eNodei(1),3))^2);
    
    cx = (xyz(eNodei(2),1)-xyz(eNodei(1),1))/Le;
    cy = (xyz(eNodei(2),2)-xyz(eNodei(1),2))/Le;
    cz = (xyz(eNodei(2),3)-xyz(eNodei(1),3))/Le;
    
    detJ = Le/2;
    invJ = 1/detJ;
    
    %%%%% Central gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% Strain-displcaement matrix: B
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %%%%% Element stiffness matrix in local coordinate system
    Ke = A(i,1)*(B'*C(i,1)*B)*detJ*gaussWt;
    
    %%%%% Transformation matrix
    Te = [cx, cy, cz,   0,  0, 0;
           0,   0, 0,  cx, cy, cz];
    %%%%% Element stiffness matrix in the global coordinate system
    KeG = Te'*Ke*Te;
    
    %%%%% Global stiffness matrix
    K(eDof,eDof) = K(eDof,eDof) + KeG;
end

%% %%% Boundary conditions
fixedP = [1;3;4];
fixedDof = [fixedP;fixedP+nP;fixedP+2*nP];

%% %%% Loading conditions
Fz = -2e9;
force = zeros(nDof,1);

force(2+2*nP,1) = Fz;%%%% force(2+2*nP,1) means Fz at node 2

%% %%%% Solution
disp = solution(nDof,fixedDof,K,force);

save('Truss_3D_Ex1.mat');

%% %%%% Post processing
u_dof = (1:nP)';
v_dof = u_dof + nP;
w_dof = u_dof + 2*nP;

scale = 1;
xyz_new = xyz + scale*[disp(u_dof,1),disp(v_dof,1),disp(w_dof,1)];

xyz_ansys = NLIST(:,1:3);
xyz_newANSYS = xyz_ansys + scale*dispANSYS;

%%%%%% Plot deformed configurations in ANSYS and FEM code
figure
for i = 1:nE
    eNodei = eNodes(i,:);
    h1 = plot3(xyz(eNodei,1),xyz(eNodei,2),xyz(eNodei,3),'k:','LineWidth',2);
    hold on
    h2 = plot3(xyz_new(eNodei,1),xyz_new(eNodei,2),xyz_new(eNodei,3),'b--','LineWidth',2);
    hold on
    eNodeiANSYS = ELIST(i,7:8);
    h3 = plot3(xyz_newANSYS(eNodeiANSYS,1),xyz_newANSYS(eNodeiANSYS,2),xyz_newANSYS(eNodeiANSYS,3),'r:','LineWidth',2);
    hold on
end
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid on
view(110,30)
legend('Undeformed', 'Deformed-FEM code','Deformed - ANSYS')
legend('boxoff')
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')
axis equal




%% %%%%%% Calculate strain and stress
strain_ele = zeros(nE,1);
stress_ele = zeros(nE,1);

strain_node = zeros(nE,2);
stress_node = zeros(nE,2);

for i = 1:nE
    %%%% nodes of element
    eNodei = eNodes(i,:);
    
    eDof = [eNodei;eNodei+nP;eNodei+2*nP];
    nnDof = length(eNodei);
    
    %%% Length of element
    Le = sqrt((xyz(eNodei(2),1)-xyz(eNodei(1),1))^2 + ...
        (xyz(eNodei(2),2)-xyz(eNodei(1),2))^2 + (xyz(eNodei(2),3)-xyz(eNodei(1),3))^2);
    cx = (xyz(eNodei(2),1)-xyz(eNodei(1),1))/Le;
    cy = (xyz(eNodei(2),2)-xyz(eNodei(1),2))/Le;
    cz = (xyz(eNodei(2),3)-xyz(eNodei(1),3))/Le;
    
    detJ = Le/2;
    invJ = 1/detJ;
    
    %%%%% Central gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% B matrix
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    %%%%% Transformation matrix
    Te = [cx, cy, cz,   0,  0, 0;
           0,   0, 0,  cx, cy, cz];
       
    %%%%% Transform glocal displacements back to local displacements
    eDof_v = [eDof(:,1);eDof(:,2)];
    %%%% see slide 9 for the order of displacement vector
    
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

