clear %%% to clear all data in the Workspace

%% %%% Material property and geometrical properties
E = 2e11;%%%% Elastic modulus
C = E;%%%% Stress-strain matrix: S = C*e

A = 1;%%%%% Cross-sectional area of the bar
Lx = 1;%%%% bar's length

%% %%%% Create mesh
dx = 0.1;%%%%% mesh size

x = (0:dx:Lx)';%%%%% node coordinates
nP = length(x);%%%% number of nodes

% % %%%%%% Check the mesh
figure
scatter(x,zeros(nP,1),20)
hold on
for i = 1:nP
    t = text(x(i,1),0.08,num2str(i));
    hold on
    t.Color = 'k';
    t.FontSize = 14;
end


%%%%%%% nodal conectivity for elements
    nE = nP-1;%%% nE = number of elements

    %%%%% eNodes will store start node and end node for each element
    eNodes = zeros(nE,2);

for i = 1:nE
   eNodes(i,1) = i;
   eNodes(i,2) = i+1;
end

%% %%%%%% Calculate stiffness matrix
K = zeros(nP,nP);%%% Stiffness K has size of nPxnP

for i = 1:nE
    %%%% nodes of element
    eDof = eNodes(i,:);%%%% eDof = [start_node, end_node]
    nnDof = length(eDof);%%%% nnDof = number of DOFs for each element
    
    %%%% length of element
    Le = abs(x(eDof(2),1)-x(eDof(1),1));%%%% Le = |x2-x1|
    
    detJ = Le/2;%%%%%% Determinant of Jacobian: det(J)
    invJ = 1/detJ;%%%%%% 1/det(J)
    
    %%%%% One point gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% Strain-displacement matrix: B
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    Ke = A*(B'*C*B)*detJ*gaussWt;


    K(eDof,eDof) = K(eDof,eDof) + Ke;
end

%% %%% Boundary conditions
fixedP = find(x == 0);
fixedDof = fixedP;
nDof = nP;

%% %%% Loading condition
force = zeros(nDof,1);
force(end) = 1e8;

%% %%%% Solution
disp = solution(nDof,fixedDof,K,force);
              

%% %%%% Post processing
%%%%%%% Compare displacement field with analytical solution
disp_Analytic = 1e8/E/A*x;

figure
plot(x,disp_Analytic,'r:','LineWidth',2);
hold on
plot(x,disp,'b--','LineWidth',2);
xlabel('x (m)');
ylabel('displacement: u (m)');
grid on
view(2)
legend('u-Analytical', 'u-FEA code')
legend('boxoff')
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%%%%%%% Compare displacement field with ANSYS
figure
plot(NLIST,uansys,'r:','LineWidth',2);
hold on
plot(x,disp,'b--','LineWidth',2);
xlabel('x (m)');
ylabel('displacement: u (m)');
grid on
view(2)
legend('u-ANSYS', 'u-FEA code')
legend('boxoff')
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%% %%%%% Strain and stress
strain_ele = zeros(nE,1); 
stress_ele = zeros(nE,1);

strain_node = zeros(nE,2);
stress_node = zeros(nE,2);
for i = 1:nE
    %%%% nodes of element
    eDof = eNodes(i,:);%%%% eDof = [start_node, end_node]
    nnDof = length(eDof);%%%% nnDof = number of DOFs for each element
    
    %%%% length of element
    Le = abs(x(eDof(2),1)-x(eDof(1),1));%%%% Le = |x2-x1|
    
    detJ = Le/2;%%%%%% Determinant of Jacobian: det(J)
    invJ = 1/detJ;%%%%%% 1/det(J)
    
    %%%%% One point gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% Strain-displacement matrix: B
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    %%%%% strain for element
    ei = B*disp(eDof);
    strain_ele(i,:) = ei;
    stress_ele(i,:) = C*ei;
    
    strain_node(i,:) = ei;
    stress_node(i,:) = C*ei;
end

%%%%%%%%% FINISH!! %%%%%%%%%%%%%%
%%%%%%%%% Please save data if you want to use it later %%%%%%%%%%
