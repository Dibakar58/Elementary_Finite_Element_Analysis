clear

%% %%% Step 1. Material properties and geometrical properties
E1 = 2e11;E2 = 1.5e11; E3 = 1e11;%%%% Elastic modulus

A1 = 1; A2 = 1/2; A3 = 1/4;%%%% Cross-sectional areas (m^2)

L1 = 1; L2 = 1; L3 = 1;%%%% elngth of three sections (m)

%%%%%%%% Stress-strain relation: S = C.E
C1 = E1; C2 = E2; C3 = E3;

%% %%%% Step 2. Create mesh
%%%%%%% 2.1. generate nodal coordinates
dx = 0.2;%%%%% mesh size

x1 = (0:dx:L1)';%%%%% node coordinates for the first section

x2 = ((L1+dx):dx:(L1+L2))';%%%%% node coordinates for the second section

x3 = ((L1+L2+dx):dx:(L1+L2+L3))';%%%%% node coordinates for the third section

x = [x1; x2; x3];%%% combine all sections

nP = length(x);%%%% number of nodes

%%%%%%%% Try to plot the nodes
figure
scatter(x,zeros(nP,1),'filled');
grid on
for i = 1:nP
    t = text(x(i,1)-0.03,0.08,num2str(i));
    hold on
    t.Color = 'k';
    t.FontSize = 14;
end

%%%%%%% 2.2. nodal conectivity for elements
nE = nP-1;%%% number of elements
eNodes = zeros(nE,2);
for i = 1:nE
   eNodes(i,1) = i;
   eNodes(i,2) = i+1;
end

%% %%%%%% Step 3. Calculate stiffness matrix: K
nDof = nP;
K = zeros(nDof,nDof);%%% Global stiffness K has size of nPxnP

for i = 1:nE
    %%%% nodes of element
    eDof = eNodes(i,:);
    nnDof = length(eDof);
    
    %%%% Length of element
    Le = abs(x(eDof(2),1)-x(eDof(1),1));%%%%% Le = |x2-x1|
    
    detJ = Le/2;
    invJ = 1/detJ;
    
    %%%%% Central gauss quadrature: (xi = 0, weight = 2)
    [shape,nDeriv] = shapeFunct_Truss(0);
    gaussWt = 2;
    
    %%%%% B matrix
    Xderiv = nDeriv*invJ;
    B = zeros(1,nnDof); B(1:nnDof) = Xderiv(:);
    
    
    if x(eDof(2),1)<= L1
        Ke = A1*(B'*C1*B)*detJ*gaussWt;
    elseif x(eDof(2),1)<= (L2+L1)
        Ke = A2*(B'*C2*B)*detJ*gaussWt;
    else
        Ke = A3*(B'*C3*B)*detJ*gaussWt;
    end

    K(eDof,eDof) = K(eDof,eDof) + Ke;
end

%% %%% Boundary conditions
fixedP = find(x == 0);
fixedDof = fixedP;

%% %%% Loading conditions
F1 = -1e6; F2 = -2e5; F3 = 3e5;
force = zeros(nP,1);

p1 = find(x(:,1) == L1);
p2 = find(x(:,1) == (L1+L2));
force(p1,1) = F1;
force(p2,1) = F2;
force(end) = F3;

%% %%%% Solution
disp = solution(nDof,fixedDof,K,force);



%% %%%% Post processing
%%%% 1: Plot displacement field and compare with ANSYS
[xANSYS,id] = sort(NLIST);
uANSYS1 = uANSYS(id);

figure
plot(x,disp,'r:','LineWidth',2);
hold on
plot(xANSYS,uANSYS1,'b--','LineWidth',2);
xlabel('x (m)');
ylabel('displacement: u (m)');
grid on
view(2)
legend('u-FEA code', 'u-ANSYS')
legend('boxoff')
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')



%%%% 2: Strain and stress field
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
    
    %%%%% strain for element and nodes
    ei = B*disp(eDof);
    
    strain_ele(i,:) = ei;
    strain_node(i,:) = ei;
    
    %%%%% stress for elements and nodes
    if x(eDof(2),1)<= L1
        Ci = C1;
    elseif x(eDof(2),1)<= (L2+L1)
        Ci = C2;
    else
        Ci = C3;
    end
    
    stress_ele(i,:) = Ci*ei;
    stress_node(i,:) = Ci*ei;
end

%%%%%%%%% FINISH!! %%%%%%%%%%%%%%
%%%%%%%%% Please save data if you want to use it later %%%%%%%%%%
