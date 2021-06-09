% clear
%% %%%%% materials
E = 2e11; pois = 0.27;dens = 7850;
h0 = 0.01;%%% Plate thickness
Lx = 1;
Ly = 1;

% % %%% Stress- strain matriz: C - Plane strain
% % C = E*(1-pois)/(1-2*pois)/(1+pois)*[1 pois/(1-pois) 0;
% %                                     pois/(1-pois) 1 0;
% %                                     0 0 (1-2*pois)/2/(1-pois)];
% % %%%% lame's constants
% % lamda = pois*E/(1-2*pois)/(1+pois);
% % muy = E/2/(1+pois);

%%%%% Stress- strain matriz: C - Plane stress
C = E/(1-pois^2)*[1 pois 0;
                  pois 1 0;
                  0 0 (1-pois)/2];
% % lamda = pois*E/(1-pois^2);
% % muy = E/2/(1+pois);%%%% Shear modulus

%% %%% Import mesh from ANSYS
%%%%% Import manually files from ansys: 

nan_ELIST = isnan(ELIST(:,1));
id = find(nan_ELIST(:,1)== 0);
ELIST = ELIST(id,:);

nan_NLIST = isnan(NLIST(:,1));
id = find(nan_NLIST(:,1)== 0);
NLIST = NLIST(id,:);


xy = NLIST(:,2:3);
nP = size(xy,1);
elNode = ELIST(:,7:10);
nE = size(elNode,1);

dx = Lx/20;%%%% mesh size
dy = dx;


% % %%%%% Try to plot mesh
% % PlotMesh(xy,elNode);

%% %%%%% GDof: global number of degrees of freedom
nDof = 2*nP;%%% degree of freedom

%% %%%% boundary conditions - the left edge is clamped
left = find(xy(:,1)==0); %% left edge
% % bottom = find(nodeCoord(:,2)==0); %%% bottom
fixedDof = [left; left + nP];%%% u = v = 0 at left edge

%% %%%%% Load
fx = 5e7;
fy = 0;%1.2e7;

force = zeros(nDof,1);
right = find(xy(:,1)== Lx);
% % top = find(nodeCoord(:,2)== Ly);

%%%% Apply Fx
force(right) = fx*dy;
force(right(1)) = fx*dy/2;
force(right(end)) = fx*dy/2;

% % %%%% Apply Fy on top
% % force(top+nnode) = fy*dx;
% % force(top(1)+nnode) = fy*dx/2;
% % force(top(end)+nnode) = fy*dx/2;

%% %%%%% Stiffness matrix
[stiff,mass] = formStiffness2D(nDof,nE,elNode,nP,xy,C,dens,h0);


%% %%%%% Solution of equilibrium
disp = solution(nDof,fixedDof,stiff,force);

save('Plate2D_ex2');

%% %%% Postprocessing
udof = (1:nP)';
vdof = udof+nP;

ux = disp(udof,1);
uy = disp(vdof,1);

scale = 1;
xy_NewFEM = xy+scale*[ux uy];
min_u_FEM = min(ux);
max_u_FEM = max(ux);
figure;
PlotFieldonMesh(xy_NewFEM,elNode,ux,min_u_FEM,max_u_FEM);
% title('u (m)-FEM');%,'Fontsize',14)
grid on


min_v_FEM = min(uy);
max_v_FEM = max(uy);
figure;
PlotFieldonMesh(xy_NewFEM,elNode,uy,min_v_FEM,max_v_FEM);
title('v (m)-FEM');%,'Fontsize',14)



%% %%%%%%%% Compare with ANSYS results
%%%%%%%%%%% Step 1: process ANSYS data
dispANSYS2 = [];
for i = 1:size(dispANSYS,1)
    if isnan(dispANSYS(i,1)) == 0
        dispANSYS2 = [dispANSYS2;dispANSYS(i,:)];
    end
end


%%%%%%%%%%% Step 2: Compare displacement fields: u
udof = (1:nP)';vdof = udof+nP;
ux = disp(udof,1);
uy = disp(vdof,1);

min_u_FEM = min(ux);
max_u_FEM = max(ux);

min_u_ANSYS = min(dispANSYS2(:,2));
max_u_ANSYS = max(dispANSYS2(:,2));

min_u = min(min_u_FEM,min_u_ANSYS);
max_u = max(max_u_FEM,max_u_ANSYS);

scale = 1;
xy_NewFEM = xy+scale*[ux uy];
xy_NewANSYS = NLIST(:,2:3) + scale*dispANSYS2(:,2:3);
% 
%%%%%%%% Plot u in FEM code
figure
PlotFieldonMesh(xy_NewFEM,elNode,ux,min_u,max_u);
title('u (m)-FEM');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')

%%%%%%%% Plot u in ANSYS
figure
PlotFieldonMesh(xy_NewANSYS,ELIST(:,7:10),dispANSYS2(:,2),min_u,max_u);
title('u (m)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%%%%%%%%%%% Step 3: Compare displacement fields: v
min_v_FEM = min(uy);
max_v_FEM = max(uy);

min_v_ANSYS = min(dispANSYS2(:,3));
max_v_ANSYS = max(dispANSYS2(:,3));

min_v = min(min_v_FEM,min_v_ANSYS);
max_v = max(max_v_FEM,max_v_ANSYS);


%%%%%%%% Plot v in FEM code
figure
PlotFieldonMesh(xy_NewFEM,elNode,uy,min_v,max_v);
title('v (m)-FEM');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%%%%%%%% Plot v in ANSYS
figure
PlotFieldonMesh(xy_NewANSYS,ELIST(:,7:10),dispANSYS2(:,3),min_v,max_v);
title('v (m)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%% %%%%%% Post-processing Stress and strain
%%%%%%%% Stress components in ANSYS
stressANSYS2 = [];
for i = 1:size(stressANSYS,1)
    if isnan(stressANSYS(i,1)) == 0
        stressANSYS2 = [stressANSYS2;stressANSYS(i,:)];
    end
end

nE = size(stressANSYS2,1)/4;
Sxx_ANSYS = zeros(nE,4);
Syy_ANSYS = zeros(nE,4);
Sxy_ANSYS = zeros(nE,4);

for i = 1:nE
   i1 = (i-1)*4 + 1;
   i2 = i1 + 3;
   Sxx_ANSYS(i,:) = stressANSYS2(i1:i2,2)';
   Syy_ANSYS(i,:) = stressANSYS2(i1:i2,3)';
   Sxy_ANSYS(i,:) = stressANSYS2(i1:i2,5)';
end

% % strainANSYS2 = [];
% % for i = 1:size(strainANSYS,1)
% %     if isnan(strainANSYS(i,1)) ==0
% %         strainANSYS2 = [strainANSYS2;strainANSYS(i,:)];
% %     end
% % end
% % 
% % exx_ANSYS = zeros(nE,4);
% % eyy_ANSYS = zeros(nE,4);
% % exy_ANSYS = zeros(nE,4);
% % 
% % for i = 1:nE
% %    i1 = (i-1)*4 + 1;
% %    i2 = i1 + 3;
% %    exx_ANSYS(i,:) = strainANSYS2(i1:i2,2)';
% %    eyy_ANSYS(i,:) = strainANSYS2(i1:i2,3)';
% %    exy_ANSYS(i,:) = strainANSYS2(i1:i2,5)';
% % end



% % %% %%%%% Method 1: Calculate strain, stress directly at each node of element
% % Sxx_FEM = zeros(nE,4);
% % Syy_FEM = zeros(nE,4);
% % Sxy_FEM = zeros(nE,4);
% % 
% % exx_FEM = zeros(nE,4);
% % eyy_FEM = zeros(nE,4);
% % exy_FEM = zeros(nE,4);
% % 
% % natural_Coord = [-1,-1; 1,-1; 1,1; -1,1];
% % %%%%%% Natural coordinates (xi, eta) for our 4-node element
% % %%%%%%           (-1,1)-------------(1,1)
% % %%%%%%              |                 |
% % %%%%%%              |                 |
% % %%%%%%              |                 |
% % %%%%%%              |                 |
% % %%%%%%           (-1,-1)-------------(1,-1)
% % 
% % for ie = 1:nE                           
% %     id = elNode(ie,:);  
% %     
% %     elDof = [id id+nP];
% %     ndof = length(id);%%% ndof = 4 for Q4
% % 
% %     %%% loop over Natural coordinates
% %     for j = 1:length(id)
% %         xi = natural_Coord(j,1);
% %         eta = natural_Coord(j,2);
% %         %%% shape functions and derivatives
% %         [shape,nDeriva] = shapeFuncQ4(xi,eta);
% %         %%% Jacobian matrix, inverse of Jacobian    
% %         [J,xyDeriva] = Jacobian(xy(id,:),nDeriva);
% %         %%%  B matrix (Linear strain - displacement matrix)
% %         B = zeros(3,2*ndof);%%% B(3x8)
% %         B(1,1:ndof)             = xyDeriva(:,1)';
% %         B(2,ndof+1:(2*ndof))    = xyDeriva(:,2)';
% %         B(3,1:ndof)             = xyDeriva(:,2)';
% %         B(3,ndof+1:(2*ndof))    = xyDeriva(:,1)';
% %         
% %         %%% strain: ei = [e_xx; e_yy; 2*e_xy];
% %         Uei = disp(elDof',1);
% %         ei = B*Uei;
% %         Si = C*ei;
% %         
% %         %%%% stress and strain matrix
% %         exx_FEM(ie,j) = ei(1,1);
% %         eyy_FEM(ie,j) = ei(2,1);
% %         exy_FEM(ie,j) = ei(3,1);
% % 
% %         Sxx_FEM(ie,j) = Si(1,1);
% %         Syy_FEM(ie,j) = Si(2,1);
% %         Sxy_FEM(ie,j) = Si(3,1);
% %     end
% % end

% %% %%%%% Method 2: Calculate strain, stress at Gauss Points and take average
% %%%%%%% values for nodes
% Sxx_FEM = zeros(nE,4);
% Syy_FEM = zeros(nE,4);
% Sxy_FEM = zeros(nE,4);
% 
% exx_FEM = zeros(nE,4);
% eyy_FEM = zeros(nE,4);
% exy_FEM = zeros(nE,4);
% 
% gauss_Points = 1/sqrt(3)*[-1,-1; 1,-1; 1,1; -1,1];
% % %%%%%% Natural coordinates (xi, eta) for our 4-node element
% % %%%%%%           (-1,1)-------------(1,1)
% % %%%%%%              |    *4      *3   |
% % %%%%%%              |                 |
% % %%%%%%              |                 |
% % %%%%%%              |    *1      *2   |
% % %%%%%%           (-1,-1)-------------(1,-1)
% 
% for ie = 1:nE                           
%     id = elNode(ie,:);  
%     
%     elDof = [id id+nP];
%     ndof = length(id);%%% ndof = 4 for Q4
% 
%     %%% loop over Natural coordinates
%     for j = 1:length(gauss_Points)
%         xi = gauss_Points(j,1);
%         eta = gauss_Points(j,2);
%         %%% shape functions and derivatives
%         [shape,nDeriva] = shapeFuncQ4(xi,eta);
%         %%% Jacobian matrix, inverse of Jacobian    
%         [J,xyDeriva] = Jacobian(xy(id,:),nDeriva);
%         %%%  BL matrix (Linear strain - displacement matrix)
%         B = zeros(3,2*ndof);%%% B(3x8)
%         B(1,1:ndof)             = xyDeriva(:,1)';
%         B(2,ndof+1:(2*ndof))    = xyDeriva(:,2)';
%         B(3,1:ndof)             = xyDeriva(:,2)';
%         B(3,ndof+1:(2*ndof))    = xyDeriva(:,1)';
%         
%         %%% Strain, stress at gauss Points: ei = [e_xx; e_yy; 2*e_xy];
%         Uei = disp(elDof',1);
%         ei = B*Uei;
%         Si = C*ei;
%         
%         %%%% stress and strain matrix at gauss points
%         exx_FEM(ie,j) = ei(1,1);
%         eyy_FEM(ie,j) = ei(2,1);
%         exy_FEM(ie,j) = ei(3,1);
% 
%         Sxx_FEM(ie,j) = Si(1,1);
%         Syy_FEM(ie,j) = Si(2,1);
%         Sxy_FEM(ie,j) = Si(3,1);
%     end
% end


%% %%%%% Method 3: Calculate Normal strain, stress (exx, eyy, Sxx, Syy) at 
%%%%%%%% Nodes, meanwhile, shear strain, shear stress (exy, Sxy) at reduced
%%%%%%%% gauss point (xi=0, eta=0)
%%%%%%% values for nodes
Sxx_FEM = zeros(nE,4);
Syy_FEM = zeros(nE,4);
Sxy_FEM = zeros(nE,4);

exx_FEM = zeros(nE,4);
eyy_FEM = zeros(nE,4);
exy_FEM = zeros(nE,4);

natural_Coord = [-1,-1; 1,-1; 1,1; -1,1];
gauss_Points = 0*[-1,-1; 1,-1; 1,1; -1,1];

% %%%%%% Natural coordinates (xi, eta) for our 4-node element
% %%%%%%           (-1,1)-------------(1,1)
% %%%%%%              |                 |
% %%%%%%              |                 |
% %%%%%%              |        *        |
% %%%%%%              |                 |
% %%%%%%              |                 |
% %%%%%%           (-1,-1)-------------(1,-1)

for ie = 1:nE                           
    id = elNode(ie,:);  
    
    elDof = [id id+nP];
    ndof = length(id);%%% ndof = 4 for Q4

    %%% loop over Natural coordinates
    for j = 1:length(gauss_Points)
        xi_norm = natural_Coord(j,1);
        eta_norm = natural_Coord(j,2);
        
        xi_shear = gauss_Points(j,1);
        eta_shear = gauss_Points(j,2);
        
        %%% shape functions and derivatives
        [shape_norm,nDeriva_norm] = shapeFuncQ4(xi_norm,eta_norm);
        
        [shape_shear,nDeriva_shear] = shapeFuncQ4(xi_shear,eta_shear);
        
        %%% Jacobian matrix, inverse of Jacobian    
        [J_norm,xyDeriva_norm] = Jacobian(xy(id,:),nDeriva_norm);
        
        [J_shear,xyDeriva_shear] = Jacobian(xy(id,:),nDeriva_shear);
        %%%  B matrix (Linear strain - displacement matrix)
        B_norm = zeros(3,2*ndof);%%% B(3x8)
        B_norm(1,1:ndof)             = xyDeriva_norm(:,1)';
        B_norm(2,ndof+1:(2*ndof))    = xyDeriva_norm(:,2)';
        B_norm(3,1:ndof)             = xyDeriva_norm(:,2)';
        B_norm(3,ndof+1:(2*ndof))    = xyDeriva_norm(:,1)';
        
        B_shear = zeros(3,2*ndof);%%% B(3x8)
        B_shear(1,1:ndof)             = xyDeriva_shear(:,1)';
        B_shear(2,ndof+1:(2*ndof))    = xyDeriva_shear(:,2)';
        B_shear(3,1:ndof)             = xyDeriva_shear(:,2)';
        B_shear(3,ndof+1:(2*ndof))    = xyDeriva_shear(:,1)';
        
        %%% Strain, stress at gauss Points: ei = [e_xx; e_yy; 2*e_xy];
        Uei = disp(elDof',1);
        %%%%%% we will take exx, eyy, Sxx, Syy from these vectors
        ei_norm = B_norm*Uei;
        Si_norm = C*ei_norm;
        
        %%%%%% we will take exy, Sxy from these vectors
        ei_shear = B_shear*Uei;
        Si_shear = C*ei_shear;
        
        %%%% stress and strain matrix at gauss points
        exx_FEM(ie,j) = ei_norm(1,1);
        eyy_FEM(ie,j) = ei_norm(2,1);
        exy_FEM(ie,j) = ei_shear(3,1);

        Sxx_FEM(ie,j) = Si_norm(1,1);
        Syy_FEM(ie,j) = Si_norm(2,1);
        Sxy_FEM(ie,j) = Si_shear(3,1);
    end
end

%%%%%% As you see, each node has multiple stress, strain values. It
%%%%%% depends on elements. ==> we will average these values for each node
Sxx_node_FEM = zeros(nP,1);
Syy_node_FEM = zeros(nP,1);
Sxy_node_FEM = zeros(nP,1);

exx_node_FEM = zeros(nP,1);
eyy_node_FEM = zeros(nP,1);
exy_node_FEM = zeros(nP,1);

Sxx_node_ANSYS = zeros(nP,1);
Syy_node_ANSYS = zeros(nP,1);
Sxy_node_ANSYS = zeros(nP,1);

exx_node_ANSYS = zeros(nP,1);
eyy_node_ANSYS = zeros(nP,1);
exy_node_ANSYS = zeros(nP,1);


for i = 1:nP
    
    ni = 0;
    
    Sxx_nodei = 0;
    Syy_nodei = 0;
    Sxy_nodei = 0;

    exx_nodei = 0;
    eyy_nodei = 0;
    exy_nodei = 0;
    
    Sxx_nodeANSi = 0;
    Syy_nodeANSi = 0;
    Sxy_nodeANSi = 0;

%     exx_nodeANSi = 0;
%     eyy_nodeANSi = 0;
%     exy_nodeANSi = 0;
    
    for j = 1:4
        idi = find(elNode(:,j)== i);
        ni = ni + length(idi);
        if ~isempty(idi)
            Sxx_nodei = Sxx_nodei + sum(Sxx_FEM(idi,j));
            Syy_nodei = Syy_nodei + sum(Syy_FEM(idi,j));
            Sxy_nodei = Sxy_nodei + sum(Sxy_FEM(idi,j));

            exx_nodei = exx_nodei + sum(exx_FEM(idi,j));
            eyy_nodei = eyy_nodei + sum(eyy_FEM(idi,j));
            exy_nodei = exy_nodei + sum(exy_FEM(idi,j));
            
            
            Sxx_nodeANSi = Sxx_nodeANSi + sum(Sxx_ANSYS(idi,j));
            Syy_nodeANSi = Syy_nodeANSi + sum(Syy_ANSYS(idi,j));
            Sxy_nodeANSi = Sxy_nodeANSi + sum(Sxy_ANSYS(idi,j));

%             exx_nodeANSi = exx_nodeANSi + sum(exx_ANSYS(idi,j));
%             eyy_nodeANSi = eyy_nodeANSi + sum(eyy_ANSYS(idi,j));
%             exy_nodeANSi = exy_nodeANSi + sum(exy_ANSYS(idi,j));
            
        end
    end
    Sxx_node_FEM(i,1) = 1/ni*Sxx_nodei;
    Syy_node_FEM(i,1) = 1/ni*Syy_nodei;
    Sxy_node_FEM(i,1) = 1/ni*Sxy_nodei;

    exx_node_FEM(i,1) = 1/ni*exx_nodei;
    eyy_node_FEM(i,1) = 1/ni*eyy_nodei;
    exy_node_FEM(i,1) = 1/ni*exy_nodei;
    
    
    Sxx_node_ANSYS(i,1) = 1/ni*Sxx_nodeANSi;
    Syy_node_ANSYS(i,1) = 1/ni*Syy_nodeANSi;
    Sxy_node_ANSYS(i,1) = 1/ni*Sxy_nodeANSi;

%     exx_node_ANSYS(i,1) = 1/ni*exx_nodeANSi;
%     eyy_node_ANSYS(i,1) = 1/ni*eyy_nodeANSi;
%     exy_node_ANSYS(i,1) = 1/ni*exy_nodeANSi;
end

%% %%%%%% plot Sxx
min_Sx = min(min(Sxx_node_ANSYS),min(Sxx_node_FEM));
max_Sx = max(max(Sxx_node_ANSYS),max(Sxx_node_FEM));

%%%%%%%% Plot S_xx in ANSYS
figure
PlotFieldonMesh(xy_NewANSYS,ELIST(:,7:10),Sxx_node_ANSYS,min_Sx,max_Sx);
title('S_{xx} (N/m^2)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%%%%%%%% Plot S_xx in FEM Code

figure
PlotFieldonMesh(xy_NewFEM,elNode,Sxx_node_FEM,min_Sx,max_Sx);
title('S_{xx} (N/m^2)-FEM');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')



%% %%%%%% plot Syy
min_Sy = min(min(Syy_node_ANSYS),min(Syy_node_FEM));
max_Sy = max(max(Syy_node_ANSYS),max(Syy_node_FEM));

%%%%%%%% Plot S_yy in ANSYS
figure
PlotFieldonMesh(xy_NewANSYS,ELIST(:,7:10),Syy_node_ANSYS,min_Sy,max_Sy);
title('S_{yy} (N/m^2)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')

%%%%%%%% Plot S_yy in FEM Code
min_Syy = min(Syy_node_FEM);
max_Syy = max(Syy_node_FEM);

figure
PlotFieldonMesh(xy_NewFEM,elNode,Syy_node_FEM,min_Sy,max_Sy);
title('S_{yy} (N/m^2)-FEM');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%% %%%%%% plot Sxy
min_Sxy = min(min(Sxy_node_ANSYS),min(Sxy_node_FEM));
max_Sxy = max(max(Sxy_node_ANSYS),max(Sxy_node_FEM));

%%%%%%%% Plot S_xy in ANSYS
figure
PlotFieldonMesh(xy_NewANSYS,ELIST(:,7:10),Sxy_node_ANSYS,min_Sxy,max_Sxy);
title('S_{xy} (N/m^2)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')

%%%%%%%% Plot S_xy in FEM Code
figure
PlotFieldonMesh(xy_NewFEM,elNode,Sxy_node_FEM,min_Sxy,max_Sxy);
title('S_{xy} (N/m^2)-FEM');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')

%%%%%%%% Plot von Mises stress in ANSYS
S_vm_ANSYS = sqrt(Sxx_node_ANSYS.^2 + Syy_node_ANSYS.^2 - Sxx_node_ANSYS.*Syy_node_ANSYS + 3*Sxy_node_ANSYS.^2);
S_vm_FEM = sqrt(Sxx_node_FEM.^2 + Syy_node_FEM.^2 - Sxx_node_FEM.*Syy_node_FEM + 3*Sxy_node_FEM.^2);

min_Svm = min(min(S_vm_ANSYS),min(S_vm_FEM));
max_Svm = max(max(S_vm_ANSYS),max(S_vm_FEM));

figure
PlotFieldonMesh(xy_NewANSYS,ELIST(:,7:10),S_vm_ANSYS,min_Svm,max_Svm);
title('S_{VM} (N/m^2)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')



figure
PlotFieldonMesh(xy_NewFEM,elNode,S_vm_FEM,min_Svm,max_Svm);
title('S_{VM} (N/m^2)-FEM');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%% %%%%%%% FINISH!! %%%%%%%%%%%%%%
%%%%%%%%% Please save data if you want to use it later %%%%%%%%%%

