% % clear
tol = 1e-6;
%% %%%%% Material properties
E = 2e11;
pois = 0.3;
dens = 7850;

%% %%%%% Geometrical properties
Lx = 1;
Ly = 0.1;
Lz = 0.1;

%% %%%%%% Stress- strain matrix: C - 3D
E1 = E/(1-2*pois)/(1+pois);
C = E1*[1-pois, pois, pois, 0, 0, 0;
        pois, 1-pois, pois, 0, 0, 0;
        pois, pois, 1-pois, 0, 0, 0;
        0, 0, 0, (1-2*pois)/2, 0, 0;
        0, 0, 0, 0, (1-2*pois)/2, 0;
        0, 0, 0, 0, 0, (1-2*pois)/2];
    
% %% %%%%% Mesh generation: Import mesh from ANSYS

nan_ELIST = isnan(ELIST(:,1));
id = find(nan_ELIST(:,1)== 0);
ELIST = ELIST(id,:);

nan_NLIST = isnan(NLIST(:,1));
id = find(nan_NLIST(:,1)== 0);
NLIST = NLIST(id,:);


xyz = NLIST(:,2:4);
nP = size(xyz,1);%%%% number of nodes

elNode = ELIST(:,7:14);
nE = size(elNode,1);%%%% number of elements

dx = Ly/5;%%%% defined based on information from ANSYS mesh
dy = dx;
dz = dx;

% % %%%%%%% Quick view on the model
% % figure
% % scatter3(xyz(:,1),xyz(:,2),xyz(:,3),50);
% % axis equal


%% %%%% GDof: global number of degrees of freedom
nDof = 3*nP;%% number of degrees of freedom

%% %%%%% Stiffness matrix
[stiff,mass] = formStiffness3D(nDof,nE,elNode,nP,xyz,C,dens);


%% %%%% Boundary conditions - the left edge is clamped
left = find(xyz(:,1)==0); %% left edge
fixedDof = [left; left + nP; left + 2*nP];

%% %%% Loading conditions
fz = -5e6;%%% N/m

force = zeros(nDof,1);
right_top1 = find(xyz(:,1)== Lx & xyz(:,3)== Lz & xyz(:,2)<=tol);
right_top2 = find(xyz(:,1)== Lx & xyz(:,3)== Lz & xyz(:,2)>tol & xyz(:,2)<Ly-tol);
right_top3 = find(xyz(:,1)== Lx & xyz(:,3)== Lz & xyz(:,2)>=Ly-tol);

% % figure
% % scatter3(xyz(:,1),xyz(:,2),xyz(:,3),50);
% % hold on
% % scatter3(xyz(left,1),xyz(left,2),xyz(left,3),'r');
% % hold on
% % scatter3(xyz(right_top1,1),xyz(right_top1,2),xyz(right_top1,3),'k');
% % hold on
% % scatter3(xyz(right_top3,1),xyz(right_top3,2),xyz(right_top3,3),'k');
% % hold on
% % scatter3(xyz(right_top2,1),xyz(right_top2,2),xyz(right_top2,3),'b');
% % axis equal

%%%% Apply Fz
force(right_top2 + 2*nP) = fz*dy;
force(right_top1 + 2*nP) = fz*dy/2;
force(right_top3 + 2*nP) = fz*dy/2;

%% %%%% Solution of equilibrium
disp = solution(nDof,fixedDof,stiff,force);

save('Struct3Dex1');



%% %%% Postprocessing

dispANSYS2 = [];
for i = 1:size(dispANSYS,1)
    if isnan(dispANSYS(i,1)) ==0
        dispANSYS2 = [dispANSYS2;dispANSYS(i,:)];
    end
end

u_ANSYS = dispANSYS2(:,2);
v_ANSYS = dispANSYS2(:,3);
w_ANSYS = dispANSYS2(:,4);


% 
%%%%%%%% Plot displacements
udof = (1:nP)';
vdof = udof+nP;
wdof = udof+2*nP;

u = disp(udof,1);
v = disp(vdof,1);
w = disp(wdof,1);

u_min = min(min(u_ANSYS),min(u));
u_max = max(max(u_ANSYS),max(u));

v_min = min(min(v_ANSYS),min(v));
v_max = max(max(v_ANSYS),max(v));

w_min = min(min(w_ANSYS),min(w));
w_max = max(max(w_ANSYS),max(w));



scale = 1;
%%%% deformed configuration in our FEM code
xyzNew = xyz+scale*[u v w];

%%%%%% plot displacement u-FEM
figure
scatter3(xyzNew(:,1),xyzNew(:,2),xyzNew(:,3),40,u,'filled');
title('u (m)-FEM');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([u_min u_max]);
zstep = (u_max - u_min)/5;
set(h, 'ytick', (u_min:zstep:u_max))

%%%%%% plot displacement u-ANSYS
%%%% deformed configuration in ANSYS
xyzNew_ANS = xyz+scale*[u v w];

figure
scatter3(xyzNew_ANS(:,1),xyzNew_ANS(:,2),xyzNew_ANS(:,3),40,u_ANSYS,'filled');
title('u (m)-ANSYS');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([u_min u_max]);
zstep = (u_max - u_min)/5;
set(h, 'ytick', (u_min:zstep:u_max))



%%%%%% plot displacement v-FEM
figure
scatter3(xyzNew(:,1),xyzNew(:,2),xyzNew(:,3),40,v,'filled');
title('v (m)-FEM');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([v_min v_max]);
zstep = (v_max - v_min)/5;
set(h, 'ytick', (v_min:zstep:v_max))

%%%%%% plot displacement v-ANSYS
figure
scatter3(xyzNew_ANS(:,1),xyzNew_ANS(:,2),xyzNew_ANS(:,3),40,v_ANSYS,'filled');
title('v (m)-ANSYS');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([v_min v_max]);
zstep = (v_max - v_min)/5;
set(h, 'ytick', (v_min:zstep:v_max))



%%%%%% plot displacement w-FEM
figure
scatter3(xyzNew(:,1),xyzNew(:,2),xyzNew(:,3),40,w,'filled');
title('w (m)-FEM');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([w_min w_max]);
zstep = (w_max - w_min)/5;
set(h, 'ytick', (w_min:zstep:w_max))
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')


%%%%%% plot displacement w-ANSYS
figure
scatter3(xyzNew_ANS(:,1),xyzNew_ANS(:,2),xyzNew_ANS(:,3),40,w_ANSYS,'filled');
% title('w (m)-ANSYS');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([w_min w_max]);
zstep = (w_max - w_min)/5;
set(h, 'ytick', (w_min:zstep:w_max))
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')



%% %%%%%%%%% Post-processing: Calculate Stress, strain
%%%%%%% Stress ANSYS
stressANSYS2 = [];
for i = 1:size(stressANSYS,1)
    if isnan(stressANSYS(i,1)) ==0
        stressANSYS2 = [stressANSYS2;stressANSYS(i,:)];
    end
end

Sxx_ANSYS = zeros(nE,8);
Syy_ANSYS = zeros(nE,8);
Szz_ANSYS = zeros(nE,8);
Sxy_ANSYS = zeros(nE,8);
Syz_ANSYS = zeros(nE,8);
Sxz_ANSYS = zeros(nE,8);

for i = 1:nE
   istart = (i-1)*8 + 1;
   iend = istart + 7;
   
   Sxx_ANSYS(i,:) = stressANSYS2(istart:iend,2);
   Syy_ANSYS(i,:) = stressANSYS2(istart:iend,3);
   Szz_ANSYS(i,:) = stressANSYS2(istart:iend,4);
   Sxy_ANSYS(i,:) = stressANSYS2(istart:iend,5);
   Syz_ANSYS(i,:) = stressANSYS2(istart:iend,6);
   Sxz_ANSYS(i,:) = stressANSYS2(istart:iend,7);
end

Sxx_node_ANSYS = zeros(nP,1);
Syy_node_ANSYS = zeros(nP,1);
Szz_node_ANSYS = zeros(nP,1);
Sxy_node_ANSYS = zeros(nP,1);
Syz_node_ANSYS = zeros(nP,1);
Sxz_node_ANSYS = zeros(nP,1);


for i = 1:nP
    
    ni = 0;
    
    Sxx_nodei = 0;
    Syy_nodei = 0;
    Szz_nodei = 0;
    Sxy_nodei = 0;
    Syz_nodei = 0;
    Sxz_nodei = 0;
    
    exx_nodei = 0;
    eyy_nodei = 0;
    ezz_nodei = 0;
    exy_nodei = 0;
    eyz_nodei = 0;
    exz_nodei = 0;
    
    for j = 1:8
        idi = find(elNode(:,j)==i);
        ni = ni + length(idi);
        if ~isempty(idi)
            Sxx_nodei = Sxx_nodei + sum(Sxx_ANSYS(idi,j));
            Syy_nodei = Syy_nodei + sum(Syy_ANSYS(idi,j));
            Szz_nodei = Szz_nodei + sum(Szz_ANSYS(idi,j));
            Sxy_nodei = Sxy_nodei + sum(Sxy_ANSYS(idi,j));
            Syz_nodei = Syz_nodei + sum(Syz_ANSYS(idi,j));
            Sxz_nodei = Sxz_nodei + sum(Sxz_ANSYS(idi,j));

        end
    end
    Sxx_node_ANSYS(i,1) = 1/ni*Sxx_nodei;
    Syy_node_ANSYS(i,1) = 1/ni*Syy_nodei;
    Szz_node_ANSYS(i,1) = 1/ni*Szz_nodei;
    Sxy_node_ANSYS(i,1) = 1/ni*Sxy_nodei;
    Syz_node_ANSYS(i,1) = 1/ni*Syz_nodei;
    Sxz_node_ANSYS(i,1) = 1/ni*Sxz_nodei;
end


%% %%%%%% Stress, strain calculation

%%%%%%%% First, we initialize matrices to store stress, strain components
Sxx_FEM = zeros(nE,8);
Syy_FEM = zeros(nE,8);
Szz_FEM = zeros(nE,8);
Sxy_FEM = zeros(nE,8);
Syz_FEM = zeros(nE,8);
Sxz_FEM = zeros(nE,8);

exx_FEM = zeros(nE,8);
eyy_FEM = zeros(nE,8);
ezz_FEM = zeros(nE,8);
exy_FEM = zeros(nE,8);
eyz_FEM = zeros(nE,8);
exz_FEM = zeros(nE,8);
% 
%%%%%%%%%% Second, we define natural coordinates (xi, eta,zeta) 
%%%%%%%%%% for our 8-node solid element
natural_Coord = 1*[ -1, -1, -1;
                   1, -1, -1;
                   1,  1, -1;
                  -1,  1, -1;

                  -1, -1, 1;
                   1, -1, 1;
                   1,  1, 1;
                  -1,  1, 1];
gauss_Points = 0*natural_Coord;

for ie = 1:nE                           
    id = elNode(ie,:); 
    elDof = [id, id+nP, id+2*nP];
    ndof = length(id);%%% ndof = 8

    %%% loop over Natural coordinates
    for j = 1:size(natural_Coord,1)
        xi_norm = natural_Coord(j,1);
        eta_norm = natural_Coord(j,2);
        zeta_norm = natural_Coord(j,3);
        
        xi_shear = gauss_Points(j,1);
        eta_shear = gauss_Points(j,2);
        zeta_shear = gauss_Points(j,3);
        
        %%%% shape functions and derivatives
        [shape_norm,nDeriva_norm]= shapeFuncH8(xi_norm,eta_norm,zeta_norm);
        
        [shape_shear,nDeriva_shear]= shapeFuncH8(xi_shear,eta_shear,zeta_shear);
        %%% Jacobian matrix, inverse of Jacobian
        [J_norm,xyzDeriv_norm] = Jacobian(xyz(id,:),nDeriva_norm);
        
        [J_shear,xyzDeriv_shear] = Jacobian(xyz(id,:),nDeriva_shear);
        %%%  B matrix (strain - displacement matrix)
        %%%%% this B matrix will be used to calculate normal strains,
        %%%%% stresses
        B_norm = zeros(6,3*ndof);%%% B(6x24)
        B_norm(1,1:ndof)             = xyzDeriv_norm(:,1)';
        B_norm(2,ndof+1:(2*ndof))    = xyzDeriv_norm(:,2)';
        B_norm(3,2*ndof+1:(3*ndof))  = xyzDeriv_norm(:,3)';

        B_norm(4,1:ndof)             = xyzDeriv_norm(:,2)';
        B_norm(4,ndof+1:(2*ndof))    = xyzDeriv_norm(:,1)';

        B_norm(5,ndof+1:(2*ndof))    = xyzDeriv_norm(:,3)';
        B_norm(5,2*ndof+1:(3*ndof))  = xyzDeriv_norm(:,2)';

        B_norm(6,2*ndof+1:(3*ndof))  = xyzDeriv_norm(:,1)';
        B_norm(6,1:ndof)             = xyzDeriv_norm(:,3)';
        
        
        
        %%%%% this B matrix will be used to calculate shear strains,
        %%%%% shear stresses
        B_shear = zeros(6,3*ndof);%%% B(6x24)
        B_shear(1,1:ndof)             = xyzDeriv_shear(:,1)';
        B_shear(2,ndof+1:(2*ndof))    = xyzDeriv_shear(:,2)';
        B_shear(3,2*ndof+1:(3*ndof))  = xyzDeriv_shear(:,3)';

        B_shear(4,1:ndof)             = xyzDeriv_shear(:,2)';
        B_shear(4,ndof+1:(2*ndof))    = xyzDeriv_shear(:,1)';

        B_shear(5,ndof+1:(2*ndof))    = xyzDeriv_shear(:,3)';
        B_shear(5,2*ndof+1:(3*ndof))  = xyzDeriv_shear(:,2)';

        B_shear(6,2*ndof+1:(3*ndof))  = xyzDeriv_shear(:,1)';
        B_shear(6,1:ndof)             = xyzDeriv_shear(:,3)';
        
        %%%% strain: e = [e_xx; e_yy; e_zz, 2*e_xy, 2*e_yz, 2*e_xz];
        %%%% We only need exx, eyy, ezz, Sxx, Syy, Szz from these below
        %%%% vectors
        Uei = disp(elDof',1);
        ei_norm = B_norm*Uei;
        Si_norm = C*ei_norm;
        
        %%%% We only need exy, eyz, exz, Sxy, Syz, Sxz from these below
        %%%% vectors
        ei_shear = B_shear*Uei;
        Si_shear = C*ei_shear;
        
        %%%% stress and strain matrix
        exx_FEM(ie,j) = ei_norm(1,1);
        eyy_FEM(ie,j) = ei_norm(2,1);
        ezz_FEM(ie,j) = ei_norm(3,1);
        exy_FEM(ie,j) = ei_shear(4,1);
        eyz_FEM(ie,j) = ei_shear(5,1);
        exz_FEM(ie,j) = ei_shear(6,1);
        

        Sxx_FEM(ie,j) = Si_norm(1,1);
        Syy_FEM(ie,j) = Si_norm(2,1);
        Szz_FEM(ie,j) = Si_norm(3,1);
        Sxy_FEM(ie,j) = Si_shear(4,1);
        Syz_FEM(ie,j) = Si_shear(5,1);
        Sxz_FEM(ie,j) = Si_shear(6,1);
    end
end


Sxx_node_FEM = zeros(nP,1);
Syy_node_FEM = zeros(nP,1);
Szz_node_FEM = zeros(nP,1);
Sxy_node_FEM = zeros(nP,1);
Syz_node_FEM = zeros(nP,1);
Sxz_node_FEM = zeros(nP,1);

exx_node_FEM = zeros(nP,1);
eyy_node_FEM = zeros(nP,1);
ezz_node_FEM = zeros(nP,1);
exy_node_FEM = zeros(nP,1);
eyz_node_FEM = zeros(nP,1);
exz_node_FEM = zeros(nP,1);

for i = 1:nP
    
    ni = 0;
    Sxx_nodei = 0;
    Syy_nodei = 0;
    Szz_nodei = 0;
    Sxy_nodei = 0;
    Syz_nodei = 0;
    Sxz_nodei = 0;
    
    exx_nodei = 0;
    eyy_nodei = 0;
    ezz_nodei = 0;
    exy_nodei = 0;
    eyz_nodei = 0;
    exz_nodei = 0;
    
    for j = 1:8
        idi = find(elNode(:,j)==i);
        ni = ni + length(idi);
        if ~isempty(idi)
            Sxx_nodei = Sxx_nodei + sum(Sxx_FEM(idi,j));
            Syy_nodei = Syy_nodei + sum(Syy_FEM(idi,j));
            Szz_nodei = Szz_nodei + sum(Szz_FEM(idi,j));
            Sxy_nodei = Sxy_nodei + sum(Sxy_FEM(idi,j));
            Syz_nodei = Syz_nodei + sum(Syz_FEM(idi,j));
            Sxz_nodei = Sxz_nodei + sum(Sxz_FEM(idi,j));

            exx_nodei = exx_nodei + sum(exx_FEM(idi,j));
            eyy_nodei = eyy_nodei + sum(eyy_FEM(idi,j));
            ezz_nodei = ezz_nodei + sum(ezz_FEM(idi,j));
            exy_nodei = exy_nodei + sum(exy_FEM(idi,j));
            eyz_nodei = eyz_nodei + sum(eyz_FEM(idi,j));
            exz_nodei = exz_nodei + sum(exz_FEM(idi,j));
        end
    end
    Sxx_node_FEM(i,1) = 1/ni*Sxx_nodei;
    Syy_node_FEM(i,1) = 1/ni*Syy_nodei;
    Szz_node_FEM(i,1) = 1/ni*Szz_nodei;
    Sxy_node_FEM(i,1) = 1/ni*Sxy_nodei;
    Syz_node_FEM(i,1) = 1/ni*Syz_nodei;
    Sxz_node_FEM(i,1) = 1/ni*Sxz_nodei;

    exx_node_FEM(i,1) = 1/ni*exx_nodei;
    eyy_node_FEM(i,1) = 1/ni*eyy_nodei;
    ezz_node_FEM(i,1) = 1/ni*ezz_nodei;
    exy_node_FEM(i,1) = 1/ni*exy_nodei;
    eyz_node_FEM(i,1) = 1/ni*eyz_nodei;
    exz_node_FEM(i,1) = 1/ni*exz_nodei;
end

Svm_node_FEM = sqrt(0.5*((Sxx_node_FEM-Syy_node_FEM).^2+...
    (Syy_node_FEM-Szz_node_FEM).^2 + (Szz_node_FEM-Sxx_node_FEM).^2)+...
    3*(Sxy_node_FEM.^2 + Syz_node_FEM.^2 + Sxz_node_FEM.^2));

Svm_node_ANSYS = sqrt(0.5*((Sxx_node_ANSYS-Syy_node_ANSYS).^2+...
    (Syy_node_ANSYS-Szz_node_ANSYS).^2 + (Szz_node_ANSYS-Sxx_node_ANSYS).^2)+...
    3*(Sxy_node_ANSYS.^2 + Syz_node_ANSYS.^2 + Sxz_node_ANSYS.^2));

%% %%%%% Plot S_vm
Svm_min = min(min(Svm_node_FEM),min(Svm_node_ANSYS));
Svm_max = max(max(Svm_node_FEM),max(Svm_node_ANSYS));

%%%%%% Svm_FEM
figure
scatter3(xyzNew(:,1),xyzNew(:,2),xyzNew(:,3),40,Svm_node_FEM,'filled');
title('S_{vm} (N/m)-FEM');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([Svm_min Svm_max]);
zstep = (Svm_max - Svm_min)/5;
set(h, 'ytick', (Svm_min:zstep:Svm_max))

%%%%%% Svm_ANSYS
figure
scatter3(xyzNew_ANS(:,1),xyzNew_ANS(:,2),xyzNew_ANS(:,3),40,Svm_node_ANSYS(:,1),'filled');
title('S_{vm} (N/m)-ANSYS');
axis equal
view(45,20)
h = colorbar;
colormap jet;
grid on;
caxis([Svm_min Svm_max]);
zstep = (Svm_max - Svm_min)/5;
set(h, 'ytick', (Svm_min:zstep:Svm_max))


