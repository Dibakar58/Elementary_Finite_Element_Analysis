% % clear
tol = 1e-5;

%% %%%% Material properties
E = 2e11; %%%% elastic modulus
pois = 0.27;%%% Poisson's ratio
dens = 7850;%%% mass density

% % %%% Stress- strain matrix: C - Plane strain condition
% % C = E*(1-pois)/(1-2*pois)/(1+pois)*[1 pois/(1-pois) 0;
% %                                     pois/(1-pois) 1 0;
% %                                     0 0 (1-2*pois)/2/(1-pois)];
% % %%%% Lame's constants
% % lamda = pois*E/(1-2*pois)/(1+pois);
% % muy = E/2/(1+pois);


%%%%% Stress- strain matrix: C - Plane stress condition
C = E/(1-pois^2)*[1 pois 0;
                  pois 1 0;
                  0 0 (1-pois)/2];
% % lamda = pois*E/(1-pois^2);
% % muy = E/2/(1+pois);%%%% Shear modulus

%% %%%% Geometrical properties
Lx = 1;
Ly = 1;
h = 0.01;%%% Plate thickness

nx = 20;%%%% number of division in x-direction
ny = 20;%%%% number of division in y-direction

dx = Lx/nx;%%%% mesh size in x-direction
dy = Ly/ny;%%%% mesh size in y-direction

%% %%% Create mesh by using "rectangularMesh" function
[xy,elNode,nE,nP] = MeshRectangular(Lx,Ly,nx,ny);
nDof = 2*nP;%%% Number of global degrees of freedom

% % %%%%% Try to plot mesh
% % PlotMesh(xy,elNode);

%% %%%% boundary conditions - the left edge is clamped
left = find(xy(:,1)==0); %%%% left edge at x = 0
fixedDof = [left; left + nP];%%% u = v = 0 at left edge

%% %%%% Loading conditions
fx = 5e8;%%%% Distributed force along the right edge: N/m
fy = 0;%1.2e8;%%%% Distributed force along the right edge: N/m

force = zeros(nDof,1);

right = find(xy(:,1)>= Lx-tol);

%%%% Convert the distributed load fx to nodal forces
force(right) = fx*dy;
force(right(1)) = fx*dy/2;
force(right(end)) = fx*dy/2;

%% %%%%% Calculate stiffness matrix
[stiff,mass] = formStiffness2D(nDof,nE,elNode,nP,xy,C,dens,h);

%% %%%% Solve the equilibrium equation
disp = solution(nDof,fixedDof,stiff,force);

save('PlStress0P27');%%% save data for later use

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
title('u (m)-FEM');%,'Fontsize',14)

min_v_FEM = min(uy);
max_v_FEM = max(uy);
figure;
PlotFieldonMesh(xy_NewFEM,elNode,uy,min_v_FEM,max_v_FEM);
title('v (m)-FEM');%,'Fontsize',14)

%% %%%%%%%% Compare with ANSYS results
%%%%%%%%%%% Step 1: process ANSYS data
NLIST2 = [];
for i = 1:size(NLIST,1)
    if isnan(NLIST(i,1)) ==0
        NLIST2 = [NLIST2;NLIST(i,:)];
    end
end

ELIST2 = [];
for i = 1:size(ELIST,1)
    if isnan(ELIST(i,1)) ==0
        ELIST2 = [ELIST2;ELIST(i,:)];
    end
end

dispANSYS2 = [];
for i = 1:size(dispANSYS,1)
    if isnan(dispANSYS(i,1)) ==0
        dispANSYS2 = [dispANSYS2;dispANSYS(i,:)];
    end
end


%%%%%%%%%%% Step 2: Compare displacement fields: u
udof = (1:nP)';vdof = udof+nP;
ux = disp(udof,1);uy = disp(vdof,1);

min_u_FEM = min(ux);
max_u_FEM = max(ux);

min_u_ANSYS = min(dispANSYS2(:,2));
max_u_ANSYS = max(dispANSYS2(:,2));

min_u = min(min_u_FEM,min_u_ANSYS);
max_u = max(max_u_FEM,max_u_ANSYS);

scale = 1;
xy_NewFEM = xy+scale*[ux uy];
xy_NewANSYS = NLIST2(:,1:2) + scale*dispANSYS2(:,2:3);

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
PlotFieldonMesh(xy_NewANSYS,ELIST2(:,7:10),dispANSYS2(:,2),min_u,max_u);
title('u (m)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')
view(2)


%%%%%%%%%%% Step 3: Compare displacement fields: u
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
PlotFieldonMesh(xy_NewANSYS,ELIST2(:,7:10),dispANSYS2(:,3),min_v,max_v);
title('v (m)-ANSYS');
xlabel('x (m)');
ylabel('y (m)');
grid on
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')

%%%%%%%%% FINISH!! %%%%%%%%%%%%%%
%%%%%%%%% Post-processing for strain and stress will be done in the next Lecture %%%%%%%%%%%%%%
%%%%%%%%% Please save data if you want to use it later %%%%%%%%%%
