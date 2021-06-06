clear ;


E = 2e11; %%%% elastic modulus
pois = 0.27;%%% Poisson's ratio
dens = 7850;%%% mass density

C=(E/(1-pois^2))*[1,pois,0;pois,1,0;0,0,(1-pois)/2] ;

Lx=1;
Ly=1;
h=0.01 ;
nx=20;
ny=20;

dx=Lx/nx ;
dy=Ly/ny ;

[xy,elnod,nE,nP] = MeshRectangular(Lx,Ly,nx,ny);

PlotMesh(xy,elnod);

ndof=nP*2 

K=stiff() ;



