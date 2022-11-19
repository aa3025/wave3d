% Wave equation in 3D, explicit central difference in space and time

clear; 
close all; clc;

%%
video=1; % record video? 
videofilename='wave_eq_3D.mp4';
videoformat='MPEG-4';
if video>0
%     close(writerObj);
    writerObj = VideoWriter(videofilename, videoformat);  % (see also lines 179+)
    open(writerObj); 
end

%% domain
Lx= 1;                  %  width (m)
Ly= 1;                  %  depth (m)
Lz= 1;                  %  height (m)
Nx=101;                 % nodes in x direction (odd number to have central plane to slice through)
Ny=101;                 % nodes in y direction (odd number to have central plane to slice through)
Nz=101;                  % nodes in z direction (odd number to have central plane to slice through)
dx=Lx/Nx;               % spacing along x
dy=Ly/Ny;               % spacing along y
dz=Lz/Nz;               % spacing along z

xc=(Nx+1)/2;% center plane of the domain
yc=(Ny+1)/2;
zc=(Nz+1)/2;

%dt=(C/v)*dx; %stable time step
tmax=7200; % max time, sec

%% Initial and boundary conditions (1st type, a.k.a Dirichlet)
% T(y,x,z) % -- wave function, structure of 3D array

T_0= 0;                 % Initial temperature in all nodes 
T_right  = 0;        % temperature at x=Lx
T_left   = 0;        % temperature at x=0 

T_top  = 0;             % temperature at z=Lz 
T_bott  = 0;            % temperature at z=0  

T_front=0;            % temperature at y=0 
T_back=0;            % temperature at y=Ly

%% Initial conditions for Steady State
T=zeros(Ny,Nx,Nz);
T_prev=zeros(Ny,Nx,Nz);
T_next=zeros(Ny,Nx,Nz);

T(:,1,:)    =  T_left;
T(:,Nx,:)   =  T_right;

T(1,:,:)    =  T_front;
T(Ny,:,:)   =  T_back;

T(:,:,1)    =  T_bott;
T(:,:,Nz)   =  T_top;

%% media props
C=.001; % Courant number <<1

% if needed we can compute time
v=0.0025;  % wave speed for time calculations
dt=C/v*dx;

%% Domain

x=linspace(-Lx/2,Lx/2,Nx);
y=linspace(-Ly/2,Ly/2,Ny);
z=linspace(-Lz/2,Lz/2,Nz);
[X Y Z]=meshgrid(x,y,z);

%% Initial Perturbation D (initial condition)
amplitude=20; 
radius=Nx; % size of perturbation

[x,y,z] = ndgrid(-1:(2/(radius-1)):1);
D = amplitude*exp(-150*(x.^2+y.^2+z.^2)); % Smooth Gaussian-style spherical distribution, the higher ceofficient ("150") -- the tighter (to the centre) teh distribution

w = size(D,1); % width 
dropcentre=0.5; % centre of perturbation in domain;

i1 = ceil(dropcentre*(Ny-w))+(1:w);
j1 = ceil(dropcentre*(Nx-w))+(1:w);
k1 = ceil(dropcentre*(Nz-w))+(1:w);

T(i1,j1,k1)=D;%T(i1,j1,k1)+D;
T_prev(i1,j1,k1)=T(i1,j1,k1); % "previous" T

%% Constants
iter=1;  % counter of iterations
error  = 1; % initial error for iterations
t=0;

%% plotting 
opengl hardware % or "software"
  T_max = .25; % limits for plotting isosurfaces and colors (need care for nicer pictures)
  T_min = -.25;

numisosurf=10; % number of isosurfaces
isovalues=linspace(T_min,T_max,numisosurf);

hhh=figure('units','pixels','position',[0 0 1920 1080]);
% set(hhh, 'units', 'normalized', 'outerposition', [0 0 2 2],'Visible','off');

daspect([1 1 1])
caxis([T_min T_max]);
colorbar
caxis manual
shading interp

whitebg([0 0 0]); % dark background
drawnow
pause
%% Solution
tic
i=2:Ny-1; %inner nodes along y
j=2:Nx-1; %inner nodes along x
k=2:Nz-1; %inner nodes along z

T_new=T; % T at next timestep (further updated below)
% T is now
% T_prev is at previous timestep
% T_next is at next timestep
% save(['data/' num2str(iter) '.mat'],'T');
cam=0;
while t<=tmax
filename_prev=['data/' num2str(iter) '.mat'];
filename=['data/' num2str(iter+1) '.mat'];
 % central in time and space central order finite differences (2nd order derivatives)
 
 if isfile(filename)
     load(filename,'T');     
 else 
     T_next(i,j,k) = 2*T(i,j,k)-T_prev(i,j,k)+C^2*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dy^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dx^2)+((T(i,j,k-1)-2*T(i,j,k)+T(i,j,k+1))/dz^2));
 end

    t=iter*dt; 
 
if mod(iter,1)==0 
    cam=cam+0.25;
    camposition=[2*cosd(cam),2*sind(cam),0.6];
    plot3D(X,Y,Z,Lx,Ly,Lz,dx,dy,dz,T,T_min,T_max,isovalues,numisosurf,iter,t,camposition); 
end
 save(['data/' num2str(iter) '.mat'],'T')
iter=iter+1

% shifting T's in time 
T_prev=T;
T = T_next;

if video >0 
         frame = getframe(gcf);
         writeVideo(writerObj,frame);
     end

 end
 
%% 
toc 
% plot3D(X,Y,Z,Lx,Ly,Lz,Tss,isovalues,iter); % final plot

 if video >0  close(writerObj); end
 


disp('COMPLETE!');




























