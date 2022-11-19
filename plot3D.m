
function plot3D(X,Y,Z,Lx,Ly,Lz,dx,dy,dz,T,T_min,T_max,isovalues,numisosurf,iter,t);
clf

% numisosurf=10; % number of isosurfaces

%% recalculate scales (disable this section to keep the static scale)
% T_max = max(T,[],'all');
% T_min =-T_max;%= min(T,[],'all');
% isovalues=linspace(T_min,T_max,numisosurf);

%% colorbar
% format bank
% cbh = colorbar ; %Create Colorbar
% cbh.Ticks = 1:numisosurf; %Create numisosurf ticks from 1 to numisosurf
% cbh.TickLabels = num2cell(round(linspace(T_min,T_max,numisosurf),2));    %Replace the labels of these numisosurf ticks with the numbers from T_min to T_max
 
%%
num=numel(isovalues);
 for i=1:num     
     p(i)=patch(isosurface(X,Y,Z,T,isovalues(i)));
     isonormals(X,Y,Z,T,p(i));
     set(p(i), 'CData',i); 
hold on
 end
set(p, 'CDataMapping','direct', 'FaceColor','flat', 'EdgeColor','none')
clr = jet(num);
colormap(clr);
caxis([0 num])
%colorbar

box off; grid on; axis tight; daspect([1 1 1])
h=gcf;
h.Color=[0 0 0];
axis off

camproj perspective
camlight; lighting gouraud; alpha(0.5);
axis([-Lx/2+dx Lx/2-dx -Ly/2+dx Ly/2-dy -Lz/2+dz Lz/2-dz]); % plot only half along Y to see inside
axis manual
%view(iter,20); % orbiting camera
view(330,20); %switch on to have static view

% view(3) % iso view
%title(['3D Wave equation solution in a box, iter = ' num2str(iter) ', t = ' num2str(t) ' s']);
drawnow

