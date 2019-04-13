function mm = Visualization(position,L,tau,T,videoname)
% ===          === %
%      Inputs      %
% ===          === %
%
% position is a matrix of x y z coordinates of each
%   particle
% N is the number of particles
% L is the length of the box
%
% ===          === %
%      Output      %
% ===          === %
%
% This code produces a video file of the simulation
%   video file name ===> 'points.avi'
% Frames of the movie are stored in the vector M

% ===   Start  === %
%mov = VideoWriter('output.avi');
% Initialize the figure space
close all;
figure();
mm=1;
% Set up the grid space. These values can be adjusted. 
% I wanted to have a grid with a bit of color, dashed grid lines, but
%   not too distracting. 
%taur = tau*0.5;
%taua = tau*0.25;
grid on


%A = [0 0 0;L 0 0;L L 0;0 L 0;0 0 0];
%B = [0 0 0;0 0 L;L 0 L;L 0 0];
%C = [L L 0;L L L;L 0 L];
%D = [0 L 0;0 L L;L L L];
%E = [0 L L;0 0 L];

%C0 = zeros(21,21,3);

ax=gca;
ax.GridColor = [0 0.5 0.5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
set(ax,'Xlim',[0 L]);
set(ax, 'Xtick',(0:tau:L));
set(ax,'Ylim',[0 L]);
set(ax, 'Ytick',(0:tau:L));
set(ax,'Zlim',[0 L]);
set(ax, 'Ztick',(0:tau:L));
ax.Layer = 'bottom';
axis manual;
hold(ax);
scatter3(0,0,0,'wo');
view(3);
%plot3(A(:,1),A(:,2),A(:,3),'k-');
%plot3(B(:,1),B(:,2),B(:,3),'k-');
%plot3(C(:,1),C(:,2),C(:,3),'k-');
%plot3(D(:,1),D(:,2),D(:,3),'k-');
%plot3(E(:,1),E(:,2),E(:,3),'k-');

% M is a vector that stores the frames of the movie. 
% M(1) places just figure contents without points into M

mov(1) = getframe(gcf);
frcounter = 2;
for t=1:T
    if(mod(t,200) == 0)
        hPlotData = scatter3(position(:,1,1,t),position(:,1,2,t),position(:,1,3,t), 'bo');
%    xc = position(i,1,1);
%    yc = position(i,1,2);
%    zc = position(i,1,3);
%    [x,y,z] = ellipsoid(xc,yc,zc,taua,taua,taua);
%    s = surf(x,y,z,C0);
        mov(frcounter) = getframe(gcf);
        delete(hPlotData);
        frcounter = frcounter + 1;
    end
    if(t>T-200)
        break;
    end
end
% Syntax for writing to a video: 
%   filehandle = VideoWriter('filename', Compression type)
video=VideoWriter(videoname)
open(video)
writeVideo(video,mov);
close(video);