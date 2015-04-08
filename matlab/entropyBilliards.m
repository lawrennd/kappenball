% Billiards2D is a two dimensional model of the dynamics of moving balls 
% contained in a rectangle. 'ctrl c' will kill the proces.
numBalls=10;
close all;
hold on;
boundBox=[-4 4 -4 4];
TableColour=[.4 .5 .8];
%Plothandle
axis(boundBox);
set(gcf,'color',[1 1 1]);
set(gca,'Color',TableColour,'xcolor',TableColour,'ycolor',...
TableColour,'PlotBoxAspectRatio',[1 abs((boundBox(3)-boundBox(4))/(boundBox(2)-boundBox(1))) 1],'xtick',[],'ytick',[])
%Set the radii=============================================================
r=0.125*ones(numBalls, 1); %0.125+0.125*rand(numBalls,1);
%Radii storage=============================================================

%Initial velocities========================================================
V=10*(-1+2*rand(numBalls,2));
X = initializeBilliards('rand', r, boundBox);
simulateBilliards(X, V, r, boundBox);
