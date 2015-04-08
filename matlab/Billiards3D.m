% Billiards3D is a three dimensional model of the dynamics of moving balls 
% contained in a box. 'ctrl c' will kill the proces.
updateBalls = true;
updateBar = true
frameRate = 24;
numBalls=40;
close all;
hold on;
drawflag=1;
DT=1e-2;
boundBox=[-4 4 -4 4 -4 4];
BoxColour=[.6 .7 .9];
maxSpeedColor = 60;
% Initial velocities
T = 1e-4; % temperature (below 0.95 K it would be solid :-(
boltzK = 1.3806503e-23; % boltzmann constant
m = 1.66053886e-27; % 1 atomic mass unit in kg.

%Set the radii=============================================================
r=0.25*ones(numBalls,1);

%Mass======================================================================
mass=2 * r*m;  %Give everything mass of Helium. %4/3*pi*r.^3;

varV = boltzK*T./mass;
vStd = sqrt(varV);

V=[ones(floor(numBalls/2), 3); ones(ceil(numBalls/2), 3)].*repmat(vStd, 1,3); %
                                                                                                            %V = randn(numBalls,3).*repmat(vStd, 1, 3);

% Bar chart bins
numBins = 20;

%Figure 
%set(gcf,'color',[1 1 1]);
clf
barAx = gca;
[barVals, barX] = hist(V(:), numBins);
barHandle = bar(barX, barVals);
set(barAx, 'xlim', [-3 3]*sqrt(mean(varV)));
set(barAx, 'ylim', [0 numBalls])
%set(barHandle, 'erasemode', 'xor')
%ballAx = gca;
figure
ballAx = gca; %subplot(1, 2, 1)
hold on
axis(boundBox);
set(ballAx,'Color',BoxColour,'xcolor',[0.5 0.5 0.5],'ycolor',[0.5 0.5 0.5],'zcolor',...
           [0.5 0.5 0.5],'CameraPosition',[10 6 7],'Projection','perspective',...
           'xtick',[],'ytick',[],'ztick',[],'LineWidth',2)
zoom(1.2);
box on;
axis equal;
%Set the mesh for plotting the balls=======================================
[ballx,bally,ballz]=sphere;
%Radii storage=============================================================
rmatrix = repmat(r, 1, length(r)) + repmat(r', length(r), 1);
% for j=1:numBalls;
%   for i=1:j,
%     rmatrix(i,j)=r(j)+r(i);
%   end;
% end;
rmatrix=rmatrix.*triu(abs(-1+eye(size(rmatrix))));


X=[(boundBox(2)-boundBox(1)-2*max(r))*rand(numBalls,1)+boundBox(1)+max(r),...
   (boundBox(4)-boundBox(3)-2*max(r))*rand(numBalls,1)+boundBox(3)+max(r),...
   (boundBox(6)-boundBox(5)-2*max(r))*rand(numBalls,1)+boundBox(5)+max(r)];
% Interball distances of centers.
distMatrix = triu(sqrt(dist2(X, X))); 
%Initial Collisiondetectionmatrix==========================================
collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
while find(collisionMatrix<=0);
    % Initial Locations
    X=[(boundBox(2)-boundBox(1)-2*max(r))*rand(numBalls,1)+boundBox(1)+max(r),...
       (boundBox(4)-boundBox(3)-2*max(r))*rand(numBalls,1)+boundBox(3)+max(r),...
       (boundBox(6)-boundBox(5)-2*max(r))*rand(numBalls,1)+boundBox(5)+max(r)];
    % X=[(boundBox(2)-boundBox(1)-2*max(r))*rand(numBalls,1)+boundBox(1)+max(r),...
    % (boundBox(4)-boundBox(3)-2*max(r))*rand(numBalls,1)+boundBox(3)+max(r),...
    % (boundBox(6)-boundBox(5)-2*max(r))*rand(numBalls,1)+boundBox(5)+max(r)];
    distMatrix = triu(sqrt(dist2(X, X)));
    collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
end

%Initial velocities========================================================
disp(sum(sum(V.*V)))
redVals = gammainc(sum(V.*V, 2).*mass./(boltzK*T), 3);
ballColour = [redVals zeros(numBalls, 2)];
%Plot startingpositions====================================================
for k=1:numBalls;
  xvals = ballx*r(k)+X(k,1);
  yvals = bally*r(k)+X(k,2);
  zvals = ballz*r(k)+X(k,3);
  ballHand(k) = surf(xvals, yvals, zvals, ...
                     'LineStyle','none',...
                     'FaceColor', ballColour(mod(k-1,length(ballColour))+1,:), ...
                       'AmbientStrength',0.5);
end
tic
light;
lighting gouraud
drawnow;
%Loop
while drawflag==1;
  %Edgedetecton positive==================================================
  d=X+repmat(r,1,3)-repmat([boundBox(2),boundBox(4),boundBox(6)],numBalls,1);
  dt=(d>=0).*d./V;
  X=X-V.*dt;
  V=V.*(2*(d>=0==0)-1);
  %Edgedetecton negative==================================================
  d=X-repmat(r,1,3)-repmat([boundBox(1),boundBox(3),boundBox(5)],numBalls,1);
  dt=(d<=0).*d./V;
  X=X-V.*dt;
  V=V.*(2*(d<=0==0)-1);
  %Distancematrix=========================================================
  distMatrix = triu(sqrt(dist2(X, X)));   
  %Collisiondetectionmatrix===============================================
  collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
  %=======================================================================
  if find(collisionMatrix<0);
    [I,J]=find(collisionMatrix<0);
    for i=1:length(I)
      normdist=normr([X(I(i),1)-X(J(i),1) X(I(i),2)-X(J(i),2) X(I(i),3)-X(J(i),3)]);      
      %Velocity component along the line connecting the two centres of
      %ball A and ball B:===============================================
      vaA=(V(I(i),1)*normdist(1)+V(I(i),2)*normdist(2)+V(I(i),3)*normdist(3));
      vaB=(V(J(i),1)*normdist(1)+V(J(i),2)*normdist(2)+V(J(i),3)*normdist(3));
      dt=abs(r(I(i))+r(J(i))-distMatrix(I(i),J(i)))/(abs(vaA)+abs(vaB));
      %Set back the positions of the colliding balls==================== 
      X(I(i):J(i),:)=X(I(i):J(i),:)-V(I(i):J(i),:)*dt;
    end
    for i=1:length(I)
      normdist=normr([X(I(i),1)-X(J(i),1) X(I(i),2)-X(J(i),2) X(I(i),3)-X(J(i),3)])';
      e2=[normdist(3),normdist(2)*normdist(3)/(normdist(1)-1),1+normdist(3)^2/(normdist(1)-1)]';
      %Coordinate transformation matrix=================================
      M=[normdist,e2,cross(normdist,e2)];
      v_old=[V(I(i),:)';V(J(i),:)'];
      %Calculate velocities in the new coordinate system================
      v_new=[M' zeros(3);zeros(3) M']*v_old;
      f=(1+mass(I(i))/mass(J(i)));
      g=(1+mass(J(i))/mass(I(i)));
      CollisionEffectMatrix=sparse([
          1-2/f 0 0   2/f 0 0 % 2:Inellastic collision======================
          0     1 0     0 0 0;
          0     0 1     0 0 0;
          2/g   0 0 1-2/g 0 0;
          0     0 0     0 1 0;
          0     0 0     0 0 1]);
      v_new_col=CollisionEffectMatrix*v_new;
      %Put the velocities in the old coordinate system==================
      V(I(i),:)=M*v_new_col(1:3);
      V(J(i),:)=M*v_new_col(4:end);

      %Update the positions after collision=============================
      X(I(i),:)=X(I(i),:)+V(I(i),:)*dt;
      X(J(i),:)=X(J(i),:)+V(J(i),:)*dt;         
    end
         %
    
  end
   %Propagation============================================================
   X=X+V*DT;
   if updateBar
       [barVals, barX] = hist(V(:), numBins);
       set(barHandle, 'ydata', barVals, 'xdata', barX);
       drawnow
   end
   %Plotting===============================================================
   if updateBalls
       if toc > 1/frameRate
           redVals = gammainc(sum(V.*V, 2).*mass./(boltzK*T), 3);
           ballColour = [redVals zeros(numBalls, 2)];
           for k=1:numBalls;
           
           xvals = ballx*r(k)+X(k,1);
           yvals = bally*r(k)+X(k,2);
           zvals = ballz*r(k)+X(k,3);
           set(ballHand(k), 'xdata', xvals, 'ydata', yvals, 'zdata', zvals, ...
                            'facecolor', ballColour(k, :));
           %set(ballHand(i), 'xdata', ballx*r(k)+X(k,1), 'ydata', ...
           %                 , 'zdata', );
           %       surf(,  , ,'LineStyle','none',...
           %     'FaceColor',ballColour(mod(k-1,length(ballColour))+1,:),'AmbientStrength',0.5);
           tic
       end
       light
       drawnow;
       
   end
   end
   % light;
   % lighting gouraud
end