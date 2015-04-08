% Billiards2D is a two dimensional model of the dynamics of moving balls 
% contained in a rectangle. 'ctrl c' will kill the proces.
NumberOfBalls=7;
close all;
hold on;
drawflag=1;
factor=82; %Adjust this factor when scaling================================
DT=2e-3;
Bound=[-4 4 -2 2];
BallColour=[[1 0 0];[1 0 0.5];[1 0.5 0];[0 1 0];[0 0 1];[1 1 0];[1 1 1];...
[0 0.3 0];[0 0 0];[0.65 0.65 0.65];[0 0.75 0.75];[0.3 0 0.6];[0.95 0.65 0.75];...
[0.5 0.25 0];[0 0.2 0.4];[0.9 0.4 0.7];[0.4 0.2 0.3];[0.65 0.55 0.15];[0.25 0.35 0.25];[0.5 0 0]];
TableColour=[.4 .5 .8];
%Plothandle
axis(Bound);
set(gcf,'color',[1 1 1]);
set(gca,'Color',TableColour,'xcolor',TableColour,'ycolor',...
TableColour,'PlotBoxAspectRatio',[1 abs((Bound(3)-Bound(4))/(Bound(2)-Bound(1))) 1],'xtick',[],'ytick',[])
%Set the radii=============================================================
r=0.25+0.25*rand(NumberOfBalls,1);
%Radii storage=============================================================
for j=1:NumberOfBalls;
   for i=1:j,
      rmatrix(i,j)=r(j)+r(i);
   end;
end;
rmatrix=rmatrix.*triu(abs(-1+eye(size(rmatrix))));
%Mass =====================================================================
mass=pi*r.^2;
X=[(Bound(2)-Bound(1)-2*max(r))*rand(NumberOfBalls,1)+Bound(1)+max(r),...
(Bound(4)-Bound(3)-2*max(r))*rand(NumberOfBalls,1)+Bound(3)+max(r)];
%==========================================================================
for j=1:NumberOfBalls;
   for i=1:j;
      distmatrix(i,j)=sqrt((X(j,1)-X(i,1))^2+(X(j,2)-X(i,2))^2);
   end;
end;
%Initial Collisiondetectionmatrix==========================================
CollisionMatrix=(distmatrix-rmatrix)+tril(abs(-1+eye(size(distmatrix))))+eye(size(distmatrix));
while find(CollisionMatrix<=0);
   X=[(Bound(2)-Bound(1)-2*max(r))*rand(NumberOfBalls,1)+Bound(1)+max(r),...
   (Bound(4)-Bound(3)-2*max(r))*rand(NumberOfBalls,1)+Bound(3)+max(r)];
   for j=1:NumberOfBalls;
      for i=1:j;
         distmatrix(i,j)=sqrt((X(j,1)-X(i,1))^2+(X(j,2)-X(i,2))^2);
      end;
   end;
   CollisionMatrix=(distmatrix-rmatrix)+tril(abs(-1+eye(size(distmatrix))))+eye(size(distmatrix));
end
%Initial velocities========================================================
V=10*(-1+2*rand(NumberOfBalls,2));
%Plot startingpositions====================================================
for k=1:NumberOfBalls;
   plot(X(k,1),X(k,2),'o','MarkerEdgeColor',BallColour(mod(k-1,length(BallColour))+1,:),...
   'MarkerFaceColor',BallColour(mod(k-1,length(BallColour))+1,:),'MarkerSize',factor*r(k));
end
drawnow;
%Loop======================================================================
while drawflag==1;
   cla;
   %Edgedetecton positive==================================================
   d=X+repmat(r,1,2)-repmat([Bound(2),Bound(4)],NumberOfBalls,1);
   dt=(d>=0).*d./V;
   X=X-V.*dt;
   V=V.*(2*(d>=0==0)-1);
   %Edgedetecton negative==================================================
   d=X-repmat(r,1,2)-repmat([Bound(1),Bound(3)],NumberOfBalls,1);
   dt=(d<=0).*d./V;
   X=X-V.*dt;
   V=V.*(2*(d<=0==0)-1);
   %Distancematrix=========================================================
   for j=1:NumberOfBalls;
      for i=1:j;
         distmatrix(i,j)=sqrt((X(j,1)-X(i,1))^2+(X(j,2)-X(i,2))^2);
      end;
   end;
   %Collisiondetectionmatrix===============================================
   CollisionMatrix=(distmatrix-rmatrix)+tril(abs(-1+eye(size(distmatrix))))+eye(size(distmatrix));
   %=======================================================================
   if find(CollisionMatrix<0);
      [I,J]=find(CollisionMatrix<0);
      for i=1:length(I)
         normdist=normr([X(I(i),1)-X(J(i),1) X(I(i),2)-X(J(i),2)]);
         %Velocity component along the line connecting the two centres of
         %ball A and ball B:===============================================
         vaA=(V(I(i),1)*normdist(1)+V(I(i),2)*normdist(2));
         vaB=(V(J(i),1)*normdist(1)+V(J(i),2)*normdist(2));
         dt=abs(r(I(i))+r(J(i))-distmatrix(I(i),J(i)))/(abs(vaA)+abs(vaB));
         %Set back the positions of the colliding balls==================== 
         %X(I(i):J(i),:)=X(I(i):J(i),:)-V(I(i):J(i),:)*dt;
      end
      for i=1:length(I)
         normdist=normr([X(I(i),1)-X(J(i),1) X(I(i),2)-X(J(i),2)]);
         %Coordinate transformation matrix=================================
         M=[normdist(1) -normdist(2);
         normdist(2)  normdist(1)];
         v_old=[V(I(i),:)';V(J(i),:)'];
         %Calculate velocity in the new coordinate system==================
         v_new=[M' zeros(2);zeros(2) M']*v_old;
         f=(1+mass(I(i))/mass(J(i)));
         g=(1+mass(J(i))/mass(I(i)));
         CollisionEffectMatrix=[
         1-2/f 0   2/f   0; % 2:Inellastic collision=======================
         0     1   0     0;
         2/g   0   1-2/g 0;
         0     0   0     1];
         v_new_col=CollisionEffectMatrix*v_new;
         %Put the velocities in the old coordinate system==================
         V(I(i),:)=M*v_new_col(1:2);
         V(J(i),:)=M*v_new_col(3:end);
         %Update the positions after collision=============================
         X(I(i),:)=X(I(i),:)+V(I(i),:)*dt;
         X(J(i),:)=X(J(i),:)+V(J(i),:)*dt;    
      end
   end
   %Propagation============================================================
   X=X+V*DT;
   %Plotting===============================================================
   for k=1:NumberOfBalls;
       plot(X(k,1),X(k,2),'o','MarkerEdgeColor',BallColour(mod(k-1,length(BallColour))+1,:),...
      'MarkerFaceColor',BallColour(mod(k-1,length(BallColour))+1,:),'MarkerSize',factor*r(k));
   end
   drawnow;
end