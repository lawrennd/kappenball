function balls = simulateBilliards(balls, barHandle)

global PAUSEDEMO
global SHOWARROWS
global RUNDEMO
global INSTANTVELOCITIES

RUNDEMO = true;
velFactor = 0.025;

if nargin < 2
  barHandle = [];
end
tableAx = get(balls(1).handle, 'parent');
X = [balls(:).X]';
r = [balls(:).r]';
V = [balls(:).V]';
ballColor = [balls(:).color]';

xlim = get(tableAx, 'xlim');
ylim = get(tableAx, 'ylim');
drawflag=1;
factor=800; 
DT=2e-3;
mass=pi*r.^2;

numBalls = length(balls);
distMatrix = triu(sqrt(dist2(X, X))); 
rmatrix = repmat(r, 1, numBalls) + repmat(r', numBalls, 1);
rmatrix = rmatrix.*triu(abs(-1+eye(size(rmatrix))));

count = 0;
while RUNDEMO
  for k=1:numBalls;
    set(balls(k).handle, 'xdata', X(k,1), 'ydata', X(k,2), 'markersize', ...
                      3*getAxisUnitInPts(tableAx)*balls(k).r*2);
    rt = balls(k).r*normr([V(k, 1) V(k, 2)]);
    set(balls(k).vhandle, 'xdata', X(k, 1) + [rt(1) V(k, 1)*velFactor], ...
                      'ydata', X(k, 2) + [rt(2) V(k, 2)*velFactor]);
  end
  drawnow
 if PAUSEDEMO
     pause(0.01)
 else
  %Edgedetecton positive==================================================
  d=X+repmat(r,1,2)-repmat([xlim(2), ylim(2)],numBalls,1);
  dt=(d>=0).*d./V;
  X=X-V.*dt;
  V=V.*(2*(d>=0==0)-1);
  %Edgedetecton negative==================================================
  d=X-repmat(r,1,2)-repmat([xlim(1), ylim(1)],numBalls,1);
  dt=(d<=0).*d./V;
  X=X-V.*dt;
  V=V.*(2*(d<=0==0)-1);
  %Distancematrix=========================================================
  distMatrix = triu(sqrt(dist2(X, X))); 
  %collisiondetectionmatrix===============================================
  collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
  %=======================================================================
  if find(collisionMatrix<0);
    [I,J]=find(collisionMatrix<0);
    for i=1:length(I)
      normdist=normr([X(I(i),1)-X(J(i),1) X(I(i),2)-X(J(i),2)]);
      %Velocity component along the line connecting the two centres of
      %ball A and ball B:===============================================
      vaA=(V(I(i),1)*normdist(1)+V(I(i),2)*normdist(2));
      vaB=(V(J(i),1)*normdist(1)+V(J(i),2)*normdist(2));
      dt=abs(r(I(i))+r(J(i))-distMatrix(I(i),J(i)))/(abs(vaA)+abs(vaB));
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
      collisionEffectMatrix=[
          1-2/f 0   2/f   0; % 2:Inellastic collision=======================
          0     1   0     0;
          2/g   0   1-2/g 0;
          0     0   0     1];
      v_new_col=collisionEffectMatrix*v_new;
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
  if ~isempty(barHandle)
    count = count + 1;
    centers = get(barHandle, 'xdata');
    barVals = hist(V(:), centers);
    % barValStore tries to estimate the marginal velocity probability.
    if exist('barValStore')
        % Drop first 1000
        if count > 1000
            if ~rem(count, 30) % retain every 30 
                barValStore = barVals + barValStore;
            end
        end
    else
        if count > 1000
            barValStore = barVals;
        else
            barValStore = zeros(size(barVals));
        end
    end
    if INSTANTVELOCITIES
        set(barHandle, 'ydata', barVals);
    else       
        set(barHandle, 'ydata', barValStore/sum(barValStore)* ...
                       numBalls);
    end

  end
 end
end