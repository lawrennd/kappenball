function balls = initializeBilliards(numBalls, type, tableAx)


global SHOWARROWS

xlim = get(tableAx, 'xlim');
ylim = get(tableAx, 'ylim');

if nargin < 1
    type = 'rand';
end
switch type
 case 'hotCold'
  T = 1e-2; % temperature (below 0.95 K it would be solid :-(
  boltzK = 1.3806503e-23; % boltzmann constant
  m = 1.66053886e-27; % 1 atomic mass unit in kg.

  mass = 2*m;  %Give everything mass of Helium. %4/3*pi*r.^3;
  varV = boltzK*T./mass*ones(numBalls, 1);
  vStd = sqrt(varV);

  V = [[2*ones(ceil(numBalls/2), 1); 0.5*ones(floor(numBalls/2), 1)] ...
       randn(numBalls, 1)*0.004].*[vStd vStd];
  spacing = (ylim(2) - ylim(1))/(numBalls);
  rVal = spacing*0.475;
  X = [((xlim(2)-xlim(1))*0.5+xlim(1))*ones(numBalls, 1) ...
       (ylim(1)+0.5*spacing:spacing:ylim(2)-0.5*spacing)'];
  for i = 1:numBalls
    balls(i).X = X(i, :)';
    balls(i).V = V(i, :)';
    if i <= ceil(numBalls/2)
        balls(i).color = [1 1 0];
    else
        balls(i).color = [0 1 1];
    end
    balls(i).r = rVal;
  end
  
 case 'uniform'
  T = 1e-2; % temperature (below 0.95 K it would be solid :-(
  boltzK = 1.3806503e-23; % boltzmann constant
  m = 1.66053886e-27; % 1 atomic mass unit in kg.

  mass = 2*m;  %Give everything mass of Helium. %4/3*pi*r.^3;
  varV = boltzK*T./mass*ones(numBalls, 1);
  vStd = sqrt(varV);

  V = [[ones(ceil(numBalls/2), 1); -ones(floor(numBalls/2), 1)] ...
       randn(numBalls, 1)*0.001].*[vStd vStd];
  spacing = (ylim(2) - ylim(1))/(numBalls);
  rVal = spacing*0.475;
  X = [((xlim(2)-xlim(1))*0.5+xlim(1))*ones(numBalls, 1) ...
       (ylim(1)+0.5*spacing:spacing:ylim(2)-0.5*spacing)'];
  for i = 1:numBalls
    balls(i).X = X(i, :)';
    balls(i).V = V(i, :)';
    balls(i).color = [1 1 0];
    balls(i).r = rVal;
  end
  
 case 'randn'
  T = 1e-2; % temperature (below 0.95 K it would be solid :-(
  boltzK = 1.3806503e-23; % boltzmann constant
  m = 1.66053886e-27; % 1 atomic mass unit in kg.
  r = 0.0625*ones(numBalls, 1);   

  mass = 2*m;  %Give everything mass of Helium. %4/3*pi*r.^3;
  varV = boltzK*T./mass*ones(numBalls, 1);
  vStd = sqrt(varV);

  rmatrix = repmat(r, 1, numBalls) + repmat(r', numBalls, 1);
  rmatrix = rmatrix.*triu(abs(-1+eye(size(rmatrix))));
  V = randn(numBalls, 2).*[vStd vStd];
  X=[(xlim(2)-xlim(1)-2*max(r))*rand(numBalls,1)+xlim(1)+max(r),...
     (ylim(2)-ylim(1)-2*max(r))*rand(numBalls,1)+ylim(1)+max(r)];
  
  distMatrix = triu(sqrt(dist2(X, X))); 
  
  collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
  while find(collisionMatrix<=0);
    X=[(xlim(2)-xlim(1)-2*max(r))*rand(numBalls,1)+xlim(1)+max(r),...
       (ylim(2)-ylim(1)-2*max(r))*rand(numBalls,1)+ylim(1)+max(r)];
    distMatrix = triu(sqrt(dist2(X, X))); 
    collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
  end
  for i = 1:numBalls
    balls(i).X = X(i, :)';
    balls(i).V = V(i, :)';
    balls(i).color = [1 1 0];
    balls(i).r = r(i);
  end

 case 'rand'
  r=0.125+0.125*rand(numBalls,1);
  rmatrix = repmat(r, 1, numBalls) + repmat(r', numBalls, 1);
  rmatrix = rmatrix.*triu(abs(-1+eye(size(rmatrix))));
  V=10*(-1+2*rand(numBalls,2));
  
  X=[(xlim(2)-xlim(1)-2*max(r))*rand(numBalls,1)+xlim(1)+max(r),...
     (ylim(2)-ylim(1)-2*max(r))*rand(numBalls,1)+ylim(1)+max(r)];
  
  distMatrix = triu(sqrt(dist2(X, X))); 
  
  collisionMatrix=(distMatrix-rmatrix)+trEil(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
  while find(collisionMatrix<=0);
    X=[(xlim(2)-xlim(1)-2*max(r))*rand(numBalls,1)+xlim(1)+max(r),...
       (ylim(2)-ylim(1)-2*max(r))*rand(numBalls,1)+ylim(1)+max(r)];
    distMatrix = triu(sqrt(dist2(X, X))); 
    collisionMatrix=(distMatrix-rmatrix)+tril(abs(-1+eye(size(distMatrix))))+eye(size(distMatrix));
  end
  ballColor=[1 0 0
             1 0 0.5
             1 0.5 0
             0 1 0
             0 0 1
             1 1 0
             1 1 1
             0 0.3 0
             0 0 0
             0.65 0.65 0.65
             0 0.75 0.75
             0.3 0 0.6
             0.95 0.65 0.75
             0.5 0.25 0
             0 0.2 0.4
             0.9 0.4 0.7
             0.4 0.2 0.3
             0.65 0.55 0.15
             0.25 0.35 0.25
             0.5 0 0];
  for i = 1:numBalls
    balls(i).X = X(i, :)';
    balls(i).V = V(i, :)';
    balls(i).color = ballColor(mod(i-1,length(ballColor))+1,:)';
    balls(i).r = r(i);
  end
end
for i=1:numBalls;
  balls(i).handle = plot(balls(i).X(1), balls(i).X(2), '.', ...
                        'MarkerEdgeColor', balls(i).color', ...
                        'MarkerFaceColor', balls(i).color', ...
                        'erasemode', 'xor');
  balls(i).vhandle = line(balls(i).X(1) + [0 balls(i).V(1)]*0.1, ...
                          balls(i).X(2) + [0 balls(i).V(2)]*0.1, ...
                          'MarkerEdgeColor', balls(i).color', ...
                          'MarkerFaceColor', balls(i).color', ...
                          'erasemode', 'xor', 'visible', 'on', ...
                          'linewidth', 3, 'color', balls(i).color');
  if SHOWARROWS
    set(balls(i).vhandle, 'visible', 'on')
  else
    set(balls(i).vhandle, 'visible', 'off')
  end
end
set(tableAx, 'xlim', xlim)
set(tableAx, 'ylim', ylim)
drawnow;
