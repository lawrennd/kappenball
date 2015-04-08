function balls = vectorBilliards(balls, barHandle)%Plot startingpositions====================================================

  for k=1:numBalls;
    set(balls(k).vhandle, 'xdata', balls(k).X(1,1) + [0 balls(k).V(1, 1)], ...
                      'ydata', balls(k).X(k,2) + [0 balls(k).V(1, 2)], ...
                      'visible', 'true');
    
  end
end