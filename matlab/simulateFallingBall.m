function ball = simulateFallingBall

global PAUSEDEMO
global RUNDEMO
global ball
global visStruct

xlim = get(visStruct.screenAx, 'xlim');
ylim = get(visStruct.screenAx, 'ylim');

cent = visStruct.holeCenter;
w = visStruct.holeWidth;

RUNDEMO = true;
dt = 0.01;
ball.score = 0;
count = 0;
while RUNDEMO
  set(ball.handle, 'xdata', ball.x(1), 'ydata', ball.x(2));
  pause(0.01)
  if PAUSEDEMO
  else
    count = count + 1;
    oldx = ball.x;
    ball.x = ball.x+ball.v*dt;
    if ball.x(1)+ball.r>xlim(2) || ball.x(1)-ball.r<xlim(1)
        ball.v(1) = -ball.v(1);
        ball.x = oldx;
    end
    ball.v(1) = randn(1)*ball.vVar + ball.v(1)*0.95;
    if(strcmp(get(visStruct.leftClick, 'visible'), 'on'))
        ball.v(1) = ball.v(1)-0.5;
        ball.energyCount = ball.energyCount+1;
        set(visStruct.energyCountTxt, 'string', num2str(ball.energyCount));
    end
    if(strcmp(get(visStruct.rightClick, 'visible'), 'on'))
        ball.v(1) = ball.v(1)+0.5;
        ball.energyCount = ball.energyCount+1;
        set(visStruct.energyCountTxt, 'string', num2str(ball.energyCount));
    end
    if ball.x(2)-ball.r<0 && ...
            ball.x(2)>0 && ...
            (ball.x(1)-ball.r < -cent-w/2 || ...
             (ball.x(1)+ball.r > -cent+w/2 && ...
              ball.x(1)-ball.r < cent-w/2) || ...
             ball.x(1)+ball.r > cent+w/2)
        set(visStruct.bangTxt, 'visible', 'on', 'position', [ball.x(1) ...
                          ball.x(2)]);
        set(ball.handle, 'visible', 'off')
        pause(0.5)
        set(ball.handle, 'visible', 'on')
        set(visStruct.bangTxt, 'visible', 'off');
        ball.x = [0 10];
        ball.v = [0 -1];
        continue
    end
    if ball.x(2)-ball.r<0 
        if (ball.x(1)-ball.r < -cent-w/2 || ...
                              (ball.x(1)+ball.r > -cent+w/2 && ...
                               ball.x(1)-ball.r < cent-w/2) || ...
                              ball.x(1)+ball.r > cent+w/2)
            ball.v(1) = -ball.v(1);
            continue
        end
    end
    if ball.x(2)+ball.r<ylim(1)
        ball.score = ball.score + 1;
        ball.x = [0 10];
        ball.v = [0 -1];
        set(visStruct.scoreTxt, 'string', num2str(ball.score));
        continue

    end
    avVal = ball.energyCount/ball.score;
    if isinf(avVal) || isnan(avVal)
        avTxt = '-';
    else
       avTxt = num2str(avVal, 3);
    end
    set(visStruct.averageTxt, 'string', avTxt);

  end
end