function resetButton

global visStruct
global ball

ball.score = 0;
ball.energyCount = 0;
set(visStruct.scoreTxt, 'string', num2str(ball.score));
set(visStruct.energyCountTxt, 'string', num2str(ball.energyCount));
set(visStruct.averageTxt, 'string', num2str(0));

a = gco;
set(a, 'Enable', 'off');
drawnow;
set(a, 'Enable', 'on');