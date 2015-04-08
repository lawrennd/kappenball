function changeVar

global ball

a = gco;
vStd = get(a, 'value');
ball.vVar = vStd*vStd;
set(a, 'Enable', 'off');
drawnow;
set(a, 'Enable', 'on');