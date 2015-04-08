function toggleRun

global RUNDEMO

RUNDEMO = ~RUNDEMO;
a = gco;
set(a, 'Enable', 'off');
drawnow;
set(a, 'Enable', 'on');