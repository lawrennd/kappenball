function toggleDemo

global PAUSEDEMO

PAUSEDEMO = ~PAUSEDEMO;
a = gco;
set(a, 'Enable', 'off');
drawnow;
set(a, 'Enable', 'on');