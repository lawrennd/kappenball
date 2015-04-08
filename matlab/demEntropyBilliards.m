% Set color of billiard table.

global PAUSEDEMO
global SHOWARROWS
global INSTANTVELOCITIES

INSTANTVELOCITIES = true;
PAUSEDEMO = true;
SHOWARROWS = false;
numBalls=20;
tableColour = [.4 .5 .8];

screenRatio = 5/6;
close all;
fh = figure('name', 'Entropy Billiards', 'NumberTitle', 'off');
set(fh, 'menubar', 'none')
ssize = get(0, 'screensize');
units = get(0, 'units');
ssize = [1 1 1024 768];

funits = get(gcf, 'units');
set(gcf, 'units', units)
set(gcf, 'position', ssize)
set(gcf, 'units', funits)
midPoint = 0.5*[ssize(3) ssize(4)];
ssize(4) = ssize(4)*screenRatio;
ssize(3) = ssize(3)*screenRatio;
ssize(1:2) = midPoint - 0.5*[ssize(3) ssize(4)];
set(gcf, 'position', ssize);

arrowButton = uicontrol('units', 'normalized', ...
                        'position', [0.12 0.1 0.06 0.05], ...                   
                        'string', 'Arrows', ...                        
                      'style', 'togglebutton', ...
                      'callback', 'toggleLines(balls)', ...
                        'value', SHOWARROWS);

pauseButton = uicontrol('units', 'normalized', ...
                        'position', [0.22 0.1 0.06 0.05], ...                   
                        'string', 'Pause', ...
                        'style', 'togglebutton', ...
                        'callback', 'toggleDemo', ...
                        'value', PAUSEDEMO);

velocitiesButton = uicontrol('units', 'normalized', ...                         
                             'position', [0.32 0.1 0.06 0.05], ...                   
                             'string', 'Instant v', ...
                             'style', 'togglebutton', ...
                             'callback', 'toggleVelocities', ...
                             'value', INSTANTVELOCITIES);

stopButton = uicontrol('units', 'normalized', ...
                       'position', [0.42 0.1 0.06 0.05], ...
                      'string', 'Exit', ...
                      'style', 'pushbutton', ...
                      'callback', 'toggleRun');

tableAx = axes('position', [0.05 0.2 0.5 0.7]);

hold on;
axis equal
axis([-2 2 -2 2]);
xlim = get(tableAx, 'xlim');
ylim = get(tableAx, 'ylim');
%axis off
%Plothandle
%set(gcf, 'color',[1 1 1]);
set(gca, 'color', tableColour, ...
         'xcolor', tableColour, ...
         'ycolor', tableColour, ...
         'xtick', [], ...
         'ytick', []);


balls = initializeBilliards(numBalls, 'hotCold', tableAx);
V = [balls(:).V]';
V = V(:);
centres = -20:2:20;
[barVals, centres] = hist(V, centres);
barAx = axes('position', [0.65 0.2 0.3 0.7]);
set(barAx, 'fontsize', 20, 'fontname', 'helvetica', 'xtick', [-40 ...
                    -20 -10 0 10 20 40]);

barHandle = bar(centres, barVals);
xlabel('velocity', 'fontsize', 20, 'fontname', 'helvetica');
ylabel('counts', 'fontsize', 20, 'fontname', 'helvetica');
set(gca, 'ylim', [0 numBalls]);
simulateBilliards(balls, barHandle);
close all 
clear all
