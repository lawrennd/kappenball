% Set color of billiard table.

numBalls=20;
tableColour = [.4 .5 .8];

close all;
tableAx = subplot(1, 2, 1);
hold on;
axis equal
axis([-1 1 -1 1]);
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



balls = initializeBilliards(numBalls, 'uniform', tableAx);
V = [balls(:).V]';
V = V(:);
centres = -10:10;
[barVals, centres] = hist(V, centres);
subplot(1, 2, 2);
barHandle = bar(centres, barVals);
set(gca, 'ylim', [0 numBalls]);
simulateBilliards(balls, barHandle);