function unitInPts = getAxisUnitInPts(ax);

fig = get(ax, 'parent');

units = get(fig, 'units');
set(fig, 'units', 'points');
pos = get(fig, 'position');
width = pos(3);
height = pos(4);
set(fig, 'units', units);
pos = get(ax, 'position');
width = width*pos(3);
height = height*pos(4);
xlim = get(ax, 'xlim');
ylim = get(ax, 'ylim');
unitInPts1 = width/(xlim(2)-xlim(1));
unitInPts2 = height/(ylim(2)-ylim(1));

unitInPts = min([unitInPts1 unitInPts2]);