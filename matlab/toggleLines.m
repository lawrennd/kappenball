function toggleLines(balls)

global SHOWARROWS 
SHOWARROWS = ~SHOWARROWS;

if SHOWARROWS
  for k = 1:length(balls)
    set(balls(k).vhandle, 'visible', 'on')
  end
else
  for k = 1:length(balls)
    set(balls(k).vhandle, 'visible', 'off')
  end
end
