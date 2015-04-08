function fallPressKey(src,evnt)

global ball
global visStruct
switch evnt.Key
  case 'leftarrow'
    set(visStruct.leftClick, 'visible', 'on');
  case 'rightarrow'
    set(visStruct.rightClick, 'visible', 'on');

end
