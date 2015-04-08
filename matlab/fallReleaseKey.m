function fallReleaseKey(src,evnt)

global ball
global visStruct

switch evnt.Key
  case 'leftarrow'
    set(visStruct.leftClick, 'visible', 'off');
  case 'rightarrow'
    set(visStruct.rightClick, 'visible', 'off');

end
