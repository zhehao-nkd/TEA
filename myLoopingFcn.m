function myLoopingFcn() 
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf
set(gcf, 'KeyPressFcn', @myKeyPressFcn)
while ~KEY_IS_PRESSED
      drawnow
      disp('looping...')
end
disp('loop ended')
function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
disp('key is pressed')
de