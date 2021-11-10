% info = f1.list;
% 
% 
% dir = 'forsim';
% mkdir(dir);
% for k = 1: length(info)
%     name = sprintf('%s/%s.wav',dir,info(k).stimuliname);
%     audiowrite(name,info(k).y,info(k).fs);
% end

%scratch to draw figures for the JNS poster
% figure D
f8b = Fenxi(8);
f8b.drawselect([5,48,49,40,41],[0,1/6]);
f8b.drawselect([5,48,49,40,41],[0,1]);

f4 = Fenxi(4);
f4.drawselect([67,31,40,34,35],[0,1]);
f4.drawselect([67,31,40,34,35],[1/25,1/5]);

f3 = Fenxi(3);
f3.drawselect([40,41],[0,4/5]);



% Figure E
f16 = Fenxi(16);
f16.drawselect([53,54],[0,2/3])

f16.drawselect([21,53,54],[0,1])

f5 = Fenxi(5);
f5.drawselect([61,62],[0,1])
f5.drawselect([10,61,62],[0,1])

f1 = Fenxi(1);
f1.drawselect([25,29],[0,1/2])
f1.drawselect([11,25,29],[0,1])

%Q2  Figure C single element presentation
f1 = Fenxi(1);
f1.drawsyl('Y606',[0.329,0.62],[1,8,9,10,11,12]);
f1.drawsyl('Y606',[0.329,0.9],[1,8,9,10,11,12,13,14,15,16,17,18]);
f1.drawsyl('Y606',[0,1],[1,8,9,10,11,12,13,14,15,16,17,18]);
f1.drawsyl('Y606',[0,1],[1,9,11,17,19])
f1.drawsinglenorm('Y606')
axis off
%f1.drawselect([12,35,36,37,38,39],[1/4,1/2])
