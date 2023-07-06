classdef AutoGui < handle
    % automatic gui robot
    % 自动点击器
    
    properties
        rob
        rectangle
        leftbutton
        writelist
    end
    
    methods
        function a = AutoGui
            import java.awt.*;
            import java.awt.event.*;
            a.rob = Robot;
            a.rectangle = Rectangle;
            a.leftbutton = java.awt.event.InputEvent.BUTTON1_MASK;
            a.genWritelist;
        end
        function [x,y] = findText(a,image, text, roi)
            if exist('roi','var')
                ocrResults = ocr(image,roi);
            else
                ocrResults = ocr(image);
            end
            bboxes = locateText(ocrResults,text,IgnoreCase=true);

            %Iocr = insertShape(image,"FilledRectangle",bboxes);figure;imshow(Iocr)
            %

            if isempty(bboxes)
                x = nan;
                y = nan;
                return
            end

            x = bboxes(1)+bboxes(3)/2;
            y =  bboxes(2)+ bboxes(4)/2;

        end

        function deletion(a)

            % Simulate pressing the Delete key
            a.rob.keyPress(java.awt.event.KeyEvent.VK_DELETE);
            a.rob.keyRelease(java.awt.event.KeyEvent.VK_DELETE);

        end

        function ctrlA(a)

            % Create a Java Robot object
         

            % Simulate pressing and holding the Ctrl key
            a.rob.keyPress(java.awt.event.KeyEvent.VK_CONTROL);

            % Simulate pressing the 'A' key
            a.rob.keyPress(java.awt.event.KeyEvent.VK_A);
            a.rob.keyRelease(java.awt.event.KeyEvent.VK_A);

            % Simulate releasing the 'A' key
            a.rob.keyRelease(java.awt.event.KeyEvent.VK_CONTROL);

        end

        function moveScroll(a,scrollAmount,wait)

            if exist('wait','var')
                pause(wait);
            end


            % Specify the number of wheel movements (positive for scrolling up, negative for scrolling down)
            % scrollAmount = 30; % Change this value as needed

            % Simulate mouse wheel movement
            a.rob.mouseWheel(-scrollAmount);

        end
        
        function move(a,x,y)
            %pause(wait);
            a.rob.mouseMove(x,y);
            
        end

        function ctrlV(a)
       

            % Simulate pressing and holding the Ctrl key
            a.rob.keyPress(java.awt.event.KeyEvent.VK_CONTROL);

            % Simulate pressing the 'V' key
            a.rob.keyPress(java.awt.event.KeyEvent.VK_V);
            a.rob.keyRelease(java.awt.event.KeyEvent.VK_V);

            % Simulate releasing the 'V' key
            a.rob.keyRelease(java.awt.event.KeyEvent.VK_CONTROL);

        end
        
        function a = genWritelist(a)
            
            
            alphabet = 'a':'z';
            keylist1 = struct;
            for k = 1: length(alphabet)
                keylist1(k).key = alphabet(k);
                keylist1(k).id = 64 + k;
            end
            
            Alphabet = 'A':'Z';
            keylist2 = struct;
            for k = 1: length(Alphabet)
                keylist2(k).key = Alphabet(k);
                keylist2(k).id = [16,64 + k];
            end
            
            nums = '0':'9';
            keylist3 = struct;
            for k = 1: length(nums)
                keylist3(k).key = nums(k);
                keylist3(k).id = 47 + k;
            end
            
            keylist4 = struct;
            keylist4(1).key = '-';
            keylist4(1).id = 45;
            keylist4(2).key = '_';
            keylist4(2).id = [16,45];
            %             keylist4(3).key = "shift";
            %             keylist4(3).id = 16;
            
            a.writelist = vertcat([keylist1,keylist2,keylist3,keylist4]);
            
            
        end
        
        function a = pressEnter(a)
            
            a.rob.keyPress(10);
            a.rob.keyRelease(10);
            
        end
        
        function click(a,x,y,wait)
            if exist('wait','var')
                pause(wait);
            end
            a.rob.mouseMove(x,y);
            a.rob.mousePress(a.leftbutton);
            a.rob.mouseRelease(a.leftbutton);
        end
        
        function doubleclick(a,x,y,wait)
            if exist('wait','var')
                pause(wait);
            end
            a.rob.mouseMove(x,y);
            a.rob.mousePress(a.leftbutton);
            a.rob.mouseRelease(a.leftbutton);
            a.rob.mousePress(a.leftbutton);
            a.rob.mouseRelease(a.leftbutton);
        end

        function Deprecated_Typewrite(a,textToType)

           

            % Specify the text you want to type
            %textToType = 'Hello, World!';

            % Convert the text to an array of characters
            characters = char(textToType);

            % Simulate typing each character
            for i = 1:length(characters)
                character = characters(i);

                % Simulate typing the character
                a.rob.keyPress(character);
                a.rob.keyRelease(character);
            end


        end
        
        function typewrite(a,txt)
            dbstop if error
            chars = convertStringsToChars(txt);
            
            for k = 1:length(chars)
                id = a.char2id(chars(k));
                
                if isempty(id)
                    disp('这个按键尚未收录');
                    pause;
                elseif length(id) == 1
                    a.rob.keyPress(id);
                    a.rob.keyRelease(id);
                elseif length(id) == 2
                    a.rob.keyPress(id(1));
                    a.rob.keyPress(id(2));
                    a.rob.keyRelease(id(2));
                    a.rob.keyRelease(id(1));
                end
                pause(0.05);
            end
            
        end
        
        function drag(a,start,destiney,innertime)
            a.rob.mouseMove(start(1),start(2));
            a.rob.mousePress(a.leftbutton)
            if exist('innertime','var')
                pause(innertime);
            end
            a.rob.mouseMove(destiney(1),destiney(2));
            a.rob.mouseRelease(a.leftbutton);
        end
        
        function img = capture(a,lefttop,rightbottom) % screenshot
            a.rectangle.x = lefttop(1);
            a.rectangle.y = lefttop(2);
            a.rectangle.width = rightbottom(1) -lefttop(1);
            a.rectangle.height = rightbottom(2) -lefttop(2);
            image = a.rob.createScreenCapture(a.rectangle);
            data = image.getData();
            % 获取像素信息
            temp = zeros(a.rectangle.width*a.rectangle.height*3,1); % 初始化数组用于存储RGB像素信息（长*宽*通道数）
            temp = data.getPixels(0,0,a.rectangle.width,a.rectangle.height,temp);
            temp = uint8(temp);
            % 提取三通道像素值
            R = temp(1:3:end);
            G = temp(2:3:end);
            B = temp(3:3:end);
            % 修改数据尺寸
            R = reshape(R,[a.rectangle.width,a.rectangle.height]);
            G = reshape(G,[a.rectangle.width,a.rectangle.height]);
            B = reshape(B,[a.rectangle.width,a.rectangle.height]);
            R = R'; G = G'; B = B';% 转置
            img = cat(3,R,G,B); % 合并
        end
        
        function mcolor = meancolor(~,image) %  判断平均颜色
            meanR = mean(image(:,:,1),'all');
            meanG = mean(image(:,:,2),'all');
            meanB = mean(image(:,:,3),'all');
            mcolor = [meanR,meanG,meanB];
        end
        
        function id = char2id(a,inputchar)
            id = a.writelist(inputchar==[a.writelist.key].').id;
        end
        
        function [draglength,clicky] = dragClickInBox(~,firsty,endy,endnum, totalnum,dragrange,targetnum)
            % firsty: 选择框第一项的y坐标      endy：选择框最末项的y坐标      endnum：选择框最末项的次序
            % totalnum: 所有可选项的数目      dragrange：拖动条从起始到使得最末项位于endy位置时所行使的距离
            % targetnum; 目标项的次序数
            clickgap = (endy - firsty)/ (endnum -1);
            draggap = dragrange /(totalnum - endnum);
            
            if targetnum <= endnum %当目标项位于选择框以内时
                draglength = 0;
                clicky = firsty + (targetnum - 1)*clickgap;
            elseif targetnum > endnum % 当目标项位于选择框之外时
                draglength = (targetnum -endnum)*draggap;
                clicky = endy;
                
            end
            
        end
    end
    
    methods(Static)
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

