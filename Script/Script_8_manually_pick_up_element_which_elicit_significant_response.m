% A script to manually label response-eliciting syllables;
global NOT_PRESSED
b = batch(pa.T);
candi = [5,21,26,53,57,58,61,36,39,40,53,58]';

close all

for idx = 1: length(candi)
    b.select(candi(idx));
    temp = b.getn;
    thisn = temp{1};
    %thisn.manpicksig;
   
    
    num = 0;
    syl = struct;
   
    
    for ww = 1: length(thisn.e)
        
        NOT_PRESSED = 1;
        
        %figure('Position',[1 41 1920 1083]); %draw.spec2(thisn.e{ww}.plty,fs);
        thisn.e{ww}.three;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf, 'KeyPressFcn', @myKeyPressFcn);
        
        while true
            if NOT_PRESSED == 0
                break;
            end
            %drawnow()
           % subplot(3,1,1)
           
            roi = drawline;
            
            g = groot;
            % True if there are no open graphics objects, false otherwise
            if ~isempty(g.Children)&& (NOT_PRESSED == 1)
                num = num + 1;
                temp = roi.Position(:,1);
                syl(num).initial = temp(1);
                syl(num).terminal = temp(2);
                
                syl(num).wholey =  thisn.e{ww}.plty;
                fs = thisn.e{ww}.sound.fs;
                syl(num).fragy =  thisn.e{ww}.plty(temp(1)*fs:temp(2)*fs);
                syl(num).name = thisn.e{ww}.sound.name;
                disp('MAGA')
            end
            
            
        end
        
        if ~isempty(g.Children)
            close all
        end
        
    end
    
    save('sylmanpick.mat','syl');
    figure; 
    sumy = []
    for oba = 1: length(syl)
        sumy = [sumy;syl(oba).fragy;zeros(0.1*32000,1)];
    end
    collect{idx} = syl;
    
    
end


function myKeyPressFcn(hObject, event)
    global NOT_PRESSED
    NOT_PRESSED  = 0;
    disp('key is pressed')
    
end




