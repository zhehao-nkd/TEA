% mannually segment syllable/note
dbstop if error


wavfolder = "C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\20210611selected20songs";
folders = Extract.filename(wavfolder,'*.wav');
syl = struct;
num = 0;
    
global NOT_PRESSED


for nnn = 1: length(folders)
    
    NOT_PRESSED = 1;
    [~,name,~] = fileparts(folders{nnn});
    [y,fs] = audioread(folders{nnn});
    
    figure('Position',[1 41 1920 1083]); Draw.spec2(y,fs);
    
    roi = drawline;
    temp = roi.Position(:,1);
    motify = y(temp(1)*fs:temp(2)*fs);
    close(gcf)
    figure('Position',[1 41 1920 1083]); Draw.spec2(motify,fs); colormap('jet')
    
    set(gcf, 'KeyPressFcn', @myKeyPressFcn);
    
    
    
    
 
    numinsong = 0; % number in the song
    while true
        if NOT_PRESSED == 0
            break;
        end
        %drawnow()
        roi = drawline;
        num = num + 1;
        numinsong = numinsong + 1;
        if NOT_PRESSED == 1
            temp = roi.Position(:,1);
            syl(num).initial = temp(1);
            syl(num).terminal = temp(2);
            syl(num).motify = motify;
            syl(num).segmenty = motify(temp(1)*fs:temp(2)*fs);
            syl(num).name = name;
            syl(num).numinsong =  numinsong;
        end
       
     
    end
    close(gcf)
end


function myKeyPressFcn(hObject, event)
    global NOT_PRESSED
    NOT_PRESSED  = 0;
    disp('key is pressed')
    
end