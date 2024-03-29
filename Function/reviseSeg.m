classdef reviseSeg
    %功能说明： revise and check the segmentations of songs
    
    methods
        
        function rs = reviseSeg(dir_or_path)
            
            if exist(dir_or_path,'dir')
                dirpath = fullfile(dir_or_path,'SegData');%"E:\Stimuli_Source\senatusTwoMotif\SegData";
                segfiles = Extract.filename(dirpath,'*.mat');
                for k = 1: length(segfiles)
                    reviseSeg.Core(segfiles{k});
                    uiwait(gcf);
                end
                
            elseif isfile(dir_or_path)
                
                reviseSeg.Core(dir_or_path);
                
            end
            
        end
        
    end
    
    
    methods(Static)
        
        function Core(path)
            
            
            load(path,'segdata')
            hFig = figure('units','normalized','outerposition',[0 0 1 1]);
            %set(hFig,'WindowState','fullscreen');
            % to draw the title
            hOscillo = subplot(5,1,1);
%             set(subplot(2,2,1),'Color','k')
            hpy = highpass(segdata.rawy,450,32000);
            plot((1:length(hpy)).'/32000,hpy,'Color','cyan');
            set(gca,'Color','k');
            xlim([0 length(hpy)/32000]);

            hAxes = subplot(5,1,2:5);
            linkaxes([hOscillo,hAxes],'x');
         
            fs = 32000;
            y = segdata.rawy;
            [S,F,T] = spectrogram(hpy,hamming(524),524-round(fs/1e3),524,fs);
            imagesc(T,F,log(1+abs(S)));
            set(gca,'YDir','normal');
            % imshow(segdata.I);
            
            cmap = PM.cmap1;
            colormap(cmap);
            title(segdata.birdid);

            
            if isfield(segdata,'eleedge') % draw element segmentations
                for idx = 1: length(segdata.eleedge)
                    %hold on
                    hVerticalLines(idx) = line([segdata.eleedge(idx),segdata.eleedge(idx)],get(hAxes,'Ylim')*1.27,'Color','white');
                end
                if length(segdata.eleedge) == 0
                    hVerticalLines = [];
                end
            else
                hVerticalLines = [];
            end
            
            for e = 1: length(segdata.syledge) % draw syllables segmentations
                hVerticalLines = [hVerticalLines line([segdata.syledge(e),segdata.syledge(e)],get(hAxes,'Ylim')*1.27,'Color','red')];
            end
            
            if isfield(segdata,'motedge') %draw motif segmentations
                for idx = 1: length(segdata.motedge)
                    %hold on
                    hVerticalLines = [hVerticalLines, line([segdata.motedge(idx),segdata.motedge(idx)],get(hAxes,'Ylim')*1.27,'Color','green')];
                end
            end

           set(gca,'Clipping','Off')
          
            
            % legend('element','syllables')
            xlabel([...
                '{\color[rgb]{0.5 .5 .5}white--element \color{red}red--syllable} Press Space to change color']);
            
            ylabel( ' Right Click--Make Left Click--Move Double Click--Delete ');
            set(hAxes,'Color','k');
            
            set(hFig,'WindowButtonDownFcn',  @mouseDown);
            set(hFig,'WindowButtonMotionFcn',@mouseMove);
            set(hFig,'WindowButtonUpFcn',    @mouseUp);
            set(hFig,'CloseRequestFcn',    @winClose);
            set(hFig,'KeyPressFcn',    @keyPress);
            
            hFig.UserData = [ 1 1 1] %this is the color to draw lines
            
            hLineToDrag = [];
            hLineToDelete = [];


            
            %hLineToCreate = [];
            
            
            function mouseDown(hObject,~)
                
                clickType = get(hObject, 'SelectionType')
                
                
                switch  clickType
                    case 'normal' % left to move
                        disp('Moved!')
                        % is the mouse down event within the axes?
                        if IsCursorInControl(hObject, hAxes)
                            currentPoint   = get(hAxes,'CurrentPoint');
                            xCurrentPoint  = currentPoint(2,1);
                            for k=1:length(hVerticalLines)
                                xVertLineCoord = get(hVerticalLines(k),'XData');
                                if abs(xCurrentPoint - xVertLineCoord(1)) < 0.005
                                    hLineToDrag = hVerticalLines(k);
                                    break;
                                end
                            end
                        end
                    case 'open' % double to delete
                        % is the mouse down event within the axes?
                        disp('Deleted')
                        if IsCursorInControl(hObject, hAxes)
                            
                            currentPoint   = get(hAxes,'CurrentPoint');
                            xCurrentPoint  = currentPoint(2,1);
                            
                            for k=1:length(hVerticalLines)
                                
                                xVertLineCoord = get(hVerticalLines(k),'XData');
                                
                                if abs(xCurrentPoint - xVertLineCoord(1)) < 0.006
                                    hLineToDelete = hVerticalLines(k);
                                    % is the mouse down event within the axes?
                                    if ~isempty(hLineToDelete) && IsCursorInControl(hObject, hAxes)
                                        currentPoint = get(hAxes,'CurrentPoint');
                                        delete(hVerticalLines(k)); % firstly delete the object
                                        hVerticalLines(k) = [];% then remove the object handle
                                    end
                                    break;
                                end
                            end
                            
                        end
                        
                        
                        
                        
                    case 'alt' % right to create
                        
                        disp('Created');
                        if IsCursorInControl(hObject, hAxes)
                            currentPoint   = get(hAxes,'CurrentPoint');
                            xCurrentPoint  = currentPoint(2,1);
                        end
                        if exist('xCurrentPoint','var')
                            hVerticalLines = [hVerticalLines line([xCurrentPoint,xCurrentPoint],get(hAxes,'Ylim')*1.26,'Color',hObject.UserData)];
                            disp(['颜色是',string(hObject.UserData)])
                        end
                        
                end
                
                
                
                
            end
            
            function mouseUp(~,~)
                
                hLineToDrag = [];
                
            end
            
            function mouseMove(hObject,~)
                
                
                
                % is the mouse down event within the axes?
                if ~isempty(hLineToDrag) && IsCursorInControl(hObject, hAxes)
                    currentPoint = get(hAxes,'CurrentPoint');
                    x            = currentPoint(2,1);
                    set(hLineToDrag, 'XData', [x x]);
                end
            end
            
            function [status] = IsCursorInControl(hCursor, hControl)
                
                status = false;
                
                % get the position of the mouse
                figCurrentPoint = get(hCursor, 'CurrentPoint');
                position      = get(hCursor, 'Position');
                xCursor       = figCurrentPoint(1,1)/position(1,3); % normalize
                yCursor       = figCurrentPoint(1,2)/position(1,4); % normalize
                
                % get the position of the axes within the GUI
                controlPos = get(hControl,'Position');
                minx    = controlPos(1);
                miny    = controlPos(2);
                maxx    = minx + controlPos(3);
                maxy    = miny + controlPos(4);
                
                % is the mouse down event within the axes?
                if xCursor >= minx && xCursor <= maxx && yCursor >= miny && yCursor <= maxy
                    status = true;
                    disp('Yes')
                end
                
            end
            
            function winClose(~,~)
                
                fprintf('Current Waveview %s is updated \n',segdata.birdid)
                
                % newsegdata = [];
                white = 0;
                red = 0;
                green = 0;
                newsegdata.syledge = [];
                newsegdata.eleedge = [];
                newsegdata.motedge = [];
                for m = 1:length(hVerticalLines)
                    if hVerticalLines(m).Color == [1 1 1]
                        white = white + 1;
                        newsegdata.eleedge(white) = hVerticalLines(m).XData(1);
                    elseif hVerticalLines(m).Color == [1 0 0]
                        red = red + 1;
                        newsegdata.syledge(red) = hVerticalLines(m).XData(1);
                    elseif hVerticalLines(m).Color == [0 1 0]
                        green = green + 1;
                        newsegdata.motedge(green) = hVerticalLines(m).XData(1);
                    end
                end
                
                if isfield(newsegdata,'eleedge')
                    if ~isempty(newsegdata.eleedge) % if ele is a field and is not empty
                        newsegdata.eleedge = sort(newsegdata.eleedge);
                    end
                end

                if isfield(newsegdata,'motedge') % motif
                    if ~isempty(newsegdata.motedge) % if ele is a field and is not empty
                        newsegdata.motedge =sort(newsegdata.motedge);
                    end
                end

                if isfield(segdata,'syledge')
                    newsegdata.syledge =sort(newsegdata.syledge);
                end

                if isfield(segdata,'y')
                    newsegdata.y = segdata.y;
                end
                
                if isfield(segdata,'rawy')
                    newsegdata.rawy = segdata.rawy;
                end
                
                if isfield(segdata,'birdid')
                    newsegdata.birdid = segdata.birdid;
                end
                
                if isfield(segdata,'I')
                    newsegdata.I = segdata.I;
                end
                
                segdata = newsegdata;
                if exist('path','var')
                    save(path,'segdata')
                end
                delete(gcf)
            end
            
            function keyPress(hObject,event)
                if strcmp(event.Key, 'space') % Green-Red-White-Green
                    if hObject.UserData  == [ 0 1 0]
                        hObject.UserData = [1 0 0];
                    elseif hObject.UserData == [1 0 0]
                        hObject.UserData = [1 1 1];
                    elseif hObject.UserData == [1 1 1]
                        hObject.UserData = [0 1 0];
                    end
                end
                disp(['颜色是',string(hObject.UserData)])
            end
            
            
        end
        
    end
    
    
end


