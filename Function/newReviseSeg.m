function newReviseSeg(path)

%close all;

if exist(path)== 2 % when it is a file
    load(path,'segdata');
    hFig = figure;
    
    % to draw the title
    hAxes = subplot(1,1,1);
    imagesc(segdata.I);
    title(segdata.birdid);
    
    for idx = 1: length(segdata.eleedge)
        %hold on
        hVerticalLines(idx) = line([segdata.eleedge(idx),segdata.eleedge(idx)],get(hAxes,'Ylim'),'Color','white');
    end
    
    if length(segdata.eleedge) == 0
        hVerticalLines = [];    
    end
    
    for e = 1: length(segdata.syledge)
        hVerticalLines = [hVerticalLines line([segdata.syledge(e),segdata.syledge(e)],get(hAxes,'Ylim'),'Color','red')];
    end
    
    % legend('element','syllables')
    xlabel(['{\color[rgb]{0.5 .5 .5}white is element \color{red}red is syllable} Press Space to change color']);
    
    set(hAxes,'Color','k');
    set(hFig,'WindowButtonDownFcn',  @mouseDown);
    set(hFig,'WindowButtonMotionFcn',@mouseMove);
    set(hFig,'WindowButtonUpFcn',    @mouseUp);
    set(hFig,'CloseRequestFcn',    @winClose);
    set(hFig,'KeyPressFcn',    @keyPress);
    
    hFig.UserData = [ 1 1 1] %this is the color to draw lines
    
    hLineToDrag = [];
    hLineToDelete = [];

elseif exist(path)== 7 % when it is a folder
    filenames = extract.filename(path,'*.mat')
    
    for m = 1: length(filenames)
        load(filenames{m},'segdata');
        hFig{m} = figure;
        
        % to draw the title
        hAxes = subplot(1,1,1);
        imagesc(segdata.I);
        title(segdata.birdid);
        
        for idx = 1: length(segdata.eleedge)
            %hold on
            hVerticalLines(idx) = line([segdata.eleedge(idx),segdata.eleedge(idx)],get(hAxes,'Ylim'),'Color','white');
        end
        
        if length(segdata.eleedge) == 0
            hVerticalLines = [];
        end
        
        for e = 1: length(segdata.syledge)
            hVerticalLines = [hVerticalLines line([segdata.syledge(e),segdata.syledge(e)],get(hAxes,'Ylim'),'Color','red')];
        end
        
        % legend('element','syllables')
        xlabel(['{\color[rgb]{0.5 .5 .5}white is element \color{red}red is syllable} Press Space to change color']);
        set(hAxes,'Color','k');
        set(hFig{m},'WindowButtonDownFcn',  @mouseDown);
        set(hFig{m},'WindowButtonMotionFcn',@mouseMove);
        set(hFig{m},'WindowButtonUpFcn',    @mouseUp);
        set(hFig{m},'CloseRequestFcn',    @winClose);
        set(hFig{m},'KeyPressFcn',    @keyPress);
        
        hFig{m}.UserData = [ 1 1 1] %this is the color to draw lines
        hLineToDrag = [];
        hLineToDelete = [];
        
        
    end
end

    function mouseDown(hObject,~)
       
        clickType = get(hObject, 'SelectionType')
        
        
        switch  clickType
            case 'normal' % left tro move
                disp('why???')
                % is the mouse down event within the axes?
                if IsCursorInControl(hObject, hAxes)
                    
                    currentPoint   = get(hAxes,'CurrentPoint');
                    xCurrentPoint  = currentPoint(2,1);
                    
                    for k=1:length(hVerticalLines)
                        
                        xVertLineCoord = get(hVerticalLines(k),'XData');
                        
                        if abs(xCurrentPoint - xVertLineCoord(1)) < 5
                            hLineToDrag = hVerticalLines(k);
                            break;
                        end
                    end
                end
            case 'open' % double to delete
                % is the mouse down event within the axes?
                if IsCursorInControl(hObject, hAxes)
                    
                    currentPoint   = get(hAxes,'CurrentPoint');
                    xCurrentPoint  = currentPoint(2,1);
                    
                    for k=1:length(hVerticalLines)
                        
                        xVertLineCoord = get(hVerticalLines(k),'XData');
                        
                        if abs(xCurrentPoint - xVertLineCoord(1)) < 0.1
                            hLineToDelete = hVerticalLines(k);
                            % is the mouse down event within the axes?
                            if ~isempty(hLineToDelete) && IsCursorInControl(hObject, hAxes)
                                currentPoint = get(hAxes,'CurrentPoint');
                                delete(hVerticalLines(k)) % firstly delete the object
                                hVerticalLines(k) = [] % then remove the object handle
                            end
                            break;
                        end
                    end
                    
                end
                
                
                
                
            case 'alt' % right to create
                
                if IsCursorInControl(hObject, hAxes)
                    
                    currentPoint   = get(hAxes,'CurrentPoint');
                    xCurrentPoint  = currentPoint(2,1);
                    
                    
                end
                
                hVerticalLines = [hVerticalLines line([xCurrentPoint,xCurrentPoint],get(hAxes,'Ylim'),'Color',hObject.UserData)];
                disp(['颜色是',string(hObject.UserData)])
                
                
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
            
        end
        disp(status)
    end

    function winClose(~,~)
        disp('GG')
        
       % newsegdata = [];
        white = 0;
        red = 0;
        newsegdata.syledge = [];
        newsegdata.eleedgem = [];
        for m = 1:length(hVerticalLines)
            if hVerticalLines(m).Color == [1 1 1]
                white = white + 1;
            newsegdata.eleedge(white) = hVerticalLines(m).XData(1);
            elseif hVerticalLines(m).Color == [1 0 0]
                red = red + 1;
                 newsegdata.syledge(red) = hVerticalLines(m).XData(1);
            end
        end
        
        if isfield(segdata,'eleedge')
            sort(newsegdata.eleedge);
        end
        
        if isfield(segdata,'syledge')
            sort(newsegdata.syledge);
        end
        
        if isfield(segdata,'y')
            newsegdata.y = segdata.y;
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
       if strcmp(event.Key, 'space')
           if hObject.UserData  == [ 1 1 1]
               hObject.UserData = [1 0 0];
           elseif hObject.UserData == [1 0 0]
               hObject.UserData = [1 1 1];
           end
       end
       disp(['颜色是',string(hObject.UserData)])
    end


end