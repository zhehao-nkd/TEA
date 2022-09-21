function oldReviseSeg(path)

%close all;
load(path,'segdata')

hFig = figure('units','normalized','outerposition',[0 0 1 1]);

%set(hFig,'WindowState','fullscreen');

 % to draw the title
hAxes = subplot(1,1,1);
imagesc(segdata.I);

cmap = [0,0,0;0.0896030038402175,0.177460786403949,0.772377257837302;1.02733538499169e-05,0.531906322143127,0.790204355661316;0,0.653257444386579,0.569214339167397;0,0.635337430754235,0.227075820665484;0,0.649944971405242,0.0135939618070814;0,0.711846281671014,3.23641129616309e-06;0,0.780331894859652,0;0,0.846265394673730,0;1.93417789036852e-06,0.891450356031797,0;0.000298602846034070,0.931186511640818,0;0.0161588836764975,0.968090832182706,0;0.110175676849356,0.991616318971059,0;0.308529040806655,0.998904633379897,0;0.467071748986709,0.999777741934953,0;0.550243862366399,0.999349631497273,0;0.619588362937050,0.998168008533641,0;0.679061809719645,0.995751484399177,0;0.727171433199057,0.991913669148313,0;0.762829197438895,0.986518494906395,0;0.790242462469082,0.980092163810161,0;0.813289168545624,0.973118946462708,0;0.835088612346230,0.965980510682342,0;0.856510285384010,0.958744526561342,0;0.877295944809658,0.951220196910804,0;0.900022649889677,0.941766143052440,0;0.913407557638204,0.934778640784699,0;0.923831035458572,0.927398110362299,0;0.932943174851906,0.919362402741546,0;0.940716028678719,0.910657554704271,0;0.947295987605561,0.901356706664461,0;0.953177109570859,0.891706597390989,0;0.958649974707923,0.881851926726321,0;0.963881308130102,0.871873098479008,0;0.968959886646116,0.861805705415255,0;0.974748837756767,0.850009050540441,0;0.978492923658549,0.842062854330428,0;0.982084007360490,0.833965480328865,0;0.985471618665228,0.825665271080110,0;0.988618110610656,0.817124130310414,0;0.991428566195205,0.808247004062456,0;0.993776943979006,0.798907812034241,0;0.995607086269905,0.789050386091939,0;0.996993329559934,0.778749061148768,0;0.998046012970266,0.768114158369140,0;0.998802186058257,0.757182484460732,0;0.999284474234472,0.745974929748419,0;0.999567501670414,0.734559337235648,0;0.999759305114816,0.723028040378702,0;0.999873973554784,0.711372246910081,0;0.999935254344718,0.699522755110929,0;0.999966940574306,0.687369885936895,0;0.999984137943231,0.674807502474124,0;0.999993282881041,0.661721557245968,0;0.999997583837261,0.647823004418224,0;0.999999241320172,0.632531384435626,0;0.999999755630015,0.615453376505742,0;0.999999914152097,0.596505395998072,0;0.999999975215456,0.575524375740835,0;0.999999992825045,0.552207180883159,0;0.999999985101992,0.526508197564670,0;0.999999928563746,0.498753746101448,0;0.999999761223866,0.469243586996188,0;0.999999204756809,0.438220102169859,0;0.999997317689850,0.405786541260666,0;0.999992234917945,0.372257950313365,0;0.999981069474533,0.338145980451264,0;0.999959234335268,0.303681152082064,0;0.999917152919158,0.269011116220679,0;0.999793166583269,0.234556031492901,0;0.999498852237756,0.200828019226696,0;0.998972998711474,0.168133955635982,0;0.998108989311878,0.136960434493596,0;0.996641131387946,0.108504742204792,0;0.994455445652611,0.0832800828699568,0;0.991507419001210,0.0614868167198256,0;0.987774157404586,0.0432279953666773,0;0.983203528873757,0.0287383041813776,0;0.977771350878884,0.0181267298746843,0;0.971614718210821,0.0107765332890833,0;0.964933131405141,0.00579031512733513,0;0.957824489571211,0.00272946408467142,0;0.950383376071615,0.00117481288078103,0;0.942754752994586,0.000489989657862349,0;0.935055103525107,0.000180185406281604,0;0.927337141894733,5.36617754521494e-05,0;0.919641995293775,1.30998401979646e-05,0;0.912031568949177,3.44158188748933e-06,0;0.904571070336757,1.13921463822851e-06,0;0.897294012301756,2.69016646004164e-07,0;0.890234257893037,4.79633680046001e-08,0;0.883500578334965,6.95260257170901e-09,3.57810024175010e-10;0.877165912151234,6.07954458691698e-09,5.19175329196685e-09;0.871217291169426,3.87237891461433e-08,3.85666292787921e-08;0.865620090633176,1.22647923331998e-07,1.22631235611760e-07;0.860402398614529,5.53186205928394e-07,5.53186205928394e-07;0.855575556366116,2.34495410446194e-06,2.34495410446194e-06;0.851062547084255,8.42714364288835e-06,8.42714364288835e-06;0.846771151778257,2.61285887560802e-05,2.61285887560802e-05;0.842641627788605,7.33811832591591e-05,7.33811832591591e-05;0.838637051459501,0.000195714811401956,0.000195714811401956;0.834718268414311,0.000493019009571090,0.000493019009571090;0.830853398327577,0.00113157824278559,0.00113157824278559;0.827038206505286,0.00223964007586136,0.00223964007586136;0.823295594836600,0.00414090790268332,0.00414090790268332;0.819697315131666,0.00774207783686642,0.00774207783686642;0.816279105291062,0.0134955646080183,0.0134955646080183;0.813064974304146,0.0216958279721834,0.0216958279721834;0.810105894039147,0.0329549376849062,0.0329549376849062;0.807493991402155,0.0483780425578325,0.0483780425578325;0.805306332296563,0.0688895256184690,0.0688895256184690;0.803537590437394,0.0944255033888212,0.0944255033888212;0.802145144572146,0.124474776438585,0.124474776438585;0.801142174218503,0.159195419942387,0.159195419942387;0.800568138156783,0.199060702629676,0.199060702629676;0.800275853085865,0.242305307545471,0.242305307545471;0.800119239815312,0.287177160441653,0.287177160441653;0.800042571108881,0.333007868679640,0.333007868679640;0.800011492563755,0.379385385449002,0.379385385449002;0.800001997178583,0.426021770703579,0.426021770703579;0.800000142114041,0.472749793986708,0.472749793986708;0.800000000000000,0.519498362403275,0.519498362403275;0.800000000000000,0.566248635336074,0.566248635336074;0.800000000000000,0.612998908268869,0.612998908268869;0.800000000000000,0.659749181201653,0.659749181201653;0.800000000000000,0.706499454134420,0.706499454134420;0.800000000000000,0.753249727067159,0.753249727067159;0.800000000000000,0.799999999999840,0.799999999999840];
colormap(cmap);
title(segdata.birdid);


if isfield(segdata,'eleedge')
    for idx = 1: length(segdata.eleedge)
        %hold on
        hVerticalLines(idx) = line([segdata.eleedge(idx),segdata.eleedge(idx)],get(hAxes,'Ylim'),'Color','white');
    end
    
    if length(segdata.eleedge) == 0
        hVerticalLines = [];
    end
    
else
    hVerticalLines = [];
end

for e = 1: length(segdata.syledge)
    
    hVerticalLines = [hVerticalLines line([segdata.syledge(e),segdata.syledge(e)],get(hAxes,'Ylim'),'Color','red')];
end

% legend('element','syllables')
xlabel([...
    '{\color[rgb]{0.5 .5 .5}white--element \color{red}red--syllable} Press Space to change color']);

ylabel( ' 右键--创建 左键--拖动 双击--删除 ');
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
                        if abs(xCurrentPoint - xVertLineCoord(1)) < 5
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
                
                disp('Created');
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
       % disp(status)
    end

    function winClose(~,~)
        
        fprintf('Current Waveview %s is updated \n',segdata.birdid)
        
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
        
        if isfield(newsegdata,'eleedge')
            if ~isempty(newsegdata.eleedge) % if ele is a field and is not empty
                sort(newsegdata.eleedge);
            end
        end
        
        if isfield(segdata,'syledge')
            sort(newsegdata.syledge);
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