classdef ImageCache
        % Saved all the necessary images of the neuron

    properties
        responseimage
        waveimage
        nameimage
    end

    methods
        function ic = ImageCache(sg,waveform)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
                ic.responseimage = struct;

             
                for k = 1: length(sg.normlist)

                    ic.responseimage(k).stimuliname = sg.normlist(k).stimuliname;
                    ic.responseimage(k).type =  sg.normlist(k).type;
                    ic.responseimage(k).name =  sg.normlist(k).name;

                    figure('Color','w','Position', [1975 615 960 227]);
                    Draw.spec(sg.normlist(k).plty,sg.normlist(k).fs);
                    yticks([0,8000,16000])
                    yticklabels({'0','8k','16k'});
                    xlabel(sg.normlist(k).stimuliname,'Interpreter','none');
                    
                    frame = getframe(gcf);
                    set(gca,'xticklabel',{[]});
                    close(gcf);
                    ic.responseimage(k).spec = frame.cdata;

                    figure('Color','w','Position', [1975 615 960 227]);
                    Draw.raster(sg.normlist(k).pltsptimes,sg.normlist(k).plty,sg.normlist(k).fs);
                    yticks([0,10])
                    yticklabels({'0','10'});
                    set(gca,'TickLength',[0 .01]);
                    frame = getframe(gcf);
                    close(gcf);
                    ic.responseimage(k).raster = frame.cdata;

                    figure('Color','w','Position', [1975 615 960 227]);
                    Draw.sdf(sg.normlist(k).plty,sg.normlist(k).fs,sg.normlist(k).pltsptimes);
                    box off
                    frame = getframe(gcf);
                    close(gcf);
                    ic.responseimage(k).sdf = frame.cdata;


                end
            

             
                figure('Color','w','Position', [1975 615 960 3*227]);
                waveform.drawAll;
                frame = getframe(gcf);
                close(gcf);
                ic.waveimage = frame.cdata;

                figure('Color','w','Position', [1975 615 960 3*227]);
                Draw.text(sg.formated_name);
                frame = getframe(gcf);
                close(gcf);
                ic.nameimage = frame.cdata;




        end

        function  drawThree(ic,formated_name,outdir)

          % images =  ic.responseimage;
          if ~exist('outdir','var')
              outdir = '.';
          end

       
            I = {};
            for idx = 1: length(ic.responseimage)

                first = ic.responseimage(idx).spec;
                second = ic.responseimage(idx).raster;
                third = ic.responseimage(idx).sdf;
          
                temp = {first,second,third};
                I{idx} = vertcat(temp{:});
                close(gcf);
            end

%             figure('Color','w','Position', [1933 672 673 497]);
%             neu.drawFirstWaveform;     % draw waveform
%             frame = getframe(gcf);
%             I{length(I)+ 1} = frame.cdata;
%             close(gcf);

            % draw blank white
            lieshu = 3;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));

            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end

            reshapedI = reshape(I, lieshu,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('%s\\鸟歌NormThree_%s.png',outdir,formated_name));


        end
    end
end