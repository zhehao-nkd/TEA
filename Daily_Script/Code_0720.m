% for creating figure of 2022 poster

% update BS_NS_Neurons
dbstop if error
NBsul = Sultan("E:\BS_NS_Neurons");
NBsul.How_Do_NCM_Neurons_respond_to_Songs;

load('O706_Z03_SPKC07_1.mat')
A.drawStimuliOnly(17,[1.95,3.40])

A.drawSelectedStimuliResp([17,25,26,27,43,44,45,101,105,118],[1.92,2.62]);

colormap(newrespmat)
% A.draw_SelectedStimuli([41,132,25,26,43,44,45,101,105])


load('R693_Z13_SPKC11_1.mat')
A.drawStimuliOnly(3,[1.95,3.8])
colormap(newrespmat)

A.drawSelectedStimuliResp([3,35,36,37,59,60,61,148,147,146],[2.98,3.68]);
colormap(newrespmat)


load('O709_P07_SPKC13_1.mat')
A.drawStimuliOnly(2,[1.95,3.85])

A.drawSelectedStimuliResp([2,132,133,134,60,61,62,50,42,27],[1.98,2.68]);
%A.drawSelectedStimuliResp([2],[1.9,2.6]);
colormap(newrespmat)



% select those sample neurons
load('O709_P15_SPKC14_1.mat')
 %twoForPoster
Draw.twoForPoster(A.list(12).plty,A.list(12).pltsptimes,32000);

load('Y675_Z24_SPKC15_1.mat')
 %twoForPoster
Draw.twoForPoster(A.list(1).plty,A.list(1).pltsptimes,32000);

load('R693_Z12_SPKC13_1.mat')
 %twoForPoster
Draw.twoForPoster(A.list(15).plty,A.list(15).pltsptimes,32000);

load('R677_Z03_SPKC09_1.mat')
 %twoForPoster
Draw.twoForPoster(A.list(10).plty,A.list(10).pltsptimes,32000);


% 筛选NSBS neurons
allsul = Sultan("E:\BS_NS_Neurons");
allsul.How_Could_NCM_Neurons_Be_Splitted_By_WL_And_FR
allsul.How_Do_NCM_Neurons_respond_to_Songs;

% 分析 BS neurons
bssul = Sultan("E:\BS_Neurons");
bssul.How_Do_NCM_Neurons_respond_to_Songs;
bssul.How_Could_NCM_Neurons_Be_Splitted_By_WL_And_FR

% Generate resp map for NS neurons "E:\NS_Neurons"
nssul = Sultan("E:\NS_Neurons");
nssul.How_Do_NCM_Neurons_respond_to_Songs;


nsbs = Sultan("E:\BS_NS_Neurons");
%nsbs.How_Do_NCM_Neurons_respond_to_Songs;
nsbs.How_Could_NCM_Neurons_Be_Splitted_By_WL_And_FR
bstruct = nsbs.Whether_neurons_from_same_birds_respond_biasedly_to_songs;

% Generate resp map for neurons of O706,但可能要去掉WNS，calls之流
sulo706 = Sultan("E:\NeuronsO706");
sulo706.How_Do_NCM_Neurons_respond_to_Songs; % Version: judgeConResp_FR

 sulo706.Whether_neurons_from_same_birds_respond_biasedly_to_songs

 testSul = Sultan("E:\TestNeurons");
 testSul.How_Do_NCM_Neurons_respond_to_Songs;
 
 
 
 % 3:26 CDF analysis
 atypes = Sultan("E:\18ANeurons");
 
 cdf_collect = atypes.Which_acoustic_feature_matters;
 
  figure('Position',[2015 490 393 407],'Color','w');  hold on;
            
      
            cdfplot(cdf_collect(11).dists1_good );
            cdfplot(cdf_collect(11).dists0_good);
            %set(gca,'LineWidth',3)
            %legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel('Relative Distance','FontSize',16.5);
            ylabel('Cumulative Percentage','FontSize',16.5);
            title('')
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
           % title(sprintf('%s P-value : %.8f',a.formated_imagename,p),'interpreter', 'none');
            
            
 
%  load('O709_P13_SPKC05_1.mat')
%  A.draw_SelectedStimuli([4,31,32,33,88,89,90,44,49,63],[2.7,3.5]);
%  
%  load('R695_Z01_SPKC06_1.mat')
%  A.draw_SelectedStimuli([6,36,37,38,103,104,105,221,245,251],[3.3,4]);
%  
% 
%  
%  b = Batch("C:\Users\Zhehao\Downloads\O709_P07F3.txt",...
%      "D:\Ephys-O709-W\P07\O709_P07F3.plx","D:\Ephys-O709-W\P07\Stimuli\O709-ReplasB51203-P7F3");
%  b = Batch("C:\Users\Zhehao\Downloads\O709_P07F5.txt",...
%      "D:\Ephys-O709-W\P07\O709_P07F5.plx","D:\Ephys-O709-W\P07\Stimuli\O709-DegsB512-P7F5");
%  
%   

% 6.23 scatter analysis

atypes.How_Do_Resp_Eliciting_Elements_Distributed_In_Space(con_allspe_inf)
 

 
 
 b.select;
 
 neuronlist = b.getn;
 
 for k = 1:length(neuronlist)
     neuronlist{k}.rawthree
 end
 
 
 
 % determine position of the axes
axp = get(gca,'Position');

% determine startpoint and endpoint for the arrows 
xs=axp(1);
xe=axp(1)+axp(3)+0.04;
ys=axp(2);
ye=axp(2)+axp(4)+0.05;

% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [xs xs],[ys ye]);

% remove old box and axes
box off

set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'YColor',get(gca,'Color'))
set(gca,'XColor',get(gca,'Color'))
grid on
 
 title('')
 xlabel('Dim1')
 ylabel('Dim2')
 
 
 
 
 
 Sultan.runBatch("E:\18ANeurons")
 
 
 %%% 2022.0627
 sul = Sultan("E:\BS_Neurons");
 bstruct = sul.Whether_neurons_from_same_birds_respond_biasedly_to_songs;
 
 