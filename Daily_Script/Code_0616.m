% % Test movenoisefile function
% 
% addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))
% getInf.moveNoiseFiles("C:\Users\Zhehao\Music\Blue732", "C:\Users\Zhehao\Music\B732_Noise")


% 写入并分析Neurons

addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))

[neuroster,malfunctions,cdrfroster,noSinChan_neuroster,nofeature_neuroster] = Archon.extractAnalysisInfo('D:');

 Archon.batch_createSinChFolder(noSinChan_neuroster); % merge single-channel folder
 
 unilist = AutoSap.script_autorun_scp(neuroster); % % 通过sap提取特征
 AutoSap.write_sapscreenshot(unilist);
 
 AutoSap.export_infofiles(1); % 导出sap信息
 AutoSap.export_datafiles(1);
 

error_neuroster = Archon.batch_genAnalysis(neuroster); % run Analysis for all neuroster
sul = Sultan("C:\Users\Zhehao\Downloads\Neuron0617");
sul.allSongallNeurons_FolderOrder;
sul.allSongallNeurons;

Sultan.applySameFunctionForAll("C:\Users\Zhehao\Downloads\Neuron0617");


sul2 = Sultan("E:\Neuron0617");
sul2.allSongallNeurons_FolderOrder;
sul2.allSongallNeurons;


sul3 = Sultan("C:\Users\Zhehao\Downloads\Test");
sul3.allSongallNeurons;


% 0619

sul = Sultan("E:\TestNeurons");
sul.allSongallNeurons

degexistinf = neuroninf([neuroninf.degexist] == 1)


%0619-0620
sul = Sultan("E:\GoodNeurons");
sul.What_percent_are_neurons_tested_with_different_stimuli_sets

all_exist_ids = mintersect(find([neuroninf.degexist]== 1),find([neuroninf.fragexist]== 1),find([neuroninf.replaexist]== 1));

alle_inf = neuroninf(all_exist_ids);
Sultan.How_Neurons_Resp_To_DegsFragsReplas(alle_inf);


sul.Whether_neurons_from_same_birds_respond_biasedly_to_songs
sul.How_Do_NCM_Neurons_respond_to_Songs;
sul.How_Could_NCM_Neurons_Be_Splitted_By_WL_And_FR;


nssul = Sultan("E:\NS_Neurons");
nssul.How_Do_NCM_Neurons_respond_to_Songs;

bssul = Sultan("E:\GoodNeurons");
bssul.How_Do_NCM_Neurons_respond_to_Songs;

figure;
pie([109,27+11,9,2,25],...
    {'CONs only','CONs/Degs','CONs/Degs/elements','CONs/Degs/replacements','All'});



new_rnames = {conallneurons.neuronname}.';
concat_respmap = horzcat(new_con_respmap,new_spe_respmap);
bfig = figure; imagesc(concat_respmap);
map = [ 1 1 1
    0 0.4470 0.7410
    0.8500 0.3250 0.0980];

colormap(bfig,map);
%colorbar;
xticks(1:size(concat_respmap,2));
xticklabels(horzcat(new_cnames,spekeywords));
%yticks(1:size(concat_respmap,1));
%yticklabels(new_rnames);
set(gca,'YTickLabel',[]);
set(gca,'TickLabelInterpreter','none');
set(gca,'TickLength',[0.001, 0.001]);
xtickangle(45)
            
bssul.Whether_neurons_from_same_birds_respond_biasedly_to_songs
            
toc


% to test evaluateConResp_FiringRate version

