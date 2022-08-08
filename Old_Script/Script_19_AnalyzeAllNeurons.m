%$ script 19   To Analyze all the neurons
dbstop if error

%tb = readtable("C:\Users\Zhehao\Desktop\TogenerateNeurons.xlsx");
%tb = readtable("C:\Users\Zhehao\Desktop\NeuronsUseNewStimuli.xlsx");

tb = readtable("C:\Users\Zhehao\Desktop\Rank12_Neurons.xlsx");

rmtb = rmmissing(tb);

tbstruct = table2struct(rmtb);

neurondata_dir = 'NeuronData';
mkdir(neurondata_dir);

% 61 cannot be saved
%for k = 84:length(tbstruct)
for k = 8:length(tbstruct)    
    txt = tbstruct(k).TXT;
    plx = tbstruct(k).PLX;
    stidir = tbstruct(k).STIMULI;
    
    if ~isfile( txt)
        continue
    end
    
    b = Batch(txt, plx,stidir);
    
    b.select;
    
    neuronlist = b.getn;
    
    for m = 1: length(neuronlist)
        thisn = neuronlist{m};
        %thisn.rawthree;
        thisn.threesingle;
    end

    try
        for w = 1:length(neuronlist)
            neuronname = sprintf('%s\\Neuron_%s',neurondata_dir,neuronlist{w}.neuronname);
            neurondata = neuronlist{w};
            save(neuronname,'neurondata');  % save neurons to mat files
        end
    catch ME
    end

end

send_mail_message('chengzh.nkd@gmail.com','Done!','Good Luck')
send_mail_message('zhehao.cheng@oist.jp','Done!','Good Luck')
