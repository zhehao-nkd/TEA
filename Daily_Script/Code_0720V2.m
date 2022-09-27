%sdfcollect = {};
sticollect = {};
sptimescollect = {};
len = 0;
for k = 1: length(A.list)
    
%sdfcollect{k} = A.list(k).sdf;
ycollect{k} = A.list(k).y;
sticollect{k} = A.list(k).features.entropy; %y as stimuli;
sptimescollect{k} = vertcat(A.list(k).sptimes{:}) + len;
len = len + length(A.list(k).y)/A.list(k).fs;
end
sumsptimes = vertcat(sptimescollect{:});

figure;
k = 2;
subplot(211)
Draw.spec(A.list(k).y,32000);
subplot(212)
plot(A.list(k).features.entropy);
xlim([0,length(A.list(k).features.entropy)])

%sumsdf = horzcat(sdfcollect{:});
sumsti = vertcat(sticollect{:});
sumy = vertcat(ycollect{:});
Stim = sumsti;
tsp = sumsptimes;
% latency = 0.0650;
% tsp = sumsptimes -latency ;

figure;
subplot(2,1,1)
Draw.spec(sumy,32000)
subplot(2,1,2)
plot(sumsti)

res_sumsdf = resample(sumsdf,length(sumsti),length(sumsdf));
res_sumsdf(res_sumsdf<0)=0;
sps_tr = res_sumsdf(1:ceil(length(sumsti)*4/5)).';
sps_tst = res_sumsdf(ceil(length(sumsti)*4/5)+1:end).';

Stim_tr = sumsti(1:ceil(length(sumsti)*4/5));
Stim_tst =sumsti(ceil(length(sumsti)*4/5)+1:end);

