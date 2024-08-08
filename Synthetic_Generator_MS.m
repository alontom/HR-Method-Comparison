% The authors give no warranty for the correct functioning of the software
% and cannot be held legally accountable.
% Intiate 200x2 cells for each dyad
% R-R intervals
dyadRR=cell(200,2);
% BPM (resampling window = min(R-R) or avg(R-R) or max(R-R) x 2
dyadBPMmin=cell(200,2);
dyadBPMavg=cell(200,2);
dyadBPMmax2=cell(200,2);

% Iterate 100 times
for it=1:100
    % Preparing 1/f noise and stretched time series

    % Generate 1/f noise
    ts = powernoise(1,1024); % length ~10 minutes
    % Linear transformation to reach ~1000 amplitude, ~100 std
    ts = ts/mean(ts); ts = round(27*(ts-min(ts)+25));
    % Same process for ts2 'weak coupling'
    tsA = powernoise(1,1024);
    tsA = tsA/mean(tsA); tsA = round(27*(tsA-min(tsA)+25)); 
    % Stretch the time series - set value x for x time points
    tempts = [];
    temptsA = [];
    for j = 1:length(ts)
        tempts(length(tempts)+1:length(tempts)+ts(j)) = ts(j);
    end
    for j = 1:length(tsA)
        temptsA(length(temptsA)+1:length(temptsA)+tsA(j)) = tsA(j);
    end

    % Condition 1: 'Coupled'
    % ts1 and ts2 from the same distribution    

    % ts1: Moving average of 1600 data points
    ts1=round(mean(reshape(tempts(1:end-rem(length(tempts),1600)),1600,[])));
    % ts2: Moving average of 1800 data points 
    ts2=mean(reshape(tempts(1:end-rem(length(tempts),1800)),1800,[]));
    % Adding noise to ts1: SNR = 1:10
    ts1=ts1+(randn(1,length(ts1))*std(ts1)/1.5);
    ts1=round(ts1);
    % Adding noise to ts2: SNR = 1:10
    ts2=ts2+(randn(1,length(ts2))*std(ts2)/1.5);
    % Equalize the sum of ts1 and ts2
    ts2=round(ts2*sum(ts1)/sum(ts2));
    ts2(end)=ts2(end)-(sum(ts2)-sum(ts1));

    % Check 1/f with detrended fluctuation analysis
    if dfa(ts1)>0.75 && dfa(ts1)<1.25 && dfa(ts2)>0.75 && dfa(ts2)<1.25
        % Assign RR and BPM
        dyadRR(it,1)={ts1};
        dyadRR(it,2)={ts2};
    else
        disp(['Check Iteration:',num2str(dfa(ts2))])
    end

    % Condition 2: 'Uncoupled' 
    % ts1 and ts2 from different distributions
    
    % ts1: Moving average of 1600 data points
    ts1=round(mean(reshape(tempts(1:end-rem(length(tempts),1600)),1600,[])));
    % ts2: Moving average of 1800 data points (on temptsA)
    ts2=mean(reshape(temptsA(1:end-rem(length(temptsA),1800)),1800,[]));
    % Adding noise to ts1: SNR = 1:10
    ts1=ts1+(randn(1,length(ts1))*std(ts1)/1.5);
    ts1=round(ts1);
    % Adding noise to ts2: SNR = 1:10
    ts2=ts2+(randn(1,length(ts2))*std(ts2)/1.5);
    % Equalize the sum of ts1 and ts2
    ts2=round(ts2*sum(ts1)/sum(ts2));
    ts2(end)=ts2(end)-(sum(ts2)-sum(ts1));

    % Check 1/f with detrended fluctuation analysis
    if dfa(ts1)>0.75 && dfa(ts1)<1.25 && dfa(ts2)>0.75 && dfa(ts2)<1.25
        % Assign RR and BPM
        dyadRR(it+100,1)={ts1};
        dyadRR(it+100,2)={ts2};
    else
        disp(['check FD',num2str(dfa(ts1))])
    end
end

% Find min, max and average R-R intervals
minR=100000;
maxR=-500;
avgR=0;
for it=1:200
    minR=min([minR,min(dyadRR{it,1}),min(dyadRR{it,2})]);
    avgR=avgR+(mean(dyadRR{it,1})+mean(dyadRR{it,2}))/2;
    maxR=max([maxR,max(dyadRR{it,1}),max(dyadRR{it,2})]);
end
avgR=round(avgR/it)/1000;
minR=round(minR/1000,3); maxR=round(maxR/1000,3); 

% Convert to BPM - ws = min(RR), avg(RR), max(RR)*2
for it=1:200
    disp(it)
    dyadBPMmin{it,1}=RR2BPM(dyadRR{it,1},minR,minR);
    dyadBPMmin{it,2}=RR2BPM(dyadRR{it,2},minR,minR);
    dyadBPMavg{it,1}=RR2BPM(dyadRR{it,1},avgR,avgR);
    dyadBPMavg{it,2}=RR2BPM(dyadRR{it,2},avgR,avgR);
    dyadBPMmax2{it,1}=RR2BPM(dyadRR{it,1},maxR*2,maxR*2);
    dyadBPMmax2{it,2}=RR2BPM(dyadRR{it,2},maxR*2,maxR*2);
end
save('Synthetic_R.mat',"dyadRR","dyadBPMmin","dyadBPMavg","dyadBPMmax2","minR","avgR","maxR")
