% The authors give no warranty for the correct functioning of the software
% and cannot be held legally accountable.
% Intiate 400x2 cells for each dyad
dyadRR=cell(400,2);
dyadBPM=cell(400,2);
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

    % Condition 1: 'strong coupling' without noise 
    % ts1 and ts2 from the same distribution

    % ts1: Moving average of 1600 data points
    ts1=round(mean(reshape(tempts(1:end-rem(length(tempts),1600)),1600,[])));
    % ts2: Moving average of 1800 data points 
    ts2=mean(reshape(tempts(1:end-rem(length(tempts),1800)),1800,[]));
    % Equalize the sum of ts1 and ts2
    ts2=round(ts2*sum(ts1)/sum(ts2));
    ts2(end)=ts2(end)-(sum(ts2)-sum(ts1));

    % Check 1/f with detrended fluctuation analysis
    if dfa(ts1)>0.75 && dfa(ts1)<1.25 && dfa(ts2)>0.75 && dfa(ts2)<1.25
        % Transform to BPM with 0.5 seconds window and no overlap
        ts1BPM=RR2BPM(ts1,0.5,0.5);
        ts2BPM=RR2BPM(ts2,0.5,0.5);
        % Assign RR and BPM
        dyadRR(it,1)={ts1};
        dyadRR(it,2)={ts2};
        dyadBPM(it,1)={ts1BPM};
        dyadBPM(it,2)={ts2BPM};
    else
        disp(['Check  Iteration:',num2str(dfa(ts1))])
    end

    % Condition 2: 'strong coupling' with noise 
    % ts1 and ts2 from the same distribution    

    % ts1: Moving average of 1600 data points
    ts1=round(mean(reshape(tempts(1:end-rem(length(tempts),1600)),1600,[])));
    % ts2: Moving average of 1800 data points 
    ts2=mean(reshape(tempts(1:end-rem(length(tempts),1800)),1800,[]));
    % Adding noise to ts1: SNR = 1:1.5
    ts1=ts1+(randn(1,length(ts1))*std(ts1)/1.5);
    ts1=round(ts1);
    % Adding noise to ts2: SNR = 1:1.5
    ts2=ts2+(randn(1,length(ts2))*std(ts2)/1.5);
    % Equalize the sum of ts1 and ts2
    ts2=round(ts2*sum(ts1)/sum(ts2));
    ts2(end)=ts2(end)-(sum(ts2)-sum(ts1));

    % Check 1/f with detrended fluctuation analysis
    if dfa(ts1)>0.75 && dfa(ts1)<1.25 && dfa(ts2)>0.75 && dfa(ts2)<1.25
        % Transform to BPM with 0.5 seconds window and no overlap
        ts1BPM=RR2BPM(ts1,0.5,0.5);
        ts2BPM=RR2BPM(ts2,0.5,0.5);
        % Assign RR and BPM
        dyadRR(it+100,1)={ts1};
        dyadRR(it+100,2)={ts2};
        dyadBPM(it+100,1)={ts1BPM};
        dyadBPM(it+100,2)={ts2BPM};
    else
        disp(['Check Iteration:',num2str(dfa(ts2))])
    end

    % Condition 3: 'weak coupling' without noise 
    % ts1 and ts2 from different distributions
    % ts1: Moving average of 1600 data points
    ts1=round(mean(reshape(tempts(1:end-rem(length(tempts),1600)),1600,[])));
    % ts2: Moving average of 1800 data points (on temptsA)
    ts2=mean(reshape(temptsA(1:end-rem(length(temptsA),1800)),1800,[]));
    % Equalize the sum of ts1 and ts2
    ts2=round(ts2*sum(ts1)/sum(ts2));
    ts2(end)=ts2(end)-(sum(ts2)-sum(ts1));

    % Check 1/f with detrended fluctuation analysis
    if dfa(ts1)>0.75 && dfa(ts1)<1.25 && dfa(ts2)>0.75 && dfa(ts2)<1.25
        % Transform to BPM with 0.5 seconds window and no overlap
        ts1BPM=RR2BPM(ts1,0.5,0.5);
        ts2BPM=RR2BPM(ts2,0.5,0.5);
        % Assign RR and BPM
        dyadRR(it+200,1)={ts1};
        dyadRR(it+200,2)={ts2};
        dyadBPM(it+200,1)={ts1BPM};
        dyadBPM(it+200,2)={ts2BPM};
    else
        disp(['Check Iteration:',num2str(dfa(ts1))])
    end

    % Condition 4: 'weak coupling' with noise 
    % ts1 and ts2 from different distributions
    
    % ts1: Moving average of 1600 data points
    ts1=round(mean(reshape(tempts(1:end-rem(length(tempts),1600)),1600,[])));
    % ts2: Moving average of 1800 data points (on temptsA)
    ts2=mean(reshape(temptsA(1:end-rem(length(temptsA),1800)),1800,[]));
    % Adding noise to ts1: SNR = 1:1.5
    ts1=ts1+(randn(1,length(ts1))*std(ts1)/1.5);
    ts1=round(ts1);
    % Adding noise to ts2: SNR = 1:1.5
    ts2=ts2+(randn(1,length(ts2))*std(ts2)/1.5);
    % Equalize the sum of ts1 and ts2
    ts2=round(ts2*sum(ts1)/sum(ts2));
    ts2(end)=ts2(end)-(sum(ts2)-sum(ts1));

    % Check 1/f with detrended fluctuation analysis
    if dfa(ts1)>0.8 && dfa(ts1)<1.2 && dfa(ts2)>0.8 && dfa(ts2)<1.2
        % Transform to BPM with 0.5 seconds window and no overlap
        ts1BPM=RR2BPM(ts1,0.5,0.5);
        ts2BPM=RR2BPM(ts2,0.5,0.5);
        % Assign RR and BPM
        dyadRR(it+300,1)={ts1};
        dyadRR(it+300,2)={ts2};
        dyadBPM(it+300,1)={ts1BPM};
        dyadBPM(it+300,2)={ts2BPM};
    else
        disp(['check FD',num2str(dfa(ts1))])
    end
end