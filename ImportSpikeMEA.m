%% by Diego Lopez-Pigozzi 2024 --- dlpigozzi@gmail.com ---
fino=isfinite(spiketable.Times);
e_names = unique(spiketable.Electrode(fino));
no_e = length(e_names);
e_table = table();
for i_elect = 1:1:no_e
%         i_elect
        my_index = strcmp(spiketable.Electrode, e_names(i_elect));
        
        Electrode = char(e_names(i_elect));
        eparts = split(Electrode, '_');
        Well = eparts(1);
        eName = eparts(2);
        No_Spikes = sum(my_index);
        Spike_Times = {spiketable.('Times')(my_index)'};
        Spike_Mags  = {spiketable.('AmplitudemV')(my_index)'};
        
        e_table(i_elect, :) = table(e_names(i_elect), Well, eName, No_Spikes, Spike_Times, Spike_Mags);
end

% i=330;j=335;
% MUA72_ts=e_table.Spike_Times{i,1};MUA72_amp=e_table.Spike_Mags{i,1};
% figure;plot(MUA72_ts,MUA72_amp);title(['Spikes electrode ',e_table.Var1{i}])
% MUA77_ts=e_table.Spike_Times{j,1};MUA77_amp=e_table.Spike_Mags{j,1};
% figure;plot(MUA77_ts,MUA77_amp);title(['Spikes electrode ',e_table.Var1{j}])
% figure;int=[-0.10,0.10];sep=400;
% [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(MUA72_ts, MUA77_ts, int);subplot(2,2,1);hist(tsOffsets,sep);title(['Crosscorrelation e',e_table.Var1{i},' VS e',e_table.Var1{j}]);xlim(int)
% [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(MUA77_ts, MUA77_ts, int);subplot(2,2,2);hist(tsOffsets,sep);title(['Autocorrelation e',e_table.Var1{i},' VS e',e_table.Var1{i}]);xlim(int)
% [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(MUA72_ts, MUA72_ts, int);subplot(2,2,3);hist(tsOffsets,sep);title(['Autocorrelation e',e_table.Var1{j},' VS e',e_table.Var1{j}]);xlim(int)
% 



