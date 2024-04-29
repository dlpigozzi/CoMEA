%% by Diego Lopez-Pigozzi 2024 --- dlpigozzi@gmail.com ---
%  copyrighted by GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
% here the link to the paper once published
%% Import spike list
%---------- Assembloids Kiki
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220324_Assem-NHDF-chip_SpikeList_Plate01_004DIV_xxxDIV_xxxDIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220330_Assem-NHDF-chip_SpikeList_Plate01_010DIV_006DIV_xxxDIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220406_Assem-NHDF-chip_SpikeList_Plate01_017DIV_013DIV_007DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220416_Assem-NHDF-chip_SpikeList_Plate01_027DIV_023DIV_017DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220420_Assem-NHDF-chip_SpikeList_Plate01_031DIV_027DIV_021DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220424_Assem-NHDF-chip_SpikeList_Plate01_035DIV_031DIV_025DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220429_Assem-NHDF-chip_SpikeList_Plate01_040DIV_036DIV_030DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220420_Assem-NHDF-chip_Plate01_031DIV_027DIV_021DIV(decrypted)(000)_spike_list.csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220503_Assem-NHDF-chip_SpikeList_Plate01_044DIV_040DIV_034DIV(decrypted).csv';
file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220526_Assem-NHDF-chip_SpikeList_Plate01_067DIV_063DIV_057DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220726_Assem-NHDF-chip_Plate01_128DIV_124DIV_118DIV(decrypted)(000)_spike_list.csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220701_Assem-NHDF-chip_SpikeList_Plate01_103DIV_099DIV_093DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220712_Assem-NHDF-chip_SpikeList_Plate01_114DIV_110DIV_104DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220715_Assem-NHDF-chip_SpikeList_Plate01_117DIV_113DIV_107DIV(decrypted).csv';
% file='C:\Users\PCO\Documents\DATA\MEA_recordings_assembloids_chips\Raw_data\20220526_Assem-NHDF-chip_Plate01_067DIV_063DIV_057DIV(decrypted)(002)_spike_list.csv';
%% Import preprocessed data
spiketable=importingCSV_SpikeList(file);ImportSpikeMEA;

%% Fixed parameters
well_name='B4';lag=10;STTC_THR=0.8;bin=0.0005;
activity_threshold=0.5;
confidence=5;

%% Organize the spatial activity per electrode
duration=ceil(spiketable.Times(end-10));sentence=strfind(file,'\');
well= strcmp(e_table.Well,well_name);
spikeTimes=(e_table.Spike_Times(well));
coords_active=(e_table.eName(well));
meanFR=zeros(8,8);big_meanFR=zeros(24,24);int_meanFR=zeros(24,24);
for i=1:length(coords_active)
    pair=coords_active{i,1};
    meanFR(str2num(pair(1)),str2num(pair(2)))=length(spikeTimes{i})/duration;
    int_meanFR(str2num(pair(1))*3-1,str2num(pair(2))*3-1)=length(spikeTimes{i})/duration;
end
big_meanFR=imresize(meanFR,15);

%% Calculate STTC (readapted from MEA-NAP toolbox)
[adjM, adjMci] = adjM_thr_parallel_diego(spikeTimes, lag, 0.0001, 12500,...
    duration, 1000);
figure; imagesc(adjM),title(['STTC lag=',num2str(lag),'ms, duration=',num2str(duration),'s']);colorbar

%% RASTER WELL
p=figure;hold on;xlim([0 duration]);ylim([0 length(spikeTimes)]);ylabel('Electrode');xlabel('Time s');
for i=1:length(spikeTimes)
    spikes=spikeTimes{i,:};
    plot(spikes,ones(1,length(spikes))*i,'.k','MarkerSize',1)
end
set(p,'Position',[395    95   970   702])

%% Recognize the active eletrodes and their coordinates
new_coords=zeros(length(coords_active),2);
for i=1:length(coords_active)
    pair=coords_active{i,1};
   new_coords(i,1)=str2double(pair(1));
   new_coords(i,2)=str2double(pair(2));
end

%% Calculate peaks of xcorr for electrodes with high STTC
sep=4000;edges=-0.1005:bin:0.1005;int=[-0.100,0.100];
[row,col] = find(adjM>STTC_THR);Mcounts=zeros(size(adjM));Mcenters=zeros(size(adjM));lag_chosen_bin=NaN(7,length(row));
for k=1:length(row)
    MUAk_ts=spikeTimes{row(k),1};MUAm_ts=spikeTimes{col(k),1};
    [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(MUAk_ts, MUAm_ts, int);
    [counts,centers]=hist(tsOffsets,edges);
    x=max(counts);avg=mean(counts(300:350));dev=std(counts(300:350));idx=find(counts>(avg+confidence*dev));
    idx=find(counts==x);
    if abs(counts(idx(1)))>abs(avg+confidence*dev)
        if length(MUAk_ts)/duration>activity_threshold
            if abs(centers(idx(1)))>0.00025
                lag_chosen_bin(1,k)=centers(idx(1));lag_chosen_bin(2,k)=counts(idx(1));lag_chosen_bin(3,k)=avg;lag_chosen_bin(4,k)=dev;lag_chosen_bin(5,k)=row(k);lag_chosen_bin(6,k)=col(k);lag_chosen_bin(7,k)=adjM(row(k),col(k));
                Mcounts(row(k),col(k))=counts(idx(1));
                Mcenters(row(k),col(k))=centers(idx(1));
            end
        end
    end
end

%% Plot the wind (arrows between synchronous electrodes based on STTC with direction determined by xcorr)
figure;imagesc(meanFR);colorbar;hold on;title(['Windplot lag=',num2str(lag),'ms']);
vectors=zeros(length(col)/2,5);
for i=1:(length(col)/2)
    if Mcenters(row(i),col(i))<0 
        vectors(i,1:4)=[new_coords(col(i),2),new_coords(col(i),1),(new_coords(row(i),2)-new_coords(col(i),2)),(new_coords(row(i),1)-new_coords(col(i),1))];
        vectors(i,5)=adjM(row(i),col(i));
        quiver(new_coords(col(i),2),new_coords(col(i),1),(new_coords(row(i),2)-new_coords(col(i),2)),(new_coords(row(i),1)-new_coords(col(i),1)),...
            'LineWidth',[(adjM(row(i),col(i)))/0.2],...
            'AutoScale','off',...
            'Color',[0 0 0])
    elseif Mcenters(row(i),col(i))>0 
        quiver(new_coords(row(i),2),new_coords(row(i),1),(-new_coords(row(i),2)+new_coords(col(i),2)),(-new_coords(row(i),1)+new_coords(col(i),1)),...
            'LineWidth',[(adjM(row(i),col(i)))/0.2],...
            'AutoScale','off',...
            'Color',[0 0 0])
    end
    
end

%% Plot the wind (here you can obtain the vectorial sum for each node)
ker = fspecial('gaussian',[3 3],1);
for i=1:5 %smooth the heat map
    kk=imfilter(big_meanFR,ker);
    big_meanFR=kk;
end
factor=length(kk)/length(meanFR);
forXticks=[1 2 3 4 5 6 7 8].*factor-ceil(factor/2)+1;
figure;imagesc(kk);colorbar;hold on;xticks(forXticks);xticklabels({'1','2','3','4','5','6','7','8'});yticks(forXticks);yticklabels({'1','2','3','4','5','6','7','8'});
c = colorbar;c.Label.String = 'Firing rate (Hz)';xlabel('Electrode');ylabel('Electrode');title('Vectorial sum of pairs');
colormap(jet(300));caxis([0 max(max(kk))])

for j=1:length(new_coords)
    hub=find((vectors(:,1)==new_coords(j,1))&(vectors(:,2)==new_coords(j,2)));neighbours=0;tosum=zeros(8,5);
        for m=1:length(hub)
                neighbours=neighbours+1;
                tosum(neighbours,:)=vectors(hub(m),:);
        end
        if neighbours>1
                quiver(tosum(1,1)*factor-ceil(factor/2)+1,tosum(1,2)*factor-ceil(factor/2)+1,sum(tosum(:,3))/neighbours*factor,sum(tosum(:,4))/neighbours*factor,...
                'LineWidth',2,...
                'AutoScale','on',...
                'Color',[0 0 0])
        elseif neighbours==1
            quiver(tosum(1,1)*factor-ceil(factor/2)+1,tosum(1,2)*factor-ceil(factor/2)+1,tosum(1,3)*factor,tosum(1,4)*factor,...
                'LineWidth',2,...
                'AutoScale','on',...
                'Color',[0 0 0])
        end
end
