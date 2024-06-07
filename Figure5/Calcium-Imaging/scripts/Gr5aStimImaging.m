%% Gr5aStimImaging
    


%% Script to analyze and record in struct array for Gr5a stim or other similar experiments
% This works with folder arrangement as in gr5a or behavior experiments,
% superfolder/ genotype folder / fly 1 / trial1 / 
subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 3.048; %frames per second


%Initialize preStim and postStim Frames for considering in a given trial

preStimFrames = 12;% 12 frames or 4 seconds
F0Frames = floor(1*fps); %number of prestim frames used for F0 calcuation
postStimFrames = 30;% 30 frames or 10 seconds
lightONPeriod = floor(2*fps); %2 seconds or 6 frames

%%
data={};

%Collect per genotype data in the following for loop
for i = 1:length(realSubFolders)
    currGenotype = realSubFolders(i).name;
    currPath = [pwd '\' realSubFolders(i).name];
    
    data(i).Genotype = currGenotype;
    data(i).Path = currPath;
    genotypeSubFolders = dir(currPath);
    genotypeSubFolders = genotypeSubFolders([genotypeSubFolders.isdir]');
    genotypeSubFolders = genotypeSubFolders (3:end);%remove blanck superdirectories and leave list of real subfolders
    data(i).numFlies = length(genotypeSubFolders); %assuming each fly is the immediate subfolder inside genotype folder

    
    
    for j=1:data(i).numFlies
        
        trialSubFolders = dir([genotypeSubFolders(j).folder '\', genotypeSubFolders(j).name]);
        trialSubFolders = trialSubFolders([trialSubFolders.isdir]');
        trialSubFolders = trialSubFolders (3:end);%remove blanck superdirectories and leave list of real subfolders
        data(i).perFlyData(j).numTrials = length(trialSubFolders);% assuming each trial is immediate subfolder inside the fly folder.
        data(i).perFlyData(j).trialSibFolderName = trialSubFolders(1).folder;
        for k = 1:data(i).perFlyData(j).numTrials
            
            trialStartFrames = xlsread([trialSubFolders(k).folder,'\',trialSubFolders(k).name,'\Stimulation_frames.xlsx']);%Stimulation Start Frames are written as frame numbers in a single column in an excel sheet in the trial subfolder
            tempData = xlsread([trialSubFolders(k).folder,'\',trialSubFolders(k).name,'\ROI.csv']); %read fluorescence data from ROI.csv file created by drawing ROIs in imageJ ROI editor>multimeasure
            tempData = tempData(:,2:end);%remove first column of csv which is just frame number.
            data(i).perFlyData(j).perTrialData(k).rawF = tempData;%add this to data structure as raw fluorescence data
            data(i).perFlyData(j).perTrialData(k).F = tempData - tempData(:,end); %subtract background assuming last ROI is background ROI
            data(i).perFlyData(j).perTrialData(k).F = data(i).perFlyData(j).perTrialData(k).F(:,1:end-1);%remove background ROI column from the new F data matrix
            data(i).perFlyData(j).perTrialData(k).stim1_F = data(i).perFlyData(j).perTrialData(k).F(trialStartFrames(1)-preStimFrames:trialStartFrames(1)+postStimFrames,:);%data around stim 1 from preStimFrames to postStimFrames
            %data(i).perFlyData(j).perTrialData(k).F = cat(2,data(i).perFlyData(j).perTrialData(k).F(trialStartFrames(1)-preStimFrames:trialStartFrames(1)+postStimFrames,:),data(i).perFlyData(j).perTrialData(k).F(trialStartFrames(2)-preStimFrames:trialStartFrames(2)+postStimFrames,:),data(i).perFlyData(j).perTrialData(k).F(trialStartFrames(3)-preStimFrames:trialStartFrames(3)+postStimFrames,:));%,data(i).perFlyData(j).perTrialData(k).F(trialStartFrames(4)-preStimFrames:trialStartFrames(4)+postStimFrames,:));
        
        end
        
        data(i).perFlyData(j).allTrials_stim1_F = cat(2,data(i).perFlyData(j).perTrialData.stim1_F);%concatenate all trials stim1 data into a matrix where every column is a new trial
        %data(i).perFlyData(j).avgTrial_stim1_F = mean(data(i).perFlyData(j).allTrials_stim1_F,2); %averaged trial data around stim1 from preStimFrames to postStimFrames.
        data(i).perFlyData(j).avgTrial_stim1_F = data(i).perFlyData(j).allTrials_stim1_F(:,1);
       
        data(i).perFlyData(j).Fb = mean(data(i).perFlyData(j).allTrials_stim1_F(preStimFrames-F0Frames:preStimFrames,:),1);%baseline or prestim F0 calculation. chose number of frames F0Frames
        data(i).perFlyData(j).Fd = data(i).perFlyData(j).allTrials_stim1_F-data(i).perFlyData(j).Fb;
        data(i).perFlyData(j).F_dF = data(i).perFlyData(j).Fd./data(i).perFlyData(j).Fb;%deltaF/F or (F-Fb)/Fb
        data(i).perFlyData(j).avgTrial_F_dF = mean(data(i).perFlyData(j).F_dF,2);%averaged trial deltaF/F
        data(i).perFlyData(j).ZFd = zscore(data(i).perFlyData(j).Fd,0,1);%z-scored deltaF as in Adesnik paper
        data(i).perFlyData(j).avgTrial_ZFd = mean(data(i).perFlyData(j).ZFd,2);%averaged trial zscored(deltaF)
        data(i).perFlyData(j).ZF_dF = zscore(data(i).perFlyData(j).F_dF,0,1);%z-scored deltaF/F
        data(i).perFlyData(j).avgTrial_ZF_dF = mean(data(i).perFlyData(j).ZF_dF,2);%averaged trial ZF_dF
        data(i).perFlyData(j).areaDeltaF_F = trapz(0:1/fps:lightONPeriod/fps,data(i).perFlyData(j).avgTrial_F_dF(preStimFrames+1:preStimFrames+lightONPeriod+1)); %area under deltaF/F curve during stim period
        %data(i).perFlyData(j).areaDeltaF_F_postStim = trapz(0:1/fps:(postStimFrames-lightONPeriod)/fps,data(i).perFlyData(j).avgTrial_F_dF(preStimFrames+lightONPeriod+1:end)); %area under deltaF/F curve in the post-stim period
        data(i).perFlyData(j).areaDeltaF_F_postStim = trapz(0:1/fps:lightONPeriod/fps,data(i).perFlyData(j).avgTrial_F_dF(preStimFrames+lightONPeriod+1:preStimFrames+lightONPeriod+lightONPeriod+1)); %area under deltaF/F curve in the post-stim period
        
    end
    
    
   
    %concatenate imporant perFly variables into matrix such that every
    %column is one fly or trial data
    data(i).catF = cat(2,data(i).perFlyData.allTrials_stim1_F);
    data(i).catFb = cat(2,data(i).perFlyData.Fb);
    data(i).catZF = zscore(cat(2,data(i).perFlyData.allTrials_stim1_F),1);
    data(i).catAvgTrial_F_dF = cat(2,data(i).perFlyData.avgTrial_F_dF);%every column is one fly averaged trial data
    data(i).catAvgTrial_ZFd = cat(2,data(i).perFlyData.avgTrial_ZFd);
    data(i).catAvgTrial_ZF_dF = cat(2,data(i).perFlyData.avgTrial_ZF_dF);
    data(i).catareaDeltaF_F = cat(2,data(i).perFlyData.areaDeltaF_F)';
    data(i).catareaDeltaF_F_postStim = cat(2,data(i).perFlyData.areaDeltaF_F_postStim)';
end




%% plot deltaF/F as bounded line plots

labels = {data.Genotype};
timeLine = -preStimFrames/fps:1/fps:postStimFrames/fps;
figure;

    hold on;
    p1=boundedline((timeLine)',mean(data(1).catAvgTrial_F_dF,2)',std(data(1).catAvgTrial_F_dF,0,2)'/sqrt(data(1).numFlies),'m','alpha');
    p2=boundedline((timeLine)',mean(data(2).catAvgTrial_F_dF,2)',std(data(2).catAvgTrial_F_dF,0,2)'/sqrt(data(2).numFlies),'g','alpha');
    p3=boundedline((timeLine)',mean(data(3).catAvgTrial_F_dF,2)',std(data(3).catAvgTrial_F_dF,0,2)'/sqrt(data(3).numFlies),'k','alpha');
    %p4=boundedline((timeLine)',mean(data(4).catAvgTrial_F_dF,2)',std(data(4).catAvgTrial_F_dF,0,2)'/sqrt(data(4).numFlies),'b','alpha');
    %p5=boundedline((timeLine)',mean(data(5).catAvgTrial_F_dF,2)',std(data(5).catAvgTrial_F_dF,0,2)'/sqrt(data(5).numFlies),'r','alpha');
    %p6=boundedline((timeLine)',mean(data(6).catAvgTrial_F_dF,2)',std(data(6).catAvgTrial_F_dF,0,2)'/sqrt(data(6).numFlies),'c','alpha');

    ylabel('deltaF/F');xlabel('time (s)');
    legend([p1,p2,p3],labels(1,[1,2,3]));
    %axis([-2,4, -10, 20]);
    line([0,0],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
%% plot are under curve for activation period.
labels = {data.Genotype};

aUC=cat(2,{data.catareaDeltaF_F});
maxLength = max(cellfun(@length,aUC));
aUCMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),aUC,'UniformOutput',false));
figure;errorbarjitter(aUCMat);title("aUC");xticklabels(labels);
%axis([xlim,0,300]);
%% plot are under curve for post-activation period.
labels = {data.Genotype};

aUC=cat(2,{data.catareaDeltaF_F_postStim});
maxLength = max(cellfun(@length,aUC));
aUCPostMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),aUC,'UniformOutput',false));
figure;errorbarjitter(aUCPostMat);title("aUC-postStim");xticklabels(labels);
%axis([xlim,0,300]);
%%

 %%

labels = {data.Genotype};
timeLine = -preStimFrames/fps:1/fps:postStimFrames/fps;
figure;

    hold on;
    p1=boundedline((timeLine)',mean(data(1).catAvgTrial_ZFd,2)',std(data(1).catAvgTrial_ZFd,0,2)'/sqrt(data(1).numFlies),'m','alpha');
    p2=boundedline((timeLine)',mean(data(2).catAvgTrial_ZFd,2)',std(data(2).catAvgTrial_ZFd,0,2)'/sqrt(data(2).numFlies),'g','alpha');
    p3=boundedline((timeLine)',mean(data(3).catAvgTrial_ZFd,2)',std(data(3).catAvgTrial_ZFd,0,2)'/sqrt(data(3).numFlies),'k','alpha');
    p4=boundedline((timeLine)',mean(data(4).catAvgTrial_ZFd,2)',std(data(4).catAvgTrial_ZFd,0,2)'/sqrt(data(4).numFlies),'b','alpha');
    p5=boundedline((timeLine)',mean(data(5).catAvgTrial_ZFd,2)',std(data(5).catAvgTrial_ZFd,0,2)'/sqrt(data(5).numFlies),'r','alpha');
    p6=boundedline((timeLine)',mean(data(6).catAvgTrial_ZFd,2)',std(data(6).catAvgTrial_ZFd,0,2)'/sqrt(data(6).numFlies),'c','alpha');
   
    ylabel('zscored(deltaF)');xlabel('time (s)');
    legend([p1,p2,p3,p4,p5,p6],labels);
    %axis([-2,4, -10, 20]);
    line([0,0],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

    %%

labels = {data.Genotype};
timeLine = -preStimFrames/fps:1/fps:postStimFrames/fps;
figure;

    hold on;
    p1=boundedline((timeLine)',mean(data(1).catAvgTrial_ZF_dF,2)',std(data(1).catAvgTrial_ZF_dF,0,2)'/sqrt(data(1).numFlies),'m','alpha');
    p2=boundedline((timeLine)',mean(data(2).catAvgTrial_ZF_dF,2)',std(data(2).catAvgTrial_ZF_dF,0,2)'/sqrt(data(2).numFlies),'g','alpha');
    p3=boundedline((timeLine)',mean(data(3).catAvgTrial_ZF_dF,2)',std(data(3).catAvgTrial_ZF_dF,0,2)'/sqrt(data(3).numFlies),'k','alpha');
    p4=boundedline((timeLine)',mean(data(4).catAvgTrial_ZF_dF,2)',std(data(4).catAvgTrial_ZF_dF,0,2)'/sqrt(data(4).numFlies),'b','alpha');
    p5=boundedline((timeLine)',mean(data(5).catAvgTrial_ZF_dF,2)',std(data(5).catAvgTrial_ZF_dF,0,2)'/sqrt(data(5).numFlies),'r','alpha');
    p6=boundedline((timeLine)',mean(data(6).catAvgTrial_ZF_dF,2)',std(data(6).catAvgTrial_ZF_dF,0,2)'/sqrt(data(6).numFlies),'c','alpha');


    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    ylabel('zscored(deltaF/F)');xlabel('time (s)');
    legend([p1,p2,p3,p4,p5,p6],labels);
    %axis([-2,4, -10, 20]);
    line([0,0],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
 %% plot deltaF/F as individual trials and mean trial as thick line
 figure;
 hold on;
 timeLine = -preStimFrames/fps:1/fps:postStimFrames/fps;
 plot(timeLine', data(5).catAvgTrial_F_dF','k');
 plot(timeLine', mean(data(5).catAvgTrial_F_dF,2)','LineWidth',2,'Color','k');
 plot(timeLine', data(6).catAvgTrial_F_dF','b');
 plot(timeLine', mean(data(6).catAvgTrial_F_dF,2)','LineWidth',2,'Color','b');
 plot(timeLine', data(2).catAvgTrial_F_dF','m');
 plot(timeLine', mean(data(2).catAvgTrial_F_dF,2)','LineWidth',2,'Color','m');
 line([0,0],ylim,'LineStyle','--','Color','k');
%%

 hold off;
 figure;
 hold on;
 timeLine = -preStimFrames/fps:1/fps:postStimFrames/fps;
 plot(timeLine', data(3).catAvgTrial_F_dF','k');
 plot(timeLine', mean(data(3).catAvgTrial_F_dF,2)','LineWidth',2,'Color','k');
 plot(timeLine', data(1).catAvgTrial_F_dF','b');
 plot(timeLine', mean(data(1).catAvgTrial_F_dF,2)','LineWidth',2,'Color','b');
 plot(timeLine', data(4).catAvgTrial_F_dF','m');
 plot(timeLine', mean(data(4).catAvgTrial_F_dF,2)','LineWidth',2,'Color','m');
 line([0,0],ylim,'LineStyle','--','Color','k');


 hold off;
%% heatmaps
figure;
imagesc(data(4).catAvgTrial_F_dF',[-1,1]);