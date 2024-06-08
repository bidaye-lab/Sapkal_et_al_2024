%% muscle GCaMP imaging
% This works with folder arrangement as in glue muscle imaging
subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(3:end);

% remove unwanted folders like those containing movies or plots from the
% list of real data-containing subfolders.
realSubFolders = realSubFolders(~startsWith({realSubFolders.name},{'movies'}));
realSubFolders = realSubFolders(~startsWith({realSubFolders.name},{'plots'}));

% Initialise video temporal parameters
fps = 50; %frames per second

maxLastFrame = 2500;

% Define pre and post-stimulus time window to consider for plotting
preStimFrames = 100;
postStimFrames = 200;

%%
data={};

%Collect per genotype data in the following for loop
for i = 1:length(realSubFolders)
    [currGenotype, Cx] = strtok(realSubFolders(i).name,'_');
    currPath = [pwd '\' realSubFolders(i).name];
    
    data(i).Genotype = currGenotype;
    data(i).Path = currPath;
    
    
    genotypeSubFolders = dir(currPath);
    genotypeSubFolders = genotypeSubFolders([genotypeSubFolders.isdir]');
    genotypeSubFolders = genotypeSubFolders (3:end);
    
    flexedsubFolders = genotypeSubFolders(startsWith({genotypeSubFolders.name},{'flex'})); % folders contained flexed Fe-Ti trials
    extendedsubFolders = genotypeSubFolders(startsWith({genotypeSubFolders.name},{'extend'})); % folders contained extended Fe-Ti trials
    
    data(i).numFlexedFlies = length(flexedsubFolders);
    data(i).numExtendedFlies = length(extendedsubFolders);

    %Collect per "flexed" fly data in this for loop
    for j=1:data(i).numFlexedFlies
        FeFilePath = [flexedsubFolders(j).folder,'\',flexedsubFolders(j).name,'\Fe']; % get the correct path to Fe muscle data
        currSubFolders = dir(FeFilePath);
        realCurrSubFolders = currSubFolders([currSubFolders.isdir]');
        realCurrSubFolders = realCurrSubFolders(3:end);
        data(i).perFlyData(j).numFlxTrials = length(realCurrSubFolders);
        
        %Collect per "trial" data for each fly data in this for loop. Each
        %fly has 1-2 trials and each trial has 4 optogenetic stimulations.
        for k = 1:data(i).perFlyData(j).numFlxTrials
            trialStartFrames = xlsread([FeFilePath, '\',realCurrSubFolders(k).name,'\stimFrames.xlsx']);%read stimulation start frame from xlsx file
            tempFlxData = xlsread([FeFilePath, '\',realCurrSubFolders(k).name,'\ROI.csv']);%read mean fluorescence per ROI data from csv file generated in ImageJ
            tempFlxData = tempFlxData(:,2:end); % remove serial number column
            data(i).perFlyData(j).perFlxTrialData(k).Flx_rawF = tempFlxData; 
            data(i).perFlyData(j).perFlxTrialData(k).FlxF = tempFlxData;%background ROI, last column not needed for 1P data due to no stimulus artifact.
            data(i).perFlyData(j).perFlxTrialData(k).FlxF = data(i).perFlyData(j).perFlxTrialData(k).FlxF(:,1:end-1);
            
            %reshape the data into matrix where every column corresponds to
            %one optogenetic stimluation
            data(i).perFlyData(j).perFlxTrialData(k).FlxTrialsF = cat(2,data(i).perFlyData(j).perFlxTrialData(k).FlxF(trialStartFrames(1)-preStimFrames:trialStartFrames(1)+postStimFrames,:),data(i).perFlyData(j).perFlxTrialData(k).FlxF(trialStartFrames(2)-preStimFrames:trialStartFrames(2)+postStimFrames,:),data(i).perFlyData(j).perFlxTrialData(k).FlxF(trialStartFrames(3)-preStimFrames:trialStartFrames(3)+postStimFrames,:),data(i).perFlyData(j).perFlxTrialData(k).FlxF(trialStartFrames(4)-preStimFrames:trialStartFrames(4)+postStimFrames,:));
        
        end
        %merge data across trials (4 stims per trial) into one matrix
        if data(i).perFlyData(j).numFlxTrials > 1
            data(i).perFlyData(j).allTrials_FlxTrialsF = cat(2,data(i).perFlyData(j).perFlxTrialData.FlxTrialsF);
        else
            data(i).perFlyData(j).allTrials_FlxTrialsF = data(i).perFlyData(j).perFlxTrialData(1).FlxTrialsF;
        end
        
        %define baseline (F0) as 10th percentile of pre-stimulus period
        data(i).perFlyData(j).FlxTrialsFb = prctile(data(i).perFlyData(j).allTrials_FlxTrialsF(1:preStimFrames,:),10,1);
        
        %Fd=F-F0
        data(i).perFlyData(j).FlxTrialsFd = data(i).perFlyData(j).allTrials_FlxTrialsF-data(i).perFlyData(j).FlxTrialsFb;
        
        %F_dF = Fd/F0 OR F-F0/F0 (usual definitiion of deltaF/F)
        data(i).perFlyData(j).FlxTrialsF_dF = data(i).perFlyData(j).FlxTrialsFd./data(i).perFlyData(j).FlxTrialsFb;
        
        %segregate F_dF into separate averaged trial variables for each muscle for easy
        %plotting
        data(i).perFlyData(j).meanFlxTrialsF_dF1 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,1:4:end),2);
        data(i).perFlyData(j).meanFlxTrialsF_dF2 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,2:4:end),2);
        data(i).perFlyData(j).meanFlxTrialsF_dF3 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,3:4:end),2);
        data(i).perFlyData(j).meanFlxTrialsF_dF4 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,4:4:end),2);
    end
    
    %do the same for extended (separate loop in case different number of
    %flies for flexed and extended
    for j=1:data(i).numExtendedFlies
        FeFilePath = [extendedsubFolders(j).folder,'\',extendedsubFolders(j).name,'\Fe'];
        currSubFolders = dir(FeFilePath);
        realCurrSubFolders = currSubFolders([currSubFolders.isdir]');
        realCurrSubFolders = realCurrSubFolders(3:end);
        data(i).perFlyData(j).numExtTrials = length(realCurrSubFolders);
        
        for k = 1:data(i).perFlyData(j).numExtTrials
            trialStartFrames = xlsread([FeFilePath, '\',realCurrSubFolders(k).name,'\stimFrames.xlsx']);
            tempExtData = xlsread([FeFilePath, '\',realCurrSubFolders(k).name,'\ROI.csv']);
            tempExtData = tempExtData(:,2:end);
            data(i).perFlyData(j).perExtTrialData(k).Ext_rawF = tempExtData;
            data(i).perFlyData(j).perExtTrialData(k).ExtF = tempExtData;
            data(i).perFlyData(j).perExtTrialData(k).ExtF = data(i).perFlyData(j).perExtTrialData(k).ExtF(:,1:end-1);
            data(i).perFlyData(j).perExtTrialData(k).ExtTrialsF = cat(2,data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(1)-preStimFrames:trialStartFrames(1)+postStimFrames,:),data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(2)-preStimFrames:trialStartFrames(2)+postStimFrames,:),data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(3)-preStimFrames:trialStartFrames(3)+postStimFrames,:),data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(4)-preStimFrames:trialStartFrames(4)+postStimFrames,:));
        
        end
        if data(i).perFlyData(j).numExtTrials > 1
            data(i).perFlyData(j).allTrials_ExtTrialsF = cat(2,data(i).perFlyData(j).perExtTrialData.ExtTrialsF);
        else
            data(i).perFlyData(j).allTrials_ExtTrialsF = data(i).perFlyData(j).perExtTrialData(1).ExtTrialsF;
        end
       
        data(i).perFlyData(j).ExtTrialsFb = prctile(data(i).perFlyData(j).allTrials_ExtTrialsF(1:preStimFrames,:),10,1);
        data(i).perFlyData(j).ExtTrialsFd = data(i).perFlyData(j).allTrials_ExtTrialsF-data(i).perFlyData(j).ExtTrialsFb;
        data(i).perFlyData(j).ExtTrialsF_dF = data(i).perFlyData(j).ExtTrialsFd./data(i).perFlyData(j).ExtTrialsFb;
        
        data(i).perFlyData(j).meanExtTrialsF_dF1 = mean(data(i).perFlyData(j).ExtTrialsF_dF(:,1:4:end),2);
        data(i).perFlyData(j).meanExtTrialsF_dF2 = mean(data(i).perFlyData(j).ExtTrialsF_dF(:,2:4:end),2);
        data(i).perFlyData(j).meanExtTrialsF_dF3 = mean(data(i).perFlyData(j).ExtTrialsF_dF(:,3:4:end),2);
        data(i).perFlyData(j).meanExtTrialsF_dF4 = mean(data(i).perFlyData(j).ExtTrialsF_dF(:,4:4:end),2);
        
    end   
       
    % Generate concatenated matrices for easy plotting per muscle
    data(i).catMeanFlxTrialsF_dF1 = cat(2,data(i).perFlyData.meanFlxTrialsF_dF1);
    data(i).catMeanFlxTrialsF_dF2 = cat(2,data(i).perFlyData.meanFlxTrialsF_dF2);
    data(i).catMeanFlxTrialsF_dF3 = cat(2,data(i).perFlyData.meanFlxTrialsF_dF3);
    data(i).catMeanFlxTrialsF_dF4 = cat(2,data(i).perFlyData.meanFlxTrialsF_dF4);

    data(i).catMeanExtTrialsF_dF1 = cat(2,data(i).perFlyData.meanExtTrialsF_dF1);
    data(i).catMeanExtTrialsF_dF2 = cat(2,data(i).perFlyData.meanExtTrialsF_dF2);
    data(i).catMeanExtTrialsF_dF3 = cat(2,data(i).perFlyData.meanExtTrialsF_dF3);
    data(i).catMeanExtTrialsF_dF4 = cat(2,data(i).perFlyData.meanExtTrialsF_dF4);

    % generate variance normalized data for each muscle only in case of
    % inhibition experiments. use variance normalization instead of z-score
    % because we want to keep 0 at real 0 that corresponds to "no muscle
    % activity".
    data(i).catMeanFlxTrialsZF_dF1 = normalize(data(i).catMeanFlxTrialsF_dF1,1,'scale');
    data(i).catMeanFlxTrialsZF_dF2 = normalize(data(i).catMeanFlxTrialsF_dF2,1,'scale');
    data(i).catMeanFlxTrialsZF_dF3 = normalize(data(i).catMeanFlxTrialsF_dF3,1,'scale');
    data(i).catMeanFlxTrialsZF_dF4 = normalize(data(i).catMeanFlxTrialsF_dF4,1,'scale');

    data(i).catMeanExtTrialsZF_dF1 = normalize(data(i).catMeanExtTrialsF_dF1,1,'scale');
    data(i).catMeanExtTrialsZF_dF2 = normalize(data(i).catMeanExtTrialsF_dF2,1,'scale');
    data(i).catMeanExtTrialsZF_dF3 = normalize(data(i).catMeanExtTrialsF_dF3,1,'scale');
    data(i).catMeanExtTrialsZF_dF4 = normalize(data(i).catMeanExtTrialsF_dF4,1,'scale');
    
end 


%% plot all flexed trials per genotype as mean+-SEM bounded line plots

trialLength = 301;
figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsF_dF1,2)',std(data(1).catMeanFlxTrialsF_dF1,0,2)'/sqrt(data(1).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsF_dF2,2)',std(data(1).catMeanFlxTrialsF_dF2,0,2)'/sqrt(data(1).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsF_dF3,2)',std(data(1).catMeanFlxTrialsF_dF3,0,2)'/sqrt(data(1).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsF_dF4,2)',std(data(1).catMeanFlxTrialsF_dF4,0,2)'/sqrt(data(1).numFlexedFlies),'k','alpha');
    
    ylabel('deltaF/F');xlabel('frames');
    title([data(1).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 13]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsF_dF1,2)',std(data(2).catMeanFlxTrialsF_dF1,0,2)'/sqrt(data(2).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsF_dF2,2)',std(data(2).catMeanFlxTrialsF_dF2,0,2)'/sqrt(data(2).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsF_dF3,2)',std(data(2).catMeanFlxTrialsF_dF3,0,2)'/sqrt(data(2).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsF_dF4,2)',std(data(2).catMeanFlxTrialsF_dF4,0,2)'/sqrt(data(2).numFlexedFlies),'k','alpha');
    
    ylabel('deltaF/F');xlabel('frames');
    title([data(2).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 13]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsF_dF1,2)',std(data(3).catMeanFlxTrialsF_dF1,0,2)'/sqrt(data(3).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsF_dF2,2)',std(data(3).catMeanFlxTrialsF_dF2,0,2)'/sqrt(data(3).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsF_dF3,2)',std(data(3).catMeanFlxTrialsF_dF3,0,2)'/sqrt(data(3).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsF_dF4,2)',std(data(3).catMeanFlxTrialsF_dF4,0,2)'/sqrt(data(3).numFlexedFlies),'k','alpha');
    
    ylabel('deltaF/F');xlabel('frames');
    title([data(3).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 13]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsF_dF1,2)',std(data(4).catMeanFlxTrialsF_dF1,0,2)'/sqrt(data(4).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsF_dF2,2)',std(data(4).catMeanFlxTrialsF_dF2,0,2)'/sqrt(data(4).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsF_dF3,2)',std(data(4).catMeanFlxTrialsF_dF3,0,2)'/sqrt(data(4).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsF_dF4,2)',std(data(4).catMeanFlxTrialsF_dF4,0,2)'/sqrt(data(4).numFlexedFlies),'k','alpha');
    
    ylabel('deltaF/F');xlabel('frames');
    title([data(4).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 13]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
%% plot all variance-normalized flexed trials per genotype as mean+-SEM bounded line plots


trialLength = 301;
figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsZF_dF1,2)',std(data(1).catMeanFlxTrialsZF_dF1,0,2)'/sqrt(data(1).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsZF_dF2,2)',std(data(1).catMeanFlxTrialsZF_dF2,0,2)'/sqrt(data(1).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsZF_dF3,2)',std(data(1).catMeanFlxTrialsZF_dF3,0,2)'/sqrt(data(1).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(1).catMeanFlxTrialsZF_dF4,2)',std(data(1).catMeanFlxTrialsZF_dF4,0,2)'/sqrt(data(1).numFlexedFlies),'k','alpha');
    
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(1).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF1,2)',std(data(2).catMeanFlxTrialsZF_dF1,0,2)'/sqrt(data(2).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF2,2)',std(data(2).catMeanFlxTrialsZF_dF2,0,2)'/sqrt(data(2).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF3,2)',std(data(2).catMeanFlxTrialsZF_dF3,0,2)'/sqrt(data(2).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF4,2)',std(data(2).catMeanFlxTrialsZF_dF4,0,2)'/sqrt(data(2).numFlexedFlies),'k','alpha');
    
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(2).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsZF_dF1,2)',std(data(3).catMeanFlxTrialsZF_dF1,0,2)'/sqrt(data(3).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsZF_dF2,2)',std(data(3).catMeanFlxTrialsZF_dF2,0,2)'/sqrt(data(3).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsZF_dF3,2)',std(data(3).catMeanFlxTrialsZF_dF3,0,2)'/sqrt(data(3).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(3).catMeanFlxTrialsZF_dF4,2)',std(data(3).catMeanFlxTrialsZF_dF4,0,2)'/sqrt(data(3).numFlexedFlies),'k','alpha');
    
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(3).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
figure;
%subplot(2,1,1);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsZF_dF1,2)',std(data(4).catMeanFlxTrialsZF_dF1,0,2)'/sqrt(data(4).numFlexedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsZF_dF2,2)',std(data(4).catMeanFlxTrialsZF_dF2,0,2)'/sqrt(data(4).numFlexedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsZF_dF3,2)',std(data(4).catMeanFlxTrialsZF_dF3,0,2)'/sqrt(data(4).numFlexedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(4).catMeanFlxTrialsZF_dF4,2)',std(data(4).catMeanFlxTrialsZF_dF4,0,2)'/sqrt(data(4).numFlexedFlies),'k','alpha');
    
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(4).Genotype,' Flexed']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
%% plot all extended trials per genotype as mean+-SEM bounded line plots

trialLength = 301;
    figure;
    hold on;
    p1=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsF_dF1,2)',std(data(1).catMeanExtTrialsF_dF1,0,2)'/sqrt(data(1).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsF_dF2,2)',std(data(1).catMeanExtTrialsF_dF2,0,2)'/sqrt(data(1).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsF_dF3,2)',std(data(1).catMeanExtTrialsF_dF3,0,2)'/sqrt(data(1).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsF_dF4,2)',std(data(1).catMeanExtTrialsF_dF4,0,2)'/sqrt(data(1).numExtendedFlies),'k','alpha');
   
    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    ylabel('deltaF/F');xlabel('frames');
    title([data(1).Genotype,' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 9]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');


    
hold off;
figure;

%subplot(2,1,2);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsF_dF1,2)',std(data(2).catMeanExtTrialsF_dF1,0,2)'/sqrt(data(2).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsF_dF2,2)',std(data(2).catMeanExtTrialsF_dF2,0,2)'/sqrt(data(2).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsF_dF3,2)',std(data(2).catMeanExtTrialsF_dF3,0,2)'/sqrt(data(2).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsF_dF4,2)',std(data(2).catMeanExtTrialsF_dF4,0,2)'/sqrt(data(2).numExtendedFlies),'k','alpha');
   
    ylabel('deltaF/F');xlabel('frames');
    title([data(2).Genotype, ' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 9]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    figure;
%subplot(2,1,2);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsF_dF1,2)',std(data(3).catMeanExtTrialsF_dF1,0,2)'/sqrt(data(3).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsF_dF2,2)',std(data(3).catMeanExtTrialsF_dF2,0,2)'/sqrt(data(3).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsF_dF3,2)',std(data(3).catMeanExtTrialsF_dF3,0,2)'/sqrt(data(3).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsF_dF4,2)',std(data(3).catMeanExtTrialsF_dF4,0,2)'/sqrt(data(3).numExtendedFlies),'k','alpha');
   
    ylabel('deltaF/F');xlabel('frames');
    title([data(3).Genotype, ' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 9]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    
    figure;
%subplot(2,1,2);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsF_dF1,2)',std(data(4).catMeanExtTrialsF_dF1,0,2)'/sqrt(data(4).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsF_dF2,2)',std(data(4).catMeanExtTrialsF_dF2,0,2)'/sqrt(data(4).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsF_dF3,2)',std(data(4).catMeanExtTrialsF_dF3,0,2)'/sqrt(data(4).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsF_dF4,2)',std(data(4).catMeanExtTrialsF_dF4,0,2)'/sqrt(data(4).numExtendedFlies),'k','alpha');
   
    ylabel('deltaF/F');xlabel('frames');
    title([data(4).Genotype, ' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 9]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');

%% plot all variance-normalized extended trials per genotype as mean+-SEM bounded line plots

trialLength = 301;
    figure;
    hold on;
    p1=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsZF_dF1,2)',std(data(1).catMeanExtTrialsZF_dF1,0,2)'/sqrt(data(1).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsZF_dF2,2)',std(data(1).catMeanExtTrialsZF_dF2,0,2)'/sqrt(data(1).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsZF_dF3,2)',std(data(1).catMeanExtTrialsZF_dF3,0,2)'/sqrt(data(1).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(1).catMeanExtTrialsZF_dF4,2)',std(data(1).catMeanExtTrialsZF_dF4,0,2)'/sqrt(data(1).numExtendedFlies),'k','alpha');
   
    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(1).Genotype,' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');


    
hold off;
figure;

%subplot(2,1,2);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF1,2)',std(data(2).catMeanExtTrialsZF_dF1,0,2)'/sqrt(data(2).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF2,2)',std(data(2).catMeanExtTrialsZF_dF2,0,2)'/sqrt(data(2).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF3,2)',std(data(2).catMeanExtTrialsZF_dF3,0,2)'/sqrt(data(2).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF4,2)',std(data(2).catMeanExtTrialsZF_dF4,0,2)'/sqrt(data(2).numExtendedFlies),'k','alpha');
   
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(2).Genotype, ' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
    figure;
%subplot(2,1,2);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsZF_dF1,2)',std(data(3).catMeanExtTrialsZF_dF1,0,2)'/sqrt(data(3).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsZF_dF2,2)',std(data(3).catMeanExtTrialsZF_dF2,0,2)'/sqrt(data(3).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsZF_dF3,2)',std(data(3).catMeanExtTrialsZF_dF3,0,2)'/sqrt(data(3).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(3).catMeanExtTrialsZF_dF4,2)',std(data(3).catMeanExtTrialsZF_dF4,0,2)'/sqrt(data(3).numExtendedFlies),'k','alpha');
   
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(3).Genotype, ' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;
    
    figure;
%subplot(2,1,2);
    hold on;
    p1=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsZF_dF1,2)',std(data(4).catMeanExtTrialsZF_dF1,0,2)'/sqrt(data(4).numExtendedFlies),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsZF_dF2,2)',std(data(4).catMeanExtTrialsZF_dF2,0,2)'/sqrt(data(4).numExtendedFlies),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsZF_dF3,2)',std(data(4).catMeanExtTrialsZF_dF3,0,2)'/sqrt(data(4).numExtendedFlies),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(4).catMeanExtTrialsZF_dF4,2)',std(data(4).catMeanExtTrialsZF_dF4,0,2)'/sqrt(data(4).numExtendedFlies),'k','alpha');
   
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(4).Genotype, ' Extended']);
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    axis([0,300, -0.5, 4]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

%% plot single fly traces.
trialLength = 301;
    figure;
    hold on;
    p1=plot((1:trialLength)',data(2).catMeanFlxTrialsZF_dF1','m');
    %pm1 = plot((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF1,2)','m','LineWidth',2);
    p2=plot((1:trialLength)',data(2).catMeanFlxTrialsZF_dF2','g');
    %pm2 = plot((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF2,2)','g','LineWidth',2);
    p3=plot((1:trialLength)',data(2).catMeanFlxTrialsZF_dF3','b');
    %pm3 = plot((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF3,2)','b','LineWidth',2);
    p4=plot((1:trialLength)',data(2).catMeanFlxTrialsZF_dF4','k');
    %pm4 = plot((1:trialLength)',mean(data(2).catMeanFlxTrialsZF_dF4,2)','k','LineWidth',2);

    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(2).Genotype,' Flexed']);
    %legend([pm1,pm2,pm3,pm4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    %axis([0,300 -1, 10]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');


    
hold off;
%% plot single fly traces.

trialLength = 301;
    figure;
    hold on;
    p1=plot((1:trialLength)',data(2).catMeanExtTrialsZF_dF1','m');
    %pm1 = plot((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF1,2)','m','LineWidth',2);
    p2=plot((1:trialLength)',data(2).catMeanExtTrialsZF_dF2','g');
    %pm2 = plot((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF2,2)','g','LineWidth',2);
    p3=plot((1:trialLength)',data(2).catMeanExtTrialsZF_dF3','b');
    %pm3 = plot((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF3,2)','b','LineWidth',2);
    p4=plot((1:trialLength)',data(2).catMeanExtTrialsZF_dF4','k');
    %pm4 = plot((1:trialLength)',mean(data(2).catMeanExtTrialsZF_dF4,2)','k','LineWidth',2);

    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    ylabel('norm(deltaF/F)');xlabel('frames');
    title([data(2).Genotype,' Extended']);
    %legend([pm1,pm2,pm3,pm4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    %axis([0,300 -1, 10]);
    line([100,100],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');


    
hold off;
%%
