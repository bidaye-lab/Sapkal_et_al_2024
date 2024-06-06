%% muscle GCaMP imaging with 2P imaging and legs glued in flexed or extended positions
    
% This works with folder arrangement as in glue muscle imaging
subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');

% remove unwanted folders like those containing movies or plots from the
% list of real data-containing subfolders.
realSubFolders = realSubFolders(~startsWith({realSubFolders.name},{'movies'}));
realSubFolders = realSubFolders(~startsWith({realSubFolders.name},{'plots'}));
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 6; %frames per second

% define trial period by setting frames before and after each opto stim
preStimFrames = 2*fps;
postStimFrames = 4*fps;

%% Generate data structure containing all relevant raw and processed data for plotting
data={};

%Collect per genotype data in the following for loop
for i = 1:length(realSubFolders)
    currGenotype= realSubFolders(i).name;
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
            tempFlxData = tempFlxData(:,2:end);% remove serial number column
            data(i).perFlyData(j).perFlxTrialData(k).Flx_rawF = tempFlxData;
            data(i).perFlyData(j).perFlxTrialData(k).FlxF = tempFlxData - tempFlxData(:,end);%subtract background fluorescence (make sure last ROI is always background ROI)
            data(i).perFlyData(j).perFlxTrialData(k).FlxF = data(i).perFlyData(j).perFlxTrialData(k).FlxF(:,1:end-1); %remove the last background ROI from the data after background subtraction

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
        data(i).perFlyData(j).FlxTrialsFb = mean(data(i).perFlyData(j).allTrials_FlxTrialsF(1:preStimFrames,:),1);
        
        %Fd=F-F0
        data(i).perFlyData(j).FlxTrialsFd = data(i).perFlyData(j).allTrials_FlxTrialsF-data(i).perFlyData(j).FlxTrialsFb;
        
        %F_dF = Fd/F0 OR F-F0/F0 (usual definitiion of deltaF/F)
        data(i).perFlyData(j).FlxTrialsF_dF = data(i).perFlyData(j).FlxTrialsFd./data(i).perFlyData(j).FlxTrialsFb;
        
        %segregate F_dF into separate variables for each muscle for easy
        %plotting
        data(i).perFlyData(j).meanFlxTrialsF_dF1 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,1:4:end),2);
        data(i).perFlyData(j).meanFlxTrialsF_dF2 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,2:4:end),2);
        data(i).perFlyData(j).meanFlxTrialsF_dF3 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,3:4:end),2);
        data(i).perFlyData(j).meanFlxTrialsF_dF4 = mean(data(i).perFlyData(j).FlxTrialsF_dF(:,4:4:end),2);
    end
    
    % do the same for extended (separate loop in case different number of
    % flies for flexed and extended
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
            data(i).perFlyData(j).perExtTrialData(k).ExtF = tempExtData - tempExtData(:,end);
            data(i).perFlyData(j).perExtTrialData(k).ExtF = data(i).perFlyData(j).perExtTrialData(k).ExtF(:,1:end-1);
            data(i).perFlyData(j).perExtTrialData(k).ExtTrialsF = cat(2,data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(1)-preStimFrames:trialStartFrames(1)+postStimFrames,:),data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(2)-preStimFrames:trialStartFrames(2)+postStimFrames,:),data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(3)-preStimFrames:trialStartFrames(3)+postStimFrames,:),data(i).perFlyData(j).perExtTrialData(k).ExtF(trialStartFrames(4)-preStimFrames:trialStartFrames(4)+postStimFrames,:));
        
        end
        if data(i).perFlyData(j).numExtTrials>1
            data(i).perFlyData(j).allTrials_ExtTrialsF = cat(2,data(i).perFlyData(j).perExtTrialData.ExtTrialsF);
        else
            data(i).perFlyData(j).allTrials_ExtTrialsF = data(i).perFlyData(j).perExtTrialData(1).ExtTrialsF;
        end
       
        data(i).perFlyData(j).ExtTrialsFb = mean(data(i).perFlyData(j).allTrials_ExtTrialsF(1:preStimFrames,:),1);
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
    
end 


%% plot all flexed trials per genotype as mean+-SEM bounded line plots

fps=6;
timeLine = -2*fps/fps:1/fps:4*fps/fps;
figure;

    hold on;
    p1=boundedline((timeLine)',mean(data(1).catMeanFlxTrialsF_dF1,2)',std(data(1).catMeanFlxTrialsF_dF1,0,2)'/sqrt(data(1).numFlexedFlies),'m','alpha');
    p2=boundedline((timeLine)',mean(data(1).catMeanFlxTrialsF_dF2,2)',std(data(1).catMeanFlxTrialsF_dF2,0,2)'/sqrt(data(1).numFlexedFlies),'g','alpha');
    p3=boundedline((timeLine)',mean(data(1).catMeanFlxTrialsF_dF3,2)',std(data(1).catMeanFlxTrialsF_dF3,0,2)'/sqrt(data(1).numFlexedFlies),'b','alpha');
    p4=boundedline((timeLine)',mean(data(1).catMeanFlxTrialsF_dF4,2)',std(data(1).catMeanFlxTrialsF_dF4,0,2)'/sqrt(data(1).numFlexedFlies),'k','alpha');
    
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-a'});
    ylabel('deltaF/F)');xlabel('time');
    title([data(1).Genotype,' Flexed']);
    axis([-2,4, -1, 12]);
    line([0,0],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

%%  plot all extended trials per genotype as mean+-SEM bounded line plots

fps=6;
timeLine = -2*fps/fps:1/fps:4*fps/fps;
figure;

    hold on;
    p1=boundedline((timeLine)',mean(data(1).catMeanExtTrialsF_dF1,2)',std(data(1).catMeanExtTrialsF_dF1,0,2)'/sqrt(data(1).numExtendedFlies),'m','alpha');
    p2=boundedline((timeLine)',mean(data(1).catMeanExtTrialsF_dF2,2)',std(data(1).catMeanExtTrialsF_dF2,0,2)'/sqrt(data(1).numExtendedFlies),'g','alpha');
    p3=boundedline((timeLine)',mean(data(1).catMeanExtTrialsF_dF3,2)',std(data(1).catMeanExtTrialsF_dF3,0,2)'/sqrt(data(1).numExtendedFlies),'b','alpha');
    p4=boundedline((timeLine)',mean(data(1).catMeanExtTrialsF_dF4,2)',std(data(1).catMeanExtTrialsF_dF4,0,2)'/sqrt(data(1).numExtendedFlies),'k','alpha');
    
    legend([p1,p2,p3,p4],{'Extensor','Flexor','AccFlexor-b','AccFlexor-b'});
    ylabel('deltaF/F)');xlabel('time');
    title([data(1).Genotype,' Extended']);
    axis([-2,4, -1, 12]);
    line([0,0],ylim,'LineStyle','--','Color','k');
    line(xlim,[0,0],'LineStyle','-','Color','k');
    hold off;

%%