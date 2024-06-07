%% muscle GCaMP imaging

% This only works with speicific folder,subfolder arrangement as in
% selected_60-30_powdered_ACR1
subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(~startsWith({realSubFolders.name},{'movies'}));
realSubFolders = realSubFolders(~startsWith({realSubFolders.name},{'plots'}));
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 50; %frames per second

preMovFrames = 20;
postMovFrames = 50;%1s


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
    data(i).numFlies = length(genotypeSubFolders);
    
    
    MovFiles = getAllExtFiles(currPath,'lightON_ROI.csv');
    MovFramesFiles = getAllExtFiles(currPath,'lightON_frames.xlsx');
    %FlxtFiles = getAllExtFiles(currPath,'flext_lightON.csv');
    data(i).numLightONTrials = length(MovFiles);
    
    for j=1:data(i).numLightONTrials
       
        tempMovData = xlsread(char(MovFiles(j)));
        tempMovData = tempMovData(:,2:end);
        tempMovFrames = xlsread(char(MovFramesFiles(j)));
        
        data(i).perFlyData(j).Mov_rawF=tempMovData;
        data(i).perFlyData(j).MovFrames=tempMovFrames;
         
        
    %stimSeq = xlsread('Stimulation_frames.xlsx');
    %data(i).rawF =data(i).rawF(1:120,2:end);%special handling for movement in one 20180704male3t1
    %data(i).rawF =data(i).rawF(:,2:end);%remove serial number column
    %data(i).rawF =data(i).rawF(1:45,2:end);%special case for P9 male 3
        
        
        data(i).perFlyData(j).MovF = data(i).perFlyData(j).Mov_rawF;% - data(i).perFlyData(j).Mov_rawF(:,end);
        data(i).perFlyData(j).MovF = data(i).perFlyData(j).MovF(:,1:end-1);
        
        %reshape ExtMovF and FlxMovF into 4 extension and 4 flexion events
        %respectively using movement initiaion frame numbers from
        %tempMovFrames
        data(i).perFlyData(j).ExtMovF = cat(2,data(i).perFlyData(j).MovF(tempMovFrames(1,1)-preMovFrames:tempMovFrames(1,1)+postMovFrames,:),data(i).perFlyData(j).MovF(tempMovFrames(2,1)-preMovFrames:tempMovFrames(2,1)+postMovFrames,:),data(i).perFlyData(j).MovF(tempMovFrames(3,1)-preMovFrames:tempMovFrames(3,1)+postMovFrames,:),data(i).perFlyData(j).MovF(tempMovFrames(4,1)-preMovFrames:tempMovFrames(4,1)+postMovFrames,:));
        data(i).perFlyData(j).FlxMovF = cat(2,data(i).perFlyData(j).MovF(tempMovFrames(1,2)-preMovFrames:tempMovFrames(1,2)+postMovFrames,:),data(i).perFlyData(j).MovF(tempMovFrames(2,2)-preMovFrames:tempMovFrames(2,2)+postMovFrames,:),data(i).perFlyData(j).MovF(tempMovFrames(3,2)-preMovFrames:tempMovFrames(3,2)+postMovFrames,:),data(i).perFlyData(j).MovF(tempMovFrames(4,2)-preMovFrames:tempMovFrames(4,2)+postMovFrames,:));
       
       data(i).perFlyData(j).ExtMovFb = median(data(i).perFlyData(j).ExtMovF(1:15,:)); % median of pre-movement 15 frames was chosen as F0 based on empirical observation. using other paramters like prctl only make qualitative changes
       data(i).perFlyData(j).ExtMovFd = data(i).perFlyData(j).ExtMovF-data(i).perFlyData(j).ExtMovFb;
       data(i).perFlyData(j).ExtMovF_dF = data(i).perFlyData(j).ExtMovFd./data(i).perFlyData(j).ExtMovFb;
       data(i).perFlyData(j).ExtMovZFd = zscore(data(i).perFlyData(j).ExtMovFd,0,2);
       data(i).perFlyData(j).ExtMovZF_dF = zscore(data(i).perFlyData(j).ExtMovF_dF,0,2);
       
       %data(i).perFlyData(j).FlxMovFb = prctile(data(i).perFlyData(j).FlxMovF(:,:),10,1);
       data(i).perFlyData(j).FlxMovFb = median(data(i).perFlyData(j).FlxMovF(1:15,:));
       data(i).perFlyData(j).FlxMovFd = data(i).perFlyData(j).FlxMovF-data(i).perFlyData(j).FlxMovFb;
       data(i).perFlyData(j).FlxMovF_dF = data(i).perFlyData(j).FlxMovFd./data(i).perFlyData(j).FlxMovFb;
       data(i).perFlyData(j).FlxMovZFd = zscore(data(i).perFlyData(j).FlxMovFd,0,2);
       data(i).perFlyData(j).FlxMovZF_dF = zscore(data(i).perFlyData(j).FlxMovF_dF,0,2);
        
        
        
    end
    
   data(i).catExtMovF = cat(2,data(i).perFlyData.ExtMovF);
   data(i).catExtMovZF = zscore(cat(2,data(i).perFlyData.ExtMovF),1);
   data(i).catExtMovF_dF = cat(2,data(i).perFlyData.ExtMovF_dF);
   data(i).catExtMovZFd = cat(2,data(i).perFlyData.ExtMovZFd);
   data(i).catExtMovFd = cat(2,data(i).perFlyData.ExtMovFd);
   data(i).catExtMovZF_dF = cat(2,data(i).perFlyData.ExtMovZF_dF);
   
   
   data(i).catFlxMovF = cat(2,data(i).perFlyData.FlxMovF);
   data(i).catFlxMovZF = zscore(cat(2,data(i).perFlyData.FlxMovF),1);
   data(i).catFlxMovF_dF = cat(2,data(i).perFlyData.FlxMovF_dF);
   data(i).catFlxMovZFd = cat(2,data(i).perFlyData.FlxMovZFd);
   data(i).catFlxMovFd = cat(2,data(i).perFlyData.FlxMovFd);
   data(i).catFlxMovZF_dF = cat(2,data(i).perFlyData.FlxMovZF_dF);
    
end 




%% plot pooled trial data as mean+-SEM plots for each genotype for flexion movement

figure;hold on;
data(2).sampleSize = size(data(2).catExtMovF_dF(:,1:4:end),2); % number of pooled trials
trialLength = preMovFrames + postMovFrames +1;

    hold on;
    p1=boundedline((1:trialLength)',mean(data(2).catFlxMovF_dF(:,1:4:end),2)',std(data(2).catFlxMovF_dF(:,1:4:end),0,2)'/sqrt(data(2).sampleSize),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(2).catFlxMovF_dF(:,2:4:end),2)',std(data(2).catFlxMovF_dF(:,2:4:end),0,2)'/sqrt(data(2).sampleSize),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(2).catFlxMovF_dF(:,3:4:end),2)',std(data(2).catFlxMovF_dF(:,3:4:end),0,2)'/sqrt(data(2).sampleSize),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(2).catFlxMovF_dF(:,4:4:end),2)',std(data(2).catFlxMovF_dF(:,4:4:end),0,2)'/sqrt(data(2).sampleSize),'k','alpha');
    
    line(xlim,[0,0],'LineStyle','-','Color','k');
    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    %ylabel('Avg Trial Velocity +-SEM (mm/s)');xlabel('time (s)');
    legend([p1,p2,p3,p4],{'ROI1','ROI2','ROI3','ROI4'});
    axis([0,71, -0.5, 1.2]);
    line([20,20],ylim,'LineStyle','--','Color','k');
figure;
data(1).sampleSize = size(data(1).catExtMovF_dF(:,1:4:end),2); % number of pooled trials

    hold on;
    p1=boundedline((1:trialLength)',mean(data(1).catFlxMovF_dF(:,1:4:end),2)',std(data(1).catFlxMovF_dF(:,1:4:end),0,2)'/sqrt(data(1).sampleSize),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(1).catFlxMovF_dF(:,2:4:end),2)',std(data(1).catFlxMovF_dF(:,2:4:end),0,2)'/sqrt(data(1).sampleSize),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(1).catFlxMovF_dF(:,3:4:end),2)',std(data(1).catFlxMovF_dF(:,3:4:end),0,2)'/sqrt(data(1).sampleSize),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(1).catFlxMovF_dF(:,4:4:end),2)',std(data(1).catFlxMovF_dF(:,4:4:end),0,2)'/sqrt(data(1).sampleSize),'k','alpha');
    
    line(xlim,[0,0],'LineStyle','-','Color','k');
    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    %ylabel('Avg Trial Velocity +-SEM (mm/s)');xlabel('time (s)');
    legend([p1,p2,p3,p4],{'ROI1','ROI2','ROI3','ROI4'});
    axis([0,71, -0.5, 1.2]);
    line([20,20],ylim,'LineStyle','--','Color','k');

%% plot pooled trial data as mean+-SEM plots for each genotype for extension movement

figure;hold on;

trialLength = preMovFrames + postMovFrames +1;
data(2).sampleSize = size(data(2).catExtMovF_dF(:,1:4:end),2); % number of pooled trials

    hold on;
    p1=boundedline((1:trialLength)',mean(data(2).catExtMovF_dF(:,1:4:end),2)',std(data(2).catExtMovF_dF(:,1:4:end),0,2)'/sqrt(data(2).sampleSize),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(2).catExtMovF_dF(:,2:4:end),2)',std(data(2).catExtMovF_dF(:,2:4:end),0,2)'/sqrt(data(2).sampleSize),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(2).catExtMovF_dF(:,3:4:end),2)',std(data(2).catExtMovF_dF(:,3:4:end),0,2)'/sqrt(data(2).sampleSize),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(2).catExtMovF_dF(:,4:4:end),2)',std(data(2).catExtMovF_dF(:,4:4:end),0,2)'/sqrt(data(2).sampleSize),'k','alpha');
    
    line(xlim,[0,0],'LineStyle','-','Color','k');
    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    %ylabel('Avg Trial Velocity +-SEM (mm/s)');xlabel('time (s)');
    legend([p1,p2,p3,p4],{'ROI1','ROI2','ROI3','ROI4'});
    axis([0,71, -0.5, 1.2]);
    line([20,20],ylim,'LineStyle','--','Color','k');
    %axis([0,90*fps, 0, 40]);
figure;
data(1).sampleSize = size(data(1).catExtMovF_dF(:,1:4:end),2); % number of pooled trials
    hold on;
    p1=boundedline((1:trialLength)',mean(data(1).catExtMovF_dF(:,1:4:end),2)',std(data(1).catExtMovF_dF(:,1:4:end),0,2)'/sqrt(data(1).sampleSize),'m','alpha');
    p2=boundedline((1:trialLength)',mean(data(1).catExtMovF_dF(:,2:4:end),2)',std(data(1).catExtMovF_dF(:,2:4:end),0,2)'/sqrt(data(1).sampleSize),'g','alpha');
    p3=boundedline((1:trialLength)',mean(data(1).catExtMovF_dF(:,3:4:end),2)',std(data(1).catExtMovF_dF(:,3:4:end),0,2)'/sqrt(data(1).sampleSize),'b','alpha');
    p4=boundedline((1:trialLength)',mean(data(1).catExtMovF_dF(:,4:4:end),2)',std(data(1).catExtMovF_dF(:,4:4:end),0,2)'/sqrt(data(1).sampleSize),'k','alpha');
    
    line(xlim,[0,0],'LineStyle','-','Color','k');
    %set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:15*fps:max(xlim),'XTickLabel',0:15:max(xlim)/fps);
    %ylabel('Avg Trial Velocity +-SEM (mm/s)');xlabel('time (s)');
    legend([p1,p2,p3,p4],{'ROI1','ROI2','ROI3','ROI4'});
    axis([0,71, -0.5, 1.2]);
    line([20,20],ylim,'LineStyle','--','Color','k');

%%