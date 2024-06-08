%% Script for analysis of food-blob interactions (Sapkal et al 2024 Fig. 5f-h)
% This only works with speicific folder,subfolder arrangement as in
% superfolder/genotype subfolders/fly folders/tracking data
% add path of edmond.de folder
% run each section in order from start to end, one section at a time.

subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 30; %frames per second

%maxLastFrameNumber = 53600;
maxLastFrameNumber = 52000; % all flies detect food-blob much before this frame

%Initialize smoothing and threshold parameters for velocity attributes
smoothWindow = 15; %frames
angSmoothWindow= 15; %frames
walkThreshold = 2; %mm/s

%% Generate basic data structure containing all tracking derived attributes
data={};

%Collect per genotype data in the following for loop
for i = 1:length(realSubFolders)
    currGenotype = realSubFolders(i).name;
    currPath = [pwd '\' realSubFolders(i).name];
    
    data(i).Genotype = currGenotype;
    data(i).Path = currPath;
    
    genotypeSubFolders = dir(currPath);
    genotypeSubFolders = genotypeSubFolders([genotypeSubFolders.isdir]');
    genotypeSubFolders = genotypeSubFolders (3:end);
    data(i).numFlies = length(genotypeSubFolders);

    %Collect per fly data in the following for loop
    for j=1:data(i).numFlies
        currFlyPath = [genotypeSubFolders(j).folder,'\',genotypeSubFolders(j).name];
        data(i).perFlyData(j).path = currFlyPath;
        
        %read data from flytracker output files for each fly
        tempMatFileFeat = matfile([currFlyPath, '\', genotypeSubFolders(j).name, '-feat.mat']);
        tempMatFileTrk = matfile([currFlyPath, '\', genotypeSubFolders(j).name, '-track.mat']);
        tempMatFileCal = matfile([currFlyPath, '\', genotypeSubFolders(j).name, '-calibration.mat']);
        
        tempCal = tempMatFileCal.calib;
        tempFeature = tempMatFileFeat.feat;
        tempTrk = tempMatFileTrk.trk;
        tempVel = tempFeature.data(:,:,strcmp(tempFeature.names,'vel'));
        if(isempty(tempVel))
            break
        end
        tempAngVel = tempFeature.data(:,:,strcmp(tempFeature.names,'ang_vel'));
        tempOri = tempTrk.data(:,:,strcmp(tempTrk.names,'ori'));
        tempPosX = tempTrk.data(:,:,strcmp(tempTrk.names,'pos x'));
        tempPosY = tempTrk.data(:,:,strcmp(tempTrk.names,'pos y'));
        
        %at times velocity arrays are 1 or 2 frames longer/shorter than
        %expected. In such cases the following condition makes vel and
        %angVel arrays into consistent size as manually provided in variable
        %maxLastFrameNumber. flag fff alerts user of these cases
        
        if length(tempVel)<maxLastFrameNumber
           %fff=2222 
           currFeatLastFrame = length(tempVel);
           tempVel(end+1:maxLastFrameNumber)=mean(tempVel);%(currFeatLastFrame);
           tempAngVel(end+1:maxLastFrameNumber)=mean(tempAngVel);%(currFeatLastFrame);
        else
           %fff=111
           tempVel = tempVel(1:maxLastFrameNumber);
           tempAngVel = tempAngVel(1:maxLastFrameNumber);
        end
        
    % Since Trk and Feat files created seperately, do if statement once more
        if length(tempOri)<maxLastFrameNumber
           currTrkLastFrame = length(tempOri);
           tempOri(end+1:maxLastFrameNumber)=tempOri(currTrkLastFrame);
           tempPosX(end+1:maxLastFrameNumber)=tempPosX(currTrkLastFrame);
           tempPosY(end+1:maxLastFrameNumber)=tempPosY(currTrkLastFrame);
        else
           tempOri = tempOri(1:maxLastFrameNumber);
           tempPosX = tempPosX(1:maxLastFrameNumber);
           tempPosY = tempPosY(1:maxLastFrameNumber);
        end
        
        timetags = ([1:length(tempVel)])/fps;%convert frames to time in s
        
        
        data(i).perFlyData(j).feat = tempFeature;%contains flytracker output for velocities
        data(i).perFlyData(j).trk = tempTrk;%contains flytracker output for position and orientation
        data(i).perFlyData(j).ppm = tempCal.PPM;%contains flytracker calibration of pixels per mm
        
        
            
        data(i).perFlyData(j).vel = tempVel;%translational velocity
        data(i).perFlyData(j).angVel = tempAngVel;%absolute angular velocity

        data(i).perFlyData(j).ori = tempOri;%orientation
        data(i).perFlyData(j).PosX = tempPosX;%X position
        data(i).perFlyData(j).PosY = tempPosY;% Y position
        
        data(i).perFlyData(j).smoothVel = (smooth(tempVel,smoothWindow))';
        
        %walked = data(i).perFlyData(j).smoothVel>walkThreshold;
        
        %data(i).perFlyData(j).meanWalkedVel = mean(data(i).perFlyData(j).smoothVel(walked));
        
        data(i).perFlyData(j).smoothAngVel = (smooth(tempAngVel,angSmoothWindow))';
       
      
    end
    %define concatenated variables for easy plotting
    data(i).catVel = cat(1,data(i).perFlyData.vel);
    data(i).catSmoothVel = cat(1,data(i).perFlyData.smoothVel);
    
    data(i).catAngVel = cat(1,data(i).perFlyData.angVel);
    data(i).catSmoothAngVel = cat(1,data(i).perFlyData.smoothAngVel);
    
    
end



%% time spent in center of arena
distThreshold = 3;%mm
velThreshold = 2;
angVelThreshold = 1.5;
%wallDistThreshold = 20;%mm
for i = 1:length(data)
    

    %Collect per fly data in the following for loop
    for j=1:data(i).numFlies
        
        currFlyPath = data(i).perFlyData(j).path;
        
        blobLocation = xlsread([currFlyPath, '\frame.csv']);%read manually annotated location of food-blob
        data(i).perFlyData(j).blobLocation = blobLocation;
        
        %calculate distance of centroid of fly from food-blob
        data(i).perFlyData(j).CentroidDist = (sqrt((data(i).perFlyData(j).PosX-blobLocation(1)).^2+(data(i).perFlyData(j).PosY-blobLocation(2)).^2))./data(i).perFlyData(j).ppm;
        
        %define centralFrames as all video frames when fly was
        %within distThreshold to food-blob
        data(i).perFlyData(j).centralFrames = data(i).perFlyData(j).CentroidDist<distThreshold;
        tempCF=data(i).perFlyData(j).centralFrames;
        
        data(i).perFlyData(j).centralTime = nnz(data(i).perFlyData(j).centralFrames)/fps;
        data(i).perFlyData(j).centralVel = data(i).perFlyData(j).smoothVel(tempCF);
        data(i).perFlyData(j).centralAngVel = data(i).perFlyData(j).smoothAngVel(tempCF);
        data(i).perFlyData(j).centralMedianVel = mean(data(i).perFlyData(j).centralVel);
        data(i).perFlyData(j).centralMedianAngVel = mean(data(i).perFlyData(j).centralAngVel);

        TempStopFrames1 = data(i).perFlyData(j).centralVel<velThreshold & data(i).perFlyData(j).centralAngVel<angVelThreshold;
        TempStopFrames2 = [false,TempStopFrames1(1:end-1)];
        TempStopFrames3 = [false,TempStopFrames2(1:end-1)];
        TempStopFrames4 = [false,TempStopFrames3(1:end-1)];
        TempStopFrames5 = [false,TempStopFrames4(1:end-1)];
        TempStop = TempStopFrames1 & TempStopFrames2 & TempStopFrames3 & TempStopFrames4 & TempStopFrames5;
        
        %data(i).perFlyData(j).centralStopTime = nnz(data(i).perFlyData(j).centralVel<velThreshold & data(i).perFlyData(j).centralVel<angVelThreshold)/fps;
        data(i).perFlyData(j).centralStopTime = nnz(TempStop)/fps;

    end
    data(i).catCentralTime = cat(2,data(i).perFlyData.centralTime)';
    data(i).catCentralStopTime = cat(2,data(i).perFlyData.centralStopTime)';
    data(i).catCentralStopFraction = data(i).catCentralStopTime./data(i).catCentralTime;
    data(i).catCentralVel = cat(2,data(i).perFlyData.centralMedianVel)';
    data(i).catCentralAngVel = cat(2,data(i).perFlyData.centralMedianAngVel)';
    
end
%% food encounter period
lowerLimit = 1;
upperLimit = maxLastFrameNumber - 100;
encounterWindow = 150;
preEncounterWindow = 0;
velThreshold = 2;
angVelThreshold = 1.5;

for i = 1:length(data)
    [currGenotype, ACR] = strtok(realSubFolders(i).name,'_');
    currPath = [pwd '\' realSubFolders(i).name];
    for j=1:data(i).numFlies
       currFlyPath = data(i).perFlyData(j).path;
       encounterFrame = xlsread([currFlyPath, '\blob.xlsx']);
       
       data(i).perFlyData(j).encounterFrame = encounterFrame;
       data(i).perFlyData(j).postEncounterFrames = maxLastFrameNumber - encounterFrame;
       
       if encounterFrame>lowerLimit && encounterFrame<upperLimit
           data(i).perFlyData(j).encounterSmoothVel = data(i).perFlyData(j).smoothVel(1,encounterFrame-preEncounterWindow:encounterFrame+encounterWindow);
           data(i).perFlyData(j).encounterSmoothAngVel = data(i).perFlyData(j).smoothAngVel(1,encounterFrame-preEncounterWindow:encounterFrame+encounterWindow);
           data(i).perFlyData(j).encounterCentralFrames = data(i).perFlyData(j).centralFrames(1,encounterFrame-preEncounterWindow:encounterFrame+encounterWindow);
           tempECF = data(i).perFlyData(j).encounterCentralFrames;
           data(i).perFlyData(j).encounterCentralVel = data(i).perFlyData(j).encounterSmoothVel(tempECF);
           data(i).perFlyData(j).encounterCentralMedianVel = median(data(i).perFlyData(j).encounterCentralVel);
           data(i).perFlyData(j).encounterCentralAngVel = data(i).perFlyData(j).encounterSmoothAngVel(tempECF);
           data(i).perFlyData(j).encounterCentralMedianAngVel = median(data(i).perFlyData(j).encounterCentralAngVel);
           data(i).perFlyData(j).encounterCentralStopFrames = data(i).perFlyData(j).encounterCentralVel<velThreshold & data(i).perFlyData(j).encounterCentralAngVel<angVelThreshold;
           
           % Define stopping event when at least 10 continuous frames are
           % less than velocity threshold. This stringent definition is
           % necessary due to food blob causing some tracking jitter and
           % empirical observation that most feeding related stopping
           % events are longer than 10 frames or 333 ms.
           TempEncounterStop1 = data(i).perFlyData(j).encounterCentralVel<velThreshold & data(i).perFlyData(j).encounterCentralAngVel<angVelThreshold;
           TempEncounterStop2 = [false,TempEncounterStop1(1:end-1)];
           TempEncounterStop3 = [false,TempEncounterStop2(1:end-1)];
           TempEncounterStop4 = [false,TempEncounterStop3(1:end-1)];
           TempEncounterStop5 = [false,TempEncounterStop4(1:end-1)];
           TempEncounterStop6 = [TempEncounterStop1(2:end),false];
           TempEncounterStop7 = [TempEncounterStop6(2:end),false];
           TempEncounterStop8 = [TempEncounterStop7(2:end),false];
           TempEncounterStop9 = [TempEncounterStop8(2:end),false];
                      
           TempEncounterStop = TempEncounterStop1 & TempEncounterStop2 & TempEncounterStop3 & TempEncounterStop4 & TempEncounterStop5 & TempEncounterStop6 & TempEncounterStop7 & TempEncounterStop8 & TempEncounterStop9;
           
           data(i).perFlyData(j).encounterCentralStopTime = nnz(TempEncounterStop)/fps;
                      
       else
           data(i).perFlyData(j).encounterSmoothVel = [];
           data(i).perFlyData(j).encounterSmoothAngVel = [];
           data(i).perFlyData(j).encounterCentralStopTime=[];
       end
    end
    data(i).catEncounterVel = cat(1,data(i).perFlyData.encounterSmoothVel);
    data(i).catEncounterMedianVel = median(data(i).catEncounterVel,2);
    
    data(i).catEncounterAngVel = cat(1,data(i).perFlyData.encounterSmoothAngVel);
    data(i).catEncounterMedianAngVel = median(data(i).catEncounterAngVel,2);
    
    data(i).catEncounterFrame = cat(1,data(i).perFlyData.encounterFrame)';
    data(i).catPostEncounterFrames = cat(1,data(i).perFlyData.postEncounterFrames)';
    
    data(i).catEncounterCentralStopTime = cat(2,data(i).perFlyData.encounterCentralStopTime)';
    data(i).catEncounterCentralVel = cat(2,data(i).perFlyData.encounterCentralMedianVel)';
    data(i).catEncounterCentralAngVel = cat(2,data(i).perFlyData.encounterCentralMedianAngVel)';
    data(i).encounterNumFlies = size(data(i).catEncounterVel,1);
end



%% plot blob encounter vel

gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(1).catEncounterVel)',std(data(1).catEncounterVel)/sqrt(data(1).numFlies)','k','alpha');
p2=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(2).catEncounterVel)',std(data(2).catEncounterVel)/sqrt(data(2).numFlies)','b','alpha');
p3=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(3).catEncounterVel)',std(data(3).catEncounterVel)/sqrt(data(3).numFlies)','g','alpha');
p4=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(4).catEncounterVel)',std(data(4).catEncounterVel)/sqrt(data(4).numFlies)','r','alpha');
%p4=boundedline((1:maxLastFrameNumber)',mean(data(4).catSmoothVel)',std(data(4).catSmoothVel)/sqrt(data(4).numFlies)','r','alpha');
legend([p1,p2,p3,p4],gen);
hold off;
%% plot blob encounter angvel

gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(1).catEncounterAngVel)',std(data(1).catEncounterAngVel)/sqrt(data(1).numFlies)','k','alpha');
p2=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(2).catEncounterAngVel)',std(data(2).catEncounterAngVel)/sqrt(data(2).numFlies)','b','alpha');
p3=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(3).catEncounterAngVel)',std(data(3).catEncounterAngVel)/sqrt(data(3).numFlies)','g','alpha');
p4=boundedline((-preEncounterWindow/fps:1/fps:encounterWindow/fps)',mean(data(4).catEncounterAngVel)',std(data(4).catEncounterAngVel)/sqrt(data(4).numFlies)','r','alpha');
%p4=boundedline((1:maxLastFrameNumber)',mean(data(4).catSmoothVel)',std(data(4).catSmoothVel)/sqrt(data(4).numFlies)','r','alpha');
legend([p1,p2,p3,p4],gen);
hold off;


%% Generate data matrix for plotting in GraphPad

labels = {data.Genotype};

% Central-Zone astopping time during 5s after finding the blob
EncounterCentralStopTime=cat(2,{data.catEncounterCentralStopTime});
maxLength = max(cellfun(@length,EncounterCentralStopTime));
EncounterCentralStopTimeMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),EncounterCentralStopTime,'UniformOutput',false));
figure;errorbarjitter(EncounterCentralStopTimeMat);
set(gca, 'XTicklabel',labels),ylabel('encounter central stop duration');

% Central-Zone translational Velocity during 5s after finding the blob
EncounterCentralVel=cat(2,{data.catEncounterCentralVel});
maxLength = max(cellfun(@length,EncounterCentralVel));
EncounterCentralVelMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),EncounterCentralVel,'UniformOutput',false));
figure;errorbarjitter(EncounterCentralVelMat);
set(gca, 'XTicklabel',labels),ylabel('encounter trans velocity');

% Central-Zone angular Velocity during 5s after finding the blob
EncounterCentralAngVel=cat(2,{data.catEncounterCentralAngVel});
maxLength = max(cellfun(@length,EncounterCentralAngVel));
EncounterCentralAngVelMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),EncounterCentralAngVel,'UniformOutput',false));
figure;errorbarjitter(EncounterCentralAngVelMat);
set(gca, 'XTicklabel',labels),ylabel('encounter angular velocity');

%% Create blob encounter data  for heat map
%Here we increased the pre-encounter and post-encounter time window for
%showing longer time-scale data. This removes a few data points compared to
%previous plots when the fly finds the food-blob too early (<20s from start) in the
%assay or too late (<33 s from end) in the assay.

lowerLimit = 600;
upperLimit = maxLastFrameNumber - 100;
encounterWindow = 1500;
preEncounterWindow = 600;

for i = 1:length(data)
    currGenotype = realSubFolders(i).name;
    currPath = [pwd '\' realSubFolders(i).name];
    for j=1:data(i).numFlies
       currFlyPath = data(i).perFlyData(j).path;
       encounterFrame = xlsread([currFlyPath, '\blob.xlsx']);
       
       data(i).perFlyData(j).encounterFrame = encounterFrame;
              
       if encounterFrame>lowerLimit && encounterFrame<upperLimit
           data(i).perFlyData(j).encounterSmoothVelHM = data(i).perFlyData(j).smoothVel(1,encounterFrame-preEncounterWindow:encounterFrame+encounterWindow);
           data(i).perFlyData(j).encounterSmoothAngVelHM = data(i).perFlyData(j).smoothAngVel(1,encounterFrame-preEncounterWindow:encounterFrame+encounterWindow);
       else
           data(i).perFlyData(j).encounterSmoothVelHM = [];
           data(i).perFlyData(j).encounterSmoothAngVelHM = [];
           
       end

    end
    data(i).catEncounterVelHM = cat(1,data(i).perFlyData.encounterSmoothVelHM);
    data(i).catEncounterAngVelHM = cat(1,data(i).perFlyData.encounterSmoothAngVelHM);
    data(i).encounterNumFliesHM = size(data(i).catEncounterVelHM,1);
end

%% Plot heat maps

% Plot translational velocity heatmap aligned to food encounter at 20s
fullSmoothVelData = cat(1,data(1).catEncounterVelHM,data(2).catEncounterVelHM,data(3).catEncounterVelHM,data(4).catEncounterVelHM);
labels = {data(1).Genotype,data(2).Genotype,data(3).Genotype,data(4).Genotype};
minVel = 0;
maxVel = 12;
n1 = cumsum([data.encounterNumFliesHM])';
n2 = ([data.encounterNumFliesHM]/2)';
n=n1-n2;
n1StringArr = cell(1,length(n1));
n1StringArr(:) = {'----'};
newN=vertcat(n,n1);
[newN,index]=sort(newN);
newLabels = horzcat(labels,n1StringArr);
newLabels = newLabels(index);

clims=[minVel,maxVel];
figure;
imagesc(fullSmoothVelData,clims);title('Translational Velocity in food-blob assay');
set(gca,'YTick',newN','YTickLabel',newLabels,'XTick',1:20*fps:length(fullSmoothVelData),'XTickLabel',0:20:length(fullSmoothVelData)/fps);
ylabel('Genotypes');xlabel('time');colorbar;


% Plot angular velocity heatmap aligned to food encounter at 20s
fullSmoothAngVelData = cat(1,data(1).catEncounterAngVelHM,data(2).catEncounterAngVelHM,data(3).catEncounterAngVelHM,data(4).catEncounterAngVelHM);
minAngVel = 0;
maxAngVel = 8;
n1 = cumsum([data.encounterNumFliesHM])';
n2 = ([data.encounterNumFliesHM]/2)';
n=n1-n2;
n1StringArr = cell(1,length(n1));
n1StringArr(:) = {'----'};
newN=vertcat(n,n1);
[newN,index]=sort(newN);
newLabels = horzcat(labels,n1StringArr);
newLabels = newLabels(index);

clims=[minAngVel,maxAngVel];
figure;
imagesc(fullSmoothAngVelData,clims);title('Angular Velocity in food-blob assay');
set(gca,'YTick',newN','YTickLabel',newLabels,'XTick',1:20*fps:length(fullSmoothAngVelData),'XTickLabel',0:20:length(fullSmoothAngVelData)/fps);
ylabel('Genotypes');xlabel('time');colorbar;

%%