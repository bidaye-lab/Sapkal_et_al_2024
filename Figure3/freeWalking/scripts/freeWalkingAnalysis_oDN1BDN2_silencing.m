%% Script to make struct array for ACR1 tracked data
% This only works with speicific folder,subfolder arrangement as in
% selected_60-30_powdered_ACR1
subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 30; %frames per second
trialDuration = 60; %seconds
lightONduration = 30; %seconds
lightOFFduration = 30; %seconds (Except trial one which starts 30s into lightOFF, all other lightOFF are 60s
realLightON = 30;
realLightOFF = 60;
croppedTime = 10;
interTrialInterval = (realLightOFF-lightOFFduration)+(realLightON-lightONduration);
startFrame = (realLightOFF-croppedTime-lightOFFduration)*fps;

trial1 = startFrame+1:startFrame+trialDuration*fps;
trial2 = trial1(end)+interTrialInterval*fps+1:trial1(end)+(interTrialInterval+trialDuration)*fps;
trial3 = trial2(end)+interTrialInterval*fps+1:trial2(end)+(interTrialInterval+trialDuration)*fps;

OFF1 = trial1(1):trial1(1)-1+(trialDuration-lightONduration)*fps;%This is defined as lightOFF period within trial1 duration
OFF2 = trial2(1):trial2(1)-1+(trialDuration-lightONduration)*fps;
OFF3 = trial3(1):trial3(1)-1+(trialDuration-lightONduration)*fps;

ON1 = OFF1(end)+1:trial1(end);%This is defined as lightON period within trial1 duration
ON2 = OFF2(end)+1:trial2(end);
ON3 = OFF3(end)+1:trial3(end);

OFFTrial1 = trial1(1):ON1(1)+trialDuration*fps -1;
OFFTrial2 = trial2(1):ON2(1)+trialDuration*fps -1;

maxLastFrameNumber = ON3(end);

%Initialize smoothing parameters for velocity attributes
smoothWindow = 15; %frames
angSmoothWindow= 15; %frames

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
    
    featureFiles = getAllExtFiles(currPath,'feat.mat');
    trackFiles = getAllExtFiles(currPath,'track.mat');
    calFiles = getAllExtFiles(currPath,'calibration.mat');
    
    %Collect per fly data in the following for loop
    for j=1:data(i).numFlies
        currFeatureFile = char(featureFiles(j));
        tempMatFileFeat = matfile(currFeatureFile);
        tempMatFileTrk = matfile(char(trackFiles(j)));
        tempMatFileCal = matfile(char(calFiles(j)));
        tempCal = tempMatFileCal.calib;
        if tempCal.FPS==50
            tempFeature = tempMatFileFeat.newFeat;
            tempTrk = tempMatFileTrk.newTrk;
            flag=1
        else
            tempFeature = tempMatFileFeat.feat;
            tempTrk = tempMatFileTrk.trk;
        end
        
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
        %andVel arrays into consistent size as manually provided in variable
        %maxLastFrameNumber. 
        
        if length(tempVel)<maxLastFrameNumber
                    
           currFeatLastFrame = length(tempVel);
           tempVel(end+1:maxLastFrameNumber)=tempVel(currFeatLastFrame);
           tempAngVel(end+1:maxLastFrameNumber)=tempAngVel(currFeatLastFrame);
        else
           tempVel = tempVel(1:maxLastFrameNumber);
           tempAngVel = tempAngVel(1:maxLastFrameNumber);
        end
        
    % Since Trk and Feat files created seperately do if statement once more
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
        
        timetags = ([1:length(tempVel)])/fps;
        
        
        data(i).perFlyData(j).feat = tempFeature;
        data(i).perFlyData(j).trk = tempTrk;
        data(i).perFlyData(j).ppm = tempCal.PPM;
        
        %The following if statement converts manually annotated backward walking time periods as
        %negative velocity time points
        if exist(strcat(currFeatureFile(1:length(currFeatureFile)-9),'-actions.mat'),'file')
                
           tempMatFileActions = matfile(strcat(currFeatureFile(1:length(currFeatureFile)-9),'-actions.mat'));
           tempbehs = tempMatFileActions.behs;
           tempbouts = tempMatFileActions.bouts;
           backBouts = cell2mat(tempbouts(1,1));
           backBouts = backBouts(:,[1 2],:);
           backBouts = backBouts/fps;
           tempCondArr = zeros(length(backBouts),length(timetags));
           for ctrtemp=(1:length(backBouts))
               tempCondArr(ctrtemp,:)=transpose(timetags>backBouts(ctrtemp,1) &timetags<backBouts(ctrtemp,2));
           end
           condBack = sum(tempCondArr,1);
           tempVel(condBack>0) = tempVel(condBack>0)*-1;
        end
        
        
        theta2 = tempOri(2:end);
        theta1 = tempOri(1:end-1);
        smooth_kernel = [1 2 1]/4;
        ori_diffRAW = mod(theta1+pi/2 - theta2,pi)-pi/2;
        tempAngVelSigned = [ori_diffRAW(1) ori_diffRAW]*fps;
        tempAngVelSigned(2:end-1) = conv(tempAngVelSigned,smooth_kernel,'valid');
        
            
        data(i).perFlyData(j).vel = tempVel;
        data(i).perFlyData(j).angVel = tempAngVel;
        data(i).perFlyData(j).angVelSigned = tempAngVelSigned;
        data(i).perFlyData(j).ori = tempOri;
        data(i).perFlyData(j).PosX = tempPosX;
        data(i).perFlyData(j).PosY = tempPosY;
        
        data(i).perFlyData(j).smoothVel = (smooth(tempVel,smoothWindow))';
 
        
        data(i).perFlyData(j).smoothVelTrial1 = data(i).perFlyData(j).smoothVel(trial1);
        data(i).perFlyData(j).medianVelTrial1 = median(data(i).perFlyData(j).smoothVelTrial1(:,901:1800),2);
        data(i).perFlyData(j).smoothVelTrial2 = data(i).perFlyData(j).smoothVel(trial2);
        data(i).perFlyData(j).medianVelTrial2 = median(data(i).perFlyData(j).smoothVelTrial2(:,901:1800),2);
        data(i).perFlyData(j).smoothVelTrial3 = data(i).perFlyData(j).smoothVel(trial3);
        data(i).perFlyData(j).medianVelTrial3 = median(data(i).perFlyData(j).smoothVelTrial3(:,901:1800),2);
        data(i).perFlyData(j).smoothVelAllTrials = cat(1,data(i).perFlyData(j).smoothVelTrial1,data(i).perFlyData(j).smoothVelTrial2,data(i).perFlyData(j).smoothVelTrial3);
        data(i).perFlyData(j).smoothVelAvgTrial = mean(data(i).perFlyData(j).smoothVelAllTrials);
        data(i).perFlyData(j).meanOfMedianVelTrial = mean([data(i).perFlyData(j).medianVelTrial1,data(i).perFlyData(j).medianVelTrial2,data(i).perFlyData(j).medianVelTrial3]);
        
        
        
    end
    data(i).catVel = cat(1,data(i).perFlyData.vel);
    
    data(i).catSmoothVel = cat(1,data(i).perFlyData.smoothVel);
    data(i).catSmoothVelAllTrials = cat(1,data(i).perFlyData.smoothVelAllTrials);
    data(i).catSmoothVelAvgTrial = cat(1,data(i).perFlyData.smoothVelAvgTrial);
    data(i).catMeanOfMedianVelTrial = cat(1,data(i).perFlyData.meanOfMedianVelTrial);
       
    
end

%% Create full data heat map
fullSmoothVelData = cat(1,data.catSmoothVel);
labels = {data.Genotype};
minVel = 0;
maxVel = 30;
n1 = cumsum([data.numFlies])';
n2 = ([data.numFlies]/2)';
n=n1-n2;
n1StringArr = cell(1,length(n1));
n1StringArr(:) = {'----'};
newN=vertcat(n,n1);
[newN,index]=sort(newN);
newLabels = horzcat(labels,n1StringArr);
newLabels = newLabels(index);

clims=[minVel,maxVel];
figure;
imagesc(fullSmoothVelData,clims);title('CsChrimson Activation');
 line([ON1(1),ON1(1)],ylim,'LineStyle','--','Color','r');
 line([ON1(end),ON1(end)],ylim,'LineStyle','--','Color','y');
 line([ON2(1),ON2(1)],ylim,'LineStyle','--','Color','r');
 line([ON2(end),ON2(end)],ylim,'LineStyle','--','Color','y');
 line([ON3(1),ON3(1)],ylim,'LineStyle','--','Color','r');
 line([ON3(end),ON3(end)],ylim,'LineStyle','--','Color','y');
set(gca,'YTick',newN','YTickLabel',newLabels,'XTick',1:30*fps:length(fullSmoothVelData),'XTickLabel',0:30:length(fullSmoothVelData)/fps);

ylabel('Genotypes');xlabel('time');colorbar;


%% Create Mean of Median per Trial Velcotiy file for plotting in Graphpad
meanOfMedianONData=cat(2,{data.catMeanOfMedianVelTrial});
maxLength = max(cellfun(@length,meanOfMedianONData));
meanOfMedianONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),meanOfMedianONData,'UniformOutput',false));
figure;errorbarjitter(meanOfMedianONMat);
%%