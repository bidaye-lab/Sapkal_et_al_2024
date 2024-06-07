%% Analysis of free walking FlyTracker output from co-activation experiments (Fig2)

% This only works with speicific folder,subfolder arrangement as in this
% data

subFolders = dir(pwd);
realSubFolders = subFolders([subFolders.isdir]');
realSubFolders = realSubFolders(3:end);

% Initialise video temporal parameters
fps = 30; %frames per second
trialDuration = 20; %seconds
lightONduration = 10; %seconds
lightOFFduration = 10; %seconds (Except trial one which starts 30s into lightOFF, all other lightOFF are 60s
realLightON = 10;
realLightOFF = 50;
interTrialInterval = (realLightOFF-lightOFFduration)+(realLightON-lightONduration);
croppedTime = 10; % first 10s of video was cropped during pre-processing
startFrame = (realLightOFF-croppedTime-lightOFFduration)*fps;

trial1 = startFrame+1:startFrame+trialDuration*fps;
trial2 = trial1(end)+interTrialInterval*fps+1:trial1(end)+(interTrialInterval+trialDuration)*fps;
trial3 = trial2(end)+interTrialInterval*fps+1:trial2(end)+(interTrialInterval+trialDuration)*fps;
trial4 = trial3(end)+interTrialInterval*fps+1:trial3(end)+(interTrialInterval+trialDuration)*fps;
trial5 = trial4(end)+interTrialInterval*fps+1:trial4(end)+(interTrialInterval+trialDuration)*fps;

%light OFF period within each trial
OFF1 = trial1(1):trial1(1)-1+(trialDuration-lightONduration)*fps;%This is defined as lightOFF period within trial1 duration
OFF2 = trial2(1):trial2(1)-1+(trialDuration-lightONduration)*fps;
OFF3 = trial3(1):trial3(1)-1+(trialDuration-lightONduration)*fps;
OFF4 = trial4(1):trial4(1)-1+(trialDuration-lightONduration)*fps;
OFF5 = trial5(1):trial5(1)-1+(trialDuration-lightONduration)*fps;

%light ON period within each trial
ON1 = OFF1(end)+1:trial1(end);
ON2 = OFF2(end)+1:trial2(end);
ON3 = OFF3(end)+1:trial3(end);
ON4 = OFF4(end)+1:trial4(end);
ON5 = OFF5(end)+1:trial5(end);


% get frame numbers corresponding to each stimulation period
maxLastFrameNumber = ON5(end);
ONframes = zeros(1,maxLastFrameNumber);
ONframes(ON1)=1;
ONframes(ON2)=1;
ONframes(ON3)=1;
ONframes(ON4)=1;
ONframes(ON5)=1;

%Initialize smoothing parameters for velocity attributes
smoothWindow = 15; %frames or 0.5s
angSmoothWindow= 15; %frames or 0.5s


%for Analysis only during part of light on.
onStart = 1;
onEnd = 2*fps;


%% Generate data structure containing all raw tracking data and derived paramaters for plotting
data={};

% Collect per genotype data in the following for loop
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
        %andVel arrays into consistent size as manually provided in variable
        %maxLastFrameNumber. flag fff (commented out) alerts user of these cases
        
        if length(tempVel)<maxLastFrameNumber
           %fff=222         
           currFeatLastFrame = length(tempVel);
           tempVel(end+1:maxLastFrameNumber)=tempVel(currFeatLastFrame);
           tempAngVel(end+1:maxLastFrameNumber)=tempAngVel(currFeatLastFrame);
        else
           %fff=111 
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
        
        timetags = ([1:length(tempVel)])/fps;%convert frames to time in s
        
        
        data(i).perFlyData(j).feat = tempFeature;
        data(i).perFlyData(j).trk = tempTrk;
        data(i).perFlyData(j).ppm = tempCal.PPM;
        
        %The following if statement converts manually annotated backward walking time periods as
        %negative velocity time points
        if exist([currFlyPath, '\', genotypeSubFolders(j).name, '-actions.mat'],'file')% behavior annotations done using visualizer function from flytracker package are stored in file ending in -actions.mat
                
           tempMatFileActions = matfile([currFlyPath, '\', genotypeSubFolders(j).name, '-actions.mat']);
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
           tempVel(condBack>0) = abs(tempVel(condBack>0))*-1;
        end
        
        
        data(i).perFlyData(j).vel = tempVel;%translational velocity
        data(i).perFlyData(j).angVel = tempAngVel;%absolute angular velocity
        data(i).perFlyData(j).ori = tempOri;%orientation (heading)
        data(i).perFlyData(j).PosX = tempPosX;% X-position
        data(i).perFlyData(j).PosY = tempPosY;% Y- position
        
        % generate smoothed parameters to overcome tracking jitter
        data(i).perFlyData(j).smoothVel = (smooth(tempVel,smoothWindow))';
        smoothVel =data(i).perFlyData(j).smoothVel;

 
        % per trial velocity parameters
        data(i).perFlyData(j).smoothVelTrial1 = data(i).perFlyData(j).smoothVel(trial1);
        data(i).perFlyData(j).smoothVelTrial2 = data(i).perFlyData(j).smoothVel(trial2);
        data(i).perFlyData(j).smoothVelTrial3 = data(i).perFlyData(j).smoothVel(trial3);
        data(i).perFlyData(j).smoothVelTrial4 = data(i).perFlyData(j).smoothVel(trial4);
        data(i).perFlyData(j).smoothVelTrial5 = data(i).perFlyData(j).smoothVel(trial5);
        
        % all trials pooled
        data(i).perFlyData(j).smoothVelAllTrials = cat(1,data(i).perFlyData(j).smoothVelTrial1,data(i).perFlyData(j).smoothVelTrial2,data(i).perFlyData(j).smoothVelTrial3,data(i).perFlyData(j).smoothVelTrial4,data(i).perFlyData(j).smoothVelTrial5);
        
        % avergage trial parameter
        data(i).perFlyData(j).smoothVelAvgTrial = mean(data(i).perFlyData(j).smoothVelAllTrials);
        
        % do same for angular velocity
        data(i).perFlyData(j).smoothAngVel = (smooth(tempAngVel,angSmoothWindow))';
        smoothAngVel = data(i).perFlyData(j).smoothAngVel;
        data(i).perFlyData(j).smoothAngVelTrial1 = data(i).perFlyData(j).smoothAngVel(trial1);
        data(i).perFlyData(j).smoothAngVelTrial2 = data(i).perFlyData(j).smoothAngVel(trial2);
        data(i).perFlyData(j).smoothAngVelTrial3 = data(i).perFlyData(j).smoothAngVel(trial3);
        data(i).perFlyData(j).smoothAngVelTrial4 = data(i).perFlyData(j).smoothAngVel(trial4);
        data(i).perFlyData(j).smoothAngVelTrial5 = data(i).perFlyData(j).smoothAngVel(trial5);
        
        data(i).perFlyData(j).smoothAngVelAllTrials = cat(1,data(i).perFlyData(j).smoothAngVelTrial1,data(i).perFlyData(j).smoothAngVelTrial2,data(i).perFlyData(j).smoothAngVelTrial3,data(i).perFlyData(j).smoothAngVelTrial4,data(i).perFlyData(j).smoothAngVelTrial5);
        data(i).perFlyData(j).smoothAngVelAvgTrial = mean(data(i).perFlyData(j).smoothAngVelAllTrials);
        
        % Calculate Distance covered during optogenetic stimulation.
        data(i).perFlyData(j).totalDistON1 = nansum(abs(tempVel(ON1(onStart:onEnd)))/fps);
        data(i).perFlyData(j).totalDistON2 = nansum(abs(tempVel(ON2(onStart:onEnd)))/fps);
        data(i).perFlyData(j).totalDistON3 = nansum(abs(tempVel(ON3(onStart:onEnd)))/fps);   
        data(i).perFlyData(j).totalDistON4 = nansum(abs(tempVel(ON4(onStart:onEnd)))/fps);
        data(i).perFlyData(j).totalDistON5 = nansum(abs(tempVel(ON5(onStart:onEnd)))/fps);
        data(i).perFlyData(j).totalDistAllON = cat(1,data(i).perFlyData(j).totalDistON1,data(i).perFlyData(j).totalDistON2,data(i).perFlyData(j).totalDistON3,data(i).perFlyData(j).totalDistON4,data(i).perFlyData(j).totalDistON5);
        data(i).perFlyData(j).totalDistAvgON = mean(data(i).perFlyData(j).totalDistAllON);
   
        % Calculate Rotation in optostim period, as area under curve of absolute angular velocity
        data(i).perFlyData(j).areaAllTurnsON1 = trapz(timetags(ON1(onStart:onEnd)),smoothAngVel(ON1(onStart:onEnd)));
        data(i).perFlyData(j).areaAllTurnsON2 = trapz(timetags(ON2(onStart:onEnd)),smoothAngVel(ON2(onStart:onEnd)));
        data(i).perFlyData(j).areaAllTurnsON3 = trapz(timetags(ON3(onStart:onEnd)),smoothAngVel(ON3(onStart:onEnd)));
        data(i).perFlyData(j).areaAllTurnsON4 = trapz(timetags(ON4(onStart:onEnd)),smoothAngVel(ON4(onStart:onEnd)));
        data(i).perFlyData(j).areaAllTurnsON5 = trapz(timetags(ON5(onStart:onEnd)),smoothAngVel(ON5(onStart:onEnd)));
        data(i).perFlyData(j).areaAllTurnsAllON = cat(1,data(i).perFlyData(j).areaAllTurnsON1,data(i).perFlyData(j).areaAllTurnsON2,data(i).perFlyData(j).areaAllTurnsON3,data(i).perFlyData(j).areaAllTurnsON4,data(i).perFlyData(j).areaAllTurnsON5);
        data(i).perFlyData(j).areaAllTurnsAvgON = mean(data(i).perFlyData(j).areaAllTurnsAllON);

        
    end
    
    % generate concatenated matrices for easy plotting
    data(i).catVel = cat(1,data(i).perFlyData.vel);   
    data(i).catSmoothVel = cat(1,data(i).perFlyData.smoothVel);
    data(i).catSmoothVelAllTrials = cat(1,data(i).perFlyData.smoothVelAllTrials);
    data(i).catSmoothVelAvgTrial = cat(1,data(i).perFlyData.smoothVelAvgTrial);
    
    data(i).catAngVel = cat(1,data(i).perFlyData.angVel);
    data(i).catSmoothAngVel = cat(1,data(i).perFlyData.smoothAngVel);
    data(i).catSmoothAngVelAllTrials = cat(1,data(i).perFlyData.smoothAngVelAllTrials);
    data(i).catSmoothAngVelAvgTrial = cat(1,data(i).perFlyData.smoothAngVelAvgTrial);

    data(i).catTotalDistAllON = cat(1,data(i).perFlyData.totalDistAllON);
    data(i).catTotalDistAvgON = cat(1,data(i).perFlyData.totalDistAvgON);

    data(i).catAreaAllTurnsAllON = cat(1,data(i).perFlyData.areaAllTurnsAllON);
    data(i).catAreaAllTurnsAvgON = cat(1,data(i).perFlyData.areaAllTurnsAvgON);

end


%% Create Distance file for plotting in Graphpad
totalDistAvgONData=cat(2,{data.catTotalDistAvgON});
maxLength = max(cellfun(@length,totalDistAvgONData));
totalDistAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),totalDistAvgONData,'UniformOutput',false));
figure;errorbarjitter(totalDistAvgONMat);

%% Create total rotation file for plotting in Graphpad
areaAllTurnsAvgONData=cat(2,{data.catAreaAllTurnsAvgON});
maxLength = max(cellfun(@length,areaAllTurnsAvgONData));
areaAllTurnsAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),areaAllTurnsAvgONData,'UniformOutput',false));
figure;errorbarjitter(areaAllTurnsAvgONMat);

%% plot bounded line plots of translational velocity
gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(1)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(1)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(1)))).numFlies)','k','alpha');
p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(2)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(2)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(2)))).numFlies)','b','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(3)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(3)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(3)))).numFlies)','g','alpha');
p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(4)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(4)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(4)))).numFlies)','m','alpha');
p5=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(5)))).catSmoothVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(5)))).catSmoothVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(5)))).numFlies)','r','alpha');

ylim = [-10,25];
xlim = [0,trialDuration*fps];
line([lightOFFduration*fps+1,lightOFFduration*fps+1],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle','--','Color','k');
set(gca,'YTick', min(ylim):5:max(ylim),'XTick',1:5*fps:max(xlim),'XTickLabel',0:5:max(xlim)/fps);
ylabel('Avg Trial Velocity +-SEM (mm/s)');xlabel('time (s)');legend([p1,p2,p3,p4,p5],[gen(1:5)]);
axis([0,trialDuration*fps, -12, 25]);

%% plot single fly translational velocities

figure;hold on;
plotGen = 1;%chose from 1 to 5 corresponding to the genotype for plotting
p1=plot((1:trialDuration*fps)',data(plotGen).catSmoothVelAllTrials,'LineWidth',0.5)';

ylim=[-18,35];
line([lightOFFduration*fps+1,lightOFFduration*fps+1],ylim,'LineStyle',':','Color','k');
%line([12*fps+1,12*fps+1],ylim,'LineStyle','--','Color','k');
line(xlim,[0,0],'LineStyle',':','Color','k');
set(gca,'YTick', -15:5:max(ylim),'XTick',1:1*fps:max(xlim),'XTickLabel',0:1:max(xlim)/fps);
ylabel('Avg Trial Velocity (mm/s)');xlabel('time (s)');%legend([p1(1),p2(1),p3(1)],{'left','right','both'});
axis([0,trialDuration*fps, -18, 35]);


%% plot bounded line plots of angular velocity
gen={data.Genotype};
figure;hold on;%plot average of all trials per genotype and std error bounds


p1=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(1)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(1)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(1)))).numFlies)','k','alpha');
p2=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(2)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(2)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(2)))).numFlies)','b','alpha');
p3=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(3)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(3)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(3)))).numFlies)','g','alpha');
p4=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(4)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(4)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(4)))).numFlies)','m','alpha');
p5=boundedline((1:trialDuration*fps)',mean(data(strcmp({data.Genotype},char(gen(5)))).catSmoothAngVelAvgTrial)',std(data(strcmp({data.Genotype},char(gen(5)))).catSmoothAngVelAvgTrial)/sqrt(data(strcmp({data.Genotype},char(gen(5)))).numFlies)','r','alpha');

ylim = [0,6];
xlim = [0,trialDuration*fps];
line([10*fps,10*fps],ylim,'LineStyle','--','Color','k');
%line(xlim,[0,0],'LineStyle',':','Color','k');
set(gca,'YTick', min(ylim):2:max(ylim),'XTick',1:5*fps:max(xlim),'XTickLabel',0:5:max(xlim)/fps);
ylabel('Avg Trial Angular Velocity +-SEM (rad/s)');xlabel('time (s)');legend([p1,p2,p3,p4,p5],[gen(1:5)]);
axis([0,trialDuration*fps, 0, 6]);

%% for P9 co-activation experiments, define pivots as on the spot rotations
turnThreshold = 2;%r/s
pivotForwardThreshold = 5;%mm/s
pivotBackwardThreshold = -2;%mm/s

for i = 1:length(realSubFolders)
    for j = 1:data(i).numFlies
        noWalk = data(i).perFlyData(j).smoothVel < pivotForwardThreshold & data(i).perFlyData(j).smoothVel>pivotBackwardThreshold;
        turn = data(i).perFlyData(j).smoothAngVel>turnThreshold;
        pivot = turn & noWalk;
        data(i).perFlyData(j).pivot = pivot;
        data(i).perFlyData(j).pivotON1 = nnz(pivot(ON1(onStart:onEnd)))/length(ON1(onStart:onEnd));
        data(i).perFlyData(j).pivotON2 = nnz(pivot(ON2(onStart:onEnd)))/length(ON2(onStart:onEnd));
        data(i).perFlyData(j).pivotON3 = nnz(pivot(ON3(onStart:onEnd)))/length(ON3(onStart:onEnd));
        data(i).perFlyData(j).pivotON4 = nnz(pivot(ON3(onStart:onEnd)))/length(ON4(onStart:onEnd));
        data(i).perFlyData(j).pivotON5 = nnz(pivot(ON3(onStart:onEnd)))/length(ON5(onStart:onEnd));
        data(i).perFlyData(j).pivotAllON = cat(1,data(i).perFlyData(j).pivotON1,data(i).perFlyData(j).pivotON2,data(i).perFlyData(j).pivotON3,data(i).perFlyData(j).pivotON4,data(i).perFlyData(j).pivotON5);
        data(i).perFlyData(j).pivotAvgON = mean(data(i).perFlyData(j).pivotAllON);

    end

    data(i).catPivot = cat(1,data(i).perFlyData.pivot);
    data(i).catPivotAllON = cat(1,data(i).perFlyData.pivotAllON);
    data(i).catPivotAvgON = cat(1,data(i).perFlyData.pivotAvgON);

end

%% For P9 co-activation experiments, create Pivot file for plotting in Graphpad
pivotAvgONData=cat(2,{data.catPivotAvgON});
maxLength = max(cellfun(@length,pivotAvgONData));
pivotAvgONMat= cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),pivotAvgONData,'UniformOutput',false));
figure;errorbarjitter(pivotAvgONMat);
%%




