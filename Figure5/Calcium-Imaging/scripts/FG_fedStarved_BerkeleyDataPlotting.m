%% read deltaF/F summary excel files from Gabriella Sterne for plotting and AUC calculations

Gr5a_fed = xlsread('Gr5a_fed.xlsx');
Gr5a_fed = Gr5a_fed(:,2:end);
Gr66a_fed = xlsread('Gr66a_fed.xlsx');
Gr66a_fed = Gr66a_fed(:,2:end);
Gr5a_starved = xlsread('Gr5a_starved.xlsx');
Gr5a_starved = Gr5a_starved(:,2:end);
numGenotypes = 3;
fps = 0.667; %Hz
lightONPeriod = 2; %seconds
stim1 = ceil(10/fps)+1:ceil(12/fps);
%%

AUC_Gr5aFed = trapz(0:1*fps:lightONPeriod*fps,Gr5a_fed(stim1,:));
AUC_Gr66aFed = trapz(0:1*fps:lightONPeriod*fps,Gr66a_fed(stim1,:));
AUC_Gr5aStarved = trapz(0:1*fps:lightONPeriod*fps,Gr5a_starved(stim1,:));

%%
stimFrame = 16;
plotWindow = 10;
minFrame = ceil(8/fps);
figure;
hold on;
p1=boundedline((-2:1*fps:10)',mean(Gr5a_starved(stimFrame-ceil(2/fps):stimFrame+floor(10/fps),:),2)',std(Gr5a_starved(stimFrame-ceil(2/fps):stimFrame+floor(10/fps),:),0,2)'/sqrt(length(AUC_Gr5aStarved)'),'g','alpha');
p2=boundedline((-2:1*fps:10)',mean(Gr5a_fed(stimFrame-ceil(2/fps):stimFrame+floor(10/fps),:),2)',std(Gr5a_fed(stimFrame-ceil(2/fps):stimFrame+floor(10/fps),:),0,2)'/sqrt(length(AUC_Gr5aFed)'),'c','alpha');
p3=boundedline((-2:1*fps:10)',mean(Gr66a_fed(stimFrame-ceil(2/fps):stimFrame+floor(10/fps),:),2)',std(Gr66a_fed(stimFrame-ceil(2/fps):stimFrame+floor(10/fps),:),0,2)'/sqrt(length(AUC_Gr66aFed)'),'k','alpha');
line([0,0],ylim,'LineStyle','--','Color','k');
hold off;
%%
plotWindow = 60;
figure;
hold on;
p1=boundedline((1:plotWindow+1)',mean(Gr5a_starved(1:plotWindow+1,:),2)',std(Gr5a_starved(1:plotWindow+1,:),0,2)'/sqrt(length(AUC_Gr5aStarved)'),'g','alpha');
p2=boundedline((1:plotWindow+1)',mean(Gr5a_fed(1:plotWindow+1,:),2)',std(Gr5a_fed(1:plotWindow+1,:),0,2)'/sqrt(length(AUC_Gr5aFed)'),'c','alpha');
p3=boundedline((1:plotWindow+1)',mean(Gr66a_fed(1:plotWindow+1,:),2)',std(Gr66a_fed(1:plotWindow+1,:),0,2)'/sqrt(length(AUC_Gr66aFed)'),'k','alpha');
hold off;
%%