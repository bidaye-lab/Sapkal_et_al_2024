function [BOUT_ENDS_AND_BEGINNINGS_indices,BOUT_ENDS_AND_BEGINNINGS]=GET_FEEDING_BURSTS(Onsets,...
    diff_Onsets__or___IFI,Criterion_for_defining_regular_feeding,MINIMAL_NUMBER_OF_TOUCHES_IN_THE_BURSTS,Plot_raster_for_check_or_not)

%Onsets=ons1{ThisFly,Condition};
%IFI=IFIs1{ThisFly,Condition};

VARIABLE_TO_USE=diff_Onsets__or___IFI;%IFI  % ATTENTION!!!!!!!!!!!!!!


Criterion=Criterion_for_defining_regular_feeding; % ATTENTION!!!!!!!!!!!!!!
if Plot_raster_for_check_or_not==1    % CHECK ADEQUACY OF CRITERIA
    figure
    subplot(2,1,1)
    plot(Onsets,ones(size(Onsets)),'.k')
    hold all
end

long_INTERVAL_END_index=[find(VARIABLE_TO_USE>Criterion),length(Onsets)];

if Plot_raster_for_check_or_not==1    
    plot(Onsets(long_INTERVAL_END_index),ones(size(long_INTERVAL_END_index)),'.r')
end

long_INTERVAL_BEGINNING_index=[1,long_INTERVAL_END_index(1:end-1)+1];
if Plot_raster_for_check_or_not==1
    plot(Onsets(long_INTERVAL_BEGINNING_index),ones(size(long_INTERVAL_BEGINNING_index)),'.g')
    % xlim([8488 12670])
end

BOUT_ENDS_AND_BEGINNINGS_indices=[long_INTERVAL_BEGINNING_index',...
    long_INTERVAL_END_index',1+(long_INTERVAL_END_index'-long_INTERVAL_BEGINNING_index')];

BOUT_ENDS_AND_BEGINNINGS_indices=BOUT_ENDS_AND_BEGINNINGS_indices(BOUT_ENDS_AND_BEGINNINGS_indices(:,3)>=MINIMAL_NUMBER_OF_TOUCHES_IN_THE_BURSTS,:);

BOUT_ENDS_AND_BEGINNINGS=[(Onsets(BOUT_ENDS_AND_BEGINNINGS_indices(:,1)))',...
    ( Onsets(BOUT_ENDS_AND_BEGINNINGS_indices(:,2)) )',...
    BOUT_ENDS_AND_BEGINNINGS_indices(:,3)];
