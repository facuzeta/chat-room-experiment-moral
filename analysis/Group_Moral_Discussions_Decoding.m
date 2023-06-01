%% Group Moral Discussions Experiment - Decoding

% This script reads the .csv files containing data from our participants
% and uses it to create a decoding model for predicting which groups will
% reach consensus using linguistic variables.

% The required files are dataS1.csv, moralConsensus.csv, and 
% features_moral_times.csv. For the latter, a custom automatically 
% generated MATLAB function is also required to properly import the data.
% This function is importfileMoralTimes.m. Additionally, the MATLAB
% function hline.m is required.


%% Loading the Required Files

clearvars; close all; clc;

a=importdata('dataS1.csv');
dataFull=a.data;

c=importdata('moralConsensus.csv');
consensus=c.data;
clear a c

featuresMoralTimes = importfileMoralTimes('features_moral_times.csv', 2, 30460);
% This file contains the number of words, word length, number of messages,
% and other useful variables 
%% Creating the Necessary Variables for Decoding

% This is the maximum number of messages (interventions) that occurred in a
% conversation. This was obtained on preliminary analysis.
maxLength=79;

% Initializations:

Consensus={[],[],[],[],[],[]};
Ranges={[],[],[],[],[],[]};
% Values chosen from the final number of responses for each question:
allTimes={nan(172,maxLength),nan(147,maxLength),nan(174,maxLength),nan(166,maxLength),nan(177,maxLength),nan(176,maxLength)};
allWordTimes={nan(172,maxLength),nan(147,maxLength),nan(174,maxLength),nan(166,maxLength),nan(177,maxLength),nan(176,maxLength)};
allLengthTimes={nan(172,maxLength),nan(147,maxLength),nan(174,maxLength),nan(166,maxLength),nan(177,maxLength),nan(176,maxLength)};
% Real durations are determined by checking when in the conversations
% the consensus was actually sent to the platform (participants could
% have continued talking after agreeing and sending their consensus; we
% use this to ignore anything that occurrs after it has been received).
allRealDurations={[],[],[],[],[],[]};

% End of initializations.

% Loop over questions:
for q=1:6
    groups=consensus(consensus(:,2)==q+8,1);
    k=1;
    % Loop over groups:
    for group=groups'
        % Check if we have data for Stage 2 for this group, and all its
        % members have completed at least Stages 1 and 2:
        if ismember(group,featuresMoralTimes.group_id) && length(dataFull(dataFull(:,7)==group,q))==3
                        
            %Real Times for this group and question:
            TimesGroup=featuresMoralTimes.seconds(featuresMoralTimes.group_id==group & featuresMoralTimes.question_id==q+8);
            %Real duration of the conversation:
            RealDuration=min(sum(featuresMoralTimes.consensus_reached(featuresMoralTimes.group_id==group & featuresMoralTimes.question_id==q+8)=='False')+1,length(featuresMoralTimes.consensus_reached(featuresMoralTimes.group_id==group & featuresMoralTimes.question_id==q+8)));      
            %Number of Words sent in each message:
            WordsTimeEvolutions=featuresMoralTimes.n_words(featuresMoralTimes.group_id==group & featuresMoralTimes.question_id==q+8);
            %Average Length per Word sent in each message:
            AverageLengthTimeEvolutions=featuresMoralTimes.words_length__mean(featuresMoralTimes.group_id==group & featuresMoralTimes.question_id==q+8);
            
            %Consensus:
            consensusQuestions=consensus(consensus(:,1)==group,2);
            consensusAnswers=consensus(consensus(:,1)==group,3);
            actualConsensus=consensusAnswers(consensusQuestions==q+8)~=-1;
            
            %Initial Range of Opinions:
            opinions=dataFull(dataFull(:,7)==group,q);
            internalRange=max(opinions)-min(opinions);
            
            %Group Members:
            [integrantes,order]=sort(dataFull(dataFull(:,7)==group,q));
            integrantesIds=dataFull(dataFull(:,7)==group,14);
            integrantesIds=integrantesIds(order);
            
            %Save actual consensus
            Consensus{q}=[Consensus{q},actualConsensus];
            %Save internal range
            Ranges{q}=[Ranges{q},internalRange];

            %Save all real times
            allTimes{q}(k,1:length(TimesGroup))=TimesGroup;
            %save all evolutions of words
            allWordTimes{q}(k,1:length(TimesGroup))=WordsTimeEvolutions;
            %save all real durations
            allRealDurations{q}=[allRealDurations{q},RealDuration];
            %save all evolution of mean lengths
            allLengthTimes{q}(k,1:length(TimesGroup))=AverageLengthTimeEvolutions;
            k=k+1;
            
        end %if condition end
    end %end loop over groups
end %end loop over questions

% Clear all unnecesary variables:
clearvars -except allLengthTimes allRealDurations allWordTimes allTimes Consensus Ranges


%% Decoding Model

% We randomly assign 75% of the data to the training set, and 25% of the
% data to the testing set. We perform many repetitions, since this
% introduces randomness into the modeling process.

% Number of Repetitions:
maxReps=10; %beware: big numbers mean long running times; 1000 is relatively time consuming.

%Choose the sliding window size: 
windowSize=10;

%Choose the sliding step:
step=1;

% Initialization:
K=1; %step counter (timelimit)

for timelimit=0:step:300-windowSize %Sliding Windows; 300 is the maximum time.

timelimit %to be able to tell how the simulation is going.

% Repetition loop:
for repetition=1:maxReps
    
    %clear unnecesary variables. Some may be overwritten in repetitions, so it is also a safety measure.
    clearvars -except timestart timelimit windowSize repetition allRealDurations allWordTimes allLengthTimes allTimes Consensus Ranges maxReps K successesBaseRate successesFullModel successesSimpleModel
    
    %Choose a random order for the data
    order = Shuffle(1:length([Ranges{:}]));
   
    %Keep 75% of the data (change the 0.75 to other values to modify this)
    percentile = floor(0.75*length([Ranges{:}]));
    
    %save all ranges as a vector; create training set.
    rangesAll=[Ranges{:}];
    rangesTraining=rangesAll(order(1:percentile));
       
    %save all consensus as a vector; create training set.    
    ConsensusAll=[Consensus{:}];
    ConsensusTraining=ConsensusAll(order(1:percentile));  
   
    %save all durations as a vector; create training set. 
    DurationsAll=[allRealDurations{:}];
    DurationsTraining=DurationsAll(order(1:percentile));
    
    %save all number of words evolutions as a matrix; create training set. 
    WordsEvolutionAll=[allWordTimes{1};allWordTimes{2};allWordTimes{3};allWordTimes{4};allWordTimes{5};allWordTimes{6}];
    WordsEvolutionTraining=WordsEvolutionAll(order(1:percentile),:);

    %save all word length evolutions as a matrix; create training set. 
    LengthEvolutionAll=[allLengthTimes{1};allLengthTimes{2};allLengthTimes{3};allLengthTimes{4};allLengthTimes{5};allLengthTimes{6}];
    LengthEvolutionTraining=LengthEvolutionAll(order(1:percentile),:);   
    
    %save all times as a matrix; create training set. 
    TimesGroup=[allTimes{1};allTimes{2};allTimes{3};allTimes{4};allTimes{5};allTimes{6}];
    AllTimesTraining=TimesGroup(order(1:percentile),:);    
    
  
    timesQuest=AllTimesTraining; %temporary
    lengthsTimesSubjects=LengthEvolutionTraining; %temporary
    wordsTimesSubjects=WordsEvolutionTraining; %temporary
    durationsQuest=DurationsTraining; %temporary
    
    % We measure the number of words per message and mean word length 
    % inside the chosen windows ("binned"). Initialize:
    BinnedWordsTraining=[];
    BinnedLengthsTraining=[];
    
    % Times within the window:
    timesChosen=(timesQuest>=timelimit) & (timesQuest<timelimit+windowSize);
    
    % Eliminate messages after the real duration:
    for i=1:length(durationsQuest)
        wordsTimesSubjects(i,durationsQuest(i)+1:end)=nan;
        lengthsTimesSubjects(i,durationsQuest(i)+1:end)=nan;
    end
    
    % Saved binned results:
    BinnedWordsTraining=[BinnedWordsTraining,nansum(wordsTimesSubjects.*timesChosen,2)./nansum(timesChosen,2)];
    BinnedLengthsTraining=[BinnedLengthsTraining,nansum(lengthsTimesSubjects.*timesChosen,2)./nansum(timesChosen,2)];

    %Normalization (for use in logistic regression); we also eliminate nans:
    rangesTrainingNorm=zscore(rangesTraining(~isnan(BinnedWordsTraining)));
    BinnedLengthsNorm=zscore(BinnedLengthsTraining(~isnan(BinnedWordsTraining)));        
    BinnedWordsNorm=zscore(BinnedWordsTraining(~isnan(BinnedWordsTraining)));        
    
    newConsensusTraining=ConsensusTraining(~isnan(BinnedWordsTraining));
    
    % Fit models:
    
    % Full model:
    mdl = fitglm([rangesTrainingNorm',BinnedWordsNorm,BinnedLengthsNorm],newConsensusTraining','Distribution','binomial');    
    
    % Simple model:
    mdl2 = fitglm(rangesTrainingNorm',newConsensusTraining','Distribution','binomial');
    
    % Save full model coefficients:
    B=mdl.Coefficients.Estimate;
    paramIntercept=B(1);
    paramRango=B(2);
    paramDensidad=B(3);
    paramLongitud=B(4);
    
    % Save simple model coefficients:
    C=mdl2.Coefficients.Estimate;
    paramIntercept2=C(1);
    paramRango2=C(2);  
    
    % Create testing sets:
    rangesTesting=rangesAll(order(percentile+1:length(rangesAll)));
    ConsensusTesting=ConsensusAll(order(percentile+1:length(ConsensusAll)));
    DurationsTesting=DurationsAll(order(percentile+1:length(DurationsAll)));
    WordEvolutionTesting=WordsEvolutionAll(order(percentile+1:length(DurationsAll)),:);
    allTimesTesting=TimesGroup(order(percentile+1:length(rangesAll)),:);
    LengthEvolutionTesting=LengthEvolutionAll(order(percentile+1:length(DurationsAll)),:);
    
    %Temporary variables:
    timesQuest=allTimesTesting;
    wordsTimesSubjects=WordEvolutionTesting;
    lengthsTimesSubjects=LengthEvolutionTesting;
    durationsQuest=DurationsTesting;   
        
    % We measure the number of words per message and mean word length 
    % inside the chosen windows ("binned"). Initialize:    
    BinnedWordsTesting=[];
    BinnedLengthsTesting=[];
    
    % Times within the window:    
    timesChosen=(timesQuest>=timelimit) & (timesQuest<timelimit+windowSize); %Ventanas Móviles 10 
    
    % Eliminate messages after the real duration:    
    for i=1:length(durationsQuest)
        wordsTimesSubjects(i,durationsQuest(i)+1:end)=nan;  
        lengthsTimesSubjects(i,durationsQuest(i)+1:end)=nan;
    end

    % Saved binned results:
    BinnedWordsTesting=[BinnedWordsTesting,nansum(wordsTimesSubjects.*timesChosen,2)./nansum(timesChosen,2)];
    BinnedLengthsTesting=[BinnedLengthsTesting,nansum(lengthsTimesSubjects.*timesChosen,2)./nansum(timesChosen,2)];
    
    % Normalization (for use in logistic regression); we also eliminate nans:
    rangesTestingNorm=zscore(rangesTesting(~isnan(BinnedWordsTesting)));
    BinnedWordsTestingNorm=zscore(BinnedWordsTesting(~isnan(BinnedWordsTesting)));  
    BinnedLengthsTestingTestingNorm=zscore(BinnedLengthsTesting(~isnan(BinnedWordsTesting)));    
    
    newConsensusTesting=ConsensusTesting(~isnan(BinnedWordsTesting));

    % Use parameter values from the fit to estimate the variables
    x01=paramIntercept2+paramRango2.*rangesTestingNorm;
    x11=paramIntercept + paramRango.*rangesTestingNorm +paramDensidad.*BinnedWordsTestingNorm' + paramLongitud.*BinnedLengthsTestingTestingNorm';
    x00=ones(1,length(rangesTestingNorm));    
    
    % Measure successes in testing set through the previous variables
    successesSimpleModel(repetition,K)=sum(round(1./(exp(-x01)+1))==(newConsensusTesting))./length(newConsensusTesting);
    successesFullModel(repetition,K)=sum(round(1./(exp(-x11)+1))==(newConsensusTesting))./length(newConsensusTesting);       
    successesBaseRate(repetition,K)=sum(round(1./(exp(x00)+1))==(newConsensusTesting))./length(newConsensusTesting);
    
    % The Base Rate is always the best possible one from predicting always
    % "consensus" or always "disensus": 
    if successesBaseRate(repetition,K)<0.5
        successesBaseRate(repetition,K)=1-successesBaseRate(repetition,K);
    end   
    
end %end of repetitions loop
    
K=K+1; %step counter
end
'finished' %so that we know it is done.
%mean(successesFullModel-successesSimpleModel)

%% Plot: Decoding Model

points=0:300-windowSize;
pointsy=mean(successesFullModel-successesSimpleModel);
figure; 
hold on
plot(points,pointsy,'ko');
hline(0,'k--') %this function creates an horizontal line.
plot(points(sum(successesFullModel>successesSimpleModel)>950),pointsy(sum(successesFullModel>successesSimpleModel)>950),'ro');
xlabel('Seconds')
ylabel('Mean Gain (%)')
ylim([-0.02,0.130001])
yticklabels([-2:2:13])
title('Full Model vs. Range Only')