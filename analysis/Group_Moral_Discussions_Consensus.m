%% Group Moral Discussions Experiment - Probability of Consensus

% This script reads the .csv files containing data from our participants
% and uses it to explore the probability of reaching consensus in relation
% linguistic complexity of the exchanged messages. 

% The required files are dataS1.csv, and moralConsensus.csv.

%% Loading the Required Files

clearvars;close all;clc

a=importdata('dataS1.csv');
dataFull=a.data;

% load('DatosMoralS1ConPalabrasUpdated')

c=importdata('moralConsensus.csv');
consensus=c.data;
clear a c

%% Creating Plot 1 - Probability of Consensus as a Function of the Number of Words.

figure();
% We need two curves, one for low initial range of opinions (case 1), and 
% one for high initial range of opinions (case 2).
for cases=1:2
    if cases==1
        maxrange=5;
        minrange=0;
    elseif cases==2
        maxrange=10;
        minrange=6;
    end
    
    % Initializations
    consensusWordsGroup=[];
    consensusMessagesGroup=[];
    consensusWordsXMessageGroup=[];
    disensusWordsGroup=[];
    disensusMessagesGroup=[];
    disensusWordsXMessageGroup=[];
    
    % Loop over questions:
    for q=1:6
    groups=consensus(consensus(:,2)==q+8,1);
    % Loop over groups:
    for group=groups'
        opinions=dataFull(dataFull(:,7)==group,q);
        initialrange=max(opinions)-min(opinions);
        consensusQuestions=consensus(consensus(:,1)==group,2);
        consensusAnswers=consensus(consensus(:,1)==group,3);        
        
        totalWordsGroup=sum(dataFull(dataFull(:,7)==group,q+14));
        totalMessagesGroup=sum(dataFull(dataFull(:,7)==group,q+20));
        wordsXMessageGroup=totalWordsGroup./totalMessagesGroup;
               
        %If the group had 3 members, and they reached consensus:
        if length(opinions)>2 && consensusAnswers(consensusQuestions==q+8)>=0 
            if initialrange>=minrange && initialrange<=maxrange %we keep only the selected ranges
            consensusWordsGroup=[consensusWordsGroup,totalWordsGroup]; %save number of words
            consensusMessagesGroup=[consensusMessagesGroup,totalMessagesGroup]; %save number of messages
            consensusWordsXMessageGroup=[consensusWordsXMessageGroup,wordsXMessageGroup]; %save number of words per message
            end
        %If the group had 3 members, and they reached disensus:
        elseif length(opinions)>2 && consensusAnswers(consensusQuestions==q+8)==-1            
            if initialrange>=minrange && initialrange<=maxrange %we keep only the selected ranges
            disensusWordsGroup=[disensusWordsGroup,totalWordsGroup]; %save number of words
            disensusMessagesGroup=[disensusMessagesGroup,totalMessagesGroup]; %save number of messages
            disensusWordsXMessageGroup=[disensusWordsXMessageGroup,wordsXMessageGroup]; %save number of words per message
            end
        end
        
    end %ends loop over groups
    end %ends loop over questions

% Now we choose an initial value of number of words per message("density"),
% and a step, and count how many times a group had that density within that
% step, and so on.

initial=0.5;
k=1;
step=0.8;
counter=[]; %initialize
endDens=17; %last density found. Can be big for good measure.
for density=initial:step:endDens
    groupsConsensusInside=sum(consensusWordsXMessageGroup>=density & consensusWordsXMessageGroup<density+step);
    groupsDisensusInside=sum(disensusWordsXMessageGroup>=density & disensusWordsXMessageGroup<density+step);
    counter(k)=groupsConsensusInside./(groupsConsensusInside+groupsDisensusInside);
    k=k+1;
end

% We now plot the resulting curves; we separate both cases.
hold on
if cases==1
    counter
    plot(initial:step:endDens,counter,'o','color',[0, 0.5, 0],'LineWidth',1.5)
    mdl = fitglm(initial:step:endDens,counter,'Distribution','binomial'); %logistic fit
    params = mdl.Coefficients.Estimate; %save the fit parameters
    hold on
    % Plot the fit:
    plot(0:step:endDens,1./(exp(-(params(1)+params(2).*(0:step:endDens)))+1),'color',[0, 0.5, 0],'LineWidth',1.5)
elseif cases==2
    plot(0.5:step:endDens,counter,'o','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)
    mdl = fitglm(0.5:step:endDens,counter,'Distribution','binomial'); %logistic fit
    params = mdl.Coefficients.Estimate; %save the fit parameters
    hold on
    % Plot the fit:
    plot(0:step:endDens,1./(exp(-(params(1)+params(2).*(0:step:endDens)))+1),'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)        
end

xlabel('Number of Words per Message')
ylabel('P(Consensus)')
xticks(1:3:16)
yticks(0:0.2:1)
xlim([0,16.1])
legend('Range<=5','Fit','Range>5','Fit')

end %end of cases loop

%% Creating Plot 2 - Probability of Consensus as a Function of the Number of Words.

figure();
% We need two curves, one for low initial range of opinions (case 1), and 
% one for high initial range of opinions (case 2). More fine tuning of the
% initial value, step and ending density was required. 
for cases=1:2
    if cases==1
        maxrange=5;
        minrange=0;
        initial=2.7;
        endDens=5;        
        step=0.12;
    elseif cases==2
        maxrange=10;
        minrange=6;
        initial=2.7;
        endDens=5;
        step=0.11;
    end

    % Initializations
    consensusWordsLengthGroup=[];
    disensusWordsLengthGroup=[];
    
    % Loop over questions:
    for q=1:6
    groups=consensus(consensus(:,2)==q+8,1);
    % Loop over groups:
    for group=groups'
        opinions=dataFull(dataFull(:,7)==group,q);
        initialrange=max(opinions)-min(opinions);
        consensusQuestions=consensus(consensus(:,1)==group,2);
        consensusAnswers=consensus(consensus(:,1)==group,3);        
        
        totalLengthWordsGroup=nanmean(dataFull(dataFull(:,7)==group,q+26));
               
        %If the group had 3 members, and they reached consensus:
        if length(opinions)>2 && consensusAnswers(consensusQuestions==q+8)>=0 
            if initialrange>=minrange && initialrange<=maxrange %we keep only the selected ranges
            consensusWordsLengthGroup=[consensusWordsLengthGroup,totalLengthWordsGroup]; %save number of words
            end
        %If the group had 3 members, and they reached disensus:
        elseif length(opinions)>2 && consensusAnswers(consensusQuestions==q+8)==-1            
            if initialrange>=minrange && initialrange<=maxrange %we keep only the selected ranges
            disensusWordsLengthGroup=[disensusWordsLengthGroup,totalLengthWordsGroup]; %save number of words
            end
        end
        
    end %ends loop over groups
    end %ends loop over questions
    
    
% Now use the previously choosen initial value mean word length and the
% step, and count how many times a group had that length within that step, 
% and so on.
    
k=1;
counter=[];
for lengthWords=initial:step:endDens
    groupsConsensusInside=sum(consensusWordsLengthGroup>=lengthWords & consensusWordsLengthGroup<lengthWords+step);
    groupsDisensusInside=sum(disensusWordsLengthGroup>=lengthWords & disensusWordsLengthGroup<lengthWords+step);
    counter(k)=groupsConsensusInside./(groupsConsensusInside+groupsDisensusInside);
    k=k+1;
end


% We now plot the resulting curves; we separate both cases.
hold on
if cases==1
    plottingRange=initial:step:endDens-2*step; %sometimes it is useful to stop it earlier to erase averages from just one data point.
    %We also erase one case in the middle, resulting from the average of one data point.
    %plot(plottingRange([1:6,8:length(plottingRange)]),counter([1:6,8:length(plottingRange)]),'o','color',[0, 0.5, 0],'LineWidth',1.5)
    plot(plottingRange,counter(1:length(plottingRange)),'o','color',[0, 0.5, 0],'LineWidth',1.5) 
    %Logistic fit:
    mdl = fitglm(plottingRange,counter(1:length(plottingRange)),'Distribution','binomial');
    params = mdl.Coefficients.Estimate; %we save the fit parameters
    hold on
    %plot the fit:
    plot(initial-0.2:step:endDens+2*step,1./(exp(-(params(1)+params(2).*(initial-0.2:step:endDens+2*step)))+1),'color',[0, 0.5, 0],'LineWidth',1.5)
elseif cases==2
    plottingRange=initial:step:endDens-step;
    plot(plottingRange,counter(1:length(plottingRange)),'o','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)
    %Logistic fit:
    mdl = fitglm(plottingRange,counter(1:length(plottingRange)),'Distribution','binomial');
    params = mdl.Coefficients.Estimate;
    hold on
    %plot the fit
    plot(initial-0.2:step:endDens+2*step,1./(exp(-(params(1)+params(2).*(initial-0.2:step:endDens+2*step)))+1),'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)        
end

xlabel('Length per Word')
ylabel('P(Consensus)')
xlim([2.5,5])
ylim([0,1])
xticks(2.5:0.5:5)
yticks(0:0.2:1)
legend('Ranges<=5','Fit','Ranges>5','Fit','Location','SouthWest')
end    
