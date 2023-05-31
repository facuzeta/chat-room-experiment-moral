%% Moral Arguments Large Experiment

% This script reads the .csv file containing data from our participants
% and uses it to analyze the results.

% The required file is moral_Large_Experiment.csv.

%% Loading Data and Creating Useful Variables

clearvars;close all;clc

tableFullData=importdata('moral_Large_Experiment.csv',';');
data=tableFullData.textdata;
allStatements=1:5;

% We will loop over the required variables in the table, and progressively
% save them in MATLAB variables.

for i=3:size(data,1)
    j=1;
for statement=allStatements
    %Due to the coding of the data, different columns correspond to
    %different variables and different statements. So we need to do this
    %loop and go over each column to save everything appropriatly. 
    
    %We thoroughly comment each variable for the first statement; the
    %procedure is exactly the same for all other statements, and all
    %variables follow the same structure, so we do not comment them.
    
    if statement==1 
        a=data(i,3);  %Arguments
        b=data(i,1);  %preAceptability
        c=data(i,7);  %postAceptability
        
        %we create a coded variable for the conditions:
        if strcmp(data(i,6),'tendencioso_a_favor') %in favor of the original position.
            d=1;
        elseif strcmp(data(i,6),'balanceado') %one in favor and one against the original position (balanced).
            d=3;
        elseif strcmp(data(i,6),'tendencioso_en_contra') %against the original position.
            d=2;
        else
            d=4; %placebo
        end
        e=data(i,2); %preConfidence
        f=data(i,8); %postConfidence
        %we also create a coded variable for the sign of the presented arguments (with respect to the statement):
        if sum(str2num(a{1})==1 | str2num(a{1})==2)==2 %two positive arguments
            g=1;
        elseif sum(str2num(a{1})==3 | str2num(a{1})==4)==2 %two negative arguments
            g=2;
        elseif sum(str2num(a{1})==5 | str2num(a{1})==6)==2 %two placebo arguments
            g=4;
        else
            g=3; %one positive and one negative.
        end
        %which specific arguments were presented:
        h1=data(i,4);
        %the perceived validity of each argument:
        h2=data(i,5);
        
        % we create the final MATLAB variables for these two variables
        % after the loop over statements
                
    elseif statement==2
        a=data(i,11);
        b=data(i,9);
        c=data(i,15);
        if strcmp(data(i,14),'tendencioso_a_favor')
            d=1;
        elseif strcmp(data(i,14),'balanceado')
            d=3;
        elseif strcmp(data(i,14),'tendencioso_en_contra')
            d=2;
        else
            d=4;
        end       
        e=data(i,10);
        f=data(i,16);
        if sum(str2num(a{1})==7 | str2num(a{1})==8)==2
            g=1;
        elseif sum(str2num(a{1})==9 | str2num(a{1})==10)==2
            g=2;
        elseif sum(str2num(a{1})==11 | str2num(a{1})==12)==2
            g=4;
        else
            g=3;
        end
        h1=data(i,12);
        h2=data(i,13);        
    elseif statement==3
        a=data(i,19);
        b=data(i,17);
        c=data(i,23);
        if strcmp(data(i,22),'tendencioso_a_favor')
            d=1;
        elseif strcmp(data(i,22),'balanceado')
            d=3;
        elseif strcmp(data(i,22),'tendencioso_en_contra')
            d=2;
        else
            d=4;
        end       
        e=data(i,18);
        f=data(i,24);
        if sum(str2num(a{1})==13 | str2num(a{1})==14)==2
            g=1;
        elseif sum(str2num(a{1})==15 | str2num(a{1})==16)==2
            g=2;
        elseif sum(str2num(a{1})==17 | str2num(a{1})==18)==2
            g=4;
        else
            g=3;
        end
        h1=data(i,20);
        h2=data(i,21);
    elseif statement==4
        a=data(i,27);
        b=data(i,25);
        c=data(i,31);
        if strcmp(data(i,30),'tendencioso_a_favor')
            d=1;
        elseif strcmp(data(i,30),'balanceado')
            d=3;
        elseif strcmp(data(i,30),'tendencioso_en_contra')
            d=2;
        else
            d=4;
        end
        e=data(i,26);
        f=data(i,32);
        if sum(str2num(a{1})==19 | str2num(a{1})==20)==2
            g=1;
        elseif sum(str2num(a{1})==21 | str2num(a{1})==22)==2
            g=2;
        elseif sum(str2num(a{1})==23 | str2num(a{1})==24)==2
            g=4;
        else
            g=3;
        end
        h1=data(i,28);
        h2=data(i,29);        
    elseif statement==5
        a=data(i,35);
        b=data(i,33);
        c=data(i,39);
        if strcmp(data(i,38),'tendencioso_a_favor')
            d=1;
        elseif strcmp(data(i,38),'balanceado')
            d=3;
        elseif strcmp(data(i,38),'tendencioso_en_contra')
            d=2;
        else
            d=4;
        end
        e=data(i,34);
        f=data(i,40);  
        if sum(str2num(a{1})==25 | str2num(a{1})==26)==2
            g=1;
        elseif sum(str2num(a{1})==27 | str2num(a{1})==28)==2
            g=2;
        elseif sum(str2num(a{1})==29 | str2num(a{1})==30)==2
            g=4;
        else
            g=3;
        end
        h1=data(i,36);
        h2=data(i,37);        
        
    end %end of statements loop
    
    % We now create the final MATLAB variables for the received arguments 
    % and their perceived validity.
    temp1bis=str2num(h1{1});
    temp2bis=str2num(h2{1});
    
    receivedArguments(i-2,1:2,j)=str2num(a{1});
    perceivedValidity(i-2,1:2,j)=[temp1bis(2),temp2bis(2)];
    
    preAcceptability(i-2,j)=str2num(b{1});
    postAcceptability(i-2,j)=str2num(c{1});
    deltaAcceptability(i-2,j)=postAcceptability(i-2,j)-preAcceptability(i-2,j);
    
    % The "experimental" condition depends on the initial stance of the
    % participant and whether the presented arguments go in favor, against,
    % mixed or placebo with respect to that stance.
    experimentalCondition(i-2,j)=d; 
    
    % The "real" condition depends on the moral statement and whether the 
    % presented arguments are positive, negative, mixed or placebo with 
    % respect to that statement.
    realCondition(i-2,j)=g;
    
    preConfidence(i-2,j)=str2num(e{1});
    postConfidence(i-2,j)=str2num(f{1}); 
    statements(i-2,j)=statement;
    j=j+1;
    
end %end of loop over statements

end %end of loop over columns
'end'

newFullTable=[preAcceptability(:),postAcceptability(:),realCondition(:),experimentalCondition(:),preConfidence(:),postConfidence(:),deltaAcceptability(:),statements(:)];
% Col 1: Initial Acceptability
% Col 2: Final Acceptability
% Col 3: Real Condition (++, --, +-, placebo -> 1,2,3,4)
% Col 4: Experimental Condition (2 in favor, 2 against, balanced, placebo -> 1,2,3,4)
% Col 5: Initial Confidence
% Col 6: Final Confidence
% Col 7: delta Acceptability
% Col 8: Statement (1,2,3,4,5)

% We will also need the average word length of all the arguments (obtained
% outside of MATLAB):

lengthsArguments=[4.875,5.2857,6.0476,5.7059,5.6,5.9,6.375,6.6429,6.5882,5.9167,5.0909,5.15,...
    5.7143,5.875,5.4074,4.8333,5.5714,5.1613,5.2121,4.3,5.7241,5.24,5.4286,4.4706,6.1176,5,...
    5.3913,5.0833,4.7059,4.9412];

% We also consider the number of sentences, which was 1 in most arguments (
% coded as 0), and 2 sentences in a few arguments (coded as 1).
sentencesArguments=[0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0];

% Clear irrelevant variables:
clear a ans b c d e f g data statement allStatements h1 h2 i j tableFullData temp1bis temp2bis

save('VariablesMoralExp2') %we save everything, to avoid doing this step over.

%% Plot 1: Probability of Increasing Acceptability (uses "real" conditions)

% Optional:
clearvars
load('VariablesMoralExp2')

% We create the appropriate variables (probability of consensus for each condition):
increasedPlus=sum(postAcceptability(realCondition==1)>preAcceptability(realCondition==1))/length(postAcceptability(realCondition==1));
increasedMinus=sum(postAcceptability(realCondition==2)>preAcceptability(realCondition==2))/length(postAcceptability(realCondition==2));
increasedPlusMinus=sum(postAcceptability(realCondition==3)>preAcceptability(realCondition==3))/length(postAcceptability(realCondition==3));
increasedNeutral=sum(postAcceptability(realCondition==4)>preAcceptability(realCondition==4))/length(postAcceptability(realCondition==4));

% We make the plot.
figure('color','w','paperposition',[0 0 11/2 5/2]);
bar(1:4,[increasedMinus,increasedPlusMinus,increasedNeutral,increasedPlus])
xlim([0,5])
ylim([0,0.55])
xlabel('')
ylabel('P(increase acceptability)')
xticklabels({'--','+-','Placebo','++'})
hold on

% Optional:
% To perform proportionality tests (chi-square), we use: 
% [h,p]=prop_test([subieronNeutral*10549,subieronPlus*10549],[10549,10549],'false')
sigstar({[1,2], [2,3], [3,4]},[0,0.583,0]) %to plot the significance bars.

box off
set(gca,'TickDir','out')
set(gca,'FontSize',9)

%saveas(gcf,'PIncreaseAcceptability','png')
print -dpdf PIncreaseAcceptability.pdf

%% Plot 2: Probability of Depolarization (uses "experimental" conditions)

% Optional:
clearvars
load('VariablesMoralExp2')

% We create the appropriate variables (probability of consensus for each condition):
% In order to assess polarization, we need to refer the acceptabilities to
% the middle point of the scale, which was 50. 
preCenteredAcceptability=abs(preAcceptability-50);
postCenteredAcceptability=abs(postAcceptability-50);

% We now check in which cases these new centered acceptabilities are bigger
% in the initial question, the last question, or remain fixed.
PolPlus=postCenteredAcceptability>preCenteredAcceptability;
PolMinus=postCenteredAcceptability<preCenteredAcceptability;
PolFixed=postCenteredAcceptability==preCenteredAcceptability;

% We now count the previous cases to get the probabilities.
loweredPlus=sum(PolMinus(experimentalCondition==1))/length(PolMinus(experimentalCondition==1));
loweredMinus=sum(PolMinus(experimentalCondition==2))/length(PolMinus(experimentalCondition==2));
loweredPlusMinus=sum(PolMinus(experimentalCondition==3))/length(PolMinus(experimentalCondition==3));
loweredNeutral=sum(PolMinus(experimentalCondition==4))/length(PolMinus(experimentalCondition==4));

% We make the plot.
figure('color','w','paperposition',[0 0 11/2 5/2]);
bar(1:4,[loweredMinus,loweredPlusMinus,loweredNeutral,loweredPlus])
xlim([0,5])
ylim([0,0.55])
xlabel('')
ylabel('P(depolarize opinion)')
xticklabels({'2 against','1 favor 1 against','Placebo','2 favor'})
hold on

% Optional:
% To perform proportionality tests (chi-square), we use: 
% [h,p,chi2]=prop_test([loweredPlusMinus*10549,loweredNeutral*10549],[10549,10549],'false')
% [h,p,chi2]=prop_test([loweredMinus*10549,loweredPlusMinus*10549],[10549,10549],'false')
% [h,p,chi2]=prop_test([loweredPlus*10549,loweredNeutral*10549],[10549,10549],'false')
sigstar({[1,2], [2,3], [3,4]},[0,0,0])

box off
set(gca,'TickDir','out')
set(gca,'FontSize',9)

%saveas(gcf,'PDepolarizeOpinion','png')
print -dpdf PDepolarizeOpinion.pdf

%% Plot 3: Mixed Regression Coefficients with Word Length and Sentences

% Optional:
clearvars
load('VariablesMoralExp2')

% Initializations
allValidities = [];
allLengths = [];
allStatements = [];
allPosArg = [];
allNegArg = [];
allPreAccept = [];
allSentences = [];

% Positive and Negative Arguments (from the spreadsheet coding).
posArg = [1,2,7,8,13,14,19,20,25,26];
negArg = [1,2,7,8,13,14,19,20,25,26]+2;

% We create the appropriate variables for the regression, including
% one for all validity ratings for every argument and statement, one
% for the corresponding length of the presented arguments, one for the
% corresponding number of sentences, one for the corresponding statements,
% and one for the corresponding acceptability of the statements.

% Loop Over Subjects
for s=1:size(perceivedValidity,1)
    
    aux=squeeze(perceivedValidity(s,:,:));
    allValidities = [allValidities;aux(:)];
    
    aux2=squeeze(receivedArguments(s,:,:));
    indexes=aux2(:);
    
    allLengths = [allLengths;lengthsArguments(indexes)'];
    allSentences = [allSentences;sentencesArguments(indexes)'];

    allStatements = [allStatements;[1,1,2,2,3,3,4,4,5,5]'];
    allPosArg = [allPosArg;ismember(indexes,posArg)];
    allNegArg = [allNegArg;ismember(indexes,negArg)];
    
    aux3 = preAcceptability(s,:);
    aux4=[aux3;aux3];
    allPreAccept = [allPreAccept;aux4(:)];
    
end
%

% Now we create the consistency variable.
allCons = sign(allPreAccept-50).*allPosArg-sign(allPreAccept-50).*allNegArg;

% Normalize lengths for use in regression:
allLengthsNorm=zscore(allLengths);

% create regression:
tbl = table(allValidities,allCons,allSentences,allLengthsNorm,allStatements);
mdl = 'allValidities ~ allSentences + allLengthsNorm + allCons + (1|allStatements)';
na=fitlme(tbl,mdl);

% Save relevant parameters from the regression:
Params=na.Coefficients.Estimate;
PVals=na.Coefficients.pValue;
consParam=Params(2);
sentsParam=Params(3);
lengthParam=Params(4);

consPval=PVals(2);
sentsPval=PVals(3);
lengthPval=PVals(4);

% Plot the bar plot with the appropriate significance bars
figure();bar([1,2,3],[lengthParam,sentsParam,consParam])
sigstar({[1,1], [2,2], [3,3]},[lengthPval,sentsPval,consPval])

xticklabels({'Word Length','# Sentences','Consistency'})
yticks(-0.1:0.1:0.3)

ylabel('Coefficient Estimates')
box off
set(gca,'TickDir','out')
set(gca,'FontSize',9)

print -dpdf CoeffEstsMoralExp2.pdf

%% Plot 4: Length Scatter Plot 

% Here we make a scatter plot to observe the effect of length on validity.

clear Sbinned
L=allLengths;
S=allValidities;

step=0.01; %binning step
range=4.3:step:6.8-step; %appropriate interval of length values

% Binning
for i=range
    Sbinned(floor(1+(i-4.3)/step))=mean(S(i<=L & L<i+step));
    SbinnedError(floor(1+(i-4.3)/step))=nanstd(S(i<=L & L<i+step))./sqrt(length(S(i<=L & L<i+step)));    
end

% Plot
figure();
errorbar(range(Sbinned>2.3),Sbinned(Sbinned>2.3),SbinnedError(Sbinned>2.3),'b.','MarkerFace','b','MarkerSize',25);

% Linear Fit
mdl=fitlm(range(Sbinned>2.3),Sbinned(Sbinned>2.3),'RobustOpts','on');
a=mdl.Coefficients.Estimate(1);
b=mdl.Coefficients.Estimate(2);
hold on
plot(1:0.1:6.8,a+b*(1:0.1:6.8),'b','LineWidth',3)

xlabel('Mean Word Length')
ylabel('Mean Score')
yticks(2.4:0.4:3.6)
ylim([2.399,3.6])
xticks(4.6:0.4:6.6)
xlim([4.599,6.6])
box off

set(gca,'FontSize',9)
set(gca,'TickDir','out')

print -dpdf LengthExp2.pdf
