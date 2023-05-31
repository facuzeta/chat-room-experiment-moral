%% Moral Altered Arguments Experiment

% This script reads the .xlsx file containing data from our participants
% and uses it to analyze the results.

% The required file is MoralExp3Data.xlsx. A custom automatically generated
% MATLAB function is also required to properly import the data.

%% Loading Data and Creating Useful Variables

clearvars; close all; clc;

% Use the custom function to import the data:
MoralExp = importExp3Moral('MoralExp3Data.xlsx','Hoja1',1,689);

% Eliminate excluded subjects (failed attentional test)

MoralExp(MoralExp.Cereza~='Cereza',:) = [];

deleted=[];
for j=1:600
if sum(~isnan(table2array(MoralExp(j,1:294))))~=78
    deleted=[deleted,j];
end
end

MoralExp(deleted,:) = [];

% We now create useful variables for further analysis:

Acceptability=[];
for i=0:5
    Acceptability=[Acceptability;[MoralExp{:,1+49*i},ones(length(MoralExp{:,1+49*i}),1)*(i+1)],(1:size(MoralExp,1))'];
end

tableFullExp3=[];
for i=0:5
    for j=1:size(MoralExp,1)
        acceptabilities=MoralExp{j,1+49*i};
        subjects=MoralExp{j,2+49*i:3:49+49*i};
        validities=subjects(~isnan(subjects))';
        aux=MoralExp{j,3+49*i:3:49+49*i};
        difficulties=aux(~isnan(aux));
        temp=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
        conditions=temp(~isnan(subjects))'; 
        % 1 = short and fewer words
        % 2 = long and fewer words
        % 3 = short and more words
        % 4 = long and more words
        consistencies=sign((acceptabilities-5)).*[1,1,-1,-1]';
        % consistency = 1 when the argument was in favor of the original
        % stance, and consistency = -1 when the argument is against it.
        signsArgs=[1,1,-1,-1]';
        if ~isempty(conditions)
            tableFullExp3=[tableFullExp3;[ones(4,1)*acceptabilities,validities,consistencies,conditions,ones(4,1)*(i+1),ones(4,1)*j,difficulties',signsArgs]];
        end
        
    end
end

Acceptability=tableFullExp3(:,1);
Validity=tableFullExp3(:,2)/10;
Consistency=tableFullExp3(:,3);
Condition=tableFullExp3(:,4);
Statement=tableFullExp3(:,5);
Subjet=tableFullExp3(:,6);
Sign=tableFullExp3(:,8);
Disfluency=tableFullExp3(:,7)/10;

ConditionLength=Condition==2 | Condition==4; % when words are longer
ConditionDensity=Condition==3 | Condition==4; % when there are more words

% This is the main table we will be using for the analysis
tblExp3 = table(Validity,Consistency,ConditionLength,ConditionDensity,Statement,Acceptability,Subjet,Condition,Sign,Disfluency);

% Clear unnecesary variables
clear acceptabilities conditions consistencies difficulties deleted aux i j signsArgs subjects temp validities MoralExp

% Save everything
save('VariablesMoralExp3')


%% Plot 1: Perceived disfluency vs. conditions

figure('color','w','paperposition',[0 0 7/2 5/2]);

% Conditions means and standard errors:
means=[mean(Disfluency(Condition==1)*10),mean(Disfluency(Condition==2)*10),mean(Disfluency(Condition==3)*10),mean(Disfluency(Condition==4)*10)];
errors=[std(Disfluency(Condition==1)*10)./sqrt(length(Disfluency(Condition==1)*10)),std(Disfluency(Condition==2)*10)./sqrt(length(Disfluency(Condition==2)*10)),std(Disfluency(Condition==3)*10)./sqrt(length(Disfluency(Condition==3)*10)),std(Disfluency(Condition==4)*10)./sqrt(length(Disfluency(Condition==4)*10))];

bar([1,2,3,4],means)
hold on
errorbar([1,2,3,4],means,errors,'ko')
xlim([0.5,4.5])
xticks(1:4)
xticklabels({'Less and Short','Longer','More','Longer and More'})
yticks(0:4)
ylim([0,4])
ylabel('Perceived Disfluency')
box off
set(gca,'TickDir','out')
set(gca,'FontSize',9)

print -dpdf InsetAllConditionsMoralExp4.pdf

%%% Optional non-corrected statistical tests:
% [h,p]=ttest2(Disfluency(Condition==1),Disfluency(Condition==4))
% [h,p]=ttest2(Disfluency(Condition==2),Disfluency(Condition==4))
% [h,p]=ttest2(Disfluency(Condition==3),Disfluency(Condition==4))
% 
% [h,p]=ttest2(Disfluency(Condition==1),Disfluency(Condition==2))
% [h,p]=ttest2(Disfluency(Condition==1),Disfluency(Condition==3))
% [h,p]=ttest2(Disfluency(Condition==2),Disfluency(Condition==3))

[p,~,stats] = anova1(Disfluency,Condition,'off');

[results,meanVals,~,gnames]=multcompare(stats,"CType","bonferroni","Display","off");
% [results,meanVals,~,gnames]=multcompare(stats,"CType","scheffe","Display","off"); 
% all corrections yield similar results.

%% Plot 2: Perceived Disfluency vs. Validity

% Optional
clearvars
load('VariablesMoralExp3')

% We plot perceived disfluency vs. the mean validity
figure('color','w','paperposition',[0 0 7/2 5/2]);

for i=0:10
validityValues(i+1)=mean(Validity(Disfluency*10==i)*10);
errorValues(i+1)=std(Validity(Disfluency*10==i)*10)./sqrt(length(Validity(Disfluency*10==i)));
end

errorbar(1:11,validityValues,errorValues,'o')

% Linear Fit:
mdl=fitlm(1:11,validityValues,'RobustOpts','on')
c=mdl.Coefficients.Estimate(1);
d=mdl.Coefficients.Estimate(2);
hold on
plot(0:12,c+(0:12)*d,'b','LineWidth',3)

ylim([4.7,6.7]);
yticks(5:0.4:6.7);
xlim([0,12]);
xticks(1:2:11);
xticklabels({0:2:10});
xlabel('Perceived Disfluency');
ylabel('Perceived Validity');

box off
set(gca,'TickDir','out')
set(gca,'FontSize',9)

print -dpdf MoralExp4Disfluency.pdf


%% Plot 3: Mediation Analysis

% The necessary values for plot 3 are obtained here.

% In order to perform this analysis, installing these two toolboxes are
% required:

% https://github.com/canlab/MediationToolbox
% https://github.com/canlab/CanlabCore

% We add both to the MATLAB path. This should be modified as needed,
% depending on your PC.

addpath(genpath('C:\Program Files\MATLAB\R2018a\bin\MediationToolbox-master'));
addpath(genpath('C:\Program Files\MATLAB\R2018a\bin\CanlabCore-master'));

% We now do mediation analysis for the conditions with longer words:

Y=tblExp3.Validity;
M=tblExp3.Disfluency;
X=double(tblExp3.ConditionLength);

[paths, stats] = mediation(X, Y, M, 'boot', 'plots', 'verbose');
% (it can take a while to finish)

% We now do mediation analysis for the conditions with more words:

Y=tblExp3.Validity;
M=tblExp3.Disfluency;
X=double(tblExp3.ConditionDensity);

[paths, stats] = mediation(X, Y, M, 'boot', 'plots', 'verbose');
% (it can take a while to finish)
