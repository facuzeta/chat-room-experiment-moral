%% Group Moral Discussions Experiment - Replications

% This script reads the .csv files containing data from our participants
% and uses it to replicate the main results from the paper:

%Navajas, J., Heduan, F. Á., Garrido, J. M., Gonzalez, P. A., Garbulsky, 
%G., Ariely, D., & Sigman, M. (2019). Reaching consensus in polarized 
%moral debates. Current Biology, 29(23), 4124-4129.

% The required files are dataS1.csv, dataS3.csv, and 
% moralConsensus.csv. The MATLAB function prop_test is also required.

%% Loading the Required Files

clearvars;close all;clc

a=importdata('dataS1.csv');
b=importdata('dataS3.csv');

dataS1=a.data;
dataS3=b.data;

c=importdata('moralConsensus.csv');

consensus=c.data;
clear a b c

%% Creating useful variables for Plot 1.

% First, we keep only the valid groups, found in moral_consensus.csv
groups=unique(consensus(:,1));

%Now we calculate the initial range of opinions ("deltas") for each group
j=1;
deltas=nan(length(groups),6);
for group=groups'
    questions=consensus(consensus(:,1)==group,2);
    for q=questions
        answersSubjects=dataS1(dataS1(:,7)==group,q-8);
        deltas(j,q-8)=max(answersSubjects)-min(answersSubjects); 
    end
    j=j+1;
end

%We now need the probability of reaching a consensus, for each possible
%value of delta. We estimate it for each question.

probaconsensus=nan(11,6);
for q=1:6
    quest_q_consensus=consensus(consensus(:,2)==(q+8),3);
    quest_q_deltas=deltas(~isnan(deltas(:,q)),q);
    for delta=0:10
        consensosDelta=quest_q_consensus(quest_q_deltas==delta);
        probaconsensus(delta+1,q)=sum(consensosDelta>=0)/length(consensosDelta);
    end
end

%% Plot 1: Probability of Consensus vs. Range

% This is the first plot, and it shows the probability of reaching group
% consensus as a function of the initial range of opinions (delta).

figure()

plot(0:10,mean(probaconsensus,2),'.','MarkerSize',50)

%We include a logistic fit of the data:
mdl = fitglm((0:10)./10,mean(probaconsensus,2),'Distribution','binomial');

B=mdl.Coefficients.Estimate;
x0=B(1)+B(2).*(0:10)./10;
fit=1./(exp(-x0)+1);
hold on

plot(0:10,fit,'color',[0 0.4470 0.7410],'LineWidth',2)

title('Mean of All Questions')
ylim([-0.1,1.1])
xlim([-1,11])
ylabel('p(consensus)')
xlabel('Range of Ratings (\delta=rmax-rmin)')
box off

print -dpdf ProbabilityOfConsensusVsRange.pdf

%% Creating useful variables for Plot 2.

%Para borrar cuando tenga los definitivos:
% load('DatosMoralS1ConfUpd')
% load('DatosMoralS3ConfUpd')
% dataS1=datosArreglados;
% dataS3=datosArreglados2;
% clear datosArreglados datosArreglados2
%Fin del borrar

% We will need the mean confidence value for each possible rating.
for q=1:6
    for rating=0:10
        confmean(q,rating+1)=nanmean(dataS1(dataS1(:,q)==rating,q+7));        
    end
end

%We also need the distribution of confidences for the lowest, middle and
%highest rating values.
confdists={};
for q=1:6
    for rating=[0,5,10]
        confdists{q,rating/5+1}=dataS1(dataS1(:,q)==rating,q+7);
    end
end


%% Plot 2: Mean Confidence vs. Rating & Selected Distributions

% This is the second plot, which shows how the mean confidence changes with
% respect to the participant's rating. It also shows selected distributions
% for the minimum rating (0), the maximum rating (10), and the middle
% rating (5).

q=1; %question chosen for the selected distributions.

%Base plot
figure()
subplot(2,3,[1,2,3])
plot(0:10,mean(confmean,1),'.','MarkerSize',30)
xlim([-0.5,10.5])
xlabel('Rating');
ylabel('Mean Confidence')
title('How Confident Do You Feel About Your Opinion')
coefs = polyfit(0:10,mean(confmean,1),2);
hold on
plot(-0.5:0.1:10.5,(-0.5:0.1:10.5).^2.*coefs(1)+(-0.5:0.1:10.5).*coefs(2)+coefs(3),'color',[0 0.4470 0.7410],'LineWidth',1.5)
box off
ylim([3.1,5])

% Distribution plots
for i=1:3
confidence=squeeze(confdists{q,i});
subplot(2,3,i+3)
counts = histcounts(confidence,5);
bar(counts)
title(['q=',num2str(q),', rating=',num2str((i-1)*5)])
xlabel('Confidence')
end

print -dpdf ConfidenceVsRating.pdf

%% Creating useful variables for Plot 3.

%We separate the answers question by question:
questions={[],[],[],[],[],[]};

for i=1:length(consensus)
    if consensus(i,2)==9 && consensus(i,3)>=0
        questions{1}=[questions{1},consensus(i,3)];
    elseif consensus(i,2)==10 && consensus(i,3)>=0
        questions{2}=[questions{2},consensus(i,3)];
    elseif consensus(i,2)==11 && consensus(i,3)>=0
        questions{3}=[questions{3},consensus(i,3)];
    elseif consensus(i,2)==12 && consensus(i,3)>=0
        questions{4}=[questions{4},consensus(i,3)];        
    elseif consensus(i,2)==13 && consensus(i,3)>=0
        questions{5}=[questions{5},consensus(i,3)];        
    elseif consensus(i,2)==14 && consensus(i,3)>=0
        questions{6}=[questions{6},consensus(i,3)];
    end
end

% We will now separate the members with the maximum opinion in the group
% from those with the minimum opinion, and from those with an opinion in
% the middle. We will do the same for their confidence ratings. We also
% separate between groups that reached consensus and groups that reached
% disensus.

maxes={[],[],[],[],[],[]};
minimums={[],[],[],[],[],[]};
centers={[],[],[],[],[],[]};
confidenceMaxes={[],[],[],[],[],[]};
confidenceMinimums={[],[],[],[],[],[]};
confidenceCenters={[],[],[],[],[],[]};
maximumsDisensus={[],[],[],[],[],[]};
minimumsDisensus={[],[],[],[],[],[]};
centersDisensus={[],[],[],[],[],[]};
confidenceMaxesDisensus={[],[],[],[],[],[]};
confidenceMinimumsDisensus={[],[],[],[],[],[]};
confidenceCentersDisensus={[],[],[],[],[],[]};

questions2={[],[],[],[],[],[]};
questions2Disensus={[],[],[],[],[],[]};

for q=1:6
    question=questions{q};
    groups=consensus(consensus(:,2)==q+8,1);
    for group=groups'
        [members,orden]=sort(dataS1(dataS1(:,7)==group,q));
        confidences=dataS1(dataS1(:,7)==group,q+7);
        confidences=confidences(orden);
        members=members(~isnan(members));
        consensusQuestions=consensus(consensus(:,1)==group,2);
        consensusAnswers=consensus(consensus(:,1)==group,3);
        if length(members)>2 && consensusAnswers(consensusQuestions==q+8)>=0
            maxes{q}=[maxes{q},members(3)];
            minimums{q}=[minimums{q},members(1)];
            centers{q}=[centers{q},members(2)];            
            questions2{q}=[questions2{q},consensusAnswers(consensusQuestions==q+8)];
            confidenceMaxes{q}=[confidenceMaxes{q},confidences(3)];
            confidenceMinimums{q}=[confidenceMinimums{q},confidences(1)];
            confidenceCenters{q}=[confidenceCenters{q},confidences(2)];
        elseif length(members)>2 && consensusAnswers(consensusQuestions==q+8)==-1
            maximumsDisensus{q}=[maximumsDisensus{q},members(3)];
            minimumsDisensus{q}=[minimumsDisensus{q},members(1)];
            centersDisensus{q}=[centersDisensus{q},members(2)];            
            questions2Disensus{q}=[questions2Disensus{q},consensusAnswers(consensusQuestions==q+8)];
            confidenceMaxesDisensus{q}=[confidenceMaxesDisensus{q},confidences(3)];
            confidenceMinimumsDisensus{q}=[confidenceMinimumsDisensus{q},confidences(1)];
            confidenceCentersDisensus{q}=[confidenceCentersDisensus{q},confidences(2)];                        
        end     
    end
end

%% Plot 3: Weight of Each Member on the Group Rating.

% Here we evaluate the weight each group member has (either the one with
% the lowest opinion value, the one with the highest opinion value, or the
% one with an intermediate value) on the group answer (consensus).

% In order to do this, we use a linear regression with each member's
% individual rating as regressors, and the group consensus as the 
% dependent variable.
mdl = fitlm([[minimums{:}]',[centers{:}]',[maxes{:}]'],[questions2{:}]);

% We compare the previous weights to the ones that we would obtain if the
% procedure of the group to reach a consensus was merely averaging their
% individual answers from the first stage.
mdl2 = fitlm([[minimums{:}]',[centers{:}]',[maxes{:}]'],([minimums{:}]'+[centers{:}]'+[maxes{:}]')/3);

estimates=mdl.Coefficients.Estimate(2:end);
SEs=mdl.Coefficients.SE(2:end);

SAests=mdl2.Coefficients.Estimate(2);
SASEs=mdl2.Coefficients.SE(2);

% We now plot the estimates for each group member.
figure()
plot(1:3,estimates,'k-','MarkerSize',5)
hold on
plot(1:3,estimates,'k.','MarkerSize',30)
errorbar(1:3,estimates,SEs,'.k');
plot(1:3,[SAests,SAests,SAests],'rs','MarkerSize',10,'MarkerFace','r')
plot(1:3,[SAests,SAests,SAests],'-ro','MarkerSize',3)
errorbar(1:3,[SAests,SAests,SAests],[SASEs,SASEs,SASEs],'.k');
xlim([0.5,3.5])
ylim([-0.05,0.8])
yticks(0:0.2:0.8)
ylabel('Weight in Group Rating');
xticks([1,2,3])
xticklabels({'P_{min}','P_{med}','P_{max}'})
box off


%% Creating New Variables for Probability of Depolarization for Discussed vs. Undiscussed Questions

k=1;
%Discussed Questions:
for q=1:6
    groups=consensus(consensus(:,2)==q+8,1);
    for group=groups'
        opinionsS1=dataS1(dataS1(:,7)==group,q); %opinions of stage 1
        range=max(opinionsS1)-min(opinionsS1); %initial range of opinions
        memberIds=dataS1(dataS1(:,7)==group,14); 
        opinionsS3=dataS3(dataS3(:,7)==group,q); %opinions of stage 3
        consensusQuestions=consensus(consensus(:,1)==group,2); %the questions (valid only)
        consensusAnswers=consensus(consensus(:,1)==group,3); %the consensus reached (valid only)
        
        %Consensus Reached
        if length(opinionsS1)>2 && length(opinionsS3)>2 && consensusAnswers(consensusQuestions==q+8)>=0
            for integr=memberIds'
                difOp=dataS3(dataS3(:,7)==group & dataS3(:,14)==integr,q)-dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q);
                if dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q)>5
                    depol=difOp<0;
                elseif dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q)<5
                    depol=difOp>0;
                else
                    depol=0;
                end
                DifferencesOpinion(k,1:7)=[difOp,q,integr,group,1,range,depol];
                k=k+1;
            end
        end
        
        %Disensus
        if length(opinionsS1)>2 && length(opinionsS3)>2 && consensusAnswers(consensusQuestions==q+8)==-1
            for integr=memberIds'
                difOp=dataS3(dataS3(:,7)==group & dataS3(:,14)==integr,q)-dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q);                
                if dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q)>5
                    depol=difOp<0;
                elseif dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q)<5
                    depol=difOp>0;
                else
                    depol=0;
                end                
                DifferencesOpinion(k,1:7)=[difOp,q,integr,group,0,range,depol];
                k=k+1;
            end
        end
        
    end
end

%Not discussed:
groups=unique(consensus(:,1));
for group=groups'
        discutidas=consensus(consensus(:,1)==group,2)-8;
        for q=1:6
            if ~ismember(q,discutidas)
                opinionsS1=dataS1(dataS1(:,7)==group,q);
                range=max(opinionsS1)-min(opinionsS1);                
                memberIds=dataS1(dataS1(:,7)==group,14);
                opinionsS3=dataS3(dataS3(:,7)==group,q);         
                if length(opinionsS1)>2 && length(opinionsS3)>2
                    for integr=memberIds'
                      difOp=dataS3(dataS3(:,7)==group & dataS3(:,14)==integr,q)-dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q);
                      if dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q)>5
                         depol=difOp<0;
                      elseif dataS1(dataS1(:,7)==group & dataS1(:,14)==integr,q)<5
                         depol=difOp>0;
                      else
                          depol=0;
                      end                         
                      DifferencesOpinion(k,1:7)=[difOp,q,integr,group,-1,range,depol];
                      k=k+1;
                    end
                end
            end
        end
end


%% Plot 4: Probability of Depolarization for Discussed vs. Undiscussed Questions

figure()
hold on
bar(1,sum(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)<=5,7))./length(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)<=5,7)) ,'b')
bar(2,sum(DifferencesOpinion(DifferencesOpinion(:,5)==-1 & DifferencesOpinion(:,6)<=5,7))./length(DifferencesOpinion(DifferencesOpinion(:,5)==-1 & DifferencesOpinion(:,6)<=5,7)) ,'r')

bar(3.5,sum(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)>5,7))./length(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)>5,7)) ,'b')
bar(4.5,sum(DifferencesOpinion(DifferencesOpinion(:,5)==-1 & DifferencesOpinion(:,6)>5,7))./length(DifferencesOpinion(DifferencesOpinion(:,5)==-1 & DifferencesOpinion(:,6)>5,7)) ,'r')

discussedSmallRange=sum(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)<=5,7));
notDiscussedSmallRange=sum(DifferencesOpinion(DifferencesOpinion(:,5)==-1 & DifferencesOpinion(:,6)<=5,7));

discussedBigRange=sum(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)>5,7));
notDiscussedBigRange=sum(DifferencesOpinion(DifferencesOpinion(:,5)==-1 & DifferencesOpinion(:,6)>5,7));

% For using the proportion tests:
lengthDiscussedSmallRange=length(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)<=5,7));
lengthNotDiscussedSmallRange=length(DifferencesOpinion(DifferencesOpinion(:,5)==1 & DifferencesOpinion(:,6)<=5,7));

lengthDiscussedBigRange=length(DifferencesOpinion(DifferencesOpinion(:,5)~=1 & DifferencesOpinion(:,6)>5,7));
lengthNotDiscussedBigRange=length(DifferencesOpinion(DifferencesOpinion(:,5)==1 & DifferencesOpinion(:,6)>5,7));

xticks([1.5,4])
legend('Discussed','Not Discussed','Location','NorthWest')
xticklabels({'R<5','R>5'})
yticks(0:0.1:0.3)
ylabel('P(depolarization)')
ylim([0,0.3])


[h,p]=prop_test([ discussedSmallRange, notDiscussedSmallRange ] , [lengthDiscussedSmallRange, lengthNotDiscussedSmallRange], false)
%-9
[h,p]=prop_test([ discussedSmallRange, notDiscussedBigRange ] , [lengthDiscussedSmallRange, lengthNotDiscussedBigRange], false)
%0.6
[h,p]=prop_test([ notDiscussedSmallRange, notDiscussedBigRange ] , [lengthNotDiscussedSmallRange, lengthNotDiscussedBigRange], false)
%-8

[h,p]=prop_test([ discussedBigRange, notDiscussedSmallRange ] , [lengthDiscussedBigRange, lengthNotDiscussedSmallRange], false)
%0
[h,p]=prop_test([ discussedBigRange, notDiscussedBigRange ] , [lengthDiscussedBigRange, lengthNotDiscussedBigRange], false)
%-4
[h,p]=prop_test([ discussedBigRange, discussedSmallRange ] , [lengthDiscussedBigRange, lengthDiscussedSmallRange], false)
%-6