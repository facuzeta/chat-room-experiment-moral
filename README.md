# chat-room-experiment-moral

This repository contains files related to the "chat rooms experiments involving moral dilemmas" chapter. 

The raw_data folder contains the data_experiment_MORAL.json file, which has all the data collected for our main experiment.

The feature_extraction folder employs the data_experiment_MORAL.json file to create the following files:

- dataS1.csv
- dataS3.csv
- features_moral_times.csv

Additionally, the data was tagged by experimenters to assess whether groups actually reached consensus or not 
(i.e. for checking cases where consensus was sent by some of the members of the group but not others; or where
the consensus that was sent differed between members; or when no consensus was sent, but it was clearly agreed
upon by them), the results of which can be found in the moralConsensus.csv file.

All four previously mentioned csv files can be found in the analysis folder, which also contains the MATLAB files 
required for analysis and for making the presented Figures.

The Group_Moral_Discussions_Replication.m file is used for creating the plots that replicate those previously found
in the paper:

Navajas, J., Heduan, F. Á., Garrido, J. M., Gonzalez, P. A., Garbulsky, G., Ariely, D., & Sigman, M. (2019). 
Reaching consensus in polarized moral debates. Current Biology, 29(23), 4124-4129.

This file also requires the MATLAB function prop_test, which is also included in the repository. Four plots can be
made with this file, and the following data files are required: dataS1.csv, dataS3.csv, moralConsensus.csv.


The Group_Moral_Discussions_Consensus.m file is used for creating the plots that show how the Probability of Consensus
changes with respect to either the average number of words per message in the conversation, or the average word length
of the words in the conversation. Two plots can be made with this file, and the following data files are required:
dataS1.csv, and moralConsensus.csv.


The Group_Moral_Discussions_Decoding.m file is used for creating the decoding model plot. The required data files are:
dataS1.csv, moralConsensus.csv, and  features_moral_times.csv. It also requires a custom automatically generated MATLAB 
function, importfileMoralTimes.m, and a plotting MATLAB function called hline.m; both are included in the repository.


The Moral_Arguments_Large_Experiment.m file is used for analyzing and creating all the plots from our second, large-scale 
experiment. The data for this experiment can be found in the moral_Large_Experiment.csv file, which is required by the former.
It is used to make four plots, and it also requires the additional MATLAB functions prop_test.m and sigstar.m (also included).


The Moral_Linguistic_Manipulations_Experiment.m file is used for analyzing and creating all the plots from our third, linguistic
manipulations experiment. The data for this experiment can be found in the MoralExp3Data.xlsx file, which is required by the former.
It is used to make three plots. A custom automatically generated  MATLAB function is also required to properly import the data, 
called importExp3Moral.m., also included in this repository. Additionally, in order to perform the mediation analysis, installing 
these two toolboxes is required:

https://github.com/canlab/MediationToolbox
https://github.com/canlab/CanlabCore
