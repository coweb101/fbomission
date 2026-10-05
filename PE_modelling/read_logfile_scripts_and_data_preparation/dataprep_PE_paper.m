function [datastruc_stimuli,datastruc_choice,datastruc_not_choice,datastruc_feedback,datastruc_feedback_type,filenames,stim_rew_prob]=dataprep_PE_paper

% Datenvorbereitung

%clear all;
addpath(genpath('./interim_datasets/'));
load 'all_vp_behav_paper.mat'; % read file created with read_data

allmatrices = {VP01_all_data,VP02_all_data,VP03_all_data,VP04_all_data,VP05_all_data,VP06_all_data,VP07_all_data,VP08_all_data,VP09_all_data,VP10_all_data,VP11_all_data,VP12_all_data,VP13_all_data,VP14_all_data,VP15_all_data,VP16_all_data,VP17_all_data,VP18_all_data,VP19_all_data,VP20_all_data,VP21_all_data,VP22_all_data,VP23_all_data,VP24_all_data,VP25_all_data,VP26_all_data,VP27_all_data,VP28_all_data,VP29_all_data,VP30_all_data,VP31_all_data,VP32_all_data,VP33_all_data,VP34_all_data,VP35_all_data,VP36_all_data,VP37_all_data,VP38_all_data,VP39_all_data,VP40_all_data,VP41_all_data,VP42_all_data,VP43_all_data,VP44_all_data,VP45_all_data,VP46_all_data,VP47_all_data,VP48_all_data};
    
nvp=length(allmatrices);

catalldata = cat(1,allmatrices{:});

ntrials_sub=length(VP01_all_data);

% die betreffenden Spalten auswählen und in einer neuen Datei abspeichern
% Spalte 1,3, 4 und 7 für choices
choices = catalldata(:,[1 5:6 9]);
% -> choices mit 4 Spalten erstellt
% Spalte 1: VPname
% Spalte 2: Stimulus
% Spalte 3: Choice
% Spalte 4: Feedback - nur Platzhalter, wird ersetzt durch no_choice

% replace "9" in column 3 with NaN
ntrials = length(choices(:,1));
for i = 1:ntrials
    if choices{i,3}==9;
        choices{i,3} = NaN;
        choices{i,4} = NaN;
    elseif choices{i,3}==1;
        choices{i,4} = 2; %not chosen option
    elseif choices{i,3}==2;
        choices{i,4} = 1; %not chosen option
    end
end

% choices now comprises 4 columns
% Spalte 1: VPname (id)
% Spalte 2: Stimulus
% Spalte 3: Choice
% Spalte 4: Nonchoice

% feedback (Spalte 9)
feedback = catalldata(:,[1 9 10]); %letzte Spalte nur Platzhalter

% -> feedback mit 2 Spalten erstellt
% Spalte 1: file
% Spalte 2: feedback

% replace "9" in feedback with NaN
for i = 1:ntrials
    if feedback{i,2}==9;
        feedback{i,2} = NaN;
        feedback{i,3} = NaN;
   elseif feedback{i,2}==11; % presented positive FB
        feedback{i,2} = 1;
        feedback{i,3} = 1;
   elseif feedback{i,2}==1; % omitted FB (i.e. positive outcome)
        feedback{i,2} = 1;
        feedback{i,3} = 0;
   elseif feedback{i,2}==-1; % omitted FB (i.e. negative outcome)
        feedback{i,2} = 0;
        feedback{i,3} = 0;
   elseif feedback{i,2}==-11; % presented negative FB
        feedback{i,2} = 0;
        feedback{i,3} = 1;
   end
end

% feedback now comprises 3 columns
% Spalte 1: VPname (id)
% Spalte 2: valence (1 = positive/omitted negative, 0 = negative/ omitted positive)
% Spalte 3: fb "modality" (1 = present, 0 = omitted)


% To check which stimulus has which reward probabaility, I also create the
% following table
stim_rew_prob=table(catalldata(:,[1]), catalldata(:,[5]), feedback(:,[2]), feedback(:,[3]), choices(:,[3]),catalldata(:,[7]),...
    'VariableNames', {'id', 'stim', 'feedback', 'modality', 'choice', 'correct'});

% create datastrucs for choices, non-choices, stimuli, feedback, feedback
% type and filenames

% create index to separate rows by participant
i = 480:480:ntrials; % i starts at 480 (trials per participant) and adds 480 as long as the total no of trials is not reached
j = 1:480:ntrials;  % j starts at 1 and adds 480 as long as the total no of trials is not reached

% Preallocate the cell arrays
datastruc_stimuli = cell(1, length(i)); 
datastruc_choice = cell(1, length(i)); 
datastruc_not_choice = cell(1, length(i));
datastruc_feedback = cell(1, length(i)); 
datastruc_feedback_type = cell(1, length(i));

for k = 1:length(i)
    
    filenames{k,1} = choices(j(k), 1);
    
    datastruc_stimuli{k} = choices(j(k):i(k), 2);
    datastruc_choice{k} = choices(j(k):i(k), 3);
    datastruc_not_choice{k} = choices(j(k):i(k), 4);
    
    datastruc_feedback{k} = feedback(j(k):i(k), 2);
    datastruc_feedback_type{k} = feedback(j(k):i(k), 3);
    
end

datastruc_stimuli = horzcat(datastruc_stimuli{:});
datastruc_stimuli = transpose(datastruc_stimuli);

datastruc_choice = horzcat(datastruc_choice{:});
datastruc_choice = transpose(datastruc_choice);


datastruc_not_choice = horzcat(datastruc_not_choice{:});
datastruc_not_choice = transpose(datastruc_not_choice);

datastruc_feedback = horzcat(datastruc_feedback{:});
datastruc_feedback = transpose(datastruc_feedback);

datastruc_feedback_type = horzcat(datastruc_feedback_type{:});
datastruc_feedback_type = transpose(datastruc_feedback_type);


save('interim_datasets\datastruc_stimuli.mat', 'datastruc_stimuli');
save('interim_datasets\datastruc_choice.mat', 'datastruc_choice');
save('interim_datasets\datastruc_not_choice.mat', 'datastruc_not_choice');
save('interim_datasets\datastruc_feedback.mat', 'datastruc_feedback');
save('interim_datasets\datastruc_feedback_type.mat', 'datastruc_feedback_type');
save('interim_datasets\feedback.mat', 'feedback');
save('interim_datasets\choices.mat', 'choices');
save('interim_datasets\filenames.mat', 'filenames');
save('interim_datasets\all_vp_behav.mat');

catalldata = cell2table(catalldata);
catalldata.Properties.VariableNames = ["id","block","trialno_full","trialno_block","stim","choice","correct","response_time","feedback","money"];

% save table of complete data for accuracy analyses   
writetable(catalldata, 'interim_datasets\FBOmiss_behaviour_immediate.csv');

end