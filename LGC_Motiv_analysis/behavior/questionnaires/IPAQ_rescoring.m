function[IPAQ_score, IPAQ, IPAQinactivity, subject_id] = IPAQ_rescoring(study_nm, subject_id, NS)
% [IPAQ_score, IPAQ, IPAQinactivity, subject_id] = IPAQ_rescoring()
% IPAQ_rescoring attempts to fix the IPAQ score for the subjects who
% misunderstood the instructions.
%
% INPUTS
% study_nm: study name (study1 by default)
%
% subject_id: list of subjects to include
%
% NS: number of subjects included
%
% OUTPUTS
% IPAQ_score: global IPAQ score
%
% IPAQ: structure with score for each activity level (low/moderate/high)
%
% IPAQinactivity: average duration of being sitted during a regular day
%
% subject_id: list of subjects included

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
if ~exist('subject_id','var') || isempty(subject_id)
    subject_id = LGCM_subject_selection(study_nm);
end
if ~exist('NS','var') || isempty(NS)
    NS = length(subject_id);
end

%% working directory
IPAQ_folder = fullfile('P:','boulot','postdoc_CarmenSandi','results',...
    'behavior','questionnaires');

%% load IPAQ data
IPAQ_table = readtable([IPAQ_folder,filesep,'IPAQ.xlsx']);
[IPAQ_high_days, IPAQ_high_time,...
    IPAQ_mod_days, IPAQ_mod_time,...
    IPAQ_low_days, IPAQ_low_time,...
    IPAQinactivity,...
    IPAQ.MET_high,...
    IPAQ.MET_mod,...
    IPAQ.MET_low,...
    IPAQ.MET] = deal(NaN(1,NS));
n_hours_threshold = 9; % threshold to use to determine that people 
% misunderstood the questionnaire and answered for 2-weeks average instead 
% of average of 1 activity period
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = find(strcmp(IPAQ_table.Q1,sub_nm));
    if (length(sub_idx) == 1 && ~ismember(sub_nm,{'100'})) ||...
            (length(sub_idx) == 2 && ismember(sub_nm,{'100'}))
        % subject 100 filled the questionnaires twice => ignore the first values
        switch sub_nm
            case '100'
                sub_idx = sub_idx(2); % ignore first occurrence
        end
        
        % extract IPAQ values for each level of effort
        % high effort
        IPAQ_high_days(iS) = IPAQ_table.Q1_1(sub_idx);
        IPAQ_high_time(iS) = IPAQ_table.Q2_1_1(sub_idx);
        % moderate effort
        IPAQ_mod_days(iS) = IPAQ_table.Q3_1(sub_idx);
        IPAQ_mod_time(iS) = IPAQ_table.Q4_1_1(sub_idx);
        % low effort
        IPAQ_low_days(iS) = IPAQ_table.Q5_1(sub_idx);
        IPAQ_low_time(iS) = IPAQ_table.Q6_1_1(sub_idx);
        % inactivity
        IPAQinactivity(iS) = IPAQ_table.Q7_1_1(sub_idx); % no need of fixing as looks clear for everyone
        
        %% correct for those who misunderstood the questionnaire and reported their
        % average level of activity over 2 weeks instead of the average duration
        % of their activity when they exert it
        if IPAQ_high_time(iS) > n_hours_threshold ||...
                IPAQ_mod_time(iS) > n_hours_threshold ||...
                IPAQ_low_time(iS) > n_hours_threshold % in this case, people reported their total level of activity over 2 weeks instead of the average duration of their activity
            % => to fix it, we are just going to divide by two to consider
            % the average level of activity over a 1 week period
            IPAQ.MET_low(iS) = 3.3*IPAQ_low_time(iS)./2;
            IPAQ.MET_mod(iS) = 4*IPAQ_mod_time(iS)./2;
            IPAQ.MET_high(iS) = 8*IPAQ_high_time(iS)./2;
        else % regular scoring
            IPAQ.MET_low(iS) = 3.3*IPAQ_low_time(iS).*IPAQ_low_days(iS);
            IPAQ.MET_mod(iS) = 4*IPAQ_mod_time(iS).*IPAQ_mod_days(iS);
            IPAQ.MET_high(iS) = 8*IPAQ_high_time(iS).*IPAQ_high_days(iS);
        end
        IPAQ.MET(iS) = IPAQ.MET_low(iS) + IPAQ.MET_mod(iS) + IPAQ.MET_high(iS);
    else
        error(['Problem with subject ',sub_nm,' IPAQ extraction: ',...
            num2str(length(sub_idx)),' occurences for the same subject']);
    end
end % subject loop

% store global score in one single variable
IPAQ_score = IPAQ.MET;
end % function