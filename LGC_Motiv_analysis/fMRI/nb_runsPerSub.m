function[nb_runs] = nb_runsPerSub(study_nm, sub_nm)
%[nb_runs] = nb_runsPerSub(study_nm, sub_nm)
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% OUTPUTS
% nb_runs: number of runs

switch study_nm
        case 'fMRI_pilots'
            switch sub_nm
                case {'pilot_s1','pilot_s3'}
                    nb_runs = 2;
                case 'pilot_s2'
                    nb_runs = 1;
                otherwise
                    nb_runs = 4;
            end
        case 'study1'
            switch sub_nm
                case '074' % CID074: ignore first run which crashed at the beginning due to high voltage
                    nb_runs = 3;
                otherwise
                    nb_runs = 4;
            end
end
    
end % function