function[runs] = runs_definition(study_nm, sub_nm)
% [runs] = runs_definition(study_nm, sub_nm)
%runs_definition will define the number and order of the physical and
%mental effort task for the study and subject entered in input.
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% OUTPUTS
% runs: structure with number of runs for each task and also the order of
% the runs

switch study_nm
    case 'fMRI_pilots'
        switch sub_nm
            case 'pilot_s1'
                runs.nb_runs.Ep = 1;
                runs.nb_runs.Em = 1;
                runs.tasks = {'Ep','Em'};
            case 'pilot_s2'
                runs.nb_runs.Ep = 1;
                runs.nb_runs.Em = 0;
                runs.tasks = {'Ep'};
        end
    case 'study1'
        error('need to attribute runs');
    case 'study2_clinical'
        error('need to attribute runs');
    otherwise
        error('study not recognized');
end
        

end % function