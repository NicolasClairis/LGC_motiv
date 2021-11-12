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
            case 'pilot_s3'
                runs.nb_runs.Ep = 1;
                runs.nb_runs.Em = 1;
                runs.tasks = {'Em','Ep'};
        end
    case 'study1'
        switch sub_nm
            case {'CID074'}
                runs.nb_runs.Ep = 2;
                runs.nb_runs.Em = 1;
                runs.tasks = {'Ep','Em','Ep'};
            case {'CID036'}
                runs.nb_runs.Ep = 2;
                runs.nb_runs.Em = 2;
                runs.tasks = {'Ep','Em','Ep','Em'};
        end
    case 'study2_clinical'
        error('need to attribute runs');
    otherwise
        error('study not recognized');
end
        

end % function