function[task_fullName] = task_fullName_extraction(task_short_nm)
% [task_fullName] = task_fullName_extraction(task_short_nm)
% task_fullName_extraction converts Ep/Em abbreviations into
% "physical"/"mental".
%
% INPUTS
% task_short_nm: task short name as a string 'Ep' for physical effort task
% or 'Em' for mental effort task
%
% OUTPUTS
% task_fullName: task full name as a string 'physical' for physical effort
% task and 'mental' for mental effort task
%

switch task_short_nm
    case 'Ep'
        task_fullName = 'physical';
    case 'Em'
        task_fullName = 'mental';
    otherwise
        error(['Task short name ',task_short_nm,...
            ' is not Ep or Em so there might be a mistake here..']);
end

end % function