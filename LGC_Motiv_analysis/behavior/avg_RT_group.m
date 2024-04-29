function[avg_RT] = avg_RT_group(study_nm, subject_id, condition)
% [avg_RT] = avg_RT_group(study_nm, subject_id, condition)
% avg_RT_group will extract the average RT per task or across tasks for the
% list of subjects entered in inputs by using the avg_RT_perSub.m function.
%
% INPUTS
% study_nm: study name
%
% subject_id: subject list
%
% condition: condition to use
%
% OUTPUTS
% avg_RT: structure with average RT
%
% See also avg_RT_perSub.m

NS = length(subject_id);
[avg_RT.Ep, avg_RT.Em, avg_RT.EpEm] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    [avg_RT_tmp] = avg_RT_perSub(study_nm, sub_nm, condition);
    avg_RT.Ep(iS) = avg_RT_tmp.Ep;
    avg_RT.Em(iS) = avg_RT_tmp.Em;
    avg_RT.EpEm(iS) = avg_RT_tmp.EpEm;
end % subject loop

end % function