function [sub_cell] = convert_sub_id_from_num_to_cell(num_data)
%[sub_cell] = convert_sub_id_from_num_to_cell(num_data)
% convert_sub_id_from_num_to_cell will convert vector of subjects with
% numeric values to a cell with the subject ID with the correct number of
% 0 before each number to match size
%
% INPUTS
%  num_data: 1xNS vector with the subject data into numeric values
%
% OUTPUTS
% sub_cell: 1xNS cell with the subject data into cell format
%


%% extract number of subjects
NS = size(num_data, 2);
%% loop through subjects
sub_cell = cell(1,NS);
for iS = 1:NS
    sub_nber = num_data(iS);
    if sub_nber < 10
        sub_cell{iS} = ['00',num2str(sub_nber)];
    elseif sub_nber >= 10 && num_data(iS) <= 100
        sub_cell{iS} = ['0',num2str(sub_nber)];
    elseif sub_nber >= 100
        sub_cell{iS} = num2str(sub_nber);
    end % adaptation according to number of subject
end % subject loop

end % function