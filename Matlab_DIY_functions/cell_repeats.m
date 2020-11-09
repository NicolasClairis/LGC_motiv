function [ cells_repeated, idx_repetitions ] = cell_repeats( cell_list )
%[ cells_repeated, idx_repetitions ] = cell_repeats( cell_list )
% cell_repeats checks if there are some elements that are repeated inside
% cell_list. If it is the case, it returns the name of the elements in
% cell_list that are repeated and also the index inside the cell_list
% vector of strings of the repeated elements.
%
% INPUTS
% cell_list: cell containing the string-elements that you want to check
%
% OUTPUTS
% cells_repeated: cell containing the string-elements that are repeated
% more than once inside cell_list
%
% idx_repetitions: list of all the indexes inside cell_list of the repeated
% elements
%
% Written by N.Clairis - 19/12/18

% initialize var of interest
cells_repeated = {};

% extract total number of elements in cell_list
n_elements = length(cell_list);

% by default all variables set to zero (no repetitions)
repeats = zeros(n_elements,1);
% loop through elements of cell_list
for iCell = 1:n_elements
    
    
    % for each element, check if it is repeated or not in the string
    check_repeats = strcmp( cell_list{iCell}, cell_list);
    if sum(check_repeats) > 1
        repeats(check_repeats) = 1;
        
        % if element is repeated and not already included in
        % cells_repeated, add it to the list
        if sum(strcmp(cells_repeated, cell_list{iCell})) == 0
            cells_repeated = [cells_repeated, cell_list{iCell}];
        end 
        
    end
    
end

%% extract all the cells_list indexes of repetitions
idx_repetitions = find(repeats ~= 0);


end

