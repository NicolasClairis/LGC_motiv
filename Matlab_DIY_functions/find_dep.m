% Shortcut to find dependencies of a script, function, variable, etc.
% will give the list of all the scripts in the current folder (or any other folder 
% you define) which call this expression (whether used or in comments)
%
% Written by Jules and NicoC 16/01/17

expression = input(['Name of the script? \n Please add a ( if it''s a function',...
    'and you want to ignore comments. \n Example: script( \n'],'s');
where_to_check = input('Check this folder (1) or elsewhere (2)?');
if where_to_check == 1
    eval(['! findstr /R /S "',expression,'" ',pwd,'\*'])
elseif where_to_check == 2
    folder_to_check = input('Name of folder where to check?');
    eval(['! findstr /R /S "',expression,'" ',folder_to_check,'\*'])
else
    return;
end