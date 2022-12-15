function[con_name] = conNumber2conName(iCon)
% [con_name] = conNumber2conName(iCon)
% conNumber2conName will convert contrast number iCon into a string with
% the appropriate number of zeros before to match with SPM standards.
%
% INPUTS
% iCon: contrast number
%
% OUTPUTS
% con_name: string with contrast name


if iCon < 10
    con_name = ['000',num2str(iCon)];
elseif iCon >= 10 && iCon < 100
    con_name = ['00',num2str(iCon)];
elseif iCon >= 100 && iCon < 1000
    con_name = ['0',num2str(iCon)];
elseif iCon >= 1000 && iCon < 10000
    con_name = num2str(iCon);
end

end % function