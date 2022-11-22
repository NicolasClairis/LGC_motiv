function[blood_labelname] = blood_label(blood_varname)
% [blood_labelname] = blood_label(blood_varname)
% blood_label transforms the blood variable name blood_varname 
% entered in input into blood_label for clearer visualization on the
% figures
%
% INPUTS
% blood_varname: string with blood variable name for table/structure
% extraction
%
% OUTPUTS
% blood_labelname: string for figure display

switch blood_varname
    case 'NAD_div_NADH'
        blood_labelname = 'NAD/NADH (μM)';
    case 'NADP_div_NADPH'
        blood_labelname = 'NADP/NADPH (μM)';
    case 'total_NAD_precursors'
        blood_labelname = {'total NAD';'precursors (μM)'};
    case 'total_NAD_with_precursors'
        blood_labelname = {'total NAD';'with precursors (μM)'};
    case 'total_NAD_with_byproducts'
        blood_labelname = {'total NAD';' with all (μM)'};
    case 'total_NAD_byproducts'
        blood_labelname = {'total NAD';'by-products (μM)'};
    otherwise
        blood_m_nm_bis = strrep(blood_varname,'_',' ');
        blood_labelname = [blood_m_nm_bis,' (μM)'];
end

end % function