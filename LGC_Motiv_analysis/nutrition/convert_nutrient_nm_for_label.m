function[nutrient_nm2] = convert_nutrient_nm_for_label(nutrient_nm)
% [nutrient_nm2] = convert_nutrient_nm_for_label(nutrient_nm);
% convert_nutrient_nm_for_label will convert the long label from
% nutrient_nm into something nicer to read in a figure.
%
% INPUTS
% nutrient_nm: original nutrient name
%
% OUTPUTS
% nutrient_nm2: nutrient name "cleaned"

%% initialize variable
nutrient_nm2 = nutrient_nm;
%% convert portions
% micrograms
nutrient_nm2 = strrep(nutrient_nm2,'_ugPerWeek',' (μg/week)');
% milligrams
nutrient_nm2 = strrep(nutrient_nm2,'__mg_portion__PerWeek',' (mg/week)');
nutrient_nm2 = strrep(nutrient_nm2,'_mg_portion__PerWeek',' (mg/week)');
nutrient_nm2 = strrep(nutrient_nm2,'_mgPerWeek',' (mg/week)');
% grams
nutrient_nm2 = strrep(nutrient_nm2,'___g_portion__PerWeek',' (g/week)');
nutrient_nm2 = strrep(nutrient_nm2,'__g_portion__PerWeek',' (g/week)');

%% modify weird labels into sth understandable
nutrient_nm2 = strrep(nutrient_nm2,'KPON','Non-Resorbable Oligosaccharides');

% amino-acids
nutrient_nm2 = strrep(nutrient_nm2,'EEA','Essential amino-acids');
nutrient_nm2 = strrep(nutrient_nm2,'ENA','Non-essential amino-acids');
nutrient_nm2 = strrep(nutrient_nm2,'EILE','Isoleucin');
nutrient_nm2 = strrep(nutrient_nm2,'ELEU','Leucin');
nutrient_nm2 = strrep(nutrient_nm2,'ELYS','Lysin');
nutrient_nm2 = strrep(nutrient_nm2,'EMET','Methionin');
nutrient_nm2 = strrep(nutrient_nm2,'EPHE','Phenylalanin');
nutrient_nm2 = strrep(nutrient_nm2,'ETYR','Tyrosin');
nutrient_nm2 = strrep(nutrient_nm2,'ETHR','Threonin');
nutrient_nm2 = strrep(nutrient_nm2,'EVAL','Valine');
nutrient_nm2 = strrep(nutrient_nm2,'EARG','Arginin');
nutrient_nm2 = strrep(nutrient_nm2,'EHIS','Histidin');
nutrient_nm2 = strrep(nutrient_nm2,'EALA','Alanin');
nutrient_nm2 = strrep(nutrient_nm2,'EPRO','Prolin');
nutrient_nm2 = strrep(nutrient_nm2,'ESER','Serin');
nutrient_nm2 = strrep(nutrient_nm2,'EP','Purin');
nutrient_nm2 = strrep(nutrient_nm2,'Trp','Tryptophan');
% minerals
nutrient_nm2 = strrep(nutrient_nm2,'MNA','Na');
nutrient_nm2 = strrep(nutrient_nm2,'MK','Potassium K+');
nutrient_nm2 = strrep(nutrient_nm2,'MCA','Calcium');
nutrient_nm2 = strrep(nutrient_nm2,'MP','P');
nutrient_nm2 = strrep(nutrient_nm2,'MS','S');
nutrient_nm2 = strrep(nutrient_nm2,'MCL','Cl');
nutrient_nm2 = strrep(nutrient_nm2,'MFE','Fe');
nutrient_nm2 = strrep(nutrient_nm2,'MZN','Zn');
nutrient_nm2 = strrep(nutrient_nm2,'MCU','Cu');
nutrient_nm2 = strrep(nutrient_nm2,'MMN','Mn');
nutrient_nm2 = strrep(nutrient_nm2,'MF','F');
nutrient_nm2 = strrep(nutrient_nm2,'MJ','I');
% niacin-related
nutrient_nm2 = strrep(nutrient_nm2,'NA','niacin');
nutrient_nm2 = strrep(nutrient_nm2,'NE','niacin equivalents');
% Other
nutrient_nm2 = strrep(nutrient_nm2,'EH','Uric acid');
nutrient_nm2 = strrep(nutrient_nm2,'GKB','Bread');
nutrient_nm2 = strrep(nutrient_nm2,'GMKO','Salt');
nutrient_nm2 = strrep(nutrient_nm2,'FC','Cholesterol');
nutrient_nm2 = strrep(nutrient_nm2,'F06','Omegas 6');
nutrient_nm2 = strrep(nutrient_nm2,'KPG','Glycogen');
nutrient_nm2 = strrep(nutrient_nm2,'KPS','Starch (Amidon)');

%% remove all underscore and replace by a space
nutrient_nm2 = strrep(nutrient_nm2,'_',' ');

end % function