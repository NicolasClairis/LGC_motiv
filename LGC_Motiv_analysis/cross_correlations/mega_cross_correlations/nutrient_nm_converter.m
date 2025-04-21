function[nutrient_nm] = nutrient_nm_converter(short_nutrient_nm)
% [nutrient_nm] = nutrient_nm_converter(short_nutrient_nm)
% nutrient_nm_converter will convert short_nutrient_nm, which is the
% abbreviated version of the nutrient name, into nutrient_nm which
% is the full nutrient name.
%
% INPUTS
% short_nutrient_nm: short nutrient name (abbreviation)
%
% OUTPUTS
% nutrient_nm: full nutrient name with meaning of the abbreviation
%

switch short_nutrient_nm
    case {'calories','proteins','fat','glucids',...
            'Mg','Glc','Lactose','sugars',...
            'Cys','Asp','Glu','Gly',...
            'omega3'}
        nutrient_nm = short_nutrient_nm;
    case 'GJEnergie'
        nutrient_nm = 'Energy';
    case 'GCALZBEnergie'
        nutrient_nm = 'Energy';
    case 'GJZBEnergie'
        nutrient_nm = 'Energy';
    case 'ZW'
        nutrient_nm = 'water';
    case 'ZB'
        nutrient_nm = 'fibers';
    case 'ZM'
        nutrient_nm = 'minerals (ashes)';
    case 'ZO'
        nutrient_nm = 'organic acid';
    case 'ZA'
        nutrient_nm = 'alcohol (ethanol)';
    case 'VA'
        nutrient_nm = 'vitamin A - retinol equivalent';
    case 'VAR'
        nutrient_nm = 'vitamin A - retinol';
    case 'VAC'
        nutrient_nm = 'vitamin A - β-caroten';
    case 'VD'
        nutrient_nm = 'vitamin D - calciferole';
    case 'VE'
        nutrient_nm = 'vitamin E - tocopherolequivalent';
    case 'VEAT'
        nutrient_nm = 'vitamin E - α-tocopherol';
    case 'VK'
        nutrient_nm = 'vitamin K - phyllochinon';
    case 'VB1'
        nutrient_nm = 'vitamin B1 - Thiamin';
    case 'VB2'
        nutrient_nm = 'vitamin B2 - Riboflavin';
    case 'NA'
        nutrient_nm = 'Niacin';
    case 'NE'
        nutrient_nm = 'Niacin equivalents';
    case 'VB5'
        nutrient_nm = 'vitamin B5 - panthotenic acid';
    case 'VB6'
        nutrient_nm = 'vitamin B6';
    case 'VB7'
        nutrient_nm = 'vitamin B7/H - biotin';
    case 'FolicAcid'
        nutrient_nm = 'folic acid';
    case 'VB12'
        nutrient_nm = 'vitamin B12 - cobalamin';
    case 'VC'
        nutrient_nm = 'vitamin C - ascorbic acid';
    case 'MNA'
        nutrient_nm = 'Na';
    case 'MK'
        nutrient_nm = 'K';
    case 'MCA'
        nutrient_nm = 'Ca';
    case 'MP'
        nutrient_nm = 'P';
    case 'MS'
        nutrient_nm = 'S';
    case 'MCL'
        nutrient_nm = 'Cl';
    case 'MFE'
        nutrient_nm = 'Fe';
    case 'MZN'
        nutrient_nm = 'Zn';
    case 'MCU'
        nutrient_nm = 'Cu';
    case 'MMN'
        nutrient_nm = 'Mn';
    case 'MF'
        nutrient_nm = 'F';
    case 'MJ'
        nutrient_nm = 'I';
    case 'KAM'
        nutrient_nm = 'Mannitol';
    case 'KAS'
        nutrient_nm = 'Sorbitol';
    case 'KAX'
        nutrient_nm = 'Xylitol';
    case 'KA'
        nutrient_nm = 'Sugar Alcohol';
    case 'KMF'
        nutrient_nm = 'Fructose';
    case 'KMG'
        nutrient_nm = 'Galactose';
    case 'KM'
        nutrient_nm = 'Monosaccharides';
    case 'KDS'
        nutrient_nm = 'Saccharose';
    case 'KDM'
        nutrient_nm = 'Maltose';
    case 'KD'
        nutrient_nm = 'Disaccharides';
    case 'KPOR'
        nutrient_nm = 'Resorbable Oligosaccharides';
    case 'KPON'
        nutrient_nm = 'Non-Resorbable Oligosaccharides';
    case 'KPG'
        nutrient_nm = 'Glycogen';
    case 'KPS'
        nutrient_nm = 'Amidon';
    case 'KP'
        nutrient_nm = 'Polysaccharides';
    case 'KBP'
        nutrient_nm = 'Poly-pentoses';
    case 'KBH'
        nutrient_nm = 'Poly-hexoses';
    case 'KBU'
        nutrient_nm = 'Polyuronic Acid';
    case 'KBC'
        nutrient_nm = 'Cellulose';
    case 'KBL'
        nutrient_nm = 'Lignin';
    case 'KBW'
        nutrient_nm = 'Soluble Fibers';
    case 'KBN'
        nutrient_nm = 'Insoluble Fibers';
    case 'EILE'
        nutrient_nm = 'Isoleu';
    case 'ELEU'
        nutrient_nm = 'Leu';
    case 'ELYS'
        nutrient_nm = 'Lys';
    case 'EMET'
        nutrient_nm = 'Met';
    case 'EPHE'
        nutrient_nm = 'Phe';
    case 'ETYR'
        nutrient_nm = 'Tyr';
    case 'ETHR'
        nutrient_nm = 'Thr';
    case 'Trp'
        nutrient_nm = 'Trp';
    case 'EVAL'
        nutrient_nm = 'Val';
    case 'EARG'
        nutrient_nm = 'Arg';
    case 'EHIS'
        nutrient_nm = 'His';
    case 'EEA'
        nutrient_nm = 'Essential Amino-Acids';
    case 'EALA'
        nutrient_nm = 'Ala';
    case 'EPRO'
        nutrient_nm = 'Pro';
    case 'ESER'
        nutrient_nm = 'Ser';
    case 'ENA'
        nutrient_nm = 'Non-essential Amino-Acids';
    case 'EH'
        nutrient_nm = 'Uric Acid';
    case 'EP'
        nutrient_nm = 'Purin';
    case 'F40'
        nutrient_nm = 'Butanoic Acid';
    case 'F60'
        nutrient_nm = 'Hexanoic Acid';
    case 'F80'
        nutrient_nm = 'Octanoic Acid';
    case 'F100'
        nutrient_nm = 'Decanoic Acid';
    case 'F120'
        nutrient_nm = 'Dodecanoic Acid';
    case 'F140'
        nutrient_nm = 'Tetradecanoic Acid';
    case 'F150'
        nutrient_nm = 'Pentadecanoic Acid';
    case 'F160'
        nutrient_nm = 'Hexadecanoic Acid';
    case 'F170'
        nutrient_nm = 'Heptadecanoic Acid';
    case 'F180'
        nutrient_nm = 'Octadecanoic Acid';
    case 'F200'
        nutrient_nm = 'Eicosanoic Acid';
    case 'F220'
        nutrient_nm = 'Decanoic Acid';
    case 'F240'
        nutrient_nm = 'Tetracosanoic Acid';
    case 'FS'
        nutrient_nm = 'Saturated Fatty Acids';
    case 'F141'
        nutrient_nm = 'Tetradecenoic Acid';
    case 'F151'
        nutrient_nm = 'Pentadecenoic Acid';
    case 'F161'
        nutrient_nm = 'Hexadecenoic Acid';
    case 'F171'
        nutrient_nm = 'Heptadecenoic Acid';
    case 'F181'
        nutrient_nm = 'OCtadecenoic Acid';
    case 'F201'
        nutrient_nm = 'Eicosenoic Acid';
    case 'F221'
        nutrient_nm = 'Decosenoic Acid';
    case 'F241'
        nutrient_nm = 'Tetracosenoic Acid';
    case 'FU'
        nutrient_nm = 'Mono-insaturated Fatty Acids';
    case 'F162'
        nutrient_nm = 'Hexadecadienoic Acid';
    case 'F164'
        nutrient_nm = 'Hexadecatetraenoic Acid';
    case 'F182'
        nutrient_nm = 'Octadecadienoic Acid';
    case 'F183'
        nutrient_nm = 'Octadecatrienoic Acid';
    case 'F184'
        nutrient_nm = 'Octradecatetraenoic Acid';
    case 'F193'
        nutrient_nm = 'Non-adecatrienoic Acid';
    case 'F202'
        nutrient_nm = 'Eicosadienoic Acid';
    case 'F203'
        nutrient_nm = 'Eicosatrienoic Acid';
    case 'F204'
        nutrient_nm = 'Eicosatetraenoic Acid';
    case 'F205'
        nutrient_nm = 'Eicosapentaenoic Acid';
    case 'F222'
        nutrient_nm = 'Docosadienoic Acid';
    case 'F223'
        nutrient_nm = 'Docosatrienoïque Acid';
    case 'F224'
        nutrient_nm = 'Docosatetraenoic Acid';
    case 'F225'
        nutrient_nm = 'Docosapentaenoic Acid';
    case 'F226'
        nutrient_nm = 'Docosahexaenoic Acid';
    case 'FP'
        nutrient_nm = 'Polyinsaturated Fatty Acids';
    case 'FK'
        nutrient_nm = 'Short-Chain Fatty Acids';
    case 'FM'
        nutrient_nm = 'Middle-Chain Fatty Acids';
    case 'FL'
        nutrient_nm = 'Long-Chain Fatty Acids';
    case 'FO6'
        nutrient_nm = 'Omegas 6';
    case 'FG'
        nutrient_nm = 'Glycerin & Lipids';
    case 'FC'
        nutrient_nm = 'Cholesterol';
    case 'GFPSPolyins'
        nutrient_nm = 'P/S ratio';
    case 'GKB'
        nutrient_nm = 'Bread';
    case 'GMKO'
        nutrient_nm = 'Salt';
    otherwise
        error([short_nutrient_nm,' not ready for conversion']);
end % nutrient name

end % function