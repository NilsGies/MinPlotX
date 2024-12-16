function [OutputData,OutputVariables] = Core_Call_SF(Data,OxList,ExtFct,ElOxDataDef)
OutputData=[];
OutputVariables=[];

% XMapTools is a free software solution for the analysis of chemical maps
% Copyright © 2022-2023 University of Bern, Institute of Geological Sciences, Pierre Lanari
%
% XMapTools is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any 
% later version.
%
% XMapTools is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with XMapTools. If not, see https://www.gnu.org/licenses.

switch ExtFct
    
    case 'aluminosilicate_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = aluminosilicate_XMapT(InputData,InputVariables,ElOxDataDef);
    
    case 'amphiboleCa_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = amphiboleCa_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'amphiboleCa_Fe3_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = amphiboleCa_Fe3_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'biotite_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = biotite_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'brucite_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O','As2O3','Sb2O3','Cs2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = brucite_XMapT(InputData,InputVariables,ElOxDataDef);
    
    case 'chlorite_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = chlorite_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'chromite_XMapT'
        InputVariables = {'Cr2O3','TiO2','Al2O3','FeO','Fe2O3','MgO','NiO','ZnO'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = chromite_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'chloritoid_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = chloritoid_XMapT(InputData,InputVariables,ElOxDataDef);
    
    case 'cordierite_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = cordierite_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'cpx_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = cpx_XMapT(InputData,InputVariables,ElOxDataDef);
     
    case 'cpx_Fe3_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = cpx_Fe3_XMapT(InputData,InputVariables,ElOxDataDef);
    
    case 'epidote_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = epidote_XMapT(InputData,InputVariables,ElOxDataDef);
            
    case 'feldspar_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = feldspar_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'magnetite_XMapT'
        InputVariables = {'FeO','Cr2O3','MgO','MnO','NiO','TiO2','SiO2','Al2O3','CaO','ZnO','Fe2O3'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = magnetite_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'garnet_XMapT' 
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = garnet_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'garnet_Fe3_XMapT' 
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = garnet_Fe3_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'ilmenite_XMapT' 
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = Ilmenite_XMapT(InputData,InputVariables,ElOxDataDef);    
    
    case 'olivine_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3','NiO'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = olivine_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'opx_XMapT'   
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};   % add Cr and Ni???
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = opx_XMapT(InputData,InputVariables,ElOxDataDef);

    case 'rutile_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O','ZrO2'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = rutile_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'serpentine_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3','NiO'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = serpentine_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'spinel_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3','NiO','ZnO'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = spinel_XMapT(InputData,InputVariables,ElOxDataDef);
        
    case 'staurolite_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = staurolite_XMapT(InputData,InputVariables,ElOxDataDef);

    case 'whiteMica_XMapT'
        InputVariables = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O'};
        InputData = GenerateInputData(Data,OxList,InputVariables);
        [OutputData,OutputVariables] = whiteMica_XMapT(InputData,InputVariables,ElOxDataDef);
        

        
end



end


function [InputData] = GenerateInputData(Data,OxList,InputVariables)
% Internal subroutine that re-organize the data to fit the Input format
% of the selected external function
%
% Pierre Lanari (last edit 14.01.21)

InputData = zeros(size(Data,1),length(InputVariables));

[IsEl,IndEl] = ismember(InputVariables,OxList);

for i = 1:length(IsEl)
    if IsEl(i)
        InputData(:,i) = Data(:,IndEl(i));
    end
end
end


