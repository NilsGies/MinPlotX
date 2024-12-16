%garnet
function [StrctFrm,apfu,options_definition]=garnet_MinPlotX(data,headers,options)

%% define empty output variables
StrctFrm=[];
apfu=[];

%% options definition
%Parameter 1
options_definition.want_ferric_Fe3.question='Calculate Fe3+ from stoichiometry?';
options_definition.want_ferric_Fe3.Value=true; %default_value
options_definition.want_ferric_Fe3.options={true,false}; % optional
%Parameter 2
options_definition.Ti_Endmembers.question='Do you want to correct for OH=2-2Ti?';
options_definition.Ti_Endmembers.Value=false; %default_value
options_definition.Ti_Endmembers.options={true,false}; % optional

%Parameter 3
options_definition.UseKnownFe3.question='Do you want to use a known Fe3 ratio';
options_definition.UseKnownFe3.Value=false;

%Parameter 4
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 5
options_definition.RecalcTotalFe.question='Recalculate FeO and Fe2O3?';
options_definition.RecalcTotalFe.Value=false;

%% Check if the data or header is empty
if not(exist('headers','var')) || not(exist('data','var')) || isempty(data) || isempty(headers)
    return % exit the function and return options_definition
end

%% If no options
if not(exist('options','var')) || isempty(options)
    options=options_definition; % use default values generated in the options definition block
end
%% start calculation
if isfield(options,'want_ferric_Fe3') && options.want_ferric_Fe3.Value==true
    if isfield(options,'Ti_Endmembers') && options.Ti_Endmembers.Value==true
        %Fe3+ calculation for ugrandite garnets
        if isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true % known fe3 garnet_skarn_Fe3known
            [StrctFrm,apfu]=garnet_skarn_Fe3known(data,headers,options);
        else
            [StrctFrm,apfu]=garnet_skarn_Fe3unknown(data,headers,options);
        end
    else
        if isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true % known fe3 garnet_skarn_Fe3known
            [StrctFrm,apfu]=garnet_Fe3known(data,headers,options);
        else
            [StrctFrm,apfu]=garnet_Fe3unknown(data,headers,options);
        end
    end
else
    %If an FeO only calculation was selected, the following procedure
    options.Fe3_ratio.Value=0;
    options.UseKnownFe3.Value=true;
    [StrctFrm,apfu]=garnet_Fe3known(data,headers,options);
end
end