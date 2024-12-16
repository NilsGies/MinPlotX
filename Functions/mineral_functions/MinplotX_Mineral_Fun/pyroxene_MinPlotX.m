%pyroxene
function [StrctFrm,apfu,options_definition]=pyroxene_MinPlotX(data,headers,options)

%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
%Parameter 1
options_definition.want_ferric_Fe3.question='Calculate Fe3+ from stoichiometry?';
options_definition.want_ferric_Fe3.options={true,false}; % optional
options_definition.want_ferric_Fe3.Value=true; %default_value
%Parameter 2
options_definition.HighPT_Endmembers.question='Use the calculation scheme for high-pressure/temperature pyroxene?';
options_definition.HighPT_Endmembers.description='For K-, Ca-Eskola-, & Ca-Tschermaks-rich pyroxenes';
options_definition.HighPT_Endmembers.Value=false; %default_value
options_definition.HighPT_Endmembers.options={true,false}; % optional
%Parameter 3
options_definition.Na_Endmembers.question='Do you wish to include Na-Fe3+-Cr endmembers?';
options_definition.Na_Endmembers.Value=true; %default_value
options_definition.Na_Endmembers.options={true,false}; % optional

%Parameter 4 @JW
options_definition.UseKnownFe3.question='Is Fe3+/FeTotal Ratio known?';
options_definition.UseKnownFe3.Value=false;

%Parameter 5 @JW
options_definition.Fe3_ratio.question='Fe3+/FeTotal Ratio?';
options_definition.Fe3_ratio.Value=0; %default_value
options_definition.Fe3_ratio.limits=[0 1]; % optional

%Parameter 6
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

%if you do not want ferric iron, the following procedure occurs
if isfield(options,'want_ferric_Fe3') && options.want_ferric_Fe3.Value==true
    if isfield(options,'HighPT_Endmembers') && options.HighPT_Endmembers.Value==true
        if isfield(options,'UseKnownFe3') &&  options.UseKnownFe3.Value==true
            [StrctFrm,apfu]=pyroxene_fe3known_HP(data,headers,options);
        else
            [StrctFrm,apfu]=pyroxene_fe3unknown_HP(data,headers,options);
        end
    else
        if isfield(options,'UseKnownFe3') &&  options.UseKnownFe3.Value==true % known fe3 garnet_skarn_Fe3known
            [StrctFrm,apfu]=pyroxene_fe3known(data,headers,options);
        else
            [StrctFrm,apfu]=pyroxene_fe3unknown(data,headers,options);
        end
    end
else
    if isfield(options,'HighPT_Endmembers') && options.HighPT_Endmembers.Value==true
        %If an FeO only calculation was selected, the following procedure
        %occurs
        options.UseKnownFe3.Value=true;
        options.Fe3_ratio.Value=0;
        [StrctFrm,apfu]=pyroxene_fe3known_HP(data,headers,options);
    else
        %If an FeO only calculation was selected, the following procedure
        %occurs
        options.UseKnownFe3.Value=true;
        options.Fe3_ratio.Value=0;
        [StrctFrm,apfu]=pyroxene_fe3known(data,headers,options);
    end
end