%%ilmenite Structural Formula
function [StrctFrm, apfu,options_definition]=spinel_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition

% %Parameter 1
options_definition.want_ferric_Fe3.question='Calculate Fe3+ from stoichiometry?';
options_definition.want_ferric_Fe3.options={true,false}; % optional
options_definition.want_ferric_Fe3.Value=true; %default_value

% %Parameter 2
options_definition.UseKnownFe3.question='Is Fe3+/FeTotal Ratio known?';
options_definition.UseKnownFe3.Value=false;

% %Parameter 3
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 4
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

%% new
if isfield(options,'want_ferric_Fe3') && options.want_ferric_Fe3.Value==true

    if isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
        [StrctFrm, apfu]=spinel_Fe3known(data,headers,options);
    else
        [StrctFrm, apfu]=spinel_Fe3unknown(data,headers,options);
    end
else
        options.UseKnownFe3.Value=true;
    options.Fe3_ratio.Value=0;
    [StrctFrm, apfu]=spinel_Fe3known(data,headers,options);
end


