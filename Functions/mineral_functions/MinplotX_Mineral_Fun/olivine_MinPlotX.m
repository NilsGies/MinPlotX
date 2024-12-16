%olivine
function [StrctFrm,apfu,options_definition]=olivine_MinPlotX(data,headers,options)

%% define empty output variables
StrctFrm=[];
apfu=[];

%% options definition
%Parameter 1
options_definition.want_ferric_Fe3.question='Calculate Fe3+ from stoichiometry?';
options_definition.want_ferric_Fe3.Value=true; %default_value
options_definition.want_ferric_Fe3.options={true,false}; % optional

%Parameter 2
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
    %old  [StrctFrm,apfu]=olivine_fe3(data,headers,options); %calls the olivine (w/Fe3+) calculation, but only outputs apfu and no endmembers
    if isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true % known fe3 garnet_skarn_Fe3known
        [StrctFrm,apfu]=olivine_Fe3known(data,headers,options);
    else
        [StrctFrm,apfu]=olivine_Fe3unknown(data,headers,options);
    end
else
    %If an FeO only calculation was selected, the following procedure
    options.UseKnownFe3.Value=true;
    options.Fe3_ratio.Value=0;
    [StrctFrm,apfu]=olivine_Fe3known(data,headers,options);
end
end