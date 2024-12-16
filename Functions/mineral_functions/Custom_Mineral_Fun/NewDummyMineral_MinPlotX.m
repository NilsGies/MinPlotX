function [StrctFrm, apfu,options_definition]=NewDummyMineral_MinPlotX(data,headers,options)
%% this function can be modified to create individual mineral calculation routines for MinPlotX
% MinPlotX Provides the chemichal data as Matrix, the header of the matrix
% as cell array and has an additional options structure argument which can
% be send from the GUI. As an output MinPlotX requires two tables StrctFrm
% and apfu (if only one of them is generated return this table twice)

% Calculation options and can be defined in the options_definition.  These
% also allow also to activate the Fe and Mn options panels. To implement
% the new mineral function into the MinPlotX GUI the mineralname and
% function name needs to be added in the MineralCalcFunctions.csv file.

% Note: this does not work in the compiled version. Please contact us if
% you do not have access to Matlab and want to implement a new mineral calculation routine.

%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
options_definition=struct();

%% GUI options 
% Parameter to active options in GUI
options_definition.GUI_Fe3_Option.Value=true;
options_definition.GUI_Fe3fromstoichiometry.Value=true;
options_definition.GUI_tetraFe3_Option.Value=true;
options_definition.GUI_Mn3_Option.Value=true;

%% Optional Calculation Parameters to get values or settings from the user
%Parameter 1
options_definition.ThisIsADummy.question='Here we can ask for settings?';
options_definition.ThisIsADummy.options={true,false}; % optional
options_definition.ThisIsADummy.Value=true; %default_value

%Parameter 2
options_definition.WhatNumber.question='Here we can ask for a number?';
options_definition.WhatNumber.limits=[1 100]; % optional
options_definition.WhatNumber.Value=50; %default_value

%% Implemented Fe and Mn options and how to read the parameters
%Parameter Fe
options_definition.want_ferric_Fe3.question='Calculate Fe3+ from stoichiometry?';
options_definition.want_ferric_Fe3.options={true,false}; % optional
options_definition.want_ferric_Fe3.Value=true; %default_value

%Parameter Fe
options_definition.UseKnownFe3.question='Is Fe3+/FeTotal Ratio known?';
options_definition.UseKnownFe3.Value=false;

%Parameter Fe
options_definition.Fe3_ratio.question='Fe3+/FeTotal Ratio?';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter Fe 
options_definition.UseKnowntetra_Fe3.question='Use known Fe3+ Tetra Ratio?';
options_definition.UseKnowntetra_Fe3.Value=false;

%Parameter Fe
options_definition.tetra_Fe3.question='Enter tetra_Fe3'; 
options_definition.tetra_Fe3.Value=0;
options_definition.tetra_Fe3.limits=[0 1]; 

%Parameter Fe
options_definition.RecalcTotalFe.question='Recalculate FeO and Fe2O3?';
options_definition.RecalcTotalFe.Value=false;

%Parameter Mn
options_definition.UseKnownMn3.question='Do you want to use a known Mn3 ratio';
options_definition.UseKnownMn3.Value=false;

%Parameter Mn
options_definition.Mn3_ratio.question='Enter Mn3_ratio';
options_definition.Mn3_ratio.Value=1;
options_definition.Mn3_ratio.limits=[0 1];

%Parameter Mn
options_definition.RecalcTotalMn.question='Recalculate MnO and Mn2O3?';
options_definition.RecalcTotalMn.Value=false;

%% Check if the data or header is empty 
if not(exist('headers','var')) || not(exist('data','var')) || isempty(data) || isempty(headers)
    return % exit the function and return options_definition
end

%% If no options use the default options
if not(exist('options','var')) || isempty(options)
options=options_definition; % use default values generated in the options definition block
end

%% start calculation
% here you can implement your own calculation routine

% Example how to use calculation options
if isfield(options,'ThisIsADummy') && options.ThisIsADummy.Value==true
%Implement your calculation with Fe3 here
CalculationResults=data.*options_definition.WhatNumber.Value;
else
%Implement your calculation without Fe3 here
CalculationResults=data./options_definition.WhatNumber.Value;

end

%% create two tables that will displayed in the MinPlotX GUI
% This dummy function will output the input data mulltiplied or devided by the 
% numeric variable WhatNumber as output tables
StrctFrm=array2table(CalculationResults,'VariableNames',headers);
apfu=StrctFrm;
end

