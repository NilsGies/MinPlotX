%generic sulfide calculator
%last modified 17.04.2024

function [StrctFrm, apfu,options_definition]=sulfide_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
%Parameter 1
options_definition.cation_normalization.question='Do you wish to calculate on a cation basis?';
options_definition.cation_normalization.description='here text'; % optional
options_definition.cation_normalization.options={true,false}; % optional
options_definition.cation_normalization.Value=true; %default_value
%Parameter 2
options_definition.moles.question='How many cations/anions do you wish to normalize to?';
options_definition.moles.description='here text'; % optional
options_definition.moles.Value=2; %default_value
options_definition.moles.limits=[1 inf]; %optional
%Parameter 3
options_definition.As_cation.question='Do you wish to treat As as a cation?';
options_definition.As_cation.Value=true; %default_value

%% Check if the data or header is empty 
if not(exist('headers','var')) || not(exist('data','var')) || isempty(data) || isempty(headers)
    return % exit the function and return options_definition
end

%% If no options 
if not(exist('options','var')) || isempty(options)
options=options_definition; % use default values generated in the options definition block
end
%% start calculation
[m,~]=size(data); %finds the x and y size of the input data matrix

%% Molecular weights
S_mw=32.06;
Fe_mw=55.8452;
Ni_mw=58.693;
Cu_mw=63.546;
Co_mw=58.933;
Pb_mw=207.2;
Zn_mw=65.382;
As_mw=74.922; 

%% Moles of elements

MC=zeros(m,8);%create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'S'))./S_mw; % for S

%calculates for As if it is included in the analysis 
if any(strcmp(headers,'As'))
    MC(:,2)=data(:,strcmp(headers,'As'))./As_mw; %for As
end


%calculates for Fe if it is included in the analysis 
if any(strcmp(headers,'Fe'))
    MC(:,3)=data(:,strcmp(headers,'Fe'))./Fe_mw; %for Fe
end


%calculates for Ni if it is included in the analysis 
if any(strcmp(headers,'Ni'))
    MC(:,4)=data(:,strcmp(headers,'Ni'))./Ni_mw; %for Ni
end

%calculates for Cu if it is included in the analysis 
if any(strcmp(headers,'Cu'))
    MC(:,5)=data(:,strcmp(headers,'Cu'))./Cu_mw; %for Cu
end


%calculates for Co if it is included in the analysis 
if any(strcmp(headers,'Co'))
    MC(:,6)=data(:,strcmp(headers,'Co'))./Co_mw; %for Co
end

%calculates for Pb if it is included in the analysis 
if any(strcmp(headers,'Pb'))
    MC(:,7)=data(:,strcmp(headers,'Pb'))./Pb_mw; %for Pb
end

%calculates for Zn if it is included in the analysis 
if any(strcmp(headers,'Zn'))
    MC(:,8)=data(:,strcmp(headers,'Zn'))./Zn_mw; %for Zn
end

%% Normalization

%prompts the user if they wish calculate on a cation or anion basis
if options.cation_normalization.Value==true
    disp('You are calculating on a cation basis.')
  %'How many cations do you wish to normalize to?: ';
    
    CM=options.moles.Value;

    if options.As_cation.Value==true
        disp('As will be treated as a cation.')
        NF=CM./sum(MC(:,2:8),2); %normalization factor w/ As
    else
        disp('As will be treated as an anion.')
        NF=CM./sum(MC(:,3:8),2); %normalization factor w/o As
    end
else
  %'How many anions do you wish to normalize to?: ';
    AM=options.moles.Value;

    if options.As_cation.Value==true
        disp('As will be treated as a cation.')
        NF=AM./MC(:,1); %normalization factor w/o As
    else
        disp('As will be treated as an anion.')
        NF=AM./(MC(:,1)+MC(:,2)); %normalization factor w/ As
    end
end

apfu=MC.*NF; %Normalized moles = moles of elements * normalization factor

apfu=array2table(apfu,'VariableNames',{'S','As','Fe','Ni','Cu','Co','Pb','Zn'});

end






