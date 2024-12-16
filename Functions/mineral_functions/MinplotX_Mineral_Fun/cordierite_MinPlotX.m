%Cordierite structural formula 
% last modified 17.07.2024

function [StrctFrm, apfu,options_definition]=cordierite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
options_definition=struct();
%Parameter 1
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

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

[m,~]=size(data); %finds the x and y size of the input data matrix

Opfu=18; %oxygens per formula unit

%% Checks for Fe3+/Fetotal ratio
%assigns Fe3/Fetotal if it is included in the analysis

if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=ones(m,1); %Fe3+/FeTotal ratio is 1´
    options.RecalcTotalFe.Value=false;
end



%% Molecular weights

SiO2_mw=60.083;
TiO2_mw=79.865;
Al2O3_mw=101.961;
Cr2O3_mw=151.989;
Fe2O3_mw=159.6874;
Mn2O3_mw=157.873;
Y2O3_mw=225.809;
NiO_mw=74.692;
ZnO_mw=81.381;
FeO_mw=71.8442;
MnO_mw=70.937;
MgO_mw=40.304;
CaO_mw=56.0774;
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
F_mw=18.998;
Cl_mw=35.45;

%% Moles of oxides
%create columns of zeros if optional data are not included
MC=zeros(m,10);

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2

%calculates for TiO2 if it is included in the analysis 
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end 

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

%calculates for FeO and Fe2O3

if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) && isfield(options,'RecalcTotalFe') && options.RecalcTotalFe.Value==true

    %recalculate total FeO based on FeO+Fe2O3:
    FeOT=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));

    %first calculate the weight percent of Fe2O3 and FeO from FeO total
    Fe3_wtper(:,1)=(FeOT.*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
    Fe2_wtper(:,1)=FeOT-(FeOT.*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    %calculate moles of cations
    MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
    MC(:,5)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,4)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,5)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO'))>=1 %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,5)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,5)=Fe2_wtper(:,1)./FeO_mw;

    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

    %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
    %combined with MnO
    MC(:,6)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,6)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,6)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        % if MnO total or Mn2O3 total is not included, Mn = 0
    end
end

MC(:,7)=(data(:,strcmp(headers,'MgO'))./MgO_mw); %for MgO

%calculates for CaO if it is included in the analysis 
if any(strcmp(headers,'CaO'))
    MC(:,8)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end 

%calculates for Na2O if it is included in the analysis 
if any(strcmp(headers,'Na2O'))
    MC(:,9)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end 

%calculates for K2O if it is included in the analysis 
if any(strcmp(headers,'K2O'))
    MC(:,10)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
end 

%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Fe2O3
O2(:,5)=MC(:,5); %for FeO
O2(:,6)=MC(:,6); %for MnO
O2(:,7)=MC(:,7); %for MgO
O2(:,8)=MC(:,8); %for CaO
O2(:,9)=MC(:,9)./2; %for Na2O
O2(:,10)=MC(:,10)./2; %for K2O

O2total=sum(O2,2); %O2 totals
MCnormfact=Opfu./O2total; %normalization factor

apfu=MCnormfact.*MC; %creates a matrix of normalized cations
apfu(:,11)=sum(apfu,2); %calculations the total, which should be close to 10

%% structural formula

StrctFrm=zeros(m,16);

%T1 site
StrctFrm(:,1)=apfu(:,1); %Si (T1)
StrctFrm(:,2)=apfu(:,4); %Fe3 (T1)

%Al (T1)
for c=1:m
    if 6-StrctFrm(c,1)-StrctFrm(c,2) > apfu(c,3) %if 6-Si-Fe3 > Al
        StrctFrm(c,3)=apfu(c,3); %all Al goes onto T1 (probably impossible)
    else
        StrctFrm(c,3)=6-StrctFrm(c,1)-StrctFrm(c,2);
    end
end

StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3); %Sum of T1

%T2 site
StrctFrm(:,5)=apfu(:,3)-StrctFrm(:,3); %Altotal-Al(T1)=Al (T2)
StrctFrm(:,6)=apfu(:,2); %Ti (T2)
StrctFrm(:,7)=StrctFrm(:,6)+StrctFrm(:,5); % sum of T2 site

%Octahedral Site

StrctFrm(:,8)=apfu(:,5); %Fe (Oct)
StrctFrm(:,9)=apfu(:,6); %Mn (Oct)
StrctFrm(:,10)=apfu(:,7); %Mg (Oct)
StrctFrm(:,11)=StrctFrm(:,10)+StrctFrm(:,9)+StrctFrm(:,8); %Sum of Oct site

%A site
StrctFrm(:,12)=apfu(:,8); %Ca (A)
StrctFrm(:,13)=apfu(:,9); %Na (A)
StrctFrm(:,14)=apfu(:,10); %K (A)
StrctFrm(:,15)=StrctFrm(:,14)+StrctFrm(:,13)+StrctFrm(:,12); %A site sum

%XMg calculation
for c=1:m
    if StrctFrm(c,10) + StrctFrm(c,8) > 0 %checks if Mg content is > 0
        StrctFrm(c,16)=StrctFrm(c,10)./(StrctFrm(c,10)+StrctFrm(c,8)); %XMg = Mg/(Mg+Fe)
    else
        %if Fe + Mg = 0, then the XMg calculation results in division by
        %zero, prevents NaN output
        StrctFrm(c,16) = 0; 
    end
end

%limit on significant digits (eliminates rounding noise)
StrctFrm(StrctFrm<1e-6) = 0;
apfu(apfu<1e-6) = 0;

StrctFrm=array2table(StrctFrm,'VariableNames',{'Si_T1','Fe3_T1','Al_T1','Sum_T1','Al_T2','Ti_T2','Sum_T2','Fe2_B','Mn2_B','Mg_B','Sum_B','Ca_A','Na_A','K_A','A_Sum','XMg'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Fe2','Mn2','Mg','Ca','Na','K','Cation_Sum'});

end