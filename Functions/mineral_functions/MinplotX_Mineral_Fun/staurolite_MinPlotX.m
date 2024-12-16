%%staurolite Structural Formula
%last modified 31.07.2024

function [StrctFrm, apfu,options_definition]=staurolite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
%Parameter 1
options_definition.Cation_Normalization.question='Normalization by 25.53 Al + Si??';
options_definition.Cation_Normalization.Value=false; %default_value

%Parameter 2
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=0.035;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 3
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

%% Checks for Fe3+/Fetotal ratio

%assigns Fe3/Fetotal if it is included in the analysis
if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=zeros (m,1); %Fe3+/FeTotal ratio is 
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
H2O_mw=18.015;

%% Moles of Oxides

MC=zeros(m,15); %create columns of zeros if optional data are not included

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
    MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,4)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,6)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

    end
end

%adds a column of zeros if Cr is not included in the calculation
if any(strcmp(headers,'Cr2O3'))
    MC(:,5)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,7)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;

    else
        %Mn2O3 and MnO are not included, Mn2 = 0

    end
end

MC(:,8)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%calculates for ZnO if it is included in the analysis
if any(strcmp(headers,'ZnO'))
    MC(:,9)=data(:,strcmp(headers,'ZnO'))./ZnO_mw; %for ZnO
end

%calculates for CaO if it is included in the analysis 
if any(strcmp(headers,'CaO'))
    MC(:,10)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,11)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%calculates for K2O if it is included in the analysis
if any(strcmp(headers,'K2O'))
    MC(:,12)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
end

%calculates for F if it is included in the analysis 
if any(strcmp(headers,'F'))
    MC(:,13)=data(:,strcmp(headers,'F'))./F_mw; %for F
end

%calculates for Cl if it is included in the analysis 
if any(strcmp(headers,'Cl'))
    MC(:,14)=data(:,strcmp(headers,'Cl'))./Cl_mw; %for Cl
end

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,15)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
end

%% Moles of O2 units

O2(:,1)=MC(:,1).*(2); %SiO2
O2(:,2)=MC(:,2).*(2); %TiO2
O2(:,3)=MC(:,3).*(3/2); %Al2O3
O2(:,4)=MC(:,4).*(3/2); %Fe2O3
O2(:,5)=MC(:,5).*(3/2); %Cr2O3
O2(:,6)=MC(:,6); %FeO
O2(:,7)=MC(:,7); %MnO
O2(:,8)=MC(:,8); %MgO
O2(:,9)=MC(:,9); %ZnO
O2(:,10)=MC(:,10); %CaO
O2(:,11)=MC(:,11).*(1/2); %Na2O
O2(:,12)=MC(:,12).*(1/2); %K2O
O2(:,13)=MC(:,13); %F
O2(:,14)=MC(:,14); %Cl
O2(:,15)=MC(:,15).*(1/2); %H2O

O2_Sum=sum(O2(:,1:15),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

O2_N=(48)./O2_Sum; %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% Iterate OH

if not(any(strcmp(headers,'H2O')))

    %initial oxygens of OH
    Hin=4-(N_Ox(:,14)+N_Ox(:,13)); %H = 4 - (F+Cl)
    O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O

    for z=1:10
        O2(:,1)=MC(:,1).*2; %SiO2
        O2(:,2)=MC(:,2).*2; %TiO2
        O2(:,3)=MC(:,3).*(3/2); %Al2O3
        O2(:,4)=MC(:,4).*(3/2); %Fe2O3
        O2(:,5)=MC(:,5).*(3/2); %Cr2O3
        O2(:,6)=MC(:,6); %FeO
        O2(:,7)=MC(:,7); %MnO
        O2(:,8)=MC(:,8); %MgO
        O2(:,9)=MC(:,9); %ZnO
        O2(:,10)=MC(:,10); %CaO
        O2(:,11)=MC(:,11).*(1/2); %Na2O
        O2(:,12)=MC(:,12).*(1/2); %K2O
        O2(:,13)=MC(:,13); %F
        O2(:,14)=MC(:,14); %Cl
        O2(:,15)=O_OH; %H2O

        O2_Sum=sum(O2(:,1:15),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

        O2_N=(48)./O2_Sum; %normalization factor

        %normalized moles of anions
        N_Ox=O2.*O2_N;

        Hin=4-(N_Ox(:,14)+N_Ox(:,13)); %H = 4 - (F+Cl)
        O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O
    end
end

%% Moles of cations (Oxygen Normalization)

apfu_O(:,1)=N_Ox(:,1)./2; %Si
apfu_O(:,2)=N_Ox(:,2)./2; %Ti
apfu_O(:,3)=N_Ox(:,3).*(2/3); %Al
apfu_O(:,4)=N_Ox(:,4).*(2/3); %Fe3
apfu_O(:,5)=N_Ox(:,5).*(2/3); %Cr3
apfu_O(:,6)=N_Ox(:,6); %FeO
apfu_O(:,7)=N_Ox(:,7); %MnO
apfu_O(:,8)=N_Ox(:,8); %MgO
apfu_O(:,9)=N_Ox(:,9); %ZnO
apfu_O(:,10)=N_Ox(:,10); %CaO
apfu_O(:,11)=N_Ox(:,11).*2; %Na2O
apfu_O(:,12)=N_Ox(:,12).*2; %K2O
apfu_O(:,13)=N_Ox(:,13); %F
apfu_O(:,14)=N_Ox(:,14); %Cl
%% 
apfu_O(:,15)=N_Ox(:,15).*2; %OH

%% normalization to 25.53 Si + Al

%STEP 1
%Normmalization to Si + Al 
Nfact=25.53./(apfu_O(:,1)+apfu_O(:,3)); %normalization

apfu_N=apfu_O.*Nfact; %normalization of the cations

%cation charges

Chrg=zeros(m,15); %create columns of zeros if optional data are not included

Chrg(:,1)=apfu_N(:,1).*4; %Si4+
Chrg(:,2)=apfu_N(:,2).*4; %Ti4+
Chrg(:,3)=apfu_N(:,3).*3; %Al3+
Chrg(:,4)=apfu_N(:,4).*3; %Fe3+
Chrg(:,5)=apfu_N(:,5).*3; %Cr3+
Chrg(:,6)=apfu_N(:,6).*2; %Fe2+
Chrg(:,7)=apfu_N(:,7).*2; %Mn2+
Chrg(:,8)=apfu_N(:,8).*2; %Mg2+
Chrg(:,9)=apfu_N(:,9).*2; %Zn2+
Chrg(:,10)=apfu_N(:,10).*2; %Ca2+
Chrg(:,11)=apfu_N(:,11); %Na1+
Chrg(:,12)=apfu_N(:,12); %K1+
Chrg(:,13)=apfu_N(:,13); %F
Chrg(:,14)=apfu_N(:,14); %Cl

%if H2O is known, the known value is taken
if any(strcmp(headers,'H2O'))
    Chrg(:,15)=apfu_N(:,15); %OH
    Crg_Ex=sum(Chrg(:,1:15),2)-96; %excess charge (>96)
else %otherwise the OH charge is 0
    Crg_Ex=sum(Chrg(:,1:15),2)-92; %excess charge (>92)
    apfu_N(:,15)=96-sum(Chrg(:,1:14),2); %reformulates OH, for O substitution following normalization
end

apfu_N(:,16)=sum(apfu_N(:,1:12),2); %calculations the cation total
apfu_N(:,17)=30-apfu_N(:,16); %total number of vacancies

%% organization

%cation normalized data
StrctFrm_N=zeros(m,19); 
StrctFrm_N=apfu_N;
%reorganization
StrctFrm_N(:,13)=apfu_N(:,17); %vacancies
StrctFrm_N(:,14)=apfu_N(:,16); %cation sum
StrctFrm_N(:,15)=apfu_N(:,13); %F
StrctFrm_N(:,16)=apfu_N(:,14); %Cl
StrctFrm_N(:,17)=apfu_N(:,15); %OH

%O on the X site
for c=1:m
    if Crg_Ex(c,1)<0
        StrctFrm_N(c,18)=0;
    else
        StrctFrm_N(c,18)=Crg_Ex(c,1); 
    end
end

StrctFrm_N(:,19)=sum(StrctFrm_N(:,15:18),2); %X site sum

%48 O normalized data
StrctFrm_O=zeros(m,19); 
StrctFrm_O=apfu_O;

apfu_O(:,16)=sum(apfu_O(:,1:12),2); %calculations the total
apfu_O(:,17)=30-apfu_O(:,16); %total number of vacancies

%reorganization
StrctFrm_O(:,13)=apfu_O(:,17); %vacancy
StrctFrm_O(:,14)=apfu_O(:,16); %cation sum
StrctFrm_O(:,15)=apfu_O(:,13); %F
StrctFrm_O(:,16)=apfu_O(:,14); %Cl
StrctFrm_O(:,17)=apfu_O(:,15); %OH
StrctFrm_O(:,19)=apfu_O(:,15)+apfu_O(:,14)+apfu_O(:,13); %Anion Sum

%limit on significant digits (eliminates rounding noise)
StrctFrm_O(StrctFrm_O<1e-4) = 0;
apfu_O(apfu_O<1e-4) = 0;

StrctFrm_N(StrctFrm_N<1e-4) = 0;
apfu_N(apfu_N<1e-4) = 0;

%Option to output oxygen normalized vs cation normalized outputs
if options.Cation_Normalization.Value==true
    apfu=array2table(apfu_N,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn2','Mg','Zn','Ca','Na','K','F','Cl','OH','Cation_Sum','vac'});
    StrctFrm=array2table(StrctFrm_N,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn2','Mg','Zn','Ca','Na','K','vac','Cation_Sum','F','Cl','OH','O','Anion_sum'});

else
    apfu=array2table(apfu_O,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn2','Mg','Zn','Ca','Na','K','F','Cl','OH','Cation_Sum','vac'});
    StrctFrm=array2table(StrctFrm_O,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn2','Mg','Zn','Ca','Na','K','vac','Cation_Sum','F','Cl','OH','O','Anion_sum'});

end

end

