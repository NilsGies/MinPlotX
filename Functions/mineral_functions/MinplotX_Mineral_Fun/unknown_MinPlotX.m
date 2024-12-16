%% unknown mineral recalculation
%last modified 25.07.2024

function [StrctFrm, apfu, options_definition]=unknown_MinPlotX(data,headers,options)

%% define empty output variables
StrctFrm=[];
apfu=[];

%% options definition
options_definition=struct();

%Parameter 1
options_definition.UseKnownFe3.question='Do you want to use a known Fe3 ratio';
options_definition.UseKnownFe3.Value=false;

%Parameter 2
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 3
options_definition.UseKnownMn3.question='Do you want to use a known Mn3 ratio';
options_definition.UseKnownMn3.Value=false;

%Parameter 4
options_definition.Mn3_ratio.question='Enter Mn3_ratio';
options_definition.Mn3_ratio.Value=0;
options_definition.Mn3_ratio.limits=[0 1];

%Parameter 5
options_definition.cation_normalization.question='Do you wish to calculate on a cation basis?';
options_definition.cation_normalization.description='here text'; % optional
options_definition.cation_normalization.options={true,false}; % optional
options_definition.cation_normalization.Value=false; %default_value

%Parameter 6
options_definition.moles.question='How many cations/anions do you wish to normalize to?';
options_definition.moles.description='here text'; % optional
options_definition.moles.Value=12; %default_value
options_definition.moles.limits=[1 inf]; %optional

%Parameter 7
options_definition.RecalcTotalFe.question='Recalculate FeO and Fe2O3?';
options_definition.RecalcTotalFe.Value=false;

%Parameter 8
options_definition.RecalcTotalMn.question='Recalculate MnO and Mn2O3?';
options_definition.RecalcTotalMn.Value=false;


Fe2O3_mw=159.6874;
FeO_mw=71.8442;
MnO_mw=70.937;
Mn2O3_mw=157.873;

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
%load molecular weight table
if isfield(options,'T_MW')
T_MW=options.T_MW.Value;
else
T_MW=readtable(('Elements_to_oxides.txt'));
end
MW=table2array(T_MW(:,4)).'; %data converted to array
Moles_element=table2array(T_MW(:,7)).'; %data converted to array
Moles_oxygen=table2array(T_MW(:,8)).'; %data converted to array

%% Checks for Fe3+/Fetotal ratio

%assigns Fe3/Fetotal if it is included in the analysis
if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=zeros(m,1); %Fe3+/FeTotal ratio is 0
    options.RecalcTotalFe.Value=false;
end

%% Checks for Mn3+/Mntotal ratio

%assigns Mn3/Mntotal if it is included in the analysis
if isfield(options,'Mn3_ratio') && isfield(options,'UseKnownMn3') && options.UseKnownMn3.Value==true
    Mn_rat=options.Mn3_ratio.Value;
    options.RecalcTotalMn.Value=true;
else
    Mn_rat(:,1)=zeros(m,1); %Mn3+/MnTotal ratio is 0
    options.RecalcTotalMn.Value=false;
end

%% organization and assignment of weight percents
W=zeros(m,92);

if any(strcmp(headers,'H2O'))
    W(:,90)=data(:,strcmp(headers,'H2O'));
end

if any(strcmp(headers,'Li2O'))
    W(:,82)=data(:,strcmp(headers,'Li2O'));
end

if any(strcmp(headers,'BeO'))
    W(:,66)=data(:,strcmp(headers,'BeO'));
end

if any(strcmp(headers,'B2O3'))
    W(:,36)=data(:,strcmp(headers,'B2O3'));
end

if any(strcmp(headers,'CO2'))
    W(:,15)=data(:,strcmp(headers,'CO2'));
end

if any(strcmp(headers,'F'))
    W(:,91)=data(:,strcmp(headers,'F'));
end

if any(strcmp(headers,'Na2O'))
    W(:,83)=data(:,strcmp(headers,'Na2O'));
end

if any(strcmp(headers,'MgO'))
    W(:,67)=data(:,strcmp(headers,'MgO'));
end

if any(strcmp(headers,'Al2O3'))
    W(:,37)=data(:,strcmp(headers,'Al2O3'));
end

if any(strcmp(headers,'SiO2'))
    W(:,16)=data(:,strcmp(headers,'SiO2'));
end

if any(strcmp(headers,'P2O5'))
    W(:,10)=data(:,strcmp(headers,'P2O5'));
end

if any(strcmp(headers,'SO2'))
    W(:,17)=data(:,strcmp(headers,'SO2'));
end

if any(strcmp(headers,'SO3'))
    W(:,3)=data(:,strcmp(headers,'SO3'));
end

if any(strcmp(headers,'Cl'))
    W(:,92)=data(:,strcmp(headers,'Cl'));
end

if any(strcmp(headers,'K2O'))
    W(:,84)=data(:,strcmp(headers,'K2O'));
end

if any(strcmp(headers,'CaO'))
    W(:,68)=data(:,strcmp(headers,'CaO'));
end

if any(strcmp(headers,'Sc2O3'))
    W(:,38)=data(:,strcmp(headers,'Sc2O3'));
end

if any(strcmp(headers,'TiO2'))
    W(:,18)=data(:,strcmp(headers,'TiO2'));
end

if any(strcmp(headers,'V2O3'))
    W(:,39)=data(:,strcmp(headers,'V2O3'));
end

if any(strcmp(headers,'V2O5'))
    W(:,11)=data(:,strcmp(headers,'V2O5'));
end

if any(strcmp(headers,'Cr2O3'))
    W(:,40)=data(:,strcmp(headers,'Cr2O3'));
end

%calculates for MnO and Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) && isfield(options,'RecalcTotalMn') && options.RecalcTotalMn.Value==true

    %recalculate total MnO based on MnO+Mn2O3:
    MnOT=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')));

    %first calculate the weight percent of Mn2O3 and MnO from MnO total
    Mn3_wtper(:,1)=(MnOT.*Mn_rat)./((2*MnO_mw)./Mn2O3_mw); %calculates Mn2O3 from Mn3+ ratio and total MnO
    Mn2_wtper(:,1)=MnOT-(MnOT.*Mn_rat); %MnO = Total MnO - Mn3+ in Total MnO

    %calculate moles of cations
    W(:,41)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
    W(:,69)=Mn2_wtper(:,1)./MnO_mw;

elseif any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

   W(:,69)=data(:,strcmp(headers,'MnO'));
    W(:,41)=data(:,strcmp(headers,'Mn2O3'));
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        %first calculate the weight percent of Mn2O3 and MnO from MnO total
        W(:,41)=(data(:,strcmp(headers,'MnO')).*Mn_rat)./((2*MnO_mw)./Mn2O3_mw); %calculates Mn2O3 from Mn3+ ratio and total MnO
        W(:,69)=data(:,strcmp(headers,'MnO'))-(data(:,strcmp(headers,'MnO')).*Mn_rat); %MnO = Total MnO - Mn3+ in Total MnO

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 is Mn2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        W(:,41)=data(:,strcmp(headers,'Mn2O3'))-(data(:,strcmp(headers,'Mn2O3')).*(1-Mn_rat)); %Mn2O3 = Total Mn2O3 - Mn2+ in Total Mn2O3 
        W(:,69)=(data(:,strcmp(headers,'Mn2O3')).*(1-Mn_rat)).*((2*MnO_mw)./Mn2O3_mw); %calculates MnO from Mn3+ ratio and total Mn2O3
    else
        %Mn2O3 and MnO are not included, Mn2 and Mn3 = 0
    end
end

if any(strcmp(headers,'MnO2'))
    W(:,19)=data(:,strcmp(headers,'MnO2'));
end

%calculates for FeO and Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) && isfield(options,'RecalcTotalFe') && options.RecalcTotalFe.Value==true

    %recalculate total FeO based on FeO+Fe2O3:
    FeOT=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));

    %first calculate the weight percent of Fe2O3 and FeO from FeO total
    Fe3_wtper(:,1)=(FeOT.*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
    Fe2_wtper(:,1)=FeOT-(FeOT.*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    %calculate moles of cations
    W(:,42)=(Fe3_wtper(:,1));
    W(:,70)=Fe2_wtper(:,1);

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs
    W(:,70)=data(:,strcmp(headers,'FeO'));
    W(:,42)=data(:,strcmp(headers,'Fe2O3'));

else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        W(:,42)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        W(:,70)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        W(:,42)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        W(:,70)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

    end
end

if any(strcmp(headers,'CoO'))
    W(:,71)=data(:,strcmp(headers,'CoO'));
end

if any(strcmp(headers,'NiO'))
    W(:,72)=data(:,strcmp(headers,'NiO'));
end

if any(strcmp(headers,'CuO'))
    W(:,73)=data(:,strcmp(headers,'CuO'));
end

if any(strcmp(headers,'ZnO'))
    W(:,74)=data(:,strcmp(headers,'ZnO'));
end

if any(strcmp(headers,'Ga2O3'))
    W(:,43)=data(:,strcmp(headers,'Ga2O3'));
end

if any(strcmp(headers,'GeO2'))
    W(:,20)=data(:,strcmp(headers,'GeO2'));
end

if any(strcmp(headers,'As2O3'))
    W(:,44)=data(:,strcmp(headers,'As2O3'));
end

if any(strcmp(headers,'SeO2'))
    W(:,21)=data(:,strcmp(headers,'SeO2'));
end

if any(strcmp(headers,'SeO3'))
    W(:,4)=data(:,strcmp(headers,'SeO3'));
end

if any(strcmp(headers,'Rb2O'))
    W(:,85)=data(:,strcmp(headers,'Rb2O'));
end

if any(strcmp(headers,'SrO'))
    W(:,75)=data(:,strcmp(headers,'SrO'));
end

if any(strcmp(headers,'Y2O3'))
    W(:,45)=data(:,strcmp(headers,'Y2O3'));
end

if any(strcmp(headers,'ZrO2'))
    W(:,22)=data(:,strcmp(headers,'ZrO2'));
end

if any(strcmp(headers,'Nb2O5'))
    W(:,12)=data(:,strcmp(headers,'Nb2O5'));
end

if any(strcmp(headers,'MoO2'))
    W(:,23)=data(:,strcmp(headers,'MoO2'));
end

if any(strcmp(headers,'MoO3'))
    W(:,5)=data(:,strcmp(headers,'MoO3'));
end

if any(strcmp(headers,'RuO2'))
    W(:,24)=data(:,strcmp(headers,'RuO2'));
end

if any(strcmp(headers,'Rh2O3'))
    W(:,46)=data(:,strcmp(headers,'Rh2O3'));
end

if any(strcmp(headers,'PdO'))
    W(:,76)=data(:,strcmp(headers,'PdO'));
end

if any(strcmp(headers,'Ag2O'))
    W(:,86)=data(:,strcmp(headers,'Ag2O'));
end

if any(strcmp(headers,'CdO'))
    W(:,77)=data(:,strcmp(headers,'CdO'));
end

if any(strcmp(headers,'In2O3'))
    W(:,47)=data(:,strcmp(headers,'In2O3'));
end

if any(strcmp(headers,'SnO2'))
    W(:,25)=data(:,strcmp(headers,'SnO2'));
end

if any(strcmp(headers,'Sb2O3'))
    W(:,48)=data(:,strcmp(headers,'Sb2O3'));
end

if any(strcmp(headers,'TeO2'))
    W(:,26)=data(:,strcmp(headers,'TeO2'));
end

if any(strcmp(headers,'TeO3'))
    W(:,6)=data(:,strcmp(headers,'TeO3'));
end

if any(strcmp(headers,'Cs2O'))
    W(:,87)=data(:,strcmp(headers,'Cs2O'));
end

if any(strcmp(headers,'BaO'))
    W(:,78)=data(:,strcmp(headers,'BaO'));
end

if any(strcmp(headers,'La2O3'))
    W(:,49)=data(:,strcmp(headers,'La2O3'));
end

if any(strcmp(headers,'Ce2O3'))
    W(:,50)=data(:,strcmp(headers,'Ce2O3'));
end

if any(strcmp(headers,'CeO2'))
    W(:,27)=data(:,strcmp(headers,'CeO2'));
end

if any(strcmp(headers,'Pr2O3'))
    W(:,51)=data(:,strcmp(headers,'Pr2O3'));
end

if any(strcmp(headers,'Nd2O3'))
    W(:,52)=data(:,strcmp(headers,'Nd2O3'));
end

if any(strcmp(headers,'Sm2O3'))
    W(:,53)=data(:,strcmp(headers,'Sm2O3'));
end

if any(strcmp(headers,'EuO'))
    W(:,79)=data(:,strcmp(headers,'EuO'));
end

if any(strcmp(headers,'Eu2O3'))
    W(:,54)=data(:,strcmp(headers,'Eu2O3'));
end

if any(strcmp(headers,'Gd2O3'))
    W(:,55)=data(:,strcmp(headers,'Gd2O3'));
end

if any(strcmp(headers,'Tb2O3'))
    W(:,56)=data(:,strcmp(headers,'Tb2O3'));
end

if any(strcmp(headers,'Dy2O3'))
    W(:,57)=data(:,strcmp(headers,'Dy2O3'));
end

if any(strcmp(headers,'Ho2O3'))
    W(:,58)=data(:,strcmp(headers,'Ho2O3'));
end

if any(strcmp(headers,'Er2O3'))
    W(:,59)=data(:,strcmp(headers,'Er2O3'));
end

if any(strcmp(headers,'Tm2O3'))
    W(:,60)=data(:,strcmp(headers,'Tm2O3'));
end

if any(strcmp(headers,'Yb2O3'))
    W(:,61)=data(:,strcmp(headers,'Yb2O3'));
end

if any(strcmp(headers,'Lu2O3'))
    W(:,62)=data(:,strcmp(headers,'Lu2O3'));
end

if any(strcmp(headers,'HfO2'))
    W(:,28)=data(:,strcmp(headers,'HfO2'));
end

if any(strcmp(headers,'Ta2O5'))
    W(:,13)=data(:,strcmp(headers,'HfO2'));
end

if any(strcmp(headers,'WO3'))
    W(:,7)=data(:,strcmp(headers,'WO3'));
end

if any(strcmp(headers,'ReO2'))
    W(:,29)=data(:,strcmp(headers,'ReO2'));
end

if any(strcmp(headers,'ReO3'))
    W(:,8)=data(:,strcmp(headers,'ReO3'));
end

if any(strcmp(headers,'Re2O7'))
    W(:,2)=data(:,strcmp(headers,'Re2O7'));
end

if any(strcmp(headers,'OsO2'))
    W(:,30)=data(:,strcmp(headers,'OsO2'));
end

if any(strcmp(headers,'OsO4'))
    W(:,1)=data(:,strcmp(headers,'OsO4'));
end

if any(strcmp(headers,'IrO2'))
    W(:,31)=data(:,strcmp(headers,'IrO2'));
end

if any(strcmp(headers,'PtO2'))
    W(:,32)=data(:,strcmp(headers,'PtO2'));
end

if any(strcmp(headers,'Au2O'))
    W(:,88)=data(:,strcmp(headers,'Au2O'));
end

if any(strcmp(headers,'Au2O3'))
    W(:,63)=data(:,strcmp(headers,'Au2O3'));
end

if any(strcmp(headers,'HgO'))
    W(:,80)=data(:,strcmp(headers,'HgO'));
end

if any(strcmp(headers,'Tl2O'))
    W(:,89)=data(:,strcmp(headers,'Tl2O'));
end

if any(strcmp(headers,'Tl2O3'))
    W(:,64)=data(:,strcmp(headers,'Tl2O3'));
end

if any(strcmp(headers,'PbO'))
    W(:,81)=data(:,strcmp(headers,'PbO'));
end

if any(strcmp(headers,'Bi2O3'))
    W(:,65)=data(:,strcmp(headers,'Bi2O3'));
end

if any(strcmp(headers,'ThO2'))
    W(:,33)=data(:,strcmp(headers,'ThO2'));
end

if any(strcmp(headers,'PaO2'))
    W(:,34)=data(:,strcmp(headers,'PaO2'));
end

if any(strcmp(headers,'Pa2O5'))
    W(:,14)=data(:,strcmp(headers,'Pa2O5'));
end

if any(strcmp(headers,'UO2'))
    W(:,35)=data(:,strcmp(headers,'UO2'));
end

if any(strcmp(headers,'UO3'))
    W(:,9)=data(:,strcmp(headers,'UO3'));
end

%% moles of cations 
MC=(W./MW).*Moles_element; 

if options.cation_normalization.Value == true
    %Cation normalization

    MC_Sum=sum(MC(:,1:92),2)-(MC(:,6)+MC(:,14)); %sum minus F + Cl

    MC_N=(options.moles.Value)./MC_Sum; %normalization factor

    %normalized apfu
    apfu2=MC.*MC_N;

    apfu2(:,93)=sum(apfu2(:,1:92),2); %cation sum

else

    %Oxygen normalization
    O2=MC.*(Moles_oxygen./Moles_element);

    O2_Sum=sum(O2(:,1:92),2)-0.5.*(O2(:,6)+O2(:,14)); %sum of O2, including F and Cl

    O2_N=(options.moles.Value)./O2_Sum; %normalization factor

    %normalized apfu
    apfu2=(O2.*O2_N).*(Moles_element./Moles_oxygen);

    apfu2(:,93)=sum(apfu2(:,1:92),2); %cation sum
end

%limit on significant digits (eliminates rounding noise)
apfu2(apfu2<1e-6) = 0;



apfu2_vars={'Os8','Re7','S6','Se6','Mo6','Te6','W' ...
    ,'Re6','U6','P','V5','Nb','Ta','Pa5','C','Si','S4','Ti','Mn4','Ge','Se4' ...
    ,'Zr','Mo4','Ru','Sn','Te4','Ce4','Hf','Re4','Os4','Ir','Pt','Th','Pa4' ...
    ,'U4','B','Al','Sc','V3','Cr','Mn3','Fe3','Ga','As3','Y','Rh','In','Sb','La' ...
    ,'Ce3','Pr','Nd','Sm','Eu3','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Au3' ...
    ,'Tl3','Bi','Be','Mg','Ca','Mn2','Fe2','Co','Ni','Cu','Zn','Sr','Pd','Cd' ...
    ,'Ba','Eu2','Hg','Pb','Li','Na','K','Rb','Ag','Cs','Au1','Tl1','H','F','Cl','Cation_Sum'};

apfu2(isnan(apfu2))=0;

empty_cols_test=not(sum(apfu2==0,1)==size(apfu2,1));

apfu=array2table(apfu2(:,empty_cols_test),'VariableNames',apfu2_vars(empty_cols_test));

StrctFrm=apfu;



end













