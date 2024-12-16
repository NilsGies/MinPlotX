%%chloritoid Structural Formula
% last modified 01.08.2024

function [StrctFrm, apfu]=chloritoid_Fe3unknown(data,headers,options)

%% start calculation

[m,~]=size(data); %finds the x and y size of the input data matrix

cat=8.0; %cations per formula unit
Opfu=12.0; %oxygens per formula unit

%% Molecular weights

SiO2_mw=60.083;
TiO2_mw=79.865;
Al2O3_mw=101.961;
Cr2O3_mw=151.989;
Fe2O3_mw=159.6874;
Y2O3_mw=225.809;
NiO_mw=74.692;
ZnO_mw=81.381;
FeO_mw=71.8442;
MnO_mw=70.937;
Mn2O3_mw=157.873;
MgO_mw=40.304;
CaO_mw=56.0774;
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
F_mw=18.998;
Cl_mw=35.45;

%% Calculate moles of cations

MC=zeros(m,8); %create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

%calculates for FeO and/or Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
     %combined with FeO
     MC(:,4)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,4)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    else %if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        MC(:,4)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw))./FeO_mw;

    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,5)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,5)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,5)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

MC(:,6)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%calculates for CaO if it is included in the analysis
if any(strcmp(headers,'CaO'))
    MC(:,7)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,8)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4); %for FeO
O2(:,5)=MCnorm(:,5); %for MnO
O2(:,6)=MCnorm(:,6); %for MgO
O2(:,7)=MCnorm(:,7); %for CaO
O2(:,8)=MCnorm(:,8)./2; %for Na2O

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

apfu(:,1)=MCnorm(:,1); %for Si
apfu(:,2)=MCnorm(:,2); %for Ti
apfu(:,3)=MCnorm(:,3); %for Al
apfu(:,6)=MCnorm(:,5); %for Mn
apfu(:,7)=MCnorm(:,6); %for Mg
apfu(:,8)=MCnorm(:,7); %for Ca
apfu(:,9)=MCnorm(:,8); %for Na


%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 12
%if so, then there is no Fe3+
%if totalO2 < 12, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(12-totalO2) then the amount
%of Fe3+ = 2*(12-totalO2), if false then, all Fe is Fe3+


for c=1:m
    if (Opfu-O2total(c,1)) > 0
        if MCnorm(c,4) > 2.*(Opfu-O2total(c,1))
            apfu(c,4)=2.*(Opfu-O2total(c,1));
        else
            apfu(c,4)=MCnorm(c,4);
        end
    else
        apfu(c,4)=0;
    end
end

apfu(:,5)=MCnorm(:,4)-apfu(:,4); %the apfu of Fe2+ equals totalFe-Fe3+

apfu(:,10)=sum(apfu,2); %calculations the total, which should be 8


% Oxygen deficiency
apfu(:,11)=Opfu-O2total; %must be greater than zero

%% Structural Formula

StrctFrm=zeros(length(m),12);

%T SITE
%Si
for c=1:m
    if apfu(c,1)<2.000
        StrctFrm(c,1)=apfu(c,1); %If Si < 2, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=2; %If Si is in excess, then Si(T) = 2
    end
end


%Al2O3 Layer (L2)
for c=1:m
    if apfu(c,3)<=3.000
        StrctFrm(c,2)=apfu(c,3); %If Al =< 3, then Al(3) = the measured Al content
    else
        StrctFrm(c,2)=3; %If Al>3 is in excess, then Al(L1) = 3
    end
end

%(Al, Ti, Fe3+) + (Fe, Mg, Mn) layer (L1)
StrctFrm(:,3)=apfu(:,3)-StrctFrm(:,2); %Al
StrctFrm(:,4)=apfu(:,2); %Ti
StrctFrm(:,5)=apfu(:,4); %Fe3+
StrctFrm(:,6)=apfu(:,5); %Fe2+
StrctFrm(:,7)=apfu(:,6); %Mn
StrctFrm(:,8)=apfu(:,7); %Mg
StrctFrm(:,9)=apfu(:,8); %Ca
StrctFrm(:,10)=apfu(:,9); %Na
StrctFrm(:,11)=sum(StrctFrm,2); %sum


XMg=zeros(m,1);
for c=1:m
    if StrctFrm(c,8) + StrctFrm(c,6) > 0 %checks if Mg content is > 0
        XMg(c,1)=StrctFrm(c,8)./(StrctFrm(c,8)+StrctFrm(c,6)); %XMg = Mg/(Mg+Fe)
    else
        %if Fe + Mg = 0, then the XMg calculation results in division by
        %zero, prevents NaN output
        XMg(c,1) = 0; 
    end
end


%% output

%limit on significant digits (eliminates rounding noise)
StrctFrm(StrctFrm<1e-6) = 0;
StrctFrm(:,12)=apfu(:,11); %O2 deficiency
XMg(XMg<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)
apfu(apfu<1e-6) = 0;
apfu(:,11)=Opfu-O2total; %also adds O2 def to the apfu output

all=[StrctFrm XMg];

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_L2','Al_L1','Ti_L1','Fe3_L1','Fe2_L1','Mn2_L1','Mg_L1','Ca_L1','Na_L1','Cation_Sum','O2_deficiency','XMg'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Fe2','Mn2','Mg','Ca','Na','Cation_Sum','O2_deficiency'});

end


