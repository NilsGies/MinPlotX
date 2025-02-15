%%olivine Structural Formula 
% last modified 01.08.2024

function [StrctFrm, apfu]=olivine_Fe3known(data,headers,options)

[m,~]=size(data); %finds the x and y size of the input data matrix

cat=3.0; %cations per formula unit
Opfu=4.0; %oxygens per formula unit

%% Checks for Fe3+/Fetotal ratio
%assigns Fe3/Fetotal if it is included in the analysis

if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=zeros(m,1); %Fe3+/FeTotal ratio is 0
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

%% Calculate cations units

MC=zeros(m,10); %create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2

%adds a column of zeros if Ti is not included in the calculation
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

%calculates for Al2O3 if it is included in the analysis 
if any(strcmp(headers,'Al2O3'))
    MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3
end

%calculates for FeO and Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) && isfield(options,'RecalcTotalFe') && options.RecalcTotalFe.Value==true

    %recalculate total FeO based on FeO+Fe2O3:
    FeOT=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));

    %first calculate the weight percent of Fe2O3 and FeO from FeO total
    Fe3_wtper(:,1)=(FeOT.*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
    Fe2_wtper(:,1)=FeOT-(FeOT.*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    %calculate moles of cations
    MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
    MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,4)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,7)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

    end
end

%calculates for Cr2O3 if it is included in the analysis 
if any(strcmp(headers,'Cr2O3'))
    MC(:,5)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
end

%calculates for NiO if it is included in the analysis 
if any(strcmp(headers,'NiO'))
    MC(:,6)=data(:,strcmp(headers,'NiO'))./NiO_mw; %for NiO
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,8)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,8)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,8)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

MC(:,9)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%calculates for CaO if it is included in the analysis 
if any(strcmp(headers,'CaO'))
    MC(:,10)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end


%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Fe2O3
O2(:,5)=MC(:,5).*(3/2); %for Cr2O3
O2(:,6)=MC(:,6); %for NiO
O2(:,7)=MC(:,7); %for FeO
O2(:,8)=MC(:,8); %for MnO
O2(:,9)=MC(:,9); %for MgO
O2(:,10)=MC(:,10); %for CaO

O2total=sum(O2,2); %O2 totals
MCnormfact=Opfu./O2total; %normalization factor

%% Atoms pfu



apfu=MCnormfact.*MC; %creates a matrix of normalized cations
apfu(:,11)=sum(apfu,2); %calculations the total, which should be close to 3
apfu(:,12)=zeros(m,1); %create columns of zeros for O2_deficiency

%% structural formula calculation

StrctFrm=zeros(m,14); %create columns of zeros if optional data are not included

%T SITE
%Si 
for c=1:m
    if apfu(c,1)<1.000 
        StrctFrm(c,1)=apfu(c,1); %If Si < 1, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=1; %If Si is in excess, then Si(T) = 1
    end
end

%Al(T)
for c=1:m
    if 1-StrctFrm(c,1)>0 %Is 1-Si > 0? If y, then some Al goes into T
        if 1-StrctFrm(c,1)>apfu(c,3) %For low Al , 1-Si may be > Al
            StrctFrm(c,2)=apfu(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=1-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to M
        end
    else
        StrctFrm(c,2)=0; %if Si=1, then no Al goes into T
    end
end

%Fe3+(T)
for c=1:m
    if 1-StrctFrm(c,1)-StrctFrm(c,2)>0 %Is 1-(Si+Al) > 0? If y, then some Fe3+ goes into T
        if 1-StrctFrm(c,1)-StrctFrm(c,2)>apfu(c,4) %For low Fe3+ Ol, 1-(Si+Al) may be > Fe3+
            StrctFrm(c,3)=apfu(c,4); %All Fe3+ goes into T
        else
            StrctFrm(c,3)=1-StrctFrm(c,1)-StrctFrm(c,2); %if there isn't enough space in T for all Fe3+, the rest will go to M
        end
    else
        StrctFrm(c,3)=0; %if Si+Al=1, then no Fe3+ goes into T
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3);

%M SITE
%Al
StrctFrm(:,5)=apfu(:,3)-StrctFrm(:,2); %Al(M1) = Total Al - Al(T)

%Ti
StrctFrm(:,6)=apfu(:,2);

%Fe3+
StrctFrm(:,7)=apfu(:,4)-StrctFrm(:,3); %Fe3+(M1) = Total Fe3+ - Fe3+(T)

%Cr
StrctFrm(:,8)=apfu(:,5);

%Ni
StrctFrm(:,9)=apfu(:,6);

%Fe
StrctFrm(:,10)=apfu(:,7);

%Mn
StrctFrm(:,11)=apfu(:,8);

%Mg
StrctFrm(:,12)=apfu(:,9);

%Ca
StrctFrm(:,13)=apfu(:,10);

%M total
StrctFrm(:,14)=sum(StrctFrm(:,5:13),2);

%% Outputs

Endmembers(:,1)=apfu(:,9)./(apfu(:,4)+apfu(:,7)+apfu(:,8)+apfu(:,9)+apfu(:,10)); % XFo
Endmembers(:,2)=(apfu(:,4)+apfu(:,7))./(apfu(:,4)+apfu(:,7)+apfu(:,8)+apfu(:,9)+apfu(:,10)); % XFa
Endmembers(:,3)=(apfu(:,8))./(apfu(:,4)+apfu(:,7)+apfu(:,8)+apfu(:,9)+apfu(:,10)); % XTe
Endmembers(:,4)=(apfu(:,10))./(apfu(:,4)+apfu(:,7)+apfu(:,8)+apfu(:,9)+apfu(:,10)); % XCa-Ol

Endmembers(Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)

all=[StrctFrm Endmembers apfu(:,12)];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu(apfu<1e-6) = 0;
all(isnan(all))=0; %remove NaN

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M','Ti_M','Fe3_M','Cr_M','Ni_M','Fe2_M','Mn2_M','Mg_M','Ca_M','Sum_M','XFo','XFa','XTep','XLrn','O2_deficiency'});

%outputs apfu table only
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Ni','Fe2','Mn2','Mg','Ca','Cation_Sum','O2_deficiency'});

end
