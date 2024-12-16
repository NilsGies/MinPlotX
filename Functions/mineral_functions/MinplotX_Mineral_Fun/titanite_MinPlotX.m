%Calculate titanite structural formula 
%last modified 15.07.2024

function [StrctFrm, apfu,options_definition]=titanite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
options_definition=struct();

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
SrO_mw=103.619;
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
F_mw=18.998;
Cl_mw=35.45;
H2O_mw=18.015;

%% moles of cations

MC=zeros(m,11);%create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; % for SiO2
MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; % for TiO2
MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; % for Al2O3

%calculates for Y2O3 if it is included in the analysis
if any(strcmp(headers,'Y2O3'))
    MC(:,4)=(data(:,strcmp(headers,'Y2O3'))./Y2O3_mw).*2; %for Y2O3
end

%calculates for FeO and/or Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, FeO is converted to Fe2O3 and
     %combined with Fe2O3

     MC(:,5)=((data(:,strcmp(headers,'FeO')).*(Fe2O3_mw./2*FeO_mw))+data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
else
    if any(strcmp(headers,'Fe2O3')) %if Fe2O3 is Fe2O3total

        MC(:,5)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;

    elseif any(strcmp(headers,'FeO')) %if FeO total is given, converted to Fe2O3 total

        %calculate moles of cations
        MC(:,5)=((data(:,strcmp(headers,'FeO')).*(Fe2O3_mw./(2*FeO_mw)))./Fe2O3_mw).*2;
    else
        %Fe2O3 and FeO are not included, Fe = 0
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

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,6)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

%calculates for MgO if it is included in the analysis
if any(strcmp(headers,'MgO'))
    MC(:,7)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
end

MC(:,8)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,9)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%calculates for K2O if it is included in the analysis 
if any(strcmp(headers,'K2O'))
    MC(:,10)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
end

%calculates for F if it is included in the analysis 
if any(strcmp(headers,'F'))
    MC(:,11)=data(:,strcmp(headers,'F'))./F_mw; %for F
end
 

%% Normalized moles of cations 

% Normalization is calculated assuming that the tetrahedral and octahedral
% sites are fully occupied by Si, Ti, Al, Fe3+, Mn, and Mg
N_fact=2./(MC(:,1)+MC(:,2)+MC(:,3)+MC(:,5)+MC(:,6)+MC(:,7)); %normalization factor

apfu=MC.*N_fact;
apfu(:,12)=sum(apfu(:,1:10),2);

%% structural formula

%T site
StrctFrm(:,1)=apfu(:,1); %Si (T)

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

%T site sum
StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %should sum to 1

%Octahedral Site
StrctFrm(:,4)=apfu(:,3)-StrctFrm(:,2); %Al 
StrctFrm(:,5)=apfu(:,2); %Ti4+
StrctFrm(:,6)=apfu(:,5); %Fe3+ 
StrctFrm(:,7)=apfu(:,6); %Mn 
StrctFrm(:,8)=apfu(:,7); %Mg 

%Decahedral Site
StrctFrm(:,9)=apfu(:,8); %Ca 
StrctFrm(:,10)=apfu(:,4); %Y3+ 
StrctFrm(:,11)=apfu(:,9); %Na
StrctFrm(:,12)=apfu(:,10); %K

%cation sum
StrctFrm(:,13)=sum(StrctFrm(:,1:2),2)+sum(StrctFrm(:,4:12),2);

%anions
StrctFrm(:,14)=apfu(:,11); %F 
StrctFrm(:,15)=(StrctFrm(:,4)+StrctFrm(:,6))-StrctFrm(:,14); %OH = (Al + Fe) - F
Charge_sum=4.*apfu(:,1)+4.*apfu(:,2)+3.*apfu(:,3)+3.*apfu(:,4)+3.*apfu(:,5)+2.*apfu(:,6)+2.*apfu(:,7)+2.*apfu(:,8)+apfu(:,9)+apfu(:,10)-StrctFrm(:,14)-StrctFrm(:,15);
StrctFrm(:,16)=5-0.5*(10-Charge_sum); %O2
%O2 anions are calculated as the sum of the charges of cations minus 1/2
%for F + OH

%anion sum
StrctFrm(:,17)=StrctFrm(:,14)+StrctFrm(:,15)+StrctFrm(:,16);

%XTtn
StrctFrm(:,18)=StrctFrm(:,5)./(StrctFrm(:,5)+StrctFrm(:,4)+StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8));

%limit on significant digits (eliminates rounding noise)
StrctFrm(StrctFrm<1e-6) = 0;
apfu(apfu<1e-6) = 0;

StrctFrm=array2table(StrctFrm,'VariableNames',{'Si_T','Al_T','T_Sum','Al_Oct','Ti_Oct','Fe3_Oct','Mn2_Oct','Mg_Oct','Ca_Dec','Y_Dec','Na_Dec','K_Dec','Cation_Sum','F','OH','O','Anion_Sum','XTtn'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Y','Fe3','Mn2','Mg','Ca','Na','K','F','Cation_Sum'});

end



