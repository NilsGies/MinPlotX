%Clinopyroxene structural formula with Fe3+ estimation
% last modified 01.08.2024

function [StrctFrm, apfu]=pyroxene_fe3unknown(data,headers,options)

[m,~]=size(data); %finds the x and y size of the input data matrix

cat=4.0; %cations per formula unit
Opfu=6.0; %oxygens per formula unit

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

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

%adds a column of zeros if Cr is not included in the calculation
if any(strcmp(headers,'Cr2O3'))
    MC(:,4)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
end

%calculates for FeO and/or Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
     %combined with FeO
     MC(:,5)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,5)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    else %if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        MC(:,5)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw))./FeO_mw;

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

MC(:,7)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

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

%% Calculate normalized cations units

MCnormfact=cat./sum(MC,2); %normalization factor
MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5); %for FeO
O2(:,6)=MCnorm(:,6); %for MnO
O2(:,7)=MCnorm(:,7); %for MgO
O2(:,8)=MCnorm(:,8); %for CaO
O2(:,9)=MCnorm(:,9)./2; %for Na2O
O2(:,10)=MCnorm(:,10)./2; %for K2O

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

apfu=zeros(m,13); %create columns of zeros if optional data are not included

apfu(:,1)=MCnorm(:,1); %for Si
apfu(:,2)=MCnorm(:,2); %for Ti
apfu(:,3)=MCnorm(:,3); %for Al
apfu(:,5)=MCnorm(:,4); %for Cr
apfu(:,7)=MCnorm(:,6); %for Mn
apfu(:,8)=MCnorm(:,7); %for Mg
apfu(:,9)=MCnorm(:,8); %for Ca
apfu(:,10)=MCnorm(:,9); %for Na
apfu(:,11)=MCnorm(:,10); %for K

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 6
%if so, then there is no Fe3+
%if totalO2 < 6, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(6-totalO2) then the amount
%of Fe3+ = 2*(6-totalO2), if false then, all Fe is Fe3+
for c=1:m
    if (Opfu-O2total(c,1)) > 0
        if MCnorm(c,5) > 2.*(Opfu-O2total(c,1))
            apfu(c,4)=2.*(Opfu-O2total(c,1));
        else
            apfu(c,4)=MCnorm(c,5);
        end
    else
        apfu(c,4)=0;
    end
end

apfu(:,6)=MCnorm(:,5)-apfu(:,4); %the apfu of Fe2+ = totalFe-Fe3+

apfu(:,12)=sum(apfu,2); %calculations the total, which should be 4

% Oxygen deficiency
apfu(:,13)=Opfu-O2total; %must be greater than zero

%XMg
XMg=apfu(:,8)./(apfu(:,8)+apfu(:,6));

%% structural formula calculation

StrctFrm=zeros(m,13); %create columns of zeros if optional data are not included

%T SITE
%Si
for c=1:m
    if apfu(c,1)<2.000
        StrctFrm(c,1)=apfu(c,1); %If Si < 2, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=2; %If Si is in excess, then Si(T) = 2
    end
end

%Al(T)
for c=1:m
    if 2-StrctFrm(c,1)>0 %Is 2-Si > 0? If y, then some Al goes into T
        if 2-StrctFrm(c,1)>apfu(c,3) %For low Al cpx, 2-Si may be > Al
            StrctFrm(c,2)=apfu(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=2-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to M1
        end
    else
        StrctFrm(c,2)=0; %if Si=2, then no Al goes into T
    end
end

%Fe3+(T)
for c=1:m
    if 2-StrctFrm(c,1)-StrctFrm(c,2)>0 %Is 2-(Si+Al) > 0? If y, then some Fe3+ goes into T
        if 2-StrctFrm(c,1)-StrctFrm(c,2)>apfu(c,4) %For low Fe3+ cpx, 2-(Si+Al) may be > Fe3+
            StrctFrm(c,3)=apfu(c,4); %All Fe3+ goes into T
        else
            StrctFrm(c,3)=2-StrctFrm(c,1)-StrctFrm(c,2); %if there isn't enough space in T for all Fe3+, the rest will go to M1
        end
    else
        StrctFrm(c,3)=0; %if Si+Al=2, then no Fe3+ goes into T
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3); %Si + Al + Fe3+ in T

%M1 SITE

%Al(M1)
StrctFrm(:,5)=apfu(:,3)-StrctFrm(:,2); %Al(M1) = Total Al - Al(T)

%Ti (M1)
StrctFrm(:,6)=apfu(:,2);

%Fe3+ (M1)
StrctFrm(:,8)=apfu(:,4)-StrctFrm(:,3); %Fe3+(M1) = Total Fe3+ - Fe3+(T)

%Cr3+ (M1)
StrctFrm(:,7)=apfu(:,5);

%Mn (M1)
StrctFrm(:,9)=apfu(:,7);

%Mg (M1)
for c=1:m
    if XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5))<apfu(c,8) %if XMg*(1-Sum(Al to Mn) in M1) is < Mg
        StrctFrm(c,10)=XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5)); %Mg(M1)=XMg*(1-Sum(Al to Mn) in M1) and some Mg goes into M2
    else
        StrctFrm(c,10)=apfu(c,8); %if not, all Mg goes into M1
    end
end

%Fe2+ (M1)
for c=1:m
    if 1-(StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,5))>0 %Is 1-Sum(Al to Mg) in M1 >0?, if yes then:
        if 1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5)>apfu(c,6) %Is 1-Sum(Al to Mg) in M1 > Fe2+?, if yes then:
            StrctFrm(c,11)=apfu(c,6); %all Fe2+ goes into M1
        else
            StrctFrm(c,11)=1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5); %if there isn't enough space in M1 for all Fe2+, some goes into M2
        end
    else
        StrctFrm(c,11)=0; %If M1 is already filled, then no Fe2+ goes into M1
    end
end


%Sum of M1 site
StrctFrm(:,12)=StrctFrm(:,5)+StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)+StrctFrm(:,9)+StrctFrm(:,10)+StrctFrm(:,11);

%M2 SITE
%Mg (M2)
for c=1:m
    if apfu(c,8)-StrctFrm(c,10)>0 %Is Mgtotal-Mg(M1) > 0?
        StrctFrm(c,13)=apfu(c,8)-StrctFrm(c,10); %if yes, then some Mg goes into M2
    else
        StrctFrm(c,13)=0; %if no, then all Mg is in M1
    end
end

%Fe2+ (M2)
for c=1:m
    if apfu(c,6)-StrctFrm(c,11)>0 %Is Fe2+total-Fe2+(M1) > 0?
        StrctFrm(c,14)=apfu(c,6)-StrctFrm(c,11); %if yes, then some Fe2+ goes into M2
    else
        StrctFrm(c,14)=0; %if no, then all Fe2+ is in M1
    end
end

%Ca (M2)
StrctFrm(:,15)=apfu(:,9);

%Na (M2)
StrctFrm(:,16)=apfu(:,10);

%K (M2)
StrctFrm(:,17)=apfu(:,11);

%Sum of M2 site
StrctFrm(:,18)=StrctFrm(:,13)+StrctFrm(:,14)+StrctFrm(:,15)+StrctFrm(:,16)+StrctFrm(:,17);

%% end member calculations

%only do endmember calculation if structural formula with endembers is selected
Endmembers=zeros(m,7); %create columns of zeros if optional data are not included

%choose to calculate endmember fractions with Na-Fe3+-Cr
if options.Na_Endmembers.Value==true
    %Na-Ca endmembers
    A(:,1)=apfu(:,9); %Ca
    A(:,2)=apfu(:,3); %Al
    A(:,3)=apfu(:,4); %Fe3+
    A(:,4)=apfu(:,8); %Mg
    A(:,5)=apfu(:,6); %Fe2+
    A(:,6)=apfu(:,5); %Cr
    AT=transpose(A); %transpose of A

    M=[2 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 2 0 0 0; 0 2 0 0 0 0; 0 0 0 0 0 1];
    X=zeros(6,m);

    for c=1:m
        X(:,c)=inv(M)*AT(:,c); %calculates endmembers
    end

    Xtot=sum(X);%sum of endmembers
    Xnorm=X./sum(X); %normalizes the endmembers to 1
    XnormT=transpose(Xnorm); %transposes back

    Endmembers(:,1)=XnormT(:,1); % XWo
    Endmembers(:,2)=XnormT(:,2); % XFs
    Endmembers(:,3)=XnormT(:,3); % XEn
    Endmembers(:,4)=XnormT(:,4); % XJd
    Endmembers(:,5)=XnormT(:,5); % XAeg
    Endmembers(:,6)=XnormT(:,6); % XKos
    Endmembers(:,7)=(XnormT(:,1)+XnormT(:,2)+XnormT(:,3))./(XnormT(:,1)+XnormT(:,2)+XnormT(:,3)+XnormT(:,4)+XnormT(:,5)+XnormT(:,6)); % XQuad

    %limit on significant digits (eliminates rounding noise)
    StrctFrm(StrctFrm<1e-6) = 0;
    Endmembers(Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)
    apfu(apfu<1e-6) = 0;
    apfu(:,13)=Opfu-O2total; %also adds O2 def to the apfu output

    all=[StrctFrm Endmembers apfu(:,13)];
    all(isnan(all))=0; %remove NaN


    apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn2','Mg','Ca','Na','K','Cation_Sum','O2_deficiency'});
    StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M1','Ti_M1','Cr_M1','Fe3_M1','Mn_M1','Mg_M1','Fe2_M1','Sum_M1','Mg_M2','Fe2_M2','Ca_M2','Na_M2','K_M2','Sum_M2','Xwo','Xfs','Xen','Xjd','Xaeg','Xkos','Xquad','O2_deficiency'});

else
    %ortho and calcic pyroxene
    %Normalization procedure follows Morimoto et al. (1988)
    Endmembers(:,1)=(apfu(:,9)./(apfu(:,9)+apfu(:,6)+apfu(:,8)+apfu(:,7)+apfu(:,4))); % XWo
    Endmembers(:,2)=(apfu(:,6)+apfu(:,7)+apfu(:,4))./(apfu(:,9)+apfu(:,6)+apfu(:,8)+apfu(:,7)+apfu(:,4)); % XFs
    Endmembers(:,3)=(apfu(:,8)./(apfu(:,9)+apfu(:,6)+apfu(:,8)+apfu(:,7)+apfu(:,4))); % XEn
    Endmembers(:,7)=Endmembers(:,1)+Endmembers(:,2)+Endmembers(:,3);

    %limit on significant digits (eliminates rounding noise)
    StrctFrm(StrctFrm<1e-6) = 0;
    Endmembers(Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)
    apfu(apfu<1e-6) = 0;
    apfu(:,13)=Opfu-O2total; %also adds O2 def to the apfu output

   
    all=[StrctFrm Endmembers apfu(:,13)];
    all(isnan(all))=0; %remove NaN


    apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn2','Mg','Ca','Na','K','Cation_Sum','O2_deficiency'});
    StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M1','Ti_M1','Cr_M1','Fe3_M1','Mn_M1','Mg_M1','Fe2_M1','Sum_M1','Mg_M2','Fe2_M2','Ca_M2','Na_M2','K_M2','Sum_M2','Xwo','Xfs','Xen','Xjd','Xaeg','Xkos','Xquad','O2_deficiency'});

end

end