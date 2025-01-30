function [ax2plot,options_definition]=plot_tern_fun(ax2plot,options)
if not(exist('ax2plot','var'))
    ax2plot=[];
end
if not(exist('options','var'))
    options=[];
end
%% options definition
%Parameter 1
options_definition.plot_tern.Value=true;
%Parameter 2

options_definition.plot_tern.Value=true;
options_definition.background_FaceAlpha.Value=1;

options_definition.type.Value='none';
options_definition.type.options={'none';
    'tern_woenfe'; 
    'tern_qjdaeg'; 
    'tern_woenfe'; 
    'tern_woenfe'; 
    'tern_qjdaeg'; 
    'tern_feldspar_simple';
    'tern_feldspar_complex';   
    'tern_ol';
    'tern_grt'; 
    'tern_grt_half'; 
    'tern_grtoct';
    'tern_grtTi';
    'tern_mscelprl';
    'tern_mica_FOHCl';
    'tern_ctd';
    'tern_sp1';
    'tern_sp2';
    'tern_sp1';
    'tern_sp2';
    'tern_ap';
    'tern_scap';
    'tern_lws'%;
    %'tern_carbonate_Ca_Mg_Fe;
    %'tern_carbonate_Ca_Mg_Mn;
    };
options_definition.type.description={'none';
    'clinopyroxene: wollastonite-enstatite-ferrosillite ternary';
    'clinopyroxene Quad-jadeite-aegirine ternary'
    'orthopyroxene: wollastonite-enstatite-ferrosillite ternary';
    'pyroxene: wollastonite-enstatite-ferrosillite ternary';
    'pyroxene Quad-jadeite-aegirine ternary'
    'feldspar: only labels the corners';
    'feldspar: labels corners ';      %  'feldspar: labels corners rotated';
    'olivine: ternary';
    'garnet: Alm+Sps-Prp-Grs ternary';
    'garnet: Alm+Sps-Prp-Grs ternary half';
    'garnet: Al-Fe3+-Ti on the octahedral site';
    'garnet: andradite-morimotoite-schorlomite ternary';
    'mica: muscovite- aluminoceladonite-pyrophyllite ternary';
    'mica: mica hydroxyl site ternary';
    'chloritoid: Fe2-Mg-Mn';
    'oxyspinel: Cr-Fe3 + 2Ti-Al ternary';
    'oxyspinel: Mn+Zn+Ni-Fe2-Mg ternary';
    'spinel: Cr-Fe3 + 2Ti-Al ternary';
    'spinel: Mn+Zn+Ni-Fe2-Mg ternary';
    'apatite: F-OH + S-Cl ternary';
    'scapolite: CO3-Cl-SO4 ternary';
    'lawsonite: 4Ti-Cr-Fe ternary'%;
    %'carbonate: Ca-Mg-Fe;
    % carbonate: Ca-Mg-Mn
    };

options_definition.custom.Value=false; %endmember text labels
options_definition.custom.String={'X';'Y';'Z'}; %endmember text labels

options_definition.FontSize.Value=12;

if not(exist('options','var')) || isempty(options)
    if not(exist('ax2plot','var'))
        ax2plot=[];
    end
    return
end

if not(exist('ax2plot','var')) || isempty(ax2plot)
    hfig=figure;
    ax2plot=nexttile;
    try
        CenterFig_fun(hfig)
    catch ME
        disp(ME.message)
    end
end
hold on

if not(isfield(options,'FontSize')) || not(isfield(options.FontSize,'Value'))
    options.FontSize.Value=options_definition.FontSize.Value;
end
ax2plot.FontSize=options.FontSize.Value;

if options.plot_tern.Value==true && not(strcmp(options.type.Value,'tern_grt_half'))
    % read from file or define matrix
    pgon=polyshape([0 0.5 1],[0 sqrt(3)/2 0]);

    plot(pgon,'FaceColor','w')
    try
        if isfield(options,'background_FaceAlpha')
            ax2plot.Children(1).FaceAlpha=options.background_FaceAlpha.Value;
        else
            ax2plot.Children(1).FaceAlpha=1;
        end
        ax2plot.Children(1).HandleVisibility='off';

        % if background only
        %     return
        % end
        %% possible adjustment
        linecolor_1=[0.3 0.3 0.3]; %color for odd plotting intervals
        linecolor_2=[0.3 0.3 0.3]; %color for even plotting intervals
        LineStyle_1=':'; %linestyle for odd plotting intervals
        LineStyle_2='--'; %linestyle for even plotting intervals
        linewidth=0.5; %linewidth

        %matrix for boundary of triangle
        matrx_1=[0 1 0 0;
            0 0.5 0 sqrt(3)/2;
            0.5 1 sqrt(3)/2 0];

        %plot triangle outline
        plot(ax2plot,matrx_1(:,1:2)', matrx_1(:,3:4)','color','k','linewidth',1.5,'HandleVisibility','off')

        %matrix of odd plotting intervals
        matrx_2=[0.45 0.55  0.779422863 0.779422863;
            0.35 0.65 0.606217783 0.606217783;
            0.25 0.75 0.433012702 0.433012702;
            0.15 0.85 0.259807621 0.259807621;
            0.05 0.95 0.08660254 0.08660254;
            0.05 0.1 0.08660254 0;
            0.15 0.3 0.259807621 0;
            0.25 0.5 0.433012702 0;
            0.35 0.7 0.606217783 0;
            0.45 0.9 0.779422863 0;
            0.95 0.9 0.08660254 0;
            0.85 0.7 0.259807621 0;
            0.75 0.5 0.433012702 0;
            0.65 0.3 0.606217783 0;
            0.55 0.1 0.779422863 0];

        %matrix of even plotting intervals
        matrx_3=[0.4 0.6 0.692820323 0.692820323;
            0.3 0.7 0.519615242 0.519615242;
            0.2 0.8 0.346410162 0.346410162;
            0.1 0.9 0.173205081 0.173205081;
            0.1 0.2 0.173205081 0;
            0.2 0.4 0.346410162 0;
            0.3 0.6 0.519615242 0;
            0.4 0.8 0.692820323 0;
            0.9 0.8 0.173205081 0;
            0.8 0.6 0.346410162 0;
            0.7 0.4 0.519615242 0;
            0.6 0.2 0.692820323 0];

        %plot odd intervals
        plot(ax2plot,matrx_2(:,1:2)', matrx_2(:,3:4)','color',linecolor_1,'LineStyle',LineStyle_1,'linewidth',linewidth,'HandleVisibility','off')
        %plot even intervals
        plot(ax2plot,matrx_3(:,1:2)', matrx_3(:,3:4)','color',linecolor_2,'LineStyle',LineStyle_2,'linewidth',linewidth,'HandleVisibility','off')

        %contour interval labels for ternary numerical labels
        str_01={'0.0';'0.2';'0.4';'0.6';'0.8';'1.0'};

        %label positions for ternary bottom
        pos_01=[-0.02,-0.039;
            0.18,-0.039;
            0.38,-0.039;
            0.58,-0.039;
            0.78,-0.039;
            0.98,-0.039];

        %label positions for ternary right side
        pos_02=[1.02,0.01;
            0.92,0.18;
            0.82,0.35;
            0.72,0.52;
            0.62,0.70;
            0.52,0.87];

        %label positions for ternary left side
        pos_03=[0.482,0.901;
            0.382,0.735;
            0.278,0.558;
            0.178,0.386;
            0.075,0.215;
            -0.025,0.040];

        %plot ternary numerical labels
        text(ax2plot,pos_01(:,1),pos_01(:,2),str_01,'FontSize',options.FontSize.Value,'HorizontalAlignment','center','Rotation',60)
        text(ax2plot,pos_02(:,1),pos_02(:,2),str_01,'FontSize',options.FontSize.Value)
        text(ax2plot,pos_03(:,1),pos_03(:,2),str_01,'FontSize',options.FontSize.Value,'HorizontalAlignment','center','Rotation',300)
    catch ME
        
    end
elseif (options.plot_tern.Value==true && strcmp(options.type.Value,'tern_grt_half')) || (isfield(options,'plot_half_tern') && options.plot_half_tern.Value==true)
    X=[0  0.2 0.8 1];
    Y=[0 0.346410162 0.346410162 0];
    patch(X,Y,'w','HandleVisibility','off')
    %% possible adjustment
    linecolor_1=[0.3 0.3 0.3]; %color for odd plotting intervals
    linecolor_2=[0.3 0.3 0.3]; %color for even plotting intervals
    LineStyle_1=':'; %linestyle for odd plotting intervals
    LineStyle_2='--'; %linestyle for even plotting intervals
    linewidth=0.5; %linewidth

    %matrix for boundary of triangle
    matrx_1=[0 1 0 0;
        0 0.2 0 0.346410162;
        0.8 1 0.346410162 0;
        0.2 0.8  0.346410162   0.346410162  ;
        ];

    %plot triangle outline
    plot(ax2plot,matrx_1(:,1:2)', matrx_1(:,3:4)','color','k','linewidth',1.5,'HandleVisibility','off')

    %matrix of odd plotting intervals
    matrx_2=[  0.15 0.85 0.259807621 0.259807621;
        0.05 0.95 0.08660254 0.08660254;
        0.05 0.1 0.08660254 0;
        0.15 0.3 0.259807621 0;
        0.3 0.5 0.346410162 0;
        0.5 0.7 0.346410162 0;
        0.7 0.9 0.346410162 0;
        0.95 0.9 0.08660254 0;
        0.85 0.7 0.259807621 0;
        0.7 0.5 0.346410162 0;
        0.5 0.3 0.346410162 0;
        0.3 0.1 0.346410162 0
        ];

    %matrix of even plotting intervals
    matrx_3=[   0.2 0.8 0.346410162 0.346410162;
        0.1 0.9 0.173205081 0.173205081;
        0.1 0.2 0.173205081 0;
        0.2 0.4 0.346410162 0;
        0.8 0.6 0.346410162 0;
        0.9 0.8 0.173205081 0;
        0.6 0.8 0.346410162 0;
        0.4 0.6 0.346410162 0;
        0.6 0.4 0.346410162 0;
        0.4 0.2 0.346410162 0];

    %plot odd intervals
    plot(ax2plot,matrx_2(:,1:2)', matrx_2(:,3:4)','color',linecolor_1,'LineStyle',LineStyle_1,'linewidth',linewidth,'HandleVisibility','off')
    %plot even intervals
    plot(ax2plot,matrx_3(:,1:2)', matrx_3(:,3:4)','color',linecolor_2,'LineStyle',LineStyle_2,'linewidth',linewidth,'HandleVisibility','off')

    %contour interval labels for ternary numerical labels
    str_01={'0.0';'0.2';'0.4';'0.6';'0.8';'1.0'};

    %label positions for ternary bottom
    pos_01=[-0.02,-0.039;
        0.18,-0.039;
        0.38,-0.039;
        0.58,-0.039;
        0.78,-0.039;
        0.98,-0.039];

    %label positions for ternary right side
    pos_02=[1.02,0.01;
        0.92,0.18;
        0.82,0.35;
        ];

    %label positions for ternary left side
    pos_03=[ 0.178,0.386;
        0.075,0.215;
        -0.025,0.040];

    %plot ternary numerical labels
    text(ax2plot,pos_01(:,1),pos_01(:,2),str_01,'FontSize',options.FontSize.Value,'HorizontalAlignment','center','Rotation',60)
    text(ax2plot,pos_02(:,1),pos_02(:,2),str_01(1:3),'FontSize',options.FontSize.Value)
    text(ax2plot,pos_03(:,1),pos_03(:,2),str_01(end-2:end),'FontSize',options.FontSize.Value,'HorizontalAlignment','center','Rotation',300)

end

%% Mineral Specific options

if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_woenfe')
    %% clinopyroxene wollastonite-enstatite-ferrosillite ternary

    %field boundaries
    fld_cpx1=[0.25 0.75 0.4330127 0.4330127;
        0.225 0.775 0.38971143 0.38971143;
        0.1 0.9 0.17320508 0.17320508;
        0.025 0.975 0.04330127 0.04330127;
        0.5 0.5 0.38971143 0.4330127;
        0.5 0.5 0 0.04330127];

    %plot field boundaries
    plot(ax2plot,fld_cpx1(:,1:2)', fld_cpx1(:,3:4)','k','linewidth',1.5,'HandleVisibility','off')

    %text string for ternary clinopyroxene labels
    str_cpx1={'en';'wo';'fs'};
    str_cpx2={'Enstatite';'Ferrosilite';'Pigeonite';'Augite';'Diopside';'Hedenbergite'};

    %positions
    pos_cpx1=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04];

    pos_cpx2=[0.20,0.02;
        0.65,0.02;
        0.425,0.1;
        0.45,0.28;
        0.29,0.412;
        0.52,0.412];

    %text plotting
    text(ax2plot,pos_cpx1(:,1),pos_cpx1(:,2),str_cpx1,'FontSize',options.FontSize.Value*1.1667)
    text(ax2plot,pos_cpx2(:,1),pos_cpx2(:,2),str_cpx2,'FontSize',options.FontSize.Value)


elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_qjdaeg')
    %% clinopyroxene field boundaries for quad-jadeite-aegirine ternary


    %field boundaries
    fld_cpx2=[0.4 0.6 0.69282032 0.69282032;
        0.1 0.9 0.17320508 0.17320508;
        0.5 0.5 0 0.69282032];

    %plot field boundaries
    plot(ax2plot,fld_cpx2(:,1:2)', fld_cpx2(:,3:4)','k','linewidth',1.5,'HandleVisibility','off')


    %text string for ternary clinopyroxene labels
    str_cpx3={'jd';'quad';'aeg'};
    str_cpx4={'Jadeite';'Aegirine';'Omphacite';'Aegirine-';'Augite'};

    %positions
    pos_cpx3=[-0.10,0.0;
        0.52,0.92;
        1.02,-0.04];

    pos_cpx4=[0.20,0.09;
        0.66,0.09;
        0.26,0.35;
        0.58,0.37;
        0.60,0.33];

    %text plotting
    text(ax2plot,pos_cpx3(:,1),pos_cpx3(:,2),str_cpx3,'FontSize',options.FontSize.Value*1.1667)
    text(ax2plot,pos_cpx4(:,1),pos_cpx4(:,2),str_cpx4,'FontSize',options.FontSize.Value)

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_feldspar_simple')
    %% feldspar simple ternary

    %ternary endmember labels
    str_fsp1={'ab';'an';'or'};

    %positions for text labels
    pos_fsp1=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04];

    %plot endmember text labels
    text(ax2plot,pos_fsp1(:,1),pos_fsp1(:,2),str_fsp1,'FontSize',options.FontSize.Value*1.1667)

elseif not(isempty(options)) && isfield(options,'type') && (strcmp(options.type.Value,'tern_feldspar_complex')|| strcmp(options.type.Value,'tern_feldspar_complex2'))
    %% feldspar ternary with subdivision boundaries

    %field/subdivision boundaries
    fld_fsp = [0.05 0.145 0.08660254 0.077942286;
        0.1 0.145 0 0.077942286;
        0.145 0.225 0.077942286 0.129903811;
        0.15 0.235 0.259807621 0.233826859;
        0.225 0.235 0.129903811 0.233826859;
        0.25 0.325 0.433012702 0.389711432;
        0.235 0.325 0.233826859 0.389711432;
        0.35 0.415 0.606217783 0.545596004;
        0.325 0.415 0.389711432 0.545596004;
        0.45 0.505 0.779422863 0.701480577;
        0.415 0.505 0.545596004 0.701480577;
        0.55 0.505 0.779422863 0.701480577;
        0.37 0.383 0 0.08660254;
        0.383 0.225 0.08660254 0.129903811;
        0.383 0.95 0.08660254 0.08660254];

    %plot subdivision boundaries
    h1=plot(ax2plot,fld_fsp(:,1:2)', fld_fsp(:,3:4)','color','k','linewidth',1.5,'HandleVisibility','off');

    str_fsp1={'ab';'an';'or'}; %ternary endmember labels
    str_fsp2={'olig';'ands';'labr';'bytw';'ano';'mc/or/sa'}; %subdivision text labels

    %positions for endmember text labels
    pos_fsp1=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04];

    %positions for subdivision text labels
    pos_fsp2=[0.24,0.18;
        0.29,0.30;
        0.38,0.46;
        0.48,0.62;
        0.30,0.14;
        0.62,0.11];

    %plot text labels
    h2=text(ax2plot,pos_fsp1(:,1),pos_fsp1(:,2),str_fsp1,'FontSize',options.FontSize.Value*1.1667); %endmembers
    h3=text(ax2plot,pos_fsp2(:,1),pos_fsp2(:,2),str_fsp2,'FontSize',options.FontSize.Value); %subdivisions

    % if strcmp(options.type.Value,'tern_feldspar_complex2')
    %
    %     pos_fsp2=pos_fsp2.*[1 3];
    %     pos_fsp1=pos_fsp1.*[1 3];
    %     fld_fsp=fld_fsp.*[1 1 3 3];
    %       h1=plot(ax2plot,fld_fsp(:,1:2)', fld_fsp(:,3:4)','color','k','linewidth',1.5,'HandleVisibility','off');
    %
    %       h2=text(ax2plot,pos_fsp1(:,1),pos_fsp1(:,2),str_fsp1,'FontSize',options.FontSize.Value*1.1667); %endmembers
    %     h3=text(ax2plot,pos_fsp2(:,1),pos_fsp2(:,2),str_fsp2,'FontSize',options.FontSize.Value); %subdivisions
    %
    %
    %
    %     center =[0.5.*(1/3)+(1/3)  (1/3)*(cos(30*pi()/180))];
    %
    % rotate(h1,center,-60)
    % rotate(h2,center,-60)
    % rotate(h3,center,-60)
    % end

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_ol')
    %% olivine ternary

    %monticellite-kirschsteinite line
    plot([0.25 0.75], [0.433012702 0.433012702],'k','linewidth',1.5,'HandleVisibility','off')

    str_ol1={'fo';'Ca-ol';'fa + tep';'mic';'kir'}; %endmember text labels

    %positions for endmember text labels
    pos_ol1=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04;
        0.16,0.44;
        0.77,0.44];

    text(ax2plot,pos_ol1(:,1),pos_ol1(:,2),str_ol1,'FontSize',options.FontSize.Value*1.1667) %endmembers


elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_grt')
    %% garnet standard ternary

    str_grt1={'alm';'+ sps';'grs';'prp'}; %endmember text labels

    %positions for endmember text labels
    pos_grt1=[-0.14,0.02;
        -0.16,-0.02;
        0.51,0.92;
        1.02,-0.04];

    text(ax2plot,pos_grt1(:,1),pos_grt1(:,2),str_grt1,'FontSize',options.FontSize.Value*1.1667) %endmembers
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_grt_half')
    %% garnet standard ternary

    str_grt1={'X_{Alm+Sps}';'X_{Grs}';'X_{Prp}'}; %endmember text labels

    %positions for endmember text labels
    pos_grt1=[-0.016,0.22;
        %-0.16,-0.02;
        1.016 0.21;
        0.5,-0.1];

    text(ax2plot,pos_grt1(1,1),pos_grt1(1,2),str_grt1(1),'FontSize',options.FontSize.Value*1.1967,'Rotation',60,'HorizontalAlignment','center','FontWeight','bold') %endmembers
    text(ax2plot,pos_grt1(2,1),pos_grt1(2,2),str_grt1(2),'FontSize',options.FontSize.Value*1.1967,'Rotation',-60,'HorizontalAlignment','center','FontWeight','bold') %endmembers
    text(ax2plot,pos_grt1(3,1),pos_grt1(3,2),str_grt1(3),'FontSize',options.FontSize.Value*1.1967,'Rotation',0,'HorizontalAlignment','center','FontWeight','bold') %endmembers


elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_grtoct')
    %% garnet  Cr-Fe3+-Al on octahedral site ternary

    str_grt2={'Al';'Fe^{3+}';'Cr + Ti'}; %endmember text labels

    %positions for endmember text labels
    pos_grt2=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04];

    text(ax2plot,pos_grt2(:,1),pos_grt2(:,2),str_grt2,'FontSize',options.FontSize.Value*1.1667) %endmembers


elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_grtTi')
    %% garnet andradite-morimotoite-schorlomite ternary

    %text labels
    str_grt3={'R^{3+}_{2}';'R^{4+}R^{2+}';'R^{4+}_{2}'}; %endmember text labels
    str_grt4={'Andradite';'Schorlomite';'Morimotoite'}; %field text labels

    %positions for endmember text labels
    pos_grt3=[-0.10,0.00;
        0.51,0.92;
        1.02,-0.05];

    %positions for field text labels
    pos_grt4=[0.19,0.14;
        0.66,0.14;
        0.42,0.50];

    text(ax2plot,pos_grt3(:,1),pos_grt3(:,2),str_grt3,'FontSize',options.FontSize.Value*1.1667) %endmembers
    text(ax2plot,pos_grt4(:,1),pos_grt4(:,2),str_grt4,'FontSize',options.FontSize.Value*1.1667) %field boundaries

    %field boundaries
    % fld_grt1=[0.5 0.5 0 sqrt(1/3)/2;
    %     0.25 0.5 0.433012702 sqrt(1/3)/2;
    %     0.75 0.5 0.433012702 sqrt(1/3)/2];

    fld_grt1=[0.25 0.624999783 0.433012702 0.216506125;
        0.75 0.5 0.433012702 0];

    %plot field boundaries
    plot(ax2plot,fld_grt1(:,1:2)', fld_grt1(:,3:4)','k','linewidth',1.5,'HandleVisibility','off')


elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_mscelprl')
    %% mica muscovite-celadonite-pyrophyllite ternary

    str_mica1={'ms';'Alcel';'prl'}; %endmember text labels

    %positions for endmember text labels
    pos_mica1=[-0.10,0.0;
        0.50,0.92;
        1.02,-0.04];

    text(ax2plot,pos_mica1(:,1),pos_mica1(:,2),str_mica1,'FontSize',options.FontSize.Value*1.1667) %endmembers

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_mica_FOHCl')
    %% mica F-OH-Cl ternary

    str_mica2={'Cl';'F';'OH'}; %endmember text labels

    %positions for endmember text labels
    pos_mica2=[-0.10,0.0;
        0.50,0.92;
        1.02,-0.04];

    text(ax2plot,pos_mica2(:,1),pos_mica2(:,2),str_mica2,'FontSize',options.FontSize.Value*1.1667) %endmembers

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_ctd')
    %% chloritoid ternary

    str_ctd={'Mg';'Fe^{2+}';'Mn'}; %endmember text labels

    %positions for endmember text labels
    pos_ctd=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04];

    text(ax2plot,pos_ctd(:,1),pos_ctd(:,2),str_ctd,'FontSize',options.FontSize.Value*1.1667) %endmembers

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_sp1')
    %% spinel Cr-Fe3 + 2Ti-Al ternary

    str_sp1={'Cr';'Fe^{3+}';'Al'}; %endmember text labels

    %positions for endmember text labels
    pos_sp1=[-0.10,0.0;
        0.51,0.92;
        1.02,-0.04];

    text(ax2plot,pos_sp1(:,1),pos_sp1(:,2),str_sp1,'FontSize',options.FontSize.Value*1.1667) %endmembers

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_sp2')
    %% spinel Mn+Zn+Ni-Fe2-Mg ternary
    str_sp2={'Mn +';'Zn + Ni';'Fe^{2+}';'Mg'}; %endmember text labels

    %positions for endmember text labels
    pos_sp2=[-0.18,0.03;
        -0.18,-0.02;
        0.51,0.92;
        1.02,-0.04];

    text(ax2plot,pos_sp2(:,1),pos_sp2(:,2),str_sp2,'FontSize',options.FontSize.Value*1.1667) %endmembers



elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_ap')
    %% apatite F-OH-Cl-S ternary
    str_ap1={'Cl';'F';'OH + S + C'}; %endmember text labels

    %positions for endmember text labels
    pos_ap1=[-0.10,0.0;
        0.50,0.92;
        1.02,-0.04];

    text(ax2plot,pos_ap1(:,1),pos_ap1(:,2),str_ap1,'FontSize',options.FontSize.Value*1.1667) %endmembers


elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_scap')
    %% scapolite CO3-Cl-SO4 ternary
    str_scap1={'CO_{3}^{2-}';'Cl';'SO_{4}^{2-}'}; %endmember text labels

    %positions for endmember text labels
    pos_scap1=[-0.14,0.0;
        0.50,0.92;
        1.02,-0.04];

    text(ax2plot,pos_scap1(:,1),pos_scap1(:,2),str_scap1,'FontSize',options.FontSize.Value*1.1667) %endmembers

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_lws')
    %% lawsonite 4Ti-Cr-Fe ternary

    str_lws1={'4 x Ti';'Cr';'Fe'}; %endmember text labels

    %positions for endmember text labels
    pos_lws1=[-0.14,0.0;
        0.50,0.92;
        1.02,-0.04];

    text(ax2plot,pos_lws1(:,1),pos_lws1(:,2),str_lws1,'FontSize',options.FontSize.Value*1.1667) %endmembers
    
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_carbonate_Ca_Mg_Fe')
    %% Carbonate Ca-Mg-Fe ternary

    str_carb1={'Ca';'Mg';'Fe'}; %endmember text labels

    %positions for endmember text labels
    pos_carb1=[-0.14,0.0;
        0.50,0.92;
        1.02,-0.04];

    text(ax2plot,pos_carb1(:,1),pos_carb1(:,2),str_carb1,'FontSize',options.FontSize.Value*1.1667) %endmembers

 elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'tern_carbonate_Ca_Mg_Mn')
    %% Carbonate Ca-Mg-Mn ternary

    str_carb2={'Ca';'Mg';'Mn'}; %endmember text labels

    %positions for endmember text labels
    pos_carb2=[-0.14,0.0;
        0.50,0.92;
        1.02,-0.04];
    text(ax2plot,pos_carb2(:,1),pos_carb2(:,2),str_carb2,'FontSize',options.FontSize.Value*1.1667) %endmembers

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'custom_tern')

    options.custom.String={'X';'Y';'Z'}; %endmember text labels

    if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'xlabel') && not(isempty(options.custom.xlabel))
        options.custom.String{1} =options.custom.xlabel;
    end
    if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'ylabel') && not(isempty(options.custom.ylabel))
        options.custom.String{2} =options.custom.ylabel;
    end
    if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'zlabel') && not(isempty(options.custom.zlabel))
        options.custom.String{3} =options.custom.zlabel;
    end

    options.custom.Value=true;



end

if not(isempty(options)) && isfield(options,'custom') &&  isfield(options.custom,'Value') && options.custom.Value==true
    if not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'Interpreter')
        Interpreter=options.custom.Interpreter;
    else
        Interpreter='tex';
    end

    if  iscell(options.custom.String) && numel(options.custom.String)==3
        str_label=options.custom.String; %endmember text labels
    else
        str_label={'X';'Y';'Z'}; %endmember text labels
    end
    %positions for endmember text labels
    pos_label=[-0.04,-0.04;
        0.50,0.95;
        1.02,-0.04];

    %         pos_label=[-0.14,0.0;
    %         0.50,0.92;
    %         1.02,-0.04];

    text(ax2plot,pos_label(1,1),pos_label(1,2),str_label(1),'FontSize',options.FontSize.Value*1.1667,'HorizontalAlignment','right','Interpreter',Interpreter) %endmembers)
    text(ax2plot,pos_label(2,1),pos_label(2,2),str_label(2),'FontSize',options.FontSize.Value*1.1667,'HorizontalAlignment','center','Interpreter',Interpreter) %endmembers)) %endmembers
    text(ax2plot,pos_label(3,1),pos_label(3,2),str_label(3),'FontSize',options.FontSize.Value*1.1667,'HorizontalAlignment','left','Interpreter',Interpreter) %endmembers
end
%axis equal
ax2plot.DataAspectRatio=[1 1 1];
ax2plot.PlotBoxAspectRatio=[1 1 1 ];

ax2plot.XAxis.Visible=false;
ax2plot.YAxis.Visible=false;
ax2plot.Color='none';
end