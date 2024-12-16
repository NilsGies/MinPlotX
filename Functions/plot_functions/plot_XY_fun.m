function [ax2plot,options_definition]=plot_XY_fun(ax2plot,options)
if not(exist('ax2plot','var'))
    ax2plot=[];
end
if not(exist('options','var'))
    options=[];
end
%% options definition
options_definition.XYplot.Value=true;

%change all to linecolor_1 and define as the following to make it failsafe
if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'linecolor_1')) || isempty(options.plot_settings.linecolor_1)
    linecolor_1=[0.3 0.3 0.3];
    options_definition.plot_settings.linecolor_1.Value=linecolor_1;
    options_definition.linecolor_1.description='Linecolor for odd plotting intervals';
else
    linecolor_1=options.plot_settings.linecolor_1.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'linecolor_2')) || isempty(options.plot_settings.linecolor_2)
    linecolor_2=[0.3 0.3 0.3];
    options_definition.plot_settings.linecolor_2.Value=linecolor_2;
    options_definition.linecolor_2.description='Linecolor for even plotting intervals';
else
    linecolor_2=options.plot_settings.linecolor_2.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'LineStyle_1')) || isempty(options.plot_settings.LineStyle_1)
    LineStyle_1=':';
    options_definition.plot_settings.LineStyle_1.Value=LineStyle_1;
    options_definition.LineStyle_1.description='LineStyle for odd plotting intervals';
else
    LineStyle_1=options.plot_settings.LineStyle_1.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'LineStyle_2')) || isempty(options.plot_settings.LineStyle_2)
    LineStyle_2='--';
    options_definition.plot_settings.LineStyle_2.Value=LineStyle_2;
    options_definition.LineStyle_2.description='LineStyle for even plotting intervals';
else
    LineStyle_2=options.plot_settings.LineStyle_2.Value;
end


if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'linewidth_1')) || isempty(options.plot_settings.linewidth_1)
    linewidth_1=0.5;
    options_definition.plot_settings.linewidth_1.Value=linewidth_1;
else
    linewidth_1=options.plot_settings.linewidth_1.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'FontSize')) || isempty(options.FontSize)
    FontSize=12;
    options_definition.FontSize.Value=12;
else
    FontSize=options.FontSize.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'AutoXlim')) || isempty(options.AutoXlim)
    AutoXlim=false;
    options_definition.AutoXlim.Value=AutoXlim;
else
    AutoXlim=options.AutoXlim.Value;
end

options_definition.type.Value='none';
options_definition.type.options={'none';
    'cpx_QJ_XY';
    'cpx_catscaes';
    'cpx_Kcpxcaes';
    'cpx_QJ_XY';
    'cpx_catscaes';
    'cpx_Kcpxcaes';
    'ol_NiO_XY';
    'amph_XY_Ca1';
    'amph_XY_Ca2';
    'amph_XY_NaCa';
    'amph_XY_Na1';
    'amph_XY_Na2';
    'amph_XY_Fe';
    'mica_XY_ph';
    'mica_XY_SiNa';
    'mica_XY_NaKFeMg';
    'mica_XY_bio';
    'chl_XY_simp';
    'chl_XY_comp';
    'chl_XY_R2Si';
    'chl_XY_XFeSi';
    'aln_XY';
    'aln_XY';
    'ttn_XY_Ti';
    'ttn_XY_F';
    'crd_XY_NaAlSi';
    'crd_XY_AlMeNaSi';
    'crd_XY_MevsSiAl';
    'lws_XY_FeTiCrAl'
    'universal_MgO_Na2O'};

options_definition.type.description={'none';
    'clinopyroxene: Q-J diagram';
    'clinopyroxene: Ca-tschermaks vs Ca-eskolaite diagram';
    'clinopyroxene: K-Cpx vs Ca-eskolaite diagram';
    'pyroxene: Q-J diagram';
    'pyroxene: Ca-tschermaks vs Ca-eskolaite diagram';
    'pyroxene: K-Cpx vs Ca-eskolaite diagram';
    'olivine: NiO diagram';
    'amphibole: first Ca-clinoamphibole classification diagram';
    'amphibole: second Ca-clinoamphibole classification diagram';
    'amphibole: NaCa-clinoamphibole classification diagram';
    'amphibole: first Na-clinoamphibole classification diagram';
    'amphibole: second Na-clinoamphibole classification diagram';
    'amphibole: amphibole Fe3/Fe vs total Fe diagram';
    'mica: Al vs Si diagram';
    'mica: Si vs Na diagram';
    'mica: Na/(Na+K) vs Mg + Fe'
    'mica: biotite mica diagram';
    'chlorite: simple chamosite-clinochlore-sudoite diagram';
    'chlorite: complex version with wider composition space';
    'chlorite: divalent cations vs Si';
    'chlorite: XFe vs Si';
    'allanite: aln_XY';
    'epidote: aln_XY';
    'titanite: Al + Fe3 vs Ti diagram';
    'titanite: Al + Fe3 vs F diagram';
    'cordierite: Na + Al vs Si diagram';
    'cordierite: Al + Me^{2+} vs Na + Si diagram';
    'cordierite: Mg + Mn + Fe vs Si + Al diagram';
    'lawsonite: Fe + Ti + Cr vs Al diagram'
    'universal: MgO vs Na2O'
    };


if not(exist('options','var')) || isempty(options)
    if not(exist('ax2plot','var'))
        ax2plot=[];
    end
    return
end

if not(exist('ax2plot','var')) || isempty(ax2plot)
    figure
    ax2plot=nexttile;
end
hold on
       

ax2plot.FontSize=FontSize;

if (isfield(options,'default') && options.default==true)
    type=options.type.Value;
    options=options_definition;
    options.type.Value=type;
end

if options.XYplot.Value==true

    %% Mineral specific options

    %Cpx Q-J diagram
    if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'cpx_QJ_XY')

        %XY limits
        xlim([0 2])
        ylim([0 2])

        % %XY ticks and labels
        % ax2plot.YAxis.TickValues=0:0.2:2;
        % ax2plot.XAxis.TickValues=0:0.2:2;
        % ax2plot.YTickLabel=({'0.0',[],'0.4',[],'0.8',[],'1.2',[],'1.6',[],'2.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.4',[],'0.8',[],'1.2',[],'1.6',[],'2.0'});

        %field boundaries
        fld_cpxQJ=[0 2 2 0;
            0 1.5 1.5 0;
            0.3 0.4 1.2 1.6;
            1.2 1.6 0.3 0.4];

        plot(ax2plot,fld_cpxQJ(:,1:2)', fld_cpxQJ(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %text labels
        str_cpx1={'quad';'NaCa';'Na'};

        %position of text labels
        pos_cpxQJ=[0.25,1.8;
            1.05,1.00;
            1.85,0.20];

        %text plotting
        text(ax2plot,pos_cpxQJ(:,1),pos_cpxQJ(:,2),str_cpx1,'FontSize',FontSize)

        %axis labels
        xlabel('J = 2Na (apfu)')
        ylabel('Q = (Ca + Mg + Fe^{2+}) (apfu)')

        axis square
        box on

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'cpx_catscaes')

        xlabel('X_{caes}')
        ylabel('X_{ks}')

        %XY limits
        xlim([0 1])
        ylim([0 1])
        % 
        % %XY ticks and labels
        % ax2plot.YAxis.TickValues=0:0.1:1;
        % ax2plot.XAxis.TickValues=0:0.1:1;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        axis square %makes the spacing of the axes intervals equal
        box on

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'cpx_Kcpxcaes')

        xlabel('X_{caes}')
        ylabel('X_{K-Cpx}')

        %XY limits
        xlim([0 1])
        ylim([0 1])

        % %XY ticks and labels
        % ax2plot.YAxis.TickValues=0:0.1:1;
        % ax2plot.XAxis.TickValues=0:0.1:1;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        axis square %makes the spacing of the axes intervals equal
        box on
        
        %olivine Fo vs NiO plot
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_NiO_XY')

        %allows adjustable XY limits
        if AutoXlim==true

            xlim('auto'); %changes X relative to the Fo limits chosen
            options.FoLow=ax2plot.XLim(1);
            options.FoHi=ax2plot.XLim(2);

            ylim('auto'); %changes X relative to the Fo limits chosen
            options.NiLow=ax2plot.YLim(1);
            options.NiHi=ax2plot.YLim(2);

        else %otherwise choose default limits

            %default XY limits
            xlim('auto'); %changes X relative to the Fo limits chosen
            options.FoLow=50;
            options.FoHi=100;

            ylim('auto'); %changes X relative to the Fo limits chosen
            options.NiLow=0;
            options.NiHi=2;

        end

        %adjust XY limits
        ax2plot=gca;
        ylim([options.NiLow options.NiHi]);
        xlim([options.FoLow options.FoHi]);

        %plot axis labels
        xlabel('fo (mol %)')
        ylabel('NiO (wt %)')

        box on

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_XY_Ca1')

        %XY limits
        xlim([0 2])
        ylim([0 1])

        %XY ticks and labels
        ax2plot.YAxis.TickValues=0:0.1:1;
        ax2plot.XAxis.TickValues=0:0.2:2;
        ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        ax2plot.XTickLabel=({'0.0',[],'0.4',[],'0.8',[],'1.2',[],'1.6',[],'2.0'});

        %field boundaries
        fld_amphCa1=[0 2 0.5 0.5;
            0.5 0.5 0.0 1.0;
            1.5 1.5 0.0 1.0];

        plot(ax2plot,fld_amphCa1(:,1:2)', fld_amphCa1(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %text labels
        str_amphCa1={'Tremolite';'Edenite';'Magnesio-Hornblende';'Pargasite';'Sadanagaite';'Tschermakite'};

        %position of text labels
        pos_amphCa1=[0.25,0.25;
            0.25,0.75;
            1.0,0.25;
            1.0,0.75;
            1.75,0.75;
            1.75,0.25];

        %text plotting
        text(ax2plot,pos_amphCa1(:,1),pos_amphCa1(:,2),str_amphCa1,'FontSize',FontSize,'HorizontalAlignment','center')

        %axis labels
        xlabel('^{C}(Al + Fe^{3+} + 2Ti) (apfu)')
        ylabel('^{A}(Na + K + 2Ca) (apfu)')

        box on

        %second Ca clinoamphibole classification diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_XY_Ca2')

        %XY limits
        xlim([5.5 8.0])
        ylim([0 1])

        %XY ticks and labels
        ax2plot.YAxis.TickValues=0:0.1:1;
        ax2plot.XAxis.TickValues=5.5:0.25:8.0;
        ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        ax2plot.XTickLabel=({'5.5',[],'6.0',[],'6.5',[],'7.0',[],'7.5',[],'8.0'});

        %field boundaries

        fld_amphCa2=[5.5 8 0.5 0.5;
            6.5 6.5 0.0 1.0;
            7.5 7.5 0.0 1.0;
            7.5 8.0 0.9 0.9];

        plot(ax2plot,fld_amphCa2(:,1:2)', fld_amphCa2(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %text labels
        str_amphCa2={'Ferro-hornblende';'Ferro-tschermakite';'Magnesio-Hornblende';
            'Tschermakite';'Ferro-actinolite';'Actinolite';'Tremolite'};

        %position of text labels
        pos_amphCa2=[7.0,0.25;
            6.0,0.25;
            7.0,0.75;
            6.0,0.75;
            7.75,0.25;
            7.75,0.75;
            7.75,0.95];

        %text plotting
        text(ax2plot,pos_amphCa2(:,1),pos_amphCa2(:,2),str_amphCa2,'FontSize',FontSize,'HorizontalAlignment','center')

        %axis labels
        xlabel('Si (apfu)')
        ylabel('X_{Mg}')

        box on

        %first NaCa clinoamphibole classification diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_XY_NaCa')

        %XY limits
        xlim([0 2])
        ylim([0 1])

        %XY ticks and labels
        ax2plot.YAxis.TickValues=0:0.1:1;
        ax2plot.XAxis.TickValues=0:0.2:2;
        ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        ax2plot.XTickLabel=({'0.0',[],'0.4',[],'0.8',[],'1.2',[],'1.6',[],'2.0'});

        %field boundaries

        fld_amphNaCa=[0 2 0.5 0.5;
            0.5 0.0 0.0 0.5;
            0.5 0.5 0.5 1.0;
            1.5 1.5 0.0 1.0];

        plot(ax2plot,fld_amphNaCa(:,1:2)', fld_amphNaCa(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %text labels
        str_amphNaCa={'Richterite';'Winchite';'Katophorite';
            'Taramite';'Barroisite'};

        %position of text labels
        pos_amphNaCa=[0.25,0.75;
            1.0,0.25;
            1.0,0.75;
            1.75,0.75;
            1.75,0.25];


        %text plotting
        text(ax2plot,pos_amphNaCa(:,1),pos_amphNaCa(:,2),str_amphNaCa,'FontSize',FontSize,'HorizontalAlignment','center')

        %axis labels
        xlabel('^{C}(Al + Fe^{3+} + 2Ti) (apfu)')
        ylabel('^{A}(Na + K + 2Ca) (apfu)')

        box on

        % first Na clinoamphibole classification diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_XY_Na1')

        %XY limits
        xlim([0 2])
        ylim([0 1])

        %XY ticks and labels
        ax2plot.YAxis.TickValues=0:0.1:1;
        ax2plot.XAxis.TickValues=0:0.2:2;
        ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        ax2plot.XTickLabel=({'0.0',[],'0.4',[],'0.8',[],'1.2',[],'1.6',[],'2.0'});


        %field boundaries
        fld_amphNa1=[1 2 0.5 0.5;
            1.5 1.5 0.5 1.0;
            0.5 1.5 1.0 0.0];

        plot(ax2plot,fld_amphNa1(:,1:2)', fld_amphNa1(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %text labels
        str_amphNa1={'Eckermannite';'Nyboite';'Glaucophane'};

        %position of text labels
        pos_amphNa1=[1.15,0.75;
            1.75,0.75;
            1.65,0.25];


        %text plotting
        text(ax2plot,pos_amphNa1(:,1),pos_amphNa1(:,2),str_amphNa1,'FontSize',FontSize,'HorizontalAlignment','center')

        %axis labels
        xlabel('^{C}(Al + Fe^{3+} + 2Ti) (apfu)')
        ylabel('^{A}(Na + K + 2Ca) (apfu)')

        box on

        % second Na clinoamphibole classification diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_XY_Na2')

        %XY limits
        xlim([0 1])
        ylim([0 1])
        ax2plot.YAxis.TickValues=0:0.1:1;
        ax2plot.XAxis.TickValues=0:0.1:1;
        ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        %field boundaries
        fld_amphNa2=[0.0 1.0 0.5 0.5;
            0.5 0.5 0.0 1.0];


        plot(ax2plot,fld_amphNa2(:,1:2)', fld_amphNa2(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %text labels
        str_amphNa2={'Glaucophane';'Ferro-Glaucophane';'Magnesio-Riebeckite';'Riebeckite'};

        %position of text labels
        pos_amphNa2=[0.25,0.25;
            0.25,0.75;
            0.75,0.25;
            0.75,0.75];

        %text plotting
        text(ax2plot,pos_amphNa2(:,1),pos_amphNa2(:,2),str_amphNa2,'FontSize',FontSize,'HorizontalAlignment','center')

        %axis labels
        xlabel('Fe^{3+}/(Fe^{3+} + Al + Ti)')
        ylabel('Fe^{2+}/(Fe^{2+} + Mg + Mn)')

        box on

        % amphibole Fe3/Fe vs total Fe diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_XY_Fe')

        %XY limits
        xlim([0 5])
        ylim([0 1])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0:0.1:1;
        % ax2plot.XAxis.TickValues=0:0.5:5;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'1.0',[],'2.0',[],'3.0',[],'4.0',[],'5.0'});

        %axis labels
        xlabel('\SigmaFe (apfu)')
        ylabel('Fe^{3+}/\SigmaFe')

        box on

        % Mica Al vs Si diagram for muscovite-celadonite solid solution
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'mica_XY_ph')

        axis equal %makes the spacing of the axes intervals equal

        %XY limits
        xlim([1 3])
        ylim([3 4])

        %sets tick interval and labels
        ax2plot.YAxis.TickValues=3.0:0.2:4.0;
        ax2plot.XAxis.TickValues=1.0:0.2:3.0;
        ax2plot.YTickLabel=({'3.0','3.2','3.4','3.6','3.8','4.0'});
        ax2plot.XTickLabel=({'1.0',[],'1.4',[],'1.8',[],'2.2',[],'2.6',[],'3.0'});

        %muscovite-celadonite solid solution
        fld_mica_ph=[3.0 1.0 3.0 4.0;
            2.78 2.82 3.06 3.14;
            2.58 2.62 3.16 3.24;
            2.38 2.42 3.26 3.34;
            2.18 2.22 3.36 3.44;
            1.98 2.02 3.46 3.54;
            1.78 1.82 3.56 3.64;
            1.58 1.62 3.66 3.74;
            1.38 1.42 3.76 3.84;
            1.18 1.22 3.86 3.94];

        plot(ax2plot,fld_mica_ph(:,1:2)', fld_mica_ph(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot solid solution

        %numerical labels
        str_mica1={'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9'};

        %position of numerical labels
        pos_mica1=[2.85,3.2;
            2.65,3.3;
            2.45,3.4;
            2.25,3.5;
            2.05,3.6;
            1.85,3.7;
            1.65,3.8;
            1.45,3.9;
            1.15,3.8];

        %text labels
        str_mica2={'Celadonite';'Muscovite';'+ Paragonite'};

        %position of text labels
        pos_mica2=[1.0,4.07;
            3.0,2.89;
            3.0,2.82];

        %text plotting
        text(ax2plot,pos_mica1(:,1),pos_mica1(:,2),str_mica1,'FontSize',FontSize,'HorizontalAlignment','center','Rotation',60)
        text(ax2plot,pos_mica2(:,1),pos_mica2(:,2),str_mica2,'FontSize',FontSize*1.1667,'HorizontalAlignment','center')

        %axis labels
        xlabel('Al (apfu)')
        ylabel('Si (apfu)')

        box on

        %allows for not square plot
        options.custom_daspect.Value=[1 1 1];

        % Mica Na vs Si diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'mica_XY_SiNa')

        xlabel('Na (apfu)')
        ylabel('Si (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %X and Y limits
        xlim([0 1])
        ylim([3 4])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=3.0:0.1:4.0;
        % ax2plot.XAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.YTickLabel=({'3.0',[],'3.2',[],'3.4',[],'3.6',[],'3.8',[],'4.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        box on

        % Mica Na vs Si diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'mica_XY_NaKFeMg')

        xlabel('Mg + Fe (apfu)')
        ylabel('Na/(Na + K)')
        axis equal %makes the spacing of the axes intervals equal

        %X and Y limits
        xlim([0 3])
        ylim([0 1])

        box on

        % biotite compositional diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'mica_XY_bio')

        axis equal %makes the spacing of the axes intervals equal

        %XY limits
        xlim([0 1])
        ylim([0 1])

        %odd compositional line positions
        fld_mica_bio1=[0.2 0.0 0.0 0.2;
            0.6 0.0 0.0 0.6;
            1.0 0.0 0.0 1.0;
            0.4 1.0 1.0 0.4;
            0.8 1.0 1.0 0.8;
            0.8 1.0 0.0 0.2;
            0.4 1.0 0.0 0.6;
            0.0 1.0 0.0 1.0;
            0.0 0.6 0.4 1.0;
            0.0 0.2 0.8 1.0];

        %Even compositional line positions
        fld_mica_bio2=[0.4 0.0 0.0 0.4;
            0.8 0.0 0.0 0.8;
            0.2 1.0 1.0 0.2;
            0.6 1.0 1.0 0.6;
            0.6 1.0 0.0 0.4;
            0.2 1.0 0.0 0.8;
            0.0 0.8 0.2 1.0;
            0.0 0.4 0.6 1.0];


        plot(ax2plot,fld_mica_bio1(:,1:2)', fld_mica_bio1(:,3:4)','LineStyle',LineStyle_1,'Color',linecolor_1,'linewidth',linewidth_1,'HandleVisibility','off')
        plot(ax2plot,fld_mica_bio2(:,1:2)', fld_mica_bio2(:,3:4)','LineStyle',LineStyle_2,'Color',linecolor_2,'linewidth',linewidth_1,'HandleVisibility','off')

        %text labels
        str_bio={'Annite';'Siderophyllite';'Eastonite';'Phlogopite'};

        %position of text labels
        pos_bio=[0,-0.07;
            0,1.05;
            1,1.05;
            1,-0.07];

        %text plotting
        text(ax2plot,pos_bio(:,1),pos_bio(:,2),str_bio,'FontSize',FontSize*1.1667,'HorizontalAlignment','center')

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.XAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        %axis labels
        xlabel('X_{Mg}')
        ylabel('Al_{M} (apfu)')

        box on

        % chlorite simple diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'chl_XY_simp')

        ylabel('X_{Mg}')
        xlabel('R^{3+}_{VI} (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %X and Y limits
        xlim([1 3])
        ylim([0 1])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.XAxis.TickValues=1.0:0.25:3.0;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'1.0',[],'1.5',[],'2.0',[],'2.5',[],'3.0'});

        %text labels
        str_chl1={'Chamosite';'Clinochlore';'Sudoite'};

        %position of text labels
        pos_chl1=[0.8,-0.07;
            0.8,1.05;
            2.85,1.05];

        %text plotting
        text(ax2plot,pos_chl1(:,1),pos_chl1(:,2),str_chl1,'FontSize',FontSize*1.1667)

        box on

        % chlorite simple diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'chl_XY_comp')

        ylabel('X_{Mg}')
        xlabel('R^{3+}_{VI} (apfu)')
        axis square %makes the spacing of the axes intervals equal

        fld_chl_comp=[1.0 1.0 0.0 1.0;
            2.0 2.0 0.0 1.0;
            3.0 3.0 0.0 1.0;
            4.0 4.0 0.0 1.0];

        plot(ax2plot,fld_chl_comp(:,1:2)', fld_chl_comp(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([0 5])
        ylim([0 1])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.XAxis.TickValues=0.0:0.5:5.0;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'1.0',[],'2.0',[],'3.0',[],'4.0',[],'5.0'});

        %text labels
        str_chl2={'Chm';'Cln';'Mg-Ame';'Fe-Ame';'Mg-Sud';'Fe-Sud'};
        str_chl3={'Trioctahedral';'Di-trioctahehdral';'Dioctahedral'};

        %position of text labels
        pos_chl2=[1.05,0.04;
            0.85,1.04;
            1.6,1.04;
            2.05,0.04;
            2.6,1.04;
            3.05,0.04];

        %position of text labels
        pos_chl3=[1.5,0.5;
            2.5,0.5;
            3.5,0.5];

        %text plotting
        text(ax2plot,pos_chl2(:,1),pos_chl2(:,2),str_chl2,'FontSize',FontSize*1.1667)
        text(ax2plot,pos_chl3(:,1),pos_chl3(:,2),str_chl3,'FontSize',FontSize*1.1667,'HorizontalAlignment','center','Rotation',90)

        box on

        %Divalent cations vs Si diagram for chlorite
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'chl_XY_R2Si')

        ylabel('Si (apfu)')
        xlabel('R^{2+} (apfu)')
        %axis square %makes the spacing of the axes intervals equal

        %Al contours
        fld_chl_R2Si=[2.0 4.0 3.0 2.0;
            2.0 4.5 3.75 2.5;
            3.0 5.0 4.0 3.0;
            4.5 5.5 4.0 3.5];

        %R2+ contours
        fld_chl_R2Si2=[4.0 6.0 2.0 4.0;
            2.5 4.5 2.0 4.0;
            2.0 3.0 3.0 4.0];

        plot(ax2plot,fld_chl_R2Si(:,1:2)', fld_chl_R2Si(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries
        plot(ax2plot,fld_chl_R2Si2(:,1:2)', fld_chl_R2Si2(:,3:4)','k','LineStyle','--','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([2 6])
        ylim([2 4])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=2.0:0.25:4.0;
        % ax2plot.XAxis.TickValues=2.0:0.5:6.0;
        % ax2plot.YTickLabel=({'2.0',[],'2.5',[],'3.0',[],'3.5',[],'4.0'});
        % ax2plot.XTickLabel=({'2.0',[],'3.0',[],'4.0',[],'5.0',[],'6.0'});

        %text labels
        str_R2Si={'Sud';'Ame';'Clc';'Srp'};
        str_R2Si2={'Al = 1';'Al = 2';'Al = 3';'Al = 4'};
        str_R2Si3={'R^{VI} = 5.0';'R^{VI} = 5.5';'R^{VI} = 6.0'};

        %position of endmember text labels
        pos_R2Si=[2.08,3.00;
            4.1,2.05;
            5.0,2.95;
            5.7,3.95];

        %position of text labels
        pos_R2Si2=[5.02,3.8;
            3.52,3.8;
            3.02,3.3;
            2.52,2.8];

        %position of text labels
        pos_R2Si3=[2.18,3.3;
            3.18,2.8;
            4.18,2.3];

        %text plotting
        text(ax2plot,pos_R2Si(:,1),pos_R2Si(:,2),str_R2Si,'FontSize',FontSize*1.1667)
        text(ax2plot,pos_R2Si2(:,1),pos_R2Si2(:,2),str_R2Si2,'FontSize',FontSize,'HorizontalAlignment','center','Rotation',-38)
        text(ax2plot,pos_R2Si3(:,1),pos_R2Si3(:,2),str_R2Si3,'FontSize',FontSize,'HorizontalAlignment','center','Rotation',57)

        box on


        %Fe vs Si diagram for chlorite
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'chl_XY_XFeSi')

        xlabel('Si (apfu)')
        ylabel('Fe_{total} (apfu)')
        %axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_chl_FeSi=[2.0 4.0 4.0 6.0;
            2.5 3.5 3.6 4.4;
            2.5 2.5 0.0 4.5;
            2.0 2.5 2.0 2.2;
            2.80 3.75 2.32 2.90;
            2.50 3.75 0.90 1.15;
            2.80 2.80 0.00 4.80;
            3.10 3.10 0.00 4.08;
            3.50 3.50 0.00 1.10];

        plot(ax2plot,fld_chl_FeSi(:,1:2)', fld_chl_FeSi(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([2 4])
        ylim([0 6])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.5:6.0;
        % ax2plot.XAxis.TickValues=2.0:0.1:4.0;
        % ax2plot.YTickLabel=({'0',[],'1',[],'2',[],'3',[],'4',[],'5',[],'6'});
        % ax2plot.XTickLabel=({'2.0',[],'2.2',[],'2.4',[],'2.6',[],'2.8',[],'3.0',[],'3.2',[],'3.4',[],'3.6',[],'3.8',[],'4.0'});
        % ax2plot.XTickLabelRotation=0;

        %text labels
        str_FeSi={'pseudothuringite';'corundophyllite';'sheridanite';'ripidiolite';
            'daphnite';'brunsvigite';'pycno-'; 'chlorite';
            'clinochlore';'diabantite';'penninite';'talc-chlorite'};

        %position of endmember text labels
        pos_FeSi=[2.05,2.90;
            2.05,0.90;
            2.51,0.40;
            2.54,2.20;
            2.54,4.10;
            2.81,3.00;
            2.865,1.75;
            2.865,1.50;
            2.820,0.45;
            3.25,1.50;
            3.175,0.45;
            3.60,0.45];

        %text plotting
        text(ax2plot,pos_FeSi(:,1),pos_FeSi(:,2),str_FeSi,'FontSize',FontSize)

        box on

        %allanite-epidote diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'aln_XY')

        ylabel('REE + Y + Th + U (apfu)')
        xlabel('Al_{M} (apfu)')
        %axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_aln=[2.0 1.0 0.0 1.0;
            3.0 2.0 0.0 1.0];

        fld_aln2=[3.0 1.75 0.0 0.25;
            3.0 1.334 0.0 0.666;
            3.0 1.333 0.0 1.0;
            3.0 1.7498 0.0 1.0];

        fld_aln3=[3.0 1.889 0.0 0.111;
            3.0 1.57 0.0 0.43;
            3.0 1.0 0.0 1.0;
            3.0 1.571 0.0 1.0;
            3.0 1.8885 0.0 1.0];


        plot(ax2plot,fld_aln(:,1:2)', fld_aln(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries
        plot(ax2plot,fld_aln2(:,1:2)', fld_aln2(:,3:4)','color',linecolor_2,'LineStyle',LineStyle_2,'linewidth',linewidth_1,'HandleVisibility','off')
        plot(ax2plot,fld_aln3(:,1:2)', fld_aln3(:,3:4)','color',linecolor_1,'LineStyle',LineStyle_1,'linewidth',linewidth_1,'HandleVisibility','off')


        %X and Y limits
        xlim([0.5 3.0])
        ylim([0 1.0])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.XAxis.TickValues=0.5:0.25:3.0;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.5',[],'1.0',[],'1.5',[],'2.0',[],'2.5',[],'3.0'});

        %text labels
        str_aln={'czo';'ep';'aln';'Fe-aln'};
        str_aln2={'0.2';'0.4';'0.6';'0.8'};

        %position of text labels
        pos_aln=[3.05,0.02;
            1.85,0.04;
            1.95,1.04;
            0.9,1.04];

        pos_aln2=[1.7,1.02;
            1.3,1.02;
            1.2,0.65;
            1.62,0.25];

        %text plotting
        text(ax2plot,pos_aln(:,1),pos_aln(:,2),str_aln,'FontSize',FontSize*1.1667)
        text(ax2plot,pos_aln2(:,1),pos_aln2(:,2),str_aln2,'FontSize',FontSize)

        box on

        %titanite Al + Fe3 vs Ti diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ttn_XY_Ti')

        ylabel('Ti (apfu)')
        xlabel('(Al + Fe^{3+})_{VI} (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_ttn1=[0.0 1.0 1.0 0.0];

        plot(ax2plot,fld_ttn1(:,1:2)', fld_ttn1(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([0 1.0])
        ylim([0 1.0])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.XAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        box on

        %titanite Al + Fe3 vs F diagram
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ttn_XY_F')

        ylabel('F (apfu)')
        xlabel('(Al + Fe^{3+})_{VI} (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_ttn1=[0.0 1.0 0.0 1.0;
            0.0 1.0 0.0 0.8;
            0.0 1.0 0.0 0.6;
            0.0 1.0 0.0 0.4;
            0.0 1.0 0.0 0.2];

        plot(ax2plot,fld_ttn1(:,1:2)', fld_ttn1(:,3:4)','color',linecolor_2,'LineStyle',LineStyle_2,'linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([0 1.0])
        ylim([0 1.0])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.XAxis.TickValues=0.0:0.1:1.0;
        % ax2plot.YTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});
        % ax2plot.XTickLabel=({'0.0',[],'0.2',[],'0.4',[],'0.6',[],'0.8',[],'1.0'});

        box on

        %text labels
        str_ttn={'F/(F + OH) = 1.0'};

        %position of text labels
        pos_ttn=[0.5,0.54];

        %text plotting
        text(ax2plot,pos_ttn(:,1),pos_ttn(:,2),str_ttn(:,1),'FontSize',FontSize,'HorizontalAlignment','center','Rotation',45)

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'crd_XY_NaAlSi')

        ylabel('Na + Al (apfu)')
        xlabel('Si (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_crd1=[4.0 5.0 5.0 3.5];

        plot(ax2plot,fld_crd1(:,1:2)', fld_crd1(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([4.0 6.0])
        ylim([3.5 5.0])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=3.5:0.1:5.0;
        % ax2plot.XAxis.TickValues=4.0:0.1:6.0;
        % ax2plot.YTickLabel=({'3.5',[],[],[],[],'4.0',[],[],[],[],'4.5',[],[],[],[],'5.0'});
        % ax2plot.XTickLabel=({'4.0',[],[],[],[],'4.5',[],[],[],[],'5.0',[],[],[],[],'5.5',[],[],[],[],'6.0'});

        %text labels
        str_crd1={'^{Ch}Na + ^{IV}Al = ^{Ch}vac + ^{IV}Si'};

        %position of text labels
        pos_crd1=[4.48,4.2];

        %text plotting
        text(ax2plot,pos_crd1(:,1),pos_crd1(:,2),str_crd1(:,1),'FontSize',FontSize,'HorizontalAlignment','center','Rotation',297)

        box on

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'crd_XY_AlMeNaSi')

        ylabel('Al + Me^{2+} (apfu)')
        xlabel('Na + Si (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_crd2=[4.0 7.0 7.0 4.0];

        plot(ax2plot,fld_crd2(:,1:2)', fld_crd2(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([4.0 7.0])
        ylim([4.0 7.0])

        % %sets tick interval and labels
        % ax2plot.YAxis.TickValues=4.0:0.1:7.0;
        % ax2plot.XAxis.TickValues=4.0:0.1:7.0;
        % ax2plot.YTickLabel=({'4.0',[],[],[],[],'4.5',[],[],[],[],'5.0',[],[],[],[],'5.5',[],[],[],[],'6.0',[],[],[],[],'6.5',[],[],[],[],'7.0'});
        % ax2plot.XTickLabel=({'4.0',[],[],[],[],'4.5',[],[],[],[],'5.0',[],[],[],[],'5.5',[],[],[],[],'6.0',[],[],[],[],'6.5',[],[],[],[],'7.0'});

        %text labels
        str_crd2={'^{Ch}Na + ^{IV}Si = ^{Ch}vac + Me^{2+} + ^{IV}Al'};

        %position of text labels
        pos_crd2=[5.4,5.5];

        %text plotting
        text(ax2plot,pos_crd2(:,1),pos_crd2(:,2),str_crd2(:,1),'FontSize',FontSize,'HorizontalAlignment','center','Rotation',315)

        box on

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'crd_XY_MevsSiAl')
        
        ylabel('Mg + Mn^{2+} + Fe^{2+} (apfu)')
        xlabel('Si + Al (apfu)')
        axis square %makes the spacing of the axes intervals equal

        %field boundaries
        fld_crd3=[8.0 9.0 3.0 2.0];
        fld_crd4=[8.4706 9.0 2.1176 2.0];
        fld_crd5=[9.2571 9.0 1.0286 2.0];

        plot(ax2plot,fld_crd3(:,1:2)', fld_crd3(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries
        plot(ax2plot,fld_crd4(:,1:2)', fld_crd4(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries
        plot(ax2plot,fld_crd5(:,1:2)', fld_crd5(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field boundaries

        %X and Y limits
        xlim([8.6 9.2])
        ylim([1.6 2.2])

        %text labels
        str_crd3={'^{Ch}Na + ^{IV}(Mg,Fe^{2+}) = ^{Ch}vac + ^{IV}Al'};
        str_crd4={'^{Ch}Na + ^{IV}Be = ^{Ch}vac + ^{IV}Al'};
        str_crd5={'^{Ch}Na + ^{IV}Li = ^{Ch}vac + ^{VI}Mg'};

        %position of text labels
        pos_crd3=[8.92,2.1];
        pos_crd4=[8.8,2.055];
        pos_crd5=[9.066,1.8];

        %text plotting
        text(ax2plot,pos_crd3(:,1),pos_crd3(:,2),str_crd3(:,1),'FontSize',FontSize,'HorizontalAlignment','center','Rotation',314)
        text(ax2plot,pos_crd4(:,1),pos_crd4(:,2),str_crd4(:,1),'FontSize',FontSize,'HorizontalAlignment','center','Rotation',347)
        text(ax2plot,pos_crd5(:,1),pos_crd5(:,2),str_crd5(:,1),'FontSize',FontSize,'HorizontalAlignment','center','Rotation',285)

        box on

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'lws_XY_FeTiCrAl')

        xlabel('^{VI}Al (apfu)')
        ylabel('Fe^{3+} + Ti + Cr (apfu)')

        %X and Y limits
        xlim([1.5 2.0])
        ylim([0.0 0.5])

        %field boundaries
        fld_lws=[2.0 1.5 0.0 0.5];

        plot(ax2plot,fld_lws(:,1:2)', fld_lws(:,3:4)','k','linewidth',linewidth_1,'HandleVisibility','off') %plot field bound

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'custom_XY')
        xlim('auto');
        ylim('auto');

        % custom ax labels
        if  isfield(options,'custom') && isfield(options.custom,'Interpreter')
            Interpreter=options.custom.Interpreter;
        else
            Interpreter='tex';
        end

        if   isfield(options,'custom') && isfield(options.custom,'xlabel') && not(isempty(options.custom.xlabel))
            xlabel(options.custom.xlabel,'Interpreter',Interpreter)
        end

        if  isfield(options,'custom') && isfield(options.custom,'ylabel') && not(isempty(options.custom.ylabel))
            ylabel(options.custom.ylabel,'Interpreter',Interpreter)
        end

        if  isfield(options,'custom') && isfield(options.custom,'XScale') && not(isempty(options.custom.XScale))
            ax2plot.XScale=options.custom.XScale;
        end

        if  isfield(options,'custom') && isfield(options.custom,'YScale') && not(isempty(options.custom.YScale))
            ax2plot.YScale=options.custom.YScale;
        end

        % To add in the future:
        % custom ax XGrid YGrid
        % custom ax YMinorGrid XMinorGrid
        % custom ax XAxisLocation YAxisLocation

    end


if not(isempty(options)) && isfield(options,'custom_daspect') && isfield(options.custom_daspect,'Value') && not(isempty(options.custom_daspect.Value)) && numel(options.custom_daspect.Value)==3
    daspect(options.custom_daspect.Value)
else
    axis square
%     ax2plot.DataAspectRatio=[1 1 1];
% ax2plot.PlotBoxAspectRatio=[1 1 1 ];
end
box on
end




