function [ax2plot,options_definition]=plot_multi_XY_fun(ax2plot,options)
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
    'amph_multiplot';
    'ol_NiO_Mn_Ca_Mg_XXXY';
    %    'ol_large_XY';
    %     'ol_test_XY';
    %     'ol_logtest1_XY'
    'ol_NiO_CaO_Mg_XXXY';
    %    'ol_H2O_large';
    'clinopyroxene_xMg_NaO2_Cr2O3_Al2O3';
    'clinopyroxene_xMg_Endmember';
   'clinopyroxene_xMg_Endmember_HPT';
    };



options_definition.type.description={'none';
    'amphibole: Multi'
    'olivine: NiO|MnO|CaO - Mg# diagram';
    % 'olivine: large diagram';
    %     'olivine: test';
    %     'olivine: log test';
    'olivine: NiO|CaO - Mg# diagram';
    %'olivine: Fo|CaO Cr|Al Ti|Li  - H2O';
    'clinopyroxene: xMg NaO2 Cr2O3 Al2O3';
    'clinopyroxene: xMg Endmember';
    'clinopyroxene: xMg Endmember HPT';
    };
options_definition.type.n_axes=[
    3
    8
    4
    4
    2
    6];


if not(exist('options','var')) || isempty(options)
    if not(exist('ax2plot','var'))
        ax2plot=[];
    end
    return
end

%if not(exist('ax2plot','var')) || isempty(ax2plot)
    hFig=figure;
%end
%hold on


fig_resize=[1 1 1 1];
if (isfield(options,'default') && options.default==true)
    type=options.type.Value;
    options=options_definition;
    options.type.Value=type;
end

if options.XYplot.Value==true

    %% Mineral specific options
x_label_str={};
y_label_str={};
if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'amph_multiplot')
    % ax2plot=tiledlayout(3,1,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(3,2,'TileSpacing','tight');
    
    options_sub=options;
    options_sub.FontSize.Value=FontSize/1.5;
    %%%%% 1
    options_sub.type.Value='amph_XY_Ca1';
    ax2plot(1)=nexttile;
    [ax2plot(1)]=plot_XY_fun(ax2plot(1),options_sub);
    %%%%% 2
    options_sub.type.Value='amph_XY_Ca2';
    ax2plot(2)=nexttile;
    [ax2plot(2)]=plot_XY_fun(ax2plot(2),options_sub);
    %%%%% 3
    options_sub.type.Value='amph_XY_Na1';
    ax2plot(3)=nexttile;
    [ax2plot(3)]=plot_XY_fun(ax2plot(3),options_sub);
    %%%%% 4
    options_sub.type.Value='amph_XY_NaCa';
    ax2plot(4)=nexttile;
    [ax2plot(4)]=plot_XY_fun(ax2plot(4),options_sub);
    %%%%% 6
    options_sub.type.Value='amph_XY_Na2';
    ax2plot(5)=nexttile;
    [ax2plot(5)]=plot_XY_fun(ax2plot(5),options_sub);
    %%%%% 6
    options_sub.type.Value='amph_XY_Fe';
    ax2plot(6)=nexttile;
    [ax2plot(6)]=plot_XY_fun(ax2plot(6),options_sub);


    fig_resize=[1 1 0.5 1];

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_NiO_Mn_Ca_Mg_XXXY')
   % ax2plot=tiledlayout(3,1,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(3,1,'TileSpacing','tight');

    x_label_str={'','','Mg#'};
    y_label_str={'NiO [wt. %]','MnO [wt. %]','CaO [wt. %]'};
    fig_resize=[1 1 0.5 1];

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_large_XY')
%    ax2plot=tiledlayout(4,2,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(4,2,'TileSpacing','tight');
    x_label_str={'Mg#', 'Mg#', 'Mg#', 'Mg#','SiO_2 [wt. %] ','SiO_2 [wt. %]','SiO_2 [wt. %]','Al_2O_3 [wt. %]'};
    y_label_str={'SiO_2 [wt. %]','TiO_2 [wt. %]','Cr_2O_3 [wt. %]','Al_2O_3 [wt. %]','Na_2O [wt. %]','CaO [wt. %]','K_2O [wt. %]','K_2O [wt. %]'};
    fig_resize=[1 1 0.5 1];
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_test_XY')
   % ax2plot=tiledlayout(3,1,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(2,2,'TileSpacing','tight');
    y_label_str={'FeO/MgO','FeO/MgO','NiO/MnO','NiO/Cr_2O_3'};
    x_label_str={'SiO_2/MnO','SiO_2/NiO','MgO/FeO','MgO/FeO'};
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_logtest1_XY')
   % ax2plot=tiledlayout(3,1,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(2,2,'TileSpacing','tight');
    y_label_str={'log(FeO/MgO)','log(FeO/MgO)','log(NiO/MnO)','log(NiO/Cr_2O_3)'};
    x_label_str={'log(SiO_2/MnO)','log(SiO_2/NiO)','log(MgO/FeO)','log(MgO/FeO)'};

    

elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_NiO_CaO_Mg_XXXY')
   % ax2plot=tiledlayout(3,1,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(2,1,'TileSpacing','tight');

    x_label_str={'','Mg#'};
    y_label_str={'NiO [wt. %]','CaO [wt. %]'};
    fig_resize=[1 1 0.5 1];

    % elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_H2O_large')
    %     %    ax2plot=tiledlayout(4,2,'padding','compact','TileSpacing','tight');
    %     ax2plot=tiledlayout(3,2,'TileSpacing','tight');
    %     x_label_str={'','','','','H_2O [µg/g]','H_2O [µg/g]'};
    %     y_label_str={'Mg#', 'Ca [µg/g]','Cr [µg/g]','Al [µg/g]','Ti [µg/g]','Li [µg/g]'};
    %     fig_resize=[1 1 0.5 1];
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_H2O_large')
    %    ax2plot=tiledlayout(4,2,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(3,2,'TileSpacing','tight');
    x_label_str={'','','','','H_2O [µg/g]','H_2O [µg/g]'};
    y_label_str={'XFo', 'Ca [µg/g]','Cr [µg/g]','Al [µg/g]','Ti [µg/g]','Li [µg/g]'};
    fig_resize=[1 1 0.5 1];
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'clinopyroxene_xMg_NaO2_Cr2O3_Al2O3')
 ax2plot=tiledlayout(3,1,'TileSpacing','tight');
    x_label_str={'','','X_{Mg}'};
    y_label_str={'NaO_2 [wt.%]', 'Cr_2O_3 [wt.%]','Al_2O_3 [wt.%]'};
        fig_resize=[1 1 0.5 1];
        elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'clinopyroxene_xMg_Endmember_HPT')
 ax2plot=tiledlayout(2,5,'TileSpacing','tight');
    x_label_str={'X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}'};
    y_label_str={'X_{jd}','X_{aeg}','X_{dihd}','X_{Cats}','X_{kos}','X_{Kkos}','X_{Kjd}','X_{Ti cpx}','X_{Caes}','X_{opx}'};

        fig_resize=[1 1 2.75 1];
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'clinopyroxene_xMg_Endmember')
 ax2plot=tiledlayout(2,4,'TileSpacing','tight');
    x_label_str={'X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}','X_{Mg}'};
    y_label_str={'X_{wo}','X_{fs}','X_{en}','X_{jd}','X_{aeg}','X_{kos}','X_{quad}'};

        fig_resize=[1 1 2.65 1];
end

for n=1:numel(x_label_str)
    ax2plot(n)=nexttile;
    ylabel(y_label_str{n})
    xlabel(x_label_str{n})
    ax2plot(n).FontSize=FontSize;
     axis square
        box on
%    pbaspect(ax2plot,[1 1 1])
end
hFig.Position=hFig.Position.*fig_resize;

end


