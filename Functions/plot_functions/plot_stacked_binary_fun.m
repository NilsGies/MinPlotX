function [ax2plot,options_definition]=plot_stacked_binary_fun(ax2plot,options)
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
    'ol_NiO_Mn_Ca_Mg_XXXY'
    'ol_large_XY'
    'ol_test_XY'
    'ol_logtest1_XY'
    };



options_definition.type.description={'none';
    'olivine: NiO|MnO|CaO - Mg# diagram';
    'olivine: large diagram';
    'olivine: test';
    'olivine: log test';
    };
options_definition.type.n_axes=[
    3
    8
    4
    4
    ];

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

if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_NiO_Mn_Ca_Mg_XXXY')
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
    fig_resize=[1 1 1 1];
elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_logtest1_XY')
   % ax2plot=tiledlayout(3,1,'padding','compact','TileSpacing','tight');
    ax2plot=tiledlayout(2,2,'TileSpacing','tight');
    y_label_str={'log(FeO/MgO)','log(FeO/MgO)','log(NiO/MnO)','log(NiO/Cr_2O_3)'};
    x_label_str={'log(SiO_2/MnO)','log(SiO_2/NiO)','log(MgO/FeO)','log(MgO/FeO)'};
    fig_resize=[1 1 1 1];

end

for n=1:numel(x_label_str)
    ax2plot(n)=nexttile;
    ylabel(y_label_str{n})
    xlabel(x_label_str{n})
    ax2plot(n).FontSize=FontSize;
end
hFig.Position=hFig.Position.*fig_resize;

end


