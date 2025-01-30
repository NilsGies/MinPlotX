function [ax2plot,options_definition]=plot_binary_fun(ax2plot,options)
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
    'ctd_XMg'
    'Ep-Czo'
    'Ol_XFo'
    };

options_definition.type.description={'none';
    'chloritoid: XMg'
    'epidote: Epidote-Clinozoisite Binary'
    'olivine: binary'
    };

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

if (isfield(options,'default') && options.default==true)
    type=options.type.Value;
    options=options_definition;
    options.type.Value=type;
end

if options.XYplot.Value==true

    %% Mineral specific options
    x_label_str={};

    if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ctd_XMg')
        x_label_str='X_{Mg}';
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'Ep-Czo')
        x_label_str='XEp';
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'Ol_XFo')
        x_label_str='fo [mol %]';
    end


 %   ax2plot.OuterPosition=[0 0.5 1 0.35];
    %changes the location of the X axis to the origin
    ax2plot.XAxisLocation='origin';
    ax2plot.YLim=[-0.5 0.5];
    ax2plot.XLim=[0 1];
    %remove Y axis
    set(ax2plot,'ytick',[])
    ax2plot.YAxis.Visible = 'off';
    %make the x axis ticks stick out in both directions
    ax2plot.TickDir='both';
    %change the position of the X axis
        xlabel(ax2plot,x_label_str)

ax2plot.XLabel.Position=[0.5 -0.25 1];
ax2plot.XLabel.HorizontalAlignment='center';
ax2plot.XLabel.FontSize=12;
ax2plot.XLabel.FontWeight='bold';
end


