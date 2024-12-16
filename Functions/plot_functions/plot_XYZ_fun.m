function [ax2plot,options_definition]=plot_XYZ_fun(ax2plot,options)
if not(exist('ax2plot','var'))
    ax2plot=[];
end
if not(exist('options','var'))
    ax2plot=[];
end
%% options definition
options_definition.XYZplot.Value=true;

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
    'ol_CaO_Al2O3_NiO_XYZ';
    'grt_CaO_Al2O3_MgO_XYZ';
    'mica_xMg_Al_Si';
    'custom_XYZ';
    };

options_definition.type.description={'none';
    'olivine: CaO Al2O3 NiO XYZ';
    'garnet: CaO Al2O3 MgO XYZ';
    'mica: xMg Al Si';
    'custom: XYZ';
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
rotate_state=true;

if options.XYZplot.Value==true
    %% Mineral specific options
    if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'ol_CaO_Al2O3_NiO_XYZ')

        xlim('auto');
        ylim('auto');
        zlim('auto');


        %plot axis labels
        xlabel('CaO [wt %]')
        ylabel('Al2O3 [wt %]')
        zlabel('NiO [wt %]')
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'grt_CaO_Al2O3_MgO_XYZ')
        xlim('auto');
        ylim('auto');
        zlim('auto');


        %plot axis labels
        xlabel('Ca apfu')
        ylabel('Al apfu')
        zlabel('Mg apfu')
    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'mica_xMg_Al_Si')
        xlim([0 1]);
        ylim([0 1]);
        zlim('auto');

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

        %axis labels
        xlabel('X_{Mg}')
        ylabel('Al_{M} (apfu)')
        zlabel('Si (apfu)')

        box on


        %plot axis labels
         

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'custom_XYZ')
        xlim('auto');
        ylim('auto');
        zlim('auto');

        if not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'Interpreter')
            Interpreter=options.custom.Interpreter;
        else
            Interpreter='tex';
        end

        if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'xlabel') && not(isempty(options.custom.xlabel))
            xlabel(options.custom.xlabel,'Interpreter',Interpreter)
        end
        if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'ylabel') && not(isempty(options.custom.ylabel))
            ylabel(options.custom.ylabel,'Interpreter',Interpreter)
        end
        if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'zlabel') && not(isempty(options.custom.zlabel))
            zlabel(options.custom.zlabel,'Interpreter',Interpreter)
        else
            zlabel('')
            rotate_state=false;
        end

        if  isfield(options,'custom') && isfield(options.custom,'XScale') && not(isempty(options.custom.XScale))
            ax2plot.XScale=options.custom.XScale;
        end

        if  isfield(options,'custom') && isfield(options.custom,'YScale') && not(isempty(options.custom.YScale))
            ax2plot.YScale=options.custom.YScale;
        end


        if  isfield(options,'custom') && isfield(options.custom,'ZScale') && not(isempty(options.custom.ZScale))
            ax2plot.ZScale=options.custom.ZScale;
        end

    else

        xlabel('X')
        ylabel('Y')
        zlabel('Z')
    end
end
%        pbaspect(ax2plot,[1 1 1])
grid on
axis square
box on
view(3)

if rotate_state==true
    rotate3d

end
end




