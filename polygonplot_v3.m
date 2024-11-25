% ----------------------------------------------------------------------- %
function polygonplot_v3(data_bound,opt_axes,opt_lines,opt_area)

N = size(data_bound,2);
nlines = 1;
% Defaults
if(nargin<4)
    opt_area.err = 'std';
    opt_area.FaceAlpha = 0.5;
end
if(nargin<3)
    opt_lines.LineWidth = 2;
    opt_lines.LineStyle = '-';
    opt_lines.Marker    = 'none';
end
if(nargin<2)
    opt_axes = [];
    opt_axes.Background = 'none';
end
if(nargin<1)
    error('Not enough parameters');
end

% Error detection

if(~isfield(opt_lines,'Color'))     % Color properties
    line_color = [];
else
    if(size(opt_lines.Color,1)~=nlines)
        error('Number of colors must be equal to the number of lines to plot.');
    else
        if(~isfield(opt_area,'Color') && M~=0)
            error('Please, specify also the colors of the shaded areas.');
        else
            if(size(opt_area.Color,1)~=nlines)
                error('Number of shaded areas must be equal to the number of lines to plot.');
            else
                line_color = opt_lines.Color;
                opt_lines = rmfield(opt_lines,'Color');
            end
        end
    end
end
if(~isfield(opt_lines,'Labels'))    % Text labels
    labels = false;
else
    labels = opt_lines.Labels;
    opt_lines = rmfield(opt_lines,'Labels');
end
if(~isfield(opt_area,'Color'))
    opt_area.Color = 0.5.*ones(nlines,3);
end
if(~isfield(opt_area,'FaceAlpha'))
    opt_area.FaceAlpha = 0.5;
end
if(~isfield(opt_lines,'Legend'))
    leg = [];
else
    leg = opt_lines.Legend;
    opt_lines = rmfield(opt_lines,'Legend');
end
if(isfield(opt_axes,'Labels'))
    if(length(opt_axes.Labels)~=N)
        error('You must provide N axis labels.');
    end
end
set(gcf,'Color',[1 1 1]);


% Computing the mean and standard deviation of the data matrix

% Plots
m_down = [data_bound(2,:)'; data_bound(2,1)];
m_up   = [data_bound(3,:)'; data_bound(3,1)];

R  = [data_bound(1,:)'; data_bound(1,1)];

% TH = (2*pi/N)*((N:-1:0)'*ones(1,nlines));
TH = (2*pi/N)*((0:1:N)');

[X,Y] = pol2cart(TH, R);

data_min = min([m_down(:);m_up(:)]);
data_max = max([m_down(:);m_up(:)]);

    opt_axes.Ticks = linspace(data_min,data_max,N+1);

r_iso = (opt_axes.Ticks(:))*ones(1,N);

th_iso = (2*pi/N)*(ones(length(opt_axes.Ticks),1)*(N:-1:1));

[x,y] = pol2cart(th_iso, r_iso);
% h_iso = line([x,x(:,1)]',[y,y(:,1)]','LineWidth',0.5,'Color',0.85.*ones(1,3));
% for iso_id = 1:1:N   % Exclude axes from legend
%     set(get(get(h_iso(iso_id),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% end
hold on;

% Compute and plot the axes depending on N
th_jump = (2*pi/N)*(ones(2,1)*(N:-1:1));

radii   = [zeros(1,N); max(abs([data_max,data_min])).*ones(1,N)];

[x,y]   = pol2cart(th_jump, radii);
h_axes  = line(x,y,'LineWidth',1,'Color','k');
for ax_id = 1:1:N   % Exclude axes from legend
    set(get(get(h_axes(ax_id),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold on;

% Display axis Ticks
%if(~isfield(opt_axes,'Background')), opt_axes.Background = 'none'; end
%loc = (opt_axes.Ticks(:)-data_min)/(data_max-data_min);
loc = (opt_axes.Ticks(:));

% Display axis labels
if(isfield(opt_axes,'Labels'))
    th_lbl = (2*pi/N)*(0:1:N-1);
    r_lbl  = max(abs([data_max,data_min])+0.01).*ones(1,N);
    [xlbl,ylbl] = pol2cart(th_lbl,r_lbl);
    
%     for lid = 1:1:N
%         text(xlbl(lid),ylbl(lid),opt_axes.Labels{lid},'fontweight','bold');
%     end
end

% Plot the data
opt_lines_str = adapt_options(opt_lines);
          % Shaded area
        % Computing the mean and standard deviation of the data matrix
        
        % Plots
        
            
         [xa,ya] = pol2cart([TH; fliplr(TH)],[m_down; fliplr(m_up)]);

            pat = fill(xa,ya,opt_area.Color);
            hold on;
            set(pat, 'EdgeColor', 'none');
            set(pat, 'FaceAlpha', opt_area.FaceAlpha);

    if(isempty(line_color))
        h_leg = plot(X,Y,opt_lines_str{:}); hold on;
    else
            h_leg = plot(X,Y,opt_lines_str{:},'Color',line_color);
    end

    axis equal;
    axis off;
hold off;
xlim([-0.9 0.9]);
ylim([-0.9 0.9]);
end


% Function adapt_options adapts the input options struct in a dynamic
% cell-object that can directly pass Name-Value pair arguments through a function
function optList = adapt_options(optStruct)
optList = {};
for optField = fieldnames(optStruct)'
    optList{end+1} = char(optField);                % Name parameter
    optList{end+1} = optStruct.(char(optField));    % Value parameter
end
end