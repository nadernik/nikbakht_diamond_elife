% Function beautify renders plots and figures using OpenGL for high-quality graphics
% suitable for publication or presentation. It also modifies axis properties for
% consistent, professional appearance.
%
% Applications of the function: beautify();
% 
% INPUTS:
%   gca - handle to current axes (optional, defaults to current axes)
%   print_ready - boolean flag for print-ready settings (optional, defaults to 0)
%
% OUTPUTS:
%   None - modifies the current figure and axes properties
%
% AUTHOR: NADER NIKBAKHT, SISSA 2013
function beautify(gca,print_ready)

% Set default value for print_ready parameter if not provided
if nargin <2
    print_ready = 0;
end

% OpenGL rendering settings (commented out but available for debugging)
% opengl software;
% opengl('OpenGLLineSmoothingBug',1); %Solve OpenGL bug with LineSmoothing property

% Set font properties for the current axes
set( gca                       , ...
    'FontName'   , 'Arial' );

% Configure axis appearance properties for professional look
set(gca, ...
    'Box'         , 'on'     , ...    % Display box around plot
    'TickDir'     , 'out'     , ...   % Tick marks point outward
    'TickLength'  , [.02 .02] , ...   % Length of tick marks
    'XMinorTick'  , 'off'      , ...  % Disable minor x-axis ticks
    'YMinorTick'  , 'off'      , ...  % Disable minor y-axis ticks
    'YGrid'       , 'off'     , ...   % Disable y-axis grid
    'LineWidth'   , 1);               % Set line width

% Set figure background to white for clean appearance
set(gcf, 'color', [1 1 1]);

% Apply print-ready settings if requested
if print_ready
    % Set renderer to Painters for vector graphics (commented out)
    % set(gcf, 'Renderer','Painters');
    
    % Increase font size for better readability in publications
    set(gca, 'fontsize',14);
end

end