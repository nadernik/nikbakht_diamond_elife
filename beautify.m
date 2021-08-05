% Function beautify rendered the plots and figures using OpenGL as a result
% you will get good looking graphics for publication or presentation. It
% also changes the axis properties,
% Applications of the function: beautify();
% AUTHOR: NADER NIKBAKHT, SISSA 2013
function beautify(gca,print_ready)
if nargin <2
    print_ready = 0;
end
% opengl software;
% opengl('OpenGLLineSmoothingBug',1); %Solve OpenGL bug with LineSmoothing property
set( gca                       , ...
    'FontName'   , 'Arial' );
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'     ,   'LineWidth'   , 1);
set(gcf, 'color', [1 1 1]);
if print_ready
%     set(gcf, 'Renderer','Painters');
    set(gca, 'fontsize',14);
end
end