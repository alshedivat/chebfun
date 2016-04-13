function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display coefficients graphically.
%   PLOTCOEFFS(F) plots the absolute values of the coefficients underlying the
%   representation of a CHEBFUN F on a semilogy scale. If F is an array-valued
%   CHEBFUN or has breakpoints, then a curve is plotted for each FUN of each
%   component (column) of F.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, as in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale.
%
%   H = PLOTCOEFFS(F) returns a column vector of handles to lineseries objects.
%
%   As of Chebfun version 5, the coefficients plotted by 'plotcoeffs' are either
%   Chebyshev coefficients for Chebyshev-based chebfuns or Fourier coefficients
%   for periodic chebfuns ('trigfuns'). 
%
% See also CHEBFUN/PLOT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% We can only plot the coefficients of one CHEBFUN at a time:
if ( any(cellfun(@(f) isa(f, 'chebfun'), varargin)) )
    error('CHEBFUN:CHEBFUN:plotcoeffs:multipleChebfuns', ...
        'Calls of the form PLOTCOEFFS(F, ''b'',  G, ''r'') are not supported.');
end

% Store the hold state of the current axis:
holdState = ishold;

% Grab some colours:
col = [];
if ( nargin > 1 && ischar(varargin{1}) && numel(varargin{1}) < 4 )
    % (Note. LineStyle definitions have at most 4 characters. Other plotting
    % parameters have more than characters.)
    col = regexp(varargin{1}, '[bgrcmykw]', 'match');
    if ( numel(col) > 1 )
        error('CHEBFUN:plotcoeffs:color', ...
            'Error in color/linetype argument.');
    elseif ( ~isempty(col) )
        col = col{:};
%         varargin(1) = []; % Don't remove as we need marker information.
    end
end
if ( isempty(col) )
    colIdx = find(cellfun(@(arg) all(ischar(arg)) && ...
        any(strfind(arg, 'color')), varargin));
    if ( colIdx )
        col = varargin{colIdx+1};
        varargin(colIdx+(0:1)) = [];
    else
        col = get(0, 'DefaultAxesColorOrder');
        col = repmat(col, ceil(numColumns(f)/7), 1);
    end
end

% Deal with 'LogLog' and input:
doLogLog = any(cellfun(@(s) strcmpi(s, 'loglog'), varargin));

% In the case where hold is on we need to check if the
% previous plot adjusted for 0 coefficients. This can
% be done by checking the y tick marks on the the current
% axis.
if ( holdState )

    % grab the tick labels
    ylabel = get(gca,'yticklabel');
  
    % if ylabel is 'Zero' revert ydata
    if ( ~isempty(ylabel) && strcmpi(ylabel{1},'zero') )

        % Pull out all the line data from the plot
        axesObjs = get(gcf,'Children');
        dataObjs = get(axesObjs, 'Children');
        lineObjs = findobj(dataObjs, 'type', 'line');
        ydata = get(lineObjs, 'YData');
        if ( isa(ydata,'double') )
            ydata = {ydata};
        end
    
        % Get the minimum
        Ymin = min(cellfun(@min,ydata));
    
        % revert smallest entries of ydata to 0
        for k = 1:numel(ydata)
            q = ydata{k}; 
            q(q==Ymin) = 0; 
            ydata{k} = q; 
        end
    
        % update plot with new ydata
        for k = 1:numel(lineObjs)
            set(lineObjs(k), 'YData', ydata{k});
        end
        
        % reset the y ticks and labels
        set(gca,'YTickLabelMode','auto');

    end

end


% Convert to a cell array for easy handling:
f = mat2cell(f);

% Initialise the output:
h = cell(1, numel(f));

% Loop over the columns in a quasimatrix:
for k = 1:numel(f)
    % Colours for this columns:
    if ( ischar(col) )
        colk = col;
    else
        colk = col(k,:);
    end
    % Call the column version:
    h{k} = columnPlotCoeffs(f{k}, colk, varargin{:});
end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Pull out all the line data from the plot
axesObjs = get(gcf,'Children');
dataObjs = get(axesObjs, 'Children');
lineObjs = findobj(dataObjs, 'type', 'line');
ydata = get(lineObjs, 'YData');
if ( isa(ydata,'double') )
    ydata = {ydata};
end

% Get the minimum
Ymin = min(cellfun(@min,ydata));

% If Ymin == 0 we take action
if ( Ymin == 0 )
    % get the maximum
    Ymax = min(cellfun(@max,ydata));

    % get the smallest nonzero number
    Ymin = Ymax;
    for k = 1:numel(ydata)
        q = ydata{k}; qmin = min(q(q>0)); 
        if ( qmin < Ymin )
            Ymin = qmin;
        end
    end

    % set Ymin to be smallest tick value
    ytick = get(gca,'ytick'); 
    Ymin = ytick(1);
    
    % update the zero entries of ydata
    for k = 1:numel(ydata)
        q = ydata{k}; 
        q(q==0) = Ymin; 
        ydata{k} = q; 
    end

    % update plot with new ydata
    for k = 1:numel(lineObjs)
        set(lineObjs(k), 'YData', ydata{k});
    end
    
    % update ylabel
    ylabel = get(gca,'yticklabel');
    ylabel{1} = 'Zero'; 
    set(gca,'yticklabel',ylabel);

end

% Set xScale to logarithmic if requested:
if ( doLogLog )
    set(gca, 'xScale', 'Log');
end

% Give an output if one was requested:
if ( nargout > 0 )
    try 
        h = cell2mat(h);
    catch
        % shrug
    end
    varargout = {h};
end

end

function h = columnPlotCoeffs(f, col, varargin)
numFuns = numel(f.funs);

% Initialise handle storage:
h = cell(1,numFuns);

% Call plotcoeffs at the tech level:
for j = 1:numFuns
    h{j} = plotcoeffs(f.funs{j}, varargin{:}, 'color', col); hold on
end

end
