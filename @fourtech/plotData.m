function data = plotData(f, g, h)
%PLOTDATA   Useful data values for plotting a FOURTECH object.
%   DATA = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. The struct DATA contains the following fields:
%
%       xLine: x-coordinates of for plotting smooth curves.
%       yLine: Function values of F at the coordinates stored in xLine.
%       xPoints: x-coordinates of the Chebyshev points used to represent F.
%       yPoints: Function value at the Chebyshev points used to represent F.
%   
%   DATA.xLine and DATA.yLine are used for plotting smooth curves (usually
%   passed to PLOT() with the '-' option). 
%
%   DATA.xPoints and DATA.yPoints contain the (x, F(x)) data at the Chebyshev
%   grid used to represent F, and are used for plots with markes (e.g.
%   PLOT(F,'-o').
%
%   DATA = PLOTDATA(F, G) is similar but for plot calls of the form PLOT(F, G),
%   where both F and G are FOURTECH objects. 
% 
%   DATA = PLOTDATA(F, G, H) is for plots of the form PLOT3(F, G, H). In this
%   instance, DATA also contains fields zLine and zPoints for the data
%   corresponding to H.
%
% See also PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    g = [];
end
if ( nargin < 3 )
    h = [];
end

% Get the number of points: (Oversample the wavelength)
len = max([length(f), length(g), length(h)]);
npts = min(max(501, round(4*pi*len)), chebtech.techPref().maxPoints);

% Initialise the output structure:
data = struct('xLine', [], 'yLine', [], 'xPoints', [], 'yPoints', [], 'yLim', [inf -inf]);
if ( isempty(g) )       
    % PLOT(F):
    
    % Values on oversampled Fourier grid (equally spaced).
    data.xLine = fourierpts(npts);
    tmp = prolong(f, npts);
    data.yLine = tmp.values;

    % Values on the Fourier grid tied to the FOURTECH F:
    data.xPoints = f.points();
    data.yPoints = f.values;

elseif ( isa(g, 'fourtech') )   
    % PLOT(F, G)
    
    % Also return the grid points used.
    % Grid data for f:
    data.fGrid.xLine = fourierpts(size(f.values,1));
    % Use the maximum of the lenghts of f, g and h to match the number of
    % values returned:
    data.fGrid.xPoints = fourierpts(len);
    
    % Grid data for g:
    data.gGrid.xLine = fourierpts(npts);
    % Use the maximum of the lenghts of f, g and h to match the number of
    % values returned:    
    data.gGrid.xPoints = fourierpts(len);
    
    % Values on oversampled Fourier grid (equally spaced)
    data.xLine = get(prolong(f, npts), 'values');
    data.yLine = get(prolong(g, npts), 'values');

    % Values on the largest uniform grid tied to the FOURTECH objects F and G:
    data.xPoints = get(prolong(f, len), 'values');
    data.yPoints = get(prolong(g, len), 'values');
    
    if ( isa(h, 'fourtech') )
        % PLOT3(F, G, H)
        
        % Grid data for h:
        data.hGrid.xLine = fourierpts(npts);
        % Use the maximum of the lenghts of f, g and h to match the number of
        % values returned:    
        data.hGrid.xPoints = fourierpts(len);
        
         % Values on oversampled uniform grid:
        data.zLine = get(prolong(h, npts), 'values');
        data.zPoints = get(prolong(h, len), 'values');  
    end
    
else
    error('CHEBFUN:FOURTECH:plotdata:DataType', ...
        'Invalid data types.');
end

end