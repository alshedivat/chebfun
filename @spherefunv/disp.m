function disp(F)
%DISP   Display a CHEBFUN2V.
%
% See also DISPLAY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

loose = strcmp(get(0,'FormatSpacing'),'loose');

% Compact version
if ( isempty( F ) )
    fprintf('empty spherefunv\n')
    return
end

if ( F.isTransposed )
    tString = 'Row vector';
else
    tString = 'Column vector';
end

disp(['   spherefunv object ' '(' tString ') containing' ])
if ( loose )
    fprintf('\n');
end

% Display its two SPHERFUN halves.
for j = 1:3
    disp( F.components{j} );
end

end