function plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridPoints)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP3.
%   PLOTMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
xx = gridPoints{1};
yy = gridPoints{2};
zz = gridPoints{3};
N = size(xx, 1);
nVars = S.numVars;
dataToPlot = plotOptions{2};
dom = S.domain;
tt = trigpts(N, dom);
xxplot = [xx, 2*xx(:,end,:) - xx(:,end-1,:)];
xxplot =  [xxplot; xxplot(1,:,:)];
xxplot = cat(3, xxplot, xxplot(:,:,1));
yyplot = [yy; 2*yy(end,:,:) - yy(end-1,:,:)];
yyplot = [yyplot, yyplot(:,1,:)];
yyplot = cat(3, yyplot, yyplot(:,:,1));
zzplot = cat(3, zz, 2*zz(:,:,end) - zz(:,:,end-1));
zzplot = [zzplot; zzplot(1,:,:)];
zzplot = [zzplot, zzplot(:,1,:)];

% Loop over the variables:
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vplot = dataToPlot(v(idx:idx+N-1,:,:));
    vplot = [vplot, vplot(:,1,:)]; %#ok<*AGROW>
    vplot = [vplot; vplot(1,:,:)];
    vplot = cat(3, vplot, vplot(:,:,1));
    
    % Loop over the surfaces:
    subplot(nVars, 1, k)
    for l = 1:length(p{k})
        pkl = p{k}(l);
        Sx = ( pkl.XData(1,1) == pkl.XData(1,2) ) && ...
            ( pkl.XData(1,1) == pkl.XData(2,1) ) && ...
            ( pkl.XData(1,1) == pkl.XData(2,2) );
        Sy = ( pkl.YData(1,1) == pkl.YData(1,2) ) && ...
            ( pkl.YData(1,1) == pkl.YData(2,1) ) && ...
            ( pkl.YData(1,1) == pkl.YData(2,2) );
        Sz = ( pkl.ZData(1,1) == pkl.ZData(1,2) ) && ...
            ( pkl.ZData(1,1) == pkl.ZData(2,1) ) && ...
            ( pkl.ZData(1,1) == pkl.ZData(2,2) );
        if ( Sx == 1 )
            pos = p{k}(l).XData;
            pos = pos(1);
            [~, id] = min(abs(tt-pos));
            set(p{k}(l), 'xdata', squeeze(xxplot(:,id,:)))
            set(p{k}(l), 'ydata', squeeze(yyplot(:,id,:)))
            set(p{k}(l), 'zdata', squeeze(zzplot(:,id,:)))
            set(p{k}(l), 'cdata', squeeze(vplot(:,id,:)))
        elseif ( Sy == 1 )
            pos = p{k}(l).YData;
            pos = pos(1);
            [~, id] = min(abs(tt-pos));
            set(p{k}(l), 'xdata', squeeze(xxplot(id,:,:)))
            set(p{k}(l), 'ydata', squeeze(yyplot(id,:,:)))
            set(p{k}(l), 'zdata', squeeze(zzplot(id,:,:)))
            set(p{k}(l), 'cdata', squeeze(vplot(id,:,:)))
        elseif ( Sz == 1 )
            pos = p{k}(l).ZData;
            pos = pos(1);
            [~, id] = min(abs(tt-pos));
            set(p{k}(l), 'xdata', squeeze(xxplot(:,:,id)))
            set(p{k}(l), 'ydata', squeeze(yyplot(:,:,id)))
            set(p{k}(l), 'zdata', squeeze(zzplot(:,:,id)))
            set(p{k}(l), 'cdata', squeeze(vplot(:,:,id)))
        end
    end

    % Update title:
    if ( k == 1 )
        title(sprintf('N = %i, dt = %1.1e, t = %.4f', N, dt, t))
    end
    drawnow
end

end