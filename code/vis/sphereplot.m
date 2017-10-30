function sphereplot( X, type )

if ( nargin == 1 )
    type = 'coeffs';
end

switch type
    case 'coeffs'
        ff = coeffs2vals( X );
    case 'vals'
        ff = X;
    otherwise
        error('Unknown input type.');
end

if ( ~isreal(ff) )
    ff = real( ff );
end

[n1, n2, n3] = size( X );

r   = chebpts( n1 );
th  = [pi*trigpts( n2 ); pi];
lam = [pi*trigpts( n3 ); pi] - pi/2;

% Remove doubled-up data
r   = r(floor(n1/2)+1:end);
lam = lam(floor(n3/2)+1:end);
ff  = ff(floor(n1/2)+1:end,:,floor(n3/2)+1:end);
ff(:,end+1,:) = ff(:,1,:);
ff(:,:,end+1) = ff(:,:,1);

[tt, rr, ll] = meshgrid(th, r, lam);

% Slices in the sphere to plot
rslice = 0;
tslice = squeeze(tt(1,[1 floor(n2/4)+1 floor(n2/2)+1 floor(3*n2/4)+1],1));
lslice = 0;

hslicer = slice(tt,rr,ll,ff,tslice,rslice,lslice);
hold on
for j = 1:numel(hslicer)
    h = hslicer(j);
    [xs,ys,zs] = sph2cart(h.XData,h.ZData,h.YData);
    surf(xs,ys,zs,h.CData,'EdgeColor','none','FaceColor','Interp');
end
delete(hslicer)
axis([-1 1 -1 1 -1 1])
daspect([1 1 1])
hold off

set(gca, 'Position', [0 0 1 1], 'CameraViewAngleMode', 'Manual')
colorbar('FontSize', 16, 'Position', [0.84 0.18 0.04 0.64])
axis off

end

function VALS = coeffs2vals( CFS )
% Convert to Chebyshev--Fourier--Fourier values

[n1, n2, n3] = size( CFS );
VALS = CFS;
for k = 1:n3
    VALS(:,:,k) = chebtech2.coeffs2vals( VALS(:,:,k) );
    VALS(:,:,k) = trigtech.coeffs2vals( VALS(:,:,k).' ).';
end
for j = 1:n2
    vj = reshape( VALS(:,j,:), n1, n3 );
    vj = trigtech.coeffs2vals( vj.' ).';
    VALS(:,j,:) = reshape( vj, n1, 1, n3 );
end

end
