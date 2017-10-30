function cylinderplot( X, type )

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

r  = chebpts( n1 );
th = [pi*trigpts( n2 ); pi];
z  = chebpts( n3 );

% Remove doubled-up data
r  =  r(floor(n1/2)+1:end);
ff = ff(floor(n1/2)+1:end,:,:);
ff(:,end+1,:) = ff(:,1,:);

[tt, rr, zz] = meshgrid(th, r, z);

% Slices in the cylinder to plot
rslice = 0;
tslice = tt(1,[1 floor(n2/4)+1 floor(n2/2)+1 floor(3*n2/4)+1],1);
zslice = squeeze(zz(1,1,[floor(n3/4)+1 floor(n3/2)+1 floor(3*n3/4)+1]));

hslicer = slice(tt,rr,zz,ff,tslice,rslice,zslice);
hold on
for j = 1:numel(hslicer)
    h = hslicer(j);
    [xs,ys,zs] = pol2cart(h.XData,h.YData,h.ZData);
    surf(xs,ys,zs,h.CData,'EdgeColor','none','FaceColor','Interp');
end
delete(hslicer)
axis([-1 1 -1 1 -1 1])
daspect([1 1 1])
hold off

set(gca, 'Position', [0 0 1 1], 'CameraViewAngleMode', 'Manual')
colorbar('FontSize', 16, 'Position', [0.84 0.09 0.04 0.8])
axis off

end

function VALS = coeffs2vals( CFS )
% Convert to Chebyshev--Fourier--Chebyshev values

[n1, n2, n3] = size( CFS );
VALS = CFS;
for k = 1:n3
    VALS(:,:,k) = chebtech2.coeffs2vals( VALS(:,:,k) );
    VALS(:,:,k) = trigtech.coeffs2vals( VALS(:,:,k).' ).';
end
for j = 1:n2
    vj = reshape( VALS(:,j,:), n1, n3 );
    vj = chebtech2.coeffs2vals( vj.' ).';
    VALS(:,j,:) = reshape( vj, n1, 1, n3 );
end

end