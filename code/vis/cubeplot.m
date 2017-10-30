function cubeplot( X, type )

if ( nargin == 1 )
    type = 'coeffs';
end

switch type
    case 'coeffs'
        u = chebfun3( X, 'coeffs' );
    case 'vals'
        u = chebfun3( X );
    otherwise
        error('Unknown input type.');
end

slice(u, 0, 0, 0)
axis([-1 1 -1 1 -1 1])
daspect([1 1 1])
set(gca, 'Position', [0 0 1 1], 'CameraViewAngleMode', 'Manual')
colorbar('FontSize', 16, 'Position', [0.87 0.09 0.04 0.8])
axis square
axis off

end
