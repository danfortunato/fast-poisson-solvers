function pass = test_poisson_solid_sphere( )
% Test the fast Poisson solver for the solid sphere.

tol = 1e-13;
pass = [];

% Test the Fourier modes from -3 to 3
n1 = 21; n2 = 22; n3 = 24;
r = chebpts( n1 );
th = pi*trigpts( n2 );
lam = pi*trigpts( n3 );
[rr, tt, ll] = ndgrid( r, th, lam );
for k = -3:3
    exact = @(r, th, lam) (1-r.^2).*r.^abs(k).*sin(lam).^abs(k).*exp(1i*k*th);
    rhs = @(r, th, lam) -2*(2*abs(k)+3).*r.^abs(k).*sin(lam).^abs(k).*exp(1i*k*th);
    EXACT = exact(rr, tt, ll);
    EXACT = vals2coeffs( EXACT );
    F = rhs(rr, tt, ll);
    F = vals2coeffs( F );
    X = poisson_solid_sphere( F, tol );
    pass(end+1) = norm( X(:) - EXACT(:) ) < 5*tol;
end

if ( all(pass) )
    pass = 1;
end

end

function CFS = vals2coeffs( VALS )
% Convert to Chebyshev--Fourier--Fourier coefficients

[n1, n2, n3] = size( VALS );
CFS = VALS;
for k = 1:n3
    CFS(:,:,k) = chebtech2.vals2coeffs( CFS(:,:,k) );
    CFS(:,:,k) = trigtech.vals2coeffs( CFS(:,:,k).' ).';
end
for j = 1:n2
    vj = reshape( CFS(:,j,:), n1, n3 );
    vj = trigtech.vals2coeffs( vj.' ).';
    CFS(:,j,:) = reshape( vj, n1, 1, n3 );
end

end
