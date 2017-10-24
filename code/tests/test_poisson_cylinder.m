function pass = test_poisson_cylinder( )
% Test the fast Poisson solver for the solid cylinder.

tol = 1e-13;
pass = [];

% Test the Fourier modes from -3 to 3
n1 = 21; n2 = 22; n3 = 24;
r = chebpts( n1 );
th = pi*trigpts( n2 );
z = chebpts( n3 );
[rr, tt, zz] = ndgrid( r, th, z );
for k = -3:3
    exact = @(r, th, z) (1-r.^2).*(1-z.^2).*r.^abs(k).*exp(1i*k*th);
    rhs = @(r, th, z) 2*(r.^2+2*(abs(k)+1)*z.^2-2*abs(k)-3).*r.^abs(k).*exp(1i*k*th);
    EXACT = exact(rr, tt, zz);
    EXACT = vals2coeffs( EXACT );
    F = rhs(rr, tt, zz);
    F = vals2coeffs( F );
    X = poisson_cylinder( F, tol );
    pass(end+1) = norm( X(:) - EXACT(:) ) < tol;
end

% Cartesian test case
n1 = 55; n2 = 84; n3 = 42;
r = chebpts( n1 );
th = pi*trigpts( n2 );
z = chebpts( n3 );
[rr, tt, zz] = ndgrid( r, th, z );
xx = rr.*cos(tt);
yy = rr.*sin(tt);
v = @(x,y,z) (1-x.^2-y.^2).*(1-z.^2).*(z.*cos(4*pi*(x.^2))+cos(4*pi*y.*z));
f = lap(chebfun3(v));
V = vals2coeffs(v(xx,yy,zz));
F = vals2coeffs(f(xx,yy,zz));
X = poisson_cylinder( F, tol );
pass(end+1) = norm( X(:) - V(:) ) < tol;

if ( all(pass) )
    pass = 1;
end

end

function CFS = vals2coeffs( VALS )
% Convert to Chebyshev--Fourier--Chebyshev coefficients

[n1, n2, n3] = size( VALS );
CFS = VALS;
for k = 1:n3
    CFS(:,:,k) = chebtech2.vals2coeffs( CFS(:,:,k) );
    CFS(:,:,k) = trigtech.vals2coeffs( CFS(:,:,k).' ).';
end
for j = 1:n2
    vj = reshape( CFS(:,j,:), n1, n3 );
    vj = chebtech2.vals2coeffs( vj.' ).';
    CFS(:,j,:) = reshape( vj, n1, 1, n3 );
end

end
