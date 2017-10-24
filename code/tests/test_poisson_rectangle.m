function pass = test_poisson_rectangle( )
% Test the fast Poisson solver for the rectangle.

%% Tests comparing to analytical solutions
tol = 1e-13;

% Simple solution:
n = 100;
f = chebfun2( @(x,y) -2*(1-x.^2) - 2*(1-y.^2) );
F = coeffs2( f, n, n );
exact = chebfun2( @(x,y) (1-x.^2).*(1-y.^2) );
EXACT = coeffs2( exact, n, n );
X = poisson_rectangle( F );
pass(1) = norm( X(:) - EXACT(:) ) < tol;

% More complex solution:
n = 100;
f = chebfun2( @(x,y) 1/2*exp((x+y)/2) .* (2*pi*sin(pi*x).*cos(pi*y) ...
                     + (2*pi*cos(pi*x)+(1-4*pi^2)*sin(pi*x)).*sin(pi*y)) );
F = coeffs2( f, n, n );
exact = chebfun2( @(x,y) sin(pi*x).*sin(pi*y).*exp((x+y)/2) );
EXACT = coeffs2( exact, n, n );
X = poisson_rectangle( F );
pass(2) = norm( X(:) - EXACT(:) ) < tol;

% Nonzero Dirichlet data with rectangular domain and discretization:
m = 90; n = 100;
a = 0; b = 2; c = -4.3; d = 1;
f = chebfun2( @(x,y) -2*pi^2*cos(pi*x).*cos(pi*y), [a b c d] );
F = coeffs2( f, m, n );
lbc = @(y) cos(pi*a)*cos(pi*y) + a + y;
rbc = @(y) cos(pi*b)*cos(pi*y) + b + y;
dbc = @(x) cos(pi*c)*cos(pi*x) + x + c;
ubc = @(x) cos(pi*d)*cos(pi*x) + x + d;
exact = chebfun2( @(x,y) cos(pi*x).*cos(pi*y) + x + y, [a b c d] );
EXACT = coeffs2( exact, m, n );
X = poisson_rectangle( F, lbc, rbc, dbc, ubc, [a b c d] );
pass(3) = norm( X(:) - EXACT(:) ) < tol;

%% Tests comparing to chebop2's solution
tol = 1e7*chebfunpref().cheb2Prefs.chebfun2eps;

% Square discretization size:
n = 101;
N = chebop2( @(u) lap(u) );
N.bc = 0;
f = chebfun2( @(x,y) 1+0*x );
F = coeffs2( f, n, n );
exact = solvepde( N, f, n, n );
EXACT = coeffs2( exact, n, n );
X = poisson_rectangle( F );
pass(4) = norm( X(:) - EXACT(:) ) < tol;

% Rectangular discretization size:
m = 87; n = 101;
N = chebop2( @(u) lap(u) );
N.bc = 0;
f = chebfun2( @(x,y) 1+0*x );
F = coeffs2( f, m, n );
exact = solvepde( N, f, m, n );
EXACT = coeffs2( exact, m, n );
X = poisson_rectangle( F );
pass(5) = norm( X(:) - EXACT(:) ) < tol;

% Nonstandard rectangular domain:  
m = 99; n = 101;
a = -3; b = 1; c = 4; d = 4.2;
N = chebop2( @(u) lap(u), [a b c d] );
N.bc = 0;
f = chebfun2( @(x,y) 1+0*x, [a b c d] );
F = coeffs2( f, m, n );
exact = solvepde( N, f, m, n );
EXACT = coeffs2( exact, m, n );
X = poisson_rectangle( F, [a b c d] );
pass(6) = norm( X(:) - EXACT(:) ) < tol;

% Nonzero, but constant, Dirichlet data:
m = 99; n = 101;
a = -3; b = 1; c = 4; d = 4.2;
N = chebop2( @(u) lap(u), [a b c d] );
N.bc = 1;
f = chebfun2( @(x,y) 1+0*x, [a b c d] );
F = coeffs2( f, m, n );
lbc = @(y) 1+0*y;
rbc = @(y) 1+0*y;
dbc = @(x) 1+0*x;
ubc = @(x) 1+0*x;
exact = solvepde( N, f, m, n );
EXACT = coeffs2( exact, m, n );
X = poisson_rectangle( F, lbc, rbc, dbc, ubc, [a b c d] );
pass(7) = norm( X(:) - EXACT(:) ) < tol;

% General Dirichlet data, given as function handle:
n = 201;
a = -3; b = 1; c = 4; d = 4.2;
p = @(x,y) x.*y + cos(3*x.^2.*(y-.2));
N = chebop2( @(u) lap(u), [a b c d] );
lbc = @(y) p(a,y); N.lbc = lbc;
rbc = @(y) p(b,y); N.rbc = rbc;
dbc = @(x) p(x,c); N.dbc = dbc;
ubc = @(x) p(x,d); N.ubc = ubc;
f = chebfun2( @(x,y) 1+0*x, [a b c d] );
F = coeffs2( f, n, n );
exact = solvepde( N, f, n, n );
EXACT = coeffs2( exact, n, n );
X = poisson_rectangle( F, lbc, rbc, dbc, ubc, [a b c d] );
pass(8) = norm( X(:) - EXACT(:) ) < 100*tol;

if ( all(pass) )
    pass = 1;
end

end
