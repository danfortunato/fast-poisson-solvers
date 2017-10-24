function pass = test_transforms( )
% Test the fast transforms.

tol = 1e-12;
mx2 = chebfun2(@(x,y) 1-x.^2 );
my2 = chebfun2(@(x,y) 1-y.^2 );

% Test ultra1mx2cheb:

% (1-x^2)(1-y^2):
A = zeros( 10 );
A(1,1) = 1;
f = mx2.*my2;
B = coeffs2( f, size(A,1), size(A,2) );
C = ultra1mx2cheb( ultra1mx2cheb( A ).' ).';
pass(1) = norm( B - C ) < tol;

% 9x(1-x^2)y(1-y^2):
A = zeros( 12 );
A(4,6) = 1;
C3 = ultrapoly( 3, 3/2 );
C5 = ultrapoly( 5, 3/2 );
f = mx2.*my2.*(C3*C5');
B = coeffs2( f, size(A,1), size(A,2) );
C = ultra1mx2cheb( ultra1mx2cheb( A ).' ).';
pass(2) = norm( B - C ) < tol;

A = zeros( 60 );
A(26,56) = 1;
Cy = ultrapoly( 25, 3/2 );
Cx = ultrapoly( 55, 3/2 );
f = mx2.*my2.*(Cy*Cx');
B = coeffs2( f, size(A,1), size(A,2) );
C = ultra1mx2cheb( ultra1mx2cheb( A ).' ).';
pass(3) = norm( B - C ) < tol;

% Test cheb2ultra:
A = zeros( 60 );
A(1,1) = 1;
Cy = ultrapoly( 0, 3/2 );
Cx = ultrapoly( 0, 3/2 );
f = (Cy*Cx');
B = coeffs2( f, size(A,1), size(A,2) );
C = cheb2ultra( cheb2ultra( B ).' ).';
pass(4) = norm( A - C ) < tol;

A = zeros( 10 );
A(6,4) = 1;
Cy = ultrapoly( 5, 3/2 );
Cx = ultrapoly( 3, 3/2 );
f = (Cy*Cx');
B = coeffs2( f, size(A,1), size(A,2) );
C = cheb2ultra( cheb2ultra( B ).' ).';
pass(5) = norm( A - C ) < tol;

A = zeros( 100 );
A(56,23) = 1;
Cy = ultrapoly( 55, 3/2 );
Cx = ultrapoly( 22, 3/2 );
f = (Cy*Cx');
B = coeffs2( f, size(A,1), size(A,2) );
C = cheb2ultra( cheb2ultra( B ).' ).';
pass(6) = norm( A - C ) < tol;

if ( all(pass) )
    pass = 1;
end

end