function X = poisson_rectangle( F, varargin )
%POISSON_RECTANGLE   Fast Poisson solver for the rectangle.
%   POISSON_RECTANGLE( F ) solves laplacian(U) = F on [-1,1]x[-1,1] with
%   zero Dirichlet boundary conditions. That is, U satisfies
%
%     U_{x,x} + U_{y,y} = F, on [-1,1]x[-1,1]    U = 0 on boundary
%
%   F is input as an M x N matrix of Chebyshev coefficients. The equation
%   is solved using an M x N discretization. G can be a scalar, a function
%   handle, or any chebfun2 object satisfying the Dirichlet data.
%
%   POISSON_RECTANGLE( F, TOL ) solves to the specified error tolerance.
%
%   POISSON_RECTANGLE( F, [A B C D] ) solves on the domain [A,B]x[C,D].
%
%   POISSON_RECTANGLE( F, [A B C D], TOL ) solves on the domain [A,B]x[C,D]
%   to the specified error tolerance.
%
%   POISSON_RECTANGLE( F, LBC, RBC, DBC, UBC ) solves using the given
%   Dirichlet data. The data are given as function handles.
%
%   POISSON_RECTANGLE( F, LBC, RBC, DBC, UBC, [A B C D] ) solves using the
%   given Dirichlet data on the domain [A,B]x[C,D].
%
%   POISSON_RECTANGLE( F, LBC, RBC, DBC, UBC, TOL ) solves using the
%   given Dirichlet data to the specified error tolerance.
%
%   POISSON_RECTANGLE( F, LBC, RBC, DBC, UBC, [A B C D], TOL ) solves using
%   given Dirichlet data on the domain [A,B]x[C,D] to the specified error
%   tolerance.

% DEVELOPER'S NOTE:
%
% METHOD: Spectral method (in coefficient space). We use a C^{(3/2)} basis
% to discretize the equation, resulting in a discretization of the form
% AX + XA = F, where A is a symmetric tridiagonal matrix.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation is solved by the
% alternating direction implicit (ADI) method.
%
% SOLVE COMPLEXITY:  O(M*N*log(MAX(M,N))*log(1/eps)) with M*N = total
% degrees of freedom.
% 
% AUTHORS: Dan Fortunato (dan.fortunato@gmail.com)
%          Alex Townsend (townsend@cornell.edu)
%
% The fast Poisson solver is based on:
%
% D. Fortunato and A. Townsend, Fast Poisson solvers for spectral methods,
% in preparation, 2017.

[m, n] = size( F );

% Default arguments
dom = [-1, 1, -1, 1];
tol = 1e-13;
BC = zeros(m, n);
nonzeroBC = false;

% Parse the inputs
if ( nargin == 2 )
    if ( numel(varargin{1}) == 1 )
        % Call is POISSON_RECTANGLE( F, tol )
        tol = varargin{1};
    else
        % Call is POISSON_RECTANGLE( F, dom )
        dom = varargin{1};
    end
elseif ( nargin == 3 )
    % Call is POISSON_RECTANGLE( F, dom, tol )
    [dom, tol] = varargin{:};
elseif ( nargin == 5 )
    % Call is POISSON_RECTANGLE( F, lbc, rbc, dbc, ubc )
    [lbc, rbc, dbc, ubc] = varargin{:};
    nonzeroBC = true;
elseif ( nargin == 6 )
    [lbc, rbc, dbc, ubc] = varargin{1:4};
    nonzeroBC = true;
    if ( numel(varargin{5}) == 1 )
        % Call is POISSON_RECTANGLE( F, lbc, rbc, dbc, ubc, tol )
        tol = varargin{5};
    else
        % Call is POISSON_RECTANGLE( F, lbc, rbc, dbc, ubc, dom )
        dom = varargin{5};
    end
elseif ( nargin == 7 )
    % Call is POISSON_RECTANGLE( F, lbc, rbc, dbc, ubc, dom, tol )
    [lbc, rbc, dbc, ubc, dom, tol] = varargin{:};
    nonzeroBC = true;
elseif ( nargin ~= 1 )
    error('POISSON_RECTANGLE:ARG', 'Invalid number of arguments.');
end

% Check that the domain is valid
if ( numel(dom) ~= 4 )
    error('POISSON_RECTANGLE:DOMAIN', 'Invalid domain specification.');
end

% Solve for u on the given domain, adjust diffmat to including scaling:
scl_x = (2/(dom(2)-dom(1)))^2;
scl_y = (2/(dom(4)-dom(3)))^2;

% Solver only deals with zero homogeneous Dirichlet conditions. Therefore,
% if nonzero Dirichlet conditions are given, we solve lap(u) = f with
% u|bc = g as u = v + w, where v|bc = g, and lap(w) = f - lap(v), w|bc = 0:
if ( nonzeroBC )
    % Check that the boundary conditions are valid
    if ~( isa(lbc, 'function_handle') && ...
          isa(rbc, 'function_handle') && ...
          isa(dbc, 'function_handle') && ...
          isa(ubc, 'function_handle') )
        error('POISSON_RECTANGLE:BC', ...
            'Dirichlet data needs to be given as a function handle.');
    end
    
    % Make sure the Dirichlet data match at the corners
    if ~( lbc(dom(3)) == dbc(dom(1)) && ...
          lbc(dom(4)) == ubc(dom(1)) && ...
          rbc(dom(3)) == dbc(dom(2)) && ...
          rbc(dom(4)) == ubc(dom(2)) )
          error('POISSON_RECTANGLE:BC', 'Corners must match.');
    end
    
    % First convert the boundary data from function handles to coefficients
    lbc_cfs = fun2coeffs( lbc, m, dom(3:4) );
    rbc_cfs = fun2coeffs( rbc, m, dom(3:4) );
    dbc_cfs = fun2coeffs( dbc, n, dom(1:2) );
    ubc_cfs = fun2coeffs( ubc, n, dom(1:2) );

    % Now compute an interpolant of the boundary data
    BC(1,:) = (ubc_cfs + dbc_cfs)/2;
    BC(2,:) = (ubc_cfs - dbc_cfs)/2;
    BC(1:2,1) = (rbc_cfs(1:2) + lbc_cfs(1:2))/2 - sum(BC(1:2,3:2:end),2);
    BC(1:2,2) = (rbc_cfs(1:2) - lbc_cfs(1:2))/2 - sum(BC(1:2,4:2:end),2);
    BC(3:end,1) = (rbc_cfs(3:end) + lbc_cfs(3:end))/2;
    BC(3:end,2) = (rbc_cfs(3:end) - lbc_cfs(3:end))/2;
    
    % Adjust the rhs
    % TODO: Remove this call to chebfun2
    u_lap = lap( chebfun2(BC, dom, 'coeffs') );
    u_lap = coeffs2( u_lap, m, n );
    F = F - u_lap;
end

% Convert rhs to C^{(3/2)} coefficients: 
F = cheb2ultra( cheb2ultra( F ).' ).';

% Construct M, the multiplication matrix for (1-x^2) in the C^(3/2) basis
jj = (0:n-1)';
dsub = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj+2);
dsup = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj);
d = -dsub - dsup;
Mn = spdiags([dsub d dsup], [-2 0 2], n, n);
% Construct D^{-1}, which undoes the scaling from the Laplacian identity
invDn = spdiags(-1./(jj.*(jj+3)+2), 0, n, n);
Tn = scl_y * invDn * Mn;

jj = (0:m-1)';
dsub = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj+2);
dsup = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj);
d = -dsub - dsup;
Mm = spdiags([dsub d dsup], [-2 0 2], m, m);
invDm = spdiags(-1./(jj.*(jj+3)+2), 0, m, m);

% Construct T = D^{-1} * M:
Tm = scl_x * invDm * Mm;
F = invDm * F * invDn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Alternating Direction Implicit method %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve TmX + XTn' = F using ADI, which requires O(n^2log(n)log(1/eps)) 
% operations:

% Calculate ADI shifts based on bounds on the eigenvalues of Tn and Tm:
a = -4/pi^2 * scl_y;
b = -39*n^-4 * scl_y;
c = 39*m^-4 * scl_x;
d = 4/pi^2 * scl_x;
[p, q] = ADIshifts(a, b, c, d, tol);

% Run the ADI method:
A = Tm; B = -Tn;
% Extract diagonals from A and B and call ADI C code
da = diag(A);    db = diag(B);
ua = diag(A,2);  ub = diag(B,2);
la = diag(A,-2); lb = diag(B,-2);
X = adi( da, ua, la, db, ub, lb, p, q, F );

% Convert back to Chebyshev
X = ultra1mx2cheb( ultra1mx2cheb( X ).' ).';
X = X + BC;

end

function cfs = fun2coeffs( fun, n, dom )
% Convert from a function handle to Chebyshev coefficients

vals = fun( chebpts(n, dom) );
cfs = chebtech2.vals2coeffs( vals );

end
