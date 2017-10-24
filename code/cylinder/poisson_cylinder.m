function X = poisson_cylinder( F, tol )
%POISSON_CYLINDER   Fast Poisson solver for the solid cylinder.
%   POISSON_CYLINDER( F ) solves laplacian(U) = F on the solid cylinder
%   with homogeneous zero Dirichlet boundary conditions. That is, U
%   satisfies
%
%     U_{rr} + U_r/r + U_{thth}/r^2 + U_{zz} = f(r,th,z)
%
%   on the solid cylinder [0,1]x[-pi,pi]x[-1,1] in (r,theta,z), written in
%   cylindrical coordinates with zero Dirichlet conditions at r = 1 and
%   z = +/-1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DEVELOPER'S NOTE
%
%  The underlying algorithm is based on four ideas:
%
%  1) Doubling-up the cylinder so that the radial variable "r" lives on
%     [-1,1] instead of [0,1]. This removes the artificial boundary at the
%     center line of the cylinder.
%
%  2) Decoupling the Fourier modes and solving the resulting 2D partial
%     differential equation for each Fourier mode independently.
%
%  3) Using partial regularity to impose constraints on the solutions of
%     the doubled-up solution so that it is a smooth function on the
%     cylinder.
%
%  4) Using ADI to solve the decoupled 2D PDEs in observed optimal
%     complexity.
%
% The resulting algorithmic cost is observed to be O(n^3). The underlying
% matrices are very mildly non-normal and therefore, we have been unable to
% prove a theoretical complexity statement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tolerance
if ( nargin < 2 )
    tol = 1e-13;
end

% n1 = dof in r, n2 = dof in th, n3 = dof in z:
[n1, n2, n3] = size( F );

% Chebyshev-Fourier-Chebyshev coefficients of the solution:
X = zeros(n1, n2, n3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIAL REGULARITY
%
% For u(r,th,z) = sum sum sum X(i,j,k) T_{i-1}(r)*exp(1i*th*j)*T_{k-1}(z)
% to be a continuously function on the cylinder certain constraints need to
% be imposed. Let u(r,th,z) = sum_j vj(r,z)*exp(1i*j*th), then we need
%
%      vj(r,z) = r^j*yj(r,z)   for some yj.
%
% This is too much to impose because high-order monomials cause numerical
% problems. Thus, we impose "partial regularity" on the solution:
%
%      vj(r,z) = r^{min(j,2)}*wj(r,z)    for some wj.
%
% Therefore, there are two cases to consider: (1) When |j|>=2 and (2) j =
% -1, 0, or 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 1
%
% Here, |j|>=2 where j is the Fourier mode. We want to solve for wj in
% vj(r,z) = r^2*wj(r,z) and then recover vj.  From the boundary conditions
% we know that wj(r,z) = (1-r^2)(1-z^2)yj(r,z). Hence, we solve for yj and
% then recover wj as wj = (1-r^2)(1-z^2)yj. Afterwards, we recover vj as
% vj = r^2wj. The function yj satisfies the following differential
% equation:
%
%  [r^2(1-r^2)yj_{rr} - (9r^3-5r)yj_r - 4(4r^2-1)yj - j^2(1-r^2)*yj](1-z^2) +
%                          r^2(1-r^2)[yj_{zz}-4zyj_{z}-2yj] = fj,
%
% where f(r,th,z) = sum_j fj(r,z)*exp(i*j*th).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tridiagonal discretization of yj_{zz}-4zyj_{z}-2yj, using C^(3/2)
% polynomials. This is specifically designed for ADI:
j = (0:n3-1)';
dsub = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j+2);
dsup = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j);
d = -dsub - dsup;
M = spdiags([dsub d dsup], [-2 0 2], n3, n3);   % Multiplication matrix for (1-x^2) in the C^(3/2) basis
D_inv = spdiags(-1./(j.*(j+3)+2), 0, n3, n3);     % D^{-1} undoes the scaling from the Laplacian identity
T = D_inv * M;                                 % Construct T = D^{-1} * M

% Tridiagonal discretization of yj_{rr}-4ryj_{r}-2yj, using C^(3/2)
% polynomials:
j = (0:n1-1)';
dsub = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j+2);
dsup = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j);
d = -dsub - dsup;
M2 = spdiags([dsub d dsup], [-2 0 2], n1, n1);   % Multiplication matrix for (1-x^2) in the C^(3/2) basis
D_inv2 = spdiags(-1./(j.*(j+3)+2), 0, n1, n1);   % D^{-1} undoes the scaling from the Laplacian identity

% A discretizes ??
I1 = speye( n1 );
Mr2 = I1 - M2;
v1 = j.*(j+1)./(2*j+3);
v2 = (j+2).*(j+3)./(2*j+3);
MD = spdiags([-v1 v2], [-1 1], n1, n1);   % u->d/dx (1-x^2)u in C^(3/2)
v1 = (j+1)./(2*j+3);
v2 = (j+2)./(2*j+3);
Mr = spdiags([v1 v2], [-1 1], n1, n1);    % Multiplication by r in C^(3/2)
A = Mr2/D_inv2 + 5*Mr*MD + 14*M2 - 10*I1;

% Solve for vj(r,z) for |j|>=2:
I3 = speye(n3);
shift = ceil((n2-1)/2);
for k = [-ceil((n2-1)/2):-2 2:floor((n2-1)/2)]
    % We want to solve the following in this loop:
    %     S12\(A-j^2*Mr)*Yj + (Mr2*Mr)*Yj*T^T = S*Fj*C1'*C2'
    % and then set Xj = S^{-1}*(Mr2*(Mr*Xj))*C3.'*C4.' to recover the
    % jth Fourier mode in the Chebyshev-Chebyshev basis.
    
    % Note: This loop costs O(n^2log(n)log(1/tol)) for each "j".
    
    % Construct rhs of matrix equation:
    Fk = reshape(F(:, k+shift+1, :), n1, n3);
    Fk = cheb2ultra(cheb2ultra(Fk).').';   % Convert to C^(3/2)
    Fk = Fk * D_inv.';
    
    %  The following ADI algorithm, is equivalent to:
    %     Xk = chebop2.bartelsStewart( (A-k^2*M2)\(Mr2*M2), I3, I1, T, (A-k^2*M2)\Fk, 0, 0);
    
    %     % Calculate near-optimal ADI shifts:
    %     % We need the largest/smallest eigenvalues of
    
    
    % This fast iterative code for eigs works ALL the time, except for
    % extremely rare instances that cannot be explained.
%     fA = @(v) (A-k^2*M2)\(Mr2*M2*v);
%     opts.tol = 1e-2;
%     a = real(eigs(fA, n1, 1, 'SR', opts));
%     b = real(eigs(fA, n1, 1, 'LR', opts));
    
    a = -min(.05, 1/k^2);
    b = -(4e-4/n1^2)*(k<=n1/2) + (2/k^2/n1^2)*(k>n1/2);
    c = 39*n3^(-4); % An estimate for the largest eigenvalue of inv(T)
    d = 4/pi^2;     % An estimate for the smallest eigenvalue of inv(T)
    
    % Optimal shifts for AX-XB = F, where A and B are normal matrices with
    % spectrum(A) in [a,b] and spectrum(B) in [c,d]:
    [p, q] = ADIshifts(a, b, c, d, tol);
    
    % The ADI method:
    % This is an unwrapped version of the ADI method so that it avoids
    % constructing dense matrices.
    Xk = zeros(n1, n3);
    BB = -T';
    Fk = (A-k^2*M2)\Fk;
    for s = 1:numel( p )
        Xk = (Fk-((A-k^2*M2)\((Mr2*M2)-p(s)*(A-k^2*M2))*Xk)) / (BB-p(s)*I3);
        Xk = ((Mr2*M2)-q(s)*(A-k^2*M2)) \ ((A-k^2*M2)*( Fk - Xk*(BB-q(s)*I3) ));
    end
    
    % Convert result from C^(3/2) to Chebyshev and multiply by (1-z^2):
    Xk = ultra1mx2cheb( Xk.' ).';
    Xk = ultra1mx2cheb( Mr2*Xk );
    
    % Further multiply result by r^2(1-r^2) to recover vj from yj:
    X(:, k+shift+1, :) = reshape( Xk, n1, 1, n3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 2
%
% Here, j = -1 or 1 where j is the Fourier mode. We want to solve for
% vj. We don't impose partial regularity in this case. From the boundary
% conditions we know that vj(r,z) = (1-r^2)(1-z^2)yj(r,z). Hence, we solve
% for yj and then recover wj as wj = (1-r^2)(1-z^2)yj. The function yj
% satisfies the following differential equation:
%
%  [r^2(1-r^2)yj_{rr} - (r-5r^3)yj_r - 4r^2yj - j^2(1-r^2)*yj](1-z^2) +
%                          r^2(1-r^2)[yj_{zz}-4zyj_{z}-2yj] = r^2fj,
%
% where f(r,th,z) = sum_j fj(r,z)*exp(i*j*th).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discretize [r^2(1-r^2)yj_{rr} - (r-5r^3)yj_r - 4r^2yj]
A = Mr/D_inv2 + 3*MD - 6*Mr;
for k = [-1 1]
    % This loop costs O(n^3).
    
    % Solve (A-j^2*Mr2*Mr)*Yj*Mr^T + (Mr2*Mr)*Yj*B^T = Mr2*Fj:
    Fk = reshape(F(:, k+shift+1, :), n1, n3);
    Fk = cheb2ultra(cheb2ultra(Fk).').';   % Convert to C^(3/2)
    Fk = Fk * D_inv.';
    
    Xk = chebop2.bartelsStewart( Mr*M2, I3, A, T, Fk, 0, 0);
    
    Xk = ultra1mx2cheb( Xk.' ).';
    Xk = ultra1mx2cheb( Mr*Xk );
    
    % Further multiply result by r^2(1-r^2) to recover vj from yj:
    X(:, k+shift+1, :) = reshape( Xk, n1, 1, n3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 3
%
% Here, j = 0 where j is the Fourier mode. We want to solve for
% vj. We don't impose partial regularity in this case. From the boundary
% conditions we know that vj(r,z) = (1-r^2)(1-z^2)yj(r,z). Hence, we solve
% for yj and then recover wj as wj = (1-r^2)(1-z^2)yj. The function yj
% satisfies the following differential equation:
%
%  [r^2(1-r^2)yj_{rr} - (r-5r^3)yj_r - 4r^2yj - j^2(1-r^2)*yj](1-z^2) +
%                          r^2(1-r^2)[yj_{zz}-4zyj_{z}-2yj] = r^2fj,
%
% where f(r,th,z) = sum_j fj(r,z)*exp(i*j*th).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = Mr2/D_inv2 + Mr*MD - 2*Mr2;
for k = 0
    
    % Solve (A-j^2*Mr2*Mr)*Yj*Mr^T + (Mr2*Mr)*Yj*B^T = Mr2*Fj:
    Fk = reshape(F(:, k+shift+1, :), n1, n3);
    Fk = cheb2ultra(cheb2ultra(Fk).').';   % Convert to C^(3/2)
    Fk = Fk * D_inv.';
    
    Xk = chebop2.bartelsStewart( Mr2*M2, I3, A, T, Mr2*Fk, 0, 0);
    
    Xk = ultra1mx2cheb( Xk.' ).';
    Xk = ultra1mx2cheb( Xk );
    
    % Further multiply result by r^2(1-r^2) to recover vj from yj:
    X(:, k+shift+1, :) = reshape( Xk, n1, 1, n3);
    
end

end