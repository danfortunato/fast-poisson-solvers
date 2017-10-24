function X = poisson_solid_sphere( F, tol )
%POISSON_SOLID_SPHERE   Fast Poisson solver for the solid sphere.
%   POISSON_SOLID_SPHERE( F ) solves laplacian(U) = F on the solid sphere
%   with homogeneous zero Dirichlet boundary conditions. That is, U
%   satisfies
%
%     (r^2U_{r})_{r}/r^2 + (sin(phi)U_{phi})_{phi}/(r^2sin(phi)) +
%                       U_{thth}/(r^2sin(phi)^2) = f(r,th,phi)
%
%   on the solid sphere [0,1]x[-pi,pi]x[0,pi] in (r,theta,phi), written
%   in spherical coordinates with zero Dirichlet condition at r = 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DEVELOPER'S NOTE
% 
%  The underlying algorithm is based on four ideas:
%
%  1) Doubling-up the solid sphere twice so that the radial variable "r"
%     lives on [-1,1] instead of [0,1] and the polar variable "phi"
%     lives on [-pi,pi] instead of [0,pi]. This removes the artificial
%     boundary at r = 0 as well as making the phi variable a periodic
%     direction.
%
%  2) Decoupling the Fourier modes (in the theta variable) and solving the
%     resulting 2D partial differential equation for each Fourier mode
%     independently.
%
%  3) Using partial regularity to impose constraints on the doubled-up
%     solution so that it is a smooth function on the solid sphere.
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

% n1 = dof in r, n2 = dof in theta, n3 = dof in phi:
[n1, n2, n3] = size( F );

% Chebyshev-Fourier-Fourier coefficients of the solution:
X = zeros(n1, n2, n3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTIAL REGULARITY
% 
% For u(r,th,phi) = sum sum sum X(i,j,k) T_{i-1}(r)*exp(1i*th*j)*exp(1i*phi*k)
% to be a continuously differentiable function on the cylinder certain
% constraints need to be imposed. Let u(r,th,phi) = sum_j vj(r,phi)*exp(1i*j*th),
% then we need
%
%      vj(r,phi) = r^j*sin(phi)^j*yj(r,phi)   for some yj.
%
% This is too much to impose because high-order monomials cause numerical
% problems. Thus, we impose "partial regularity" on the solution:
%
%      vj(r,phi) = (r*sin(phi))^{min(j,2)}*wj(r,phi)    for some wj.
%
% Therefore, there are two cases to consider: (1) When |j|>=2 and (2) j =
% -1, 0, or 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 1
% 
% Here, |j|>=2 where j is the Fourier mode. We want to solve for wj in
% vj(r,phi) = (r*sin(phi))^2*wj(r,phi) and then recover vj.  From the
% boundary conditions we know that wj(r,phi) = (1-r^2)yj(r,phi). Hence, we
% solve for  yj and then recover wj as wj = (1-r^2)yj. Afterwards, we
% recover vj as vj = (r*sin(phi))^2wj. The function yj satisfies the
% following differential equation:
%
%  [r^2(1-r^2)yj_{rr} - 2r(5r^2-3)yj_r + (6-20r^2)yj]sin(phi)^2 +
%      + (1-r^2)[sin(phi)^2yj_{phi,phi}+5sin(phi)cos(phi)yj_{phi}+(3cos(2phi)+1)yj -j^2yj] = fj,
%
% where f(r,th,phi) = sum_j fj(r,phi)*exp(1i*j*th).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discretize [r^2(1-r^2)yj_{rr} - 2r(5r^2-3)yj_r + (6-20r^2)yj]
D1 = ultraS.diffmat(n1, 1);
D2 = ultraS.diffmat(n1, 2);
M1 = ultraS.multmat(n1, [.125;0;0;0;-.125], 2 );   %r^2(1-r^2)
M2 = ultraS.multmat(n1, [0;.75;0;1.25], 1);       %r(5r^2-3)
M3 = ultraS.multmat(n1, [-4;0;-10], 0);           %(6-20r^2)
S01 = ultraS.convertmat(n1, 0, 0);
S02 = ultraS.convertmat(n1, 0, 1);
S12 = ultraS.convertmat(n1, 1, 1);
% A discretizes 
%   r^2(1-r^2)yj_{rr} - 2r(5r^2-3)yj_r + (6-20r^2)yj
A = S12 \ ( M1*D2 - 2*S12*M2*D1 + S02*M3 ) / S01;

% Discretize
%  [sin(phi)^2yj_{phi,phi}+5sin(phi)cos(phi)yj_{phi}+(3cos(2phi)+1)yj]
DF2 = trigspec.diffmat( n3, 2);
DF1 = trigspec.diffmat( n3, 1); 
Msin2 = trigspec.multmat(n3, [-.25;0;.5;0;-.25] );    %sin(phi)^2
Mcossin = trigspec.multmat( n3, [.25i;0;0;0;-.25i]);   %cos(phi)sin(phi)
Mc = trigspec.multmat(n3, [1.5;0;1;0;1.5]);           %(3cos(2phi)+1)
% B discretizes 
%   [sin(phi)^2yj_{phi,phi}+5sin(phi)cos(phi)yj_{phi}+(3cos(2phi)+1)yj]
B = Msin2*DF2 + 5*Mcossin*DF1 + Mc; 

% Loop over the Fourier modes |j|>=2:
I3 = speye( n3 );
Mr = ultraS.multmat( n1, [.5;0;-.5], 1);  %(1-r^2)
Mr2 = ultraS.multmat(n1, [.5;0;.5], 1);   % r^2

% Solve for vj(r,phi) for |j|>=2:
shift = ceil((n2-1)/2);
for j = [-ceil((n2-1)/2):-2 2:floor((n2-1)/2)]
    % We want to solve the following in this loop: 
    %        A*Yj*Msin2^T + Mr*Yj*Bj^T = Fj   for Yj
    % and then set Xj = S01\(Mr2*Mr*Xj)*Msin2.' to recover the 
    % jth Fourier mode in the Chebyshev-Fourier basis. 

    % Note: This loop costs O(n^2log(n)log(1/tol)) for each "j". 

    % Construct rhs of matrix equation: 
    Fj = reshape(F(:, j+shift+1, :), n1, n3);
    Fj = S01*Fj;
    Fj = Mr\(Fj/Msin2.');
    
    % Calculate near-optimal ADI shifts:
    % We need the largest/smallest eigenvalues of 
    %         Msin2\(B-j^2I);
    Bj = B - j^2*I3;
    fB = @(v) Msin2\(Bj*v);
    a = -n1.^4/10; % An estimate for the smallest eigenvalue of Mr^{-1}*A
    b = 2.04;      % An estimate for the largest eigenvalue of Mr^{-1}*A
    c = -real(eigs(fB, n3, 1, 'LR')); % An estimate for the largest eigenvalue of inv(T)
    d = -real(eigs(fB, n3, 1, 'SR')); % An estimate for the smallest eigenvalue of inv(T)
    
    % Optimal shifts for AX-XB = F, where A and B are normal matrices with 
    % spectrum(A) in [a,b] and spectrum(B) in [c,d]: 
    [p, q] = ADIshifts(a, b, c, d, tol);
    
    % The ADI method
    % This is an unwrapped version of the ADI method so that it avoids
    % constructing dense matrices. 
    Xj = zeros(n1, n3);
    for s = 1:numel(p)
        Xj = -((Fj - Mr\(A-p(s)*Mr)*Xj)*Msin2) / (Bj.'+p(s)*Msin2);
        Xj = (A-q(s)*Mr) \ (Mr*(Fj+Xj*(Bj.'+q(s)*Msin2)/Msin2));
    end
    
    % Multiply result by r^2(1-r^2)sin(phi)^2 to recover vj from yj:
    X(:, j+shift+1, :) = reshape( S01\(Mr2*Mr*Xj)*Msin2.', n1, 1, n3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 2
% 
% Here, j = -1, 0, or 1, where j is the Fourier mode. We want to solve for
% vj. We don't impose partial regularity in this case. From the boundary
% conditions we know that vj(r,phi) = (1-r^2)yj(r,phi). Hence, we solve
% for yj and then recover wj as wj = (1-r^2)yj. The function yj
% satisfies the following differential equation:
%
%  [r^2(1-r^2)yj_{rr} + 2r(1-3r^2)yj_r -6r^2yj]*sin(phi)^2 +
%    (1-r^2)[sin(phi)^2yj_{phi,phi}+sin(phi)cos(phi)yj_{phi}-j^2*yj] = r^2sin(phi)^2fj,
%
% where f(r,th,phi) = sum_j fj(r,phi)*exp(1i*j*th).

% Discretize [r^2(1-r^2)yj_{rr} + 4r(1-2r^2)yj_r + 2(1-6r^2)yj]: 
D1 = ultraS.diffmat(n1, 1);
D2 = ultraS.diffmat(n1, 2);
M1 = ultraS.multmat(n1, [.125;0;0;0;.125], 2 );   %r^2(1-r^2)
M2 = ultraS.multmat(n1, [0;-1.25;0;-.75], 1 );       %r(1-2r^2)
M3 = ultraS.multmat(n1, [.5;0;.5], 0 );           %(1-6r^2)
S02 = ultraS.convertmat(n1, 0, 1);
S12 = ultraS.convertmat(n1, 1, 1);
% A discretizes 
%        [r^2(1-r^2)yj_{rr} + 4r(1-2r^2)yj_r + 2(1-6r^2)yj]
A = S02 \ ( M1*D2 + 2*S12*M2*D1 - 6*S02*M3 );

% Discretize [sin(phi)^2yj_{phi,phi}+sin(phi)cos(phi)yj_{phi}]:
DF2 = trigspec.diffmat( n3, 2);
DF1 = trigspec.diffmat( n3, 1); 
Msin2 = trigspec.multmat(n3, [-.25;0;.5;0;-.25] );    %sin(phi)^2
Mcossin = trigspec.multmat( n3, [.25i;0;0;0;-.25i]);  %cos(phi)sin(phi)
% B discretizes 
%      [sin(phi)^2yj_{phi,phi}+sin(phi)cos(phi)yj_{phi}]
B = Msin2*DF2 + Mcossin*DF1; 

% Loop over the Fourier modes |j|<2:
Mr = ultraS.multmat( n1, [.5;0;-.5], 0);  %(1-r^2)
Mr2 = ultraS.multmat(n1, [.5;0;.5], 0);   % r^2
for j = -1:1
    % This loop costs O(n^3) for each "j", but there are only three j's
    % here.
    
    % Tag on the Fourier mode to B: 
    Bj = B - j^2*I3;
    
    % Solve Aj*Xj*Msin2^T + Mr*Xj*B^T = Fj for Xj:
    Fj = reshape(F(:, j+shift+1, :), n1, n3);
    Fj = Mr2*Fj*Msin2.';
    Xj = chebop2.bartelsStewart(A, Msin2, Mr, Bj, Fj, 0, 0); 
    
    % Multiply result by (1-r^2)sin(phi)^2 to recover vj from yj:
    X(:, j+shift+1, :) = reshape( Mr*Xj, n1, 1, n3);
    
end
end