function X = FiniteDifference_ADI( F, tol )
% Solve Poisson's equation in 2D using a FD method
%
% Inputs: 
%    f = righthand side as a chebfun2 
%    n = discretization size 
% 
% Outputs: 
%    u = solution as a chebfun2
%    Te = tridiagonal matrix for even modes
%    To = tridiagonal matrix for odd modes
%
% Dan Fortunato and Alex Townsend, June 2016 

if ( nargin < 2 ) 
   tol = 1e-12; 
end

% Get discretization size: 
[m, n] = size( F );

% Second-order FD diffmat:
h_m = 2/(m-1);
h_n = 2/(n-1);
A = (1/h_m^2)*spdiags( ones(m,1)*[1 -2 1], [-1 0 1], m-2, m-2);
B = -(1/h_n^2)*spdiags( ones(n,1)*[1 -2 1], [-1 0 1], n-2, n-2); 

% Right-hand side:
F = F(2:end-1, 2:end-1); 

% Bounds on the eigenvalues of A: 
a = -4/h_m^2*sin(pi*m/2/(m-1))^2;
b = -4/h_m^2*sin(pi/2/(m-1))^2; 
c = 4/h_n^2*sin(pi/2/(n-1))^2;
d = 4/h_n^2*sin(pi*n/2/(n-1))^2;

% Calculate no. of ADI iterations and shifts: 
[p, q] = ADIshifts(a, b, c, d, tol);

% Solve A*X-X*B=F. 
% Below computes X = lyap(A, -B, -F) in O(n^2log n) operations:
% X = zeros(m-2, n-2);
% Im = speye(m-2);
% In = speye(n-2);
% for j = 1:numel(p)
%     X = (F-(A+q(j)*Im)*X) / (B+q(j)*In);
%     X = (A+p(j)*Im) \ ( F - X*(B+p(j)*In) );
% end

% It must be true that A == B and p == q for the following to work:
d = diag(A);
u = diag(A,1);
l = diag(A,-1);
X = adi_tri( d, u, l, p, F );

% Add zero boundary values: 
Zn = zeros(1, n-2); 
Zm = zeros(m-2, 1);
X = [   0  Zn  0  ; 
        Zm  X  Zm ;
        0  Zn  0  ];

end