%=====================================%
%  Numerical solution of the Blasius  %
%             equation                %
%=====================================%
%
% The Blasius nonlinear ODE of third order 
%
% f''' + 0.5*ff'' = 0 
%
% with boundary conditions: f(0)=0, f'(0)=0, f'(infty)=1
% is transformed into a system of three nonlinear ODEs of first order
% by applying the following change of variables:
%
% f' = u
% f''= xi = u'
% 
% which turns the Blasius equation into
%
% xi' + 0.5*f*u=0.
%
% This equation, together with the definitions above forms
% a system of three of three nonlinear coupled ODE 
%
% ----                          B.C.
% |                             
% | xi' = -0.5*f*xi             integral[0,infty) xi deta = 1
% | u'  = xi                    u(0)=0
% | f'  = u                     f(0)=0
% |
% -----
%
% A vector y of unknowns and a vector F for the r.h.s. are defined as:
%
%     --   --              --         --           
%     |     |              |           |         
%     |  xi |              | -0.5*f*xi |
% y = |  u  |        F(y)= |     xi    |       
%     |  f  |              |     u     |
%     |     |              |           |
%     --   --              --         --
%
% so that the system to be solved can be rewritten as
%
%   dy
%   --   == F(y)
%  deta
%

close all
clear all

% INPUT PARAMETERS
% ----------------
etamax=20;      % interval [0,infty) is discretized into [0,etamax]
n=1000;         % number of nodes
% ----------------


deta=etamax/(n-1);

%                              dy
%  Discretize the derivative   --   with forward Euler finite differences
%                             deta
%            dy
%  so that   --   = Ay, where A is a matrix of finite difference coeff.
%           deta

Axi = (-diag(ones(n,1))+diag(ones(n-1,1),1))/deta;
Au  = (-diag(ones(n,1))+diag(ones(n-1,1),1))/deta;
Af  = (-diag(ones(n,1))+diag(ones(n-1,1),1))/deta;

Axi(end,:) = deta*[0.5 ones(1,n-2) 0.5];   % integral BC 
                                           % with trapezioidal integration
Af(end,:)  = [1 zeros(1,n-1)];             % BC u(0)=0
Au(end,:)  = [1 zeros(1,n-1)];             % BC f(0)=0
                            %    ---     ---
                            %    | Axi     |
A = blkdiag(Axi, Au, Af);   % A =|    Au   |
                            %    |      Af |
                            %    ---     ---

% Initial guess for the slution vector (linear profiles which respect B.C.)
y = [0.5*linspace(1,0,n)'*deta; linspace(0,1,n)'; linspace(0,1,n)'];

% Iterative solution of the nonlinear system Ay - F(y)=0 with Newton
toll=1e-6; normy = 1e10;

% Ay - F(y) = H(y) = 0 is solved with Newton by iterative linearization
%
% H(y) ~ H(y0) + dH(y0)/dy * (y-y0) = 0 
%
% if we put dH(y0)/dy = Jt(y0) then we can write
%
% (y-y0) = -inv(Jt(y0))*H(y0)
%
% in out specific case Jt = A - dF/dy and H=Ay - F(y)
%
% In the following J = dF/dy only

while normy>toll
   % Compute F
   F = [ -0.5*y(1:n).*y(2*n+1:3*n); y(1:n); y(n+1:2*n)];
   F(n) = 1; F(2*n) = 0; F(3*n) = 0; % B.C.
   % Compute J
   J11 = -0.5*diag(y(2*n+1:3*n)); J12 = diag(zeros(1,n));  J13 = -0.5*diag(y(1:n));
   J21 = diag(ones(n,1));         J22 = diag(zeros(n,1));  J23 = diag(zeros(n,1));
   J31 = diag(zeros(n,1));        J32 = diag(ones(n,1));   J33 = diag(zeros(n,1));
   J   = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
   J(n,:) = zeros(1,3*n);    %B.C.
   J(2*n,:) = zeros(1,3*n);  %B.C.
   J(3*n,:) = zeros(1,3*n);  %B.C.
   % Solve
   dy =  - (A-J)\(A*y-F);  normy=norm(dy)
   y = y + dy;
end

% Rename the part of the solution vector with our usual name
eta     = 0:deta:etamax;
fsecond = y(1:n); 
fprime  = y(n+1:2*n);
f       = y(2*n+1:3*n);

% Plot
figure(1)
set(gcf(), 'Units', 'centimeters', 'Position', [0 0 10 10])
set(gca(), 'FontSize', 14)
plot(f, eta, fprime, eta, fsecond, eta)
xlim([0 1.01])
ylim([0 10])
h=legend('f', 'f^\prime', 'f^{\prime\prime}');
set(h, 'Fontsize', 14)
ylabel('\eta = y (u_\infty/(\nu x))^{0.5}')

figure(2)
set(gcf(), 'Units', 'centimeters', 'Position', [0 0 10 10])
set(gca(), 'FontSize', 14)
plot(eta, eta.*fprime'-f')
xlim([0 10])
h=legend('\etaf^\prime - f');
set(h, 'Fontsize', 14)
xlabel('\eta','FontSize', 14)
