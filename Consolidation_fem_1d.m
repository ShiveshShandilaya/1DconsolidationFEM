% Date of creation: 28 May 2022
close all;
clear all;
clc;
%-------------------------------------------------------------------------%
% 1] Geometry and material specifications
%-------------------------------------------------------------------------%
E = 1800;                % Elasticity modulus of the material used (kN/m3)
k = 6.28E-5;             % Permeability (m/s)
z = 7;                   % Height of the sample(m)
delt = 0.1;              
mv = 5.6E-4;                  % coefficient of volume consolidation (m2/kN)
Yw = 9.81;
Alpha = 0.5;
Load = -100;        % Load (kN) applied on the top of the sample
%%
% H = 1;
% Tv = 0.01;               % Time factor
% Yw = 9.81;               % KN/m3
% v = 0;                   % Poisson's Ratio
% k=0.001;                 % m/day
% Eoed=((1-v)*E)/((1+v)*(1-2*v));
% Cv=k*Eoed/Yw;
% delt = Tv*H^2/Cv;              

numMesh = 8;             % Number of mesh 
%-------------------------------------------------------------------------%
% 2] Create a mesh and connectivity matrix for the specified domain
%-------------------------------------------------------------------------%
Nnodes = numMesh+1;
x0 = 0;
xs = linspace(x0, x0+z, Nnodes)';
ncoord = zeros(Nnodes,2);
le =7/8;
m = 0;
for i=2:Nnodes
    ncoord(1,:) = [0 xs(1,1)];
    ncoord(i,:) = [0 xs(i,1)];
end
eConn = [[1:Nnodes-1]', [2:Nnodes]'];
Plot_Initial(ncoord,eConn,le)
%-------------------------------------------------------------------------%
% 3] Specification of loads, boundary conditions
%-------------------------------------------------------------------------%
u0_given = 1; u0 = 0;
uL_given = 1; uL = 0;
Q0_given = 0; Q0 = 100;
QL_given = 1; QL = 100;
%-------------------------------------------------------------------------%
% 4] Compute K and F for the system in assembled form
%-------------------------------------------------------------------------%
Nnodes = length(xs);
Ndof = Nnodes;
Nelems = Nnodes-1;
%initializations
Kg = zeros(Ndof, Ndof);
Mg = zeros(Ndof, Ndof);
Qg = zeros(Ndof, 1);
for e = 1:Nelems
   n1 = eConn(e,1); %left node of the element
   n2 = eConn(e,2); %right node of the element
   
   %initialization
   Ke = zeros(2,2);
   Me = zeros(2,2);
   Qe = zeros(2,1);
   h = z/numMesh;
   syms p
   ksi_1 = 1-(p/h);
   ksi_2 = (p/h);
   B = [diff(ksi_1,p) diff(ksi_2,p)];
   ksi = [ksi_1 ksi_2];
   Ke = Ke + k.*[int(B(1,1).'.*B(1,1),p,0,h) int(B(1,1).'.*B(1,2),p,0,h).......
         ;int(B(1,2).'.*B(1,1),p,0,h) int(B(1,2).'.*B(1,2),p,0,h)];
   Me = Me + (mv*Yw).*[int(ksi_1.*ksi_1,p,0,h) int(ksi_1.*ksi_2,p,0,h);.....
       int(ksi_2.*ksi_1,p,0,h) int(ksi_2.*ksi_2,p,0,h)];
   Qe = Qe + k.*[(1.*subs(ksi_1,0)-0.*subs(ksi_1,h));(0.*subs(ksi_2,0)-1.*.....
       subs(ksi_2,h))];
   %Assembly
   Kg([n1 n2], [n1 n2]) = Kg([n1 n2], [n1 n2]) + Ke;
   Mg([n1 n2], [n1 n2]) = Mg([n1 n2], [n1 n2]) + Me;
   Qg([n1 n2], 1) = Qg([n1 n2], 1) + Qe;
end
%-------------------------------------------------------------------------%
% 5] Application of boundary conditions
%-------------------------------------------------------------------------%
% Neumann conditions
if Q0_given == 1
   Qg(1) = Qg(1) + Q0;
end
if QL_given == 1
   Qg(end) = Qg(end) + QL;
end
% Dirichlet conditions
if u0_given == 1
   Kg(1,:) = 0;
   Kg(1,1) = 1;
   Mg(1,:) = 0;
   Mg(1,1) = 1;
   Qg(1) = u0;
end
if uL_given == 1
   Kg(end,:) = 0;
   Kg(end,end) = 1;
   Mg(end,:) = 0;
   Mg(end,end) = 1;
   Qg(end) = uL;
end
%-------------------------------------------------------------------------%
% 6] Solve for 'u'
%-------------------------------------------------------------------------%
u1 = [0 100 100 100 100 100 100 100 0];     % u(z,0) = 100 kpa all initial stress will be taken by pore water
% for s =1:delt:300
%     u_2 = u1*(Mg-(1-Alpha)*delt*Kg)/(Mg+(Alpha)*delt*Kg);
%     u1 = u_2;
% end
% u_fem = u_2';
% u_fem(1,1) = 0;
% u_fem(end,1) = 0;
u_f = zeros(9,1);
t=0;
while t<90
    u_2 = u1*(Mg-(1-Alpha)*delt*Kg)/(Mg+(Alpha)*delt*Kg);
    t=t+delt;
    u1 = u_2;
    u_f = [u_2];
end
u_fem = u_f';
u_fem(1,1) = 0;
u_fem(end,1) = 0;
%-------------------------------------------------------------------------%
% 7] Post Processing
%-------------------------------------------------------------------------%
f=1:9;
figure;hold on;grid on;
plot(f, u_fem); title('Numerical Solution');
ylabel('Pore Pressure'); xlabel('Nodes along Depth');
legend('Numerical Solution');