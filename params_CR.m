%------------------------------------------------------------------------%
%   [Ir, Kc, Lc, Rc, Rr] = params_CR(J_sol, L_sol, R1_sol, R2_sol)
%   solves for the cylindrical shell parameters and ring parameters
%   using the equations defined by case CR
%   case CR equations [b1; b3; b5; b7; Bz];
%
%   INPUTS
%   J_sol: J value for solenoid
%   L_sol: length of solenoid
%   R_1: inner radius of solenoid
%   R_2: outer radius of solenoid
%
%   OUTPUTS
%   Ir: current for ring
%   Kc: current for cylinder
%   Lc: length of cylinder
%   Rc: radius of cylinder
%   Rr: radius of ring
%
%   REMARKS
%   case CR = [b1; b3; b5; b7; Bz];
%   The optimal paramters of solved by using an LM algorithum
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS:
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function [Ir, Kc, Lc, Rc, Rr, soln_found] = params_CR(J_sol, L_sol, R1_sol, R2_sol)
format long;
%% solenoid parameters
realS_J = J_sol;
realS_L = L_sol;
realS_R1 = R1_sol;
realS_R2 = R2_sol;

delta_R = realS_R2 - realS_R1;
real_z = 0;
real_mu = 1;%(4*pi)*1e-7;

c = c_vector_CR(realS_J, realS_L, realS_R1, realS_R2, real_mu, real_z);

%% inputs to LM function
ravg = (realS_R2+realS_R1)/2;
rdiff = (realS_R2-realS_R1)/2;

inGuess = realS_R1 + rdiff/2;
outGuess = realS_R2 - rdiff/2;

%intial guess for the parameters
cur_scale_ring = c(1)/(pi*((inGuess)^2+(outGuess)^2));
cur_scale_shell = cur_scale_ring/realS_L;

gammaGuess = [1; 1; realS_L; inGuess; outGuess];

%same results as above
% gammaGuess = [realS_J*delta_R*realS_L/2; realS_J*delta_R/2; realS_L; inGuess; outGuess];

tol = 1e-5;

%basically going to run forever unless solution is found or manually
%stopped
maxIterations = 500;
maxTime = 60*5;
dispVals = 0;
maxChange = [];

%% LM method
[ currentEstimate, cost, soln_found ] = LM_LeastSquares( gammaGuess, @(x)costF(x,c,cur_scale_ring,cur_scale_shell), tol, maxIterations, maxTime, dispVals, maxChange );

% if( any(currentEstimate <0 ))
%     warning('Solution Found is non-physical. Resolving...');
%     [ currentEstimate, cost, soln_found ] = LM_LeastSquares( abs(currentEstimate), @(x)costF(x,c,cur_scale_ring,cur_scale_shell), tol, maxIterations, maxTime, dispVals, maxChange );
%     currentEstimate = abs(currentEstimate);
%   
% end

% currentEstimate
% cost
currentEstimate = abs(currentEstimate);
soln_found = cost < 1e-3 && ~any(currentEstimate <0 );

%% final values
Ir = currentEstimate(1)*cur_scale_ring;
Kc = currentEstimate(2)*cur_scale_shell;
Lc = currentEstimate(3);
Rc = currentEstimate(4);
Rr = currentEstimate(5);




end

%% cost function
function [error, Jacobian] = costF(gammaGuess,c,cur_scale_ring, cur_scale_shell)
real_z = 0;
real_mu = 1;%(4*pi)*1e-7;

f =              f_vector_CR(gammaGuess(1)*cur_scale_ring, gammaGuess(2)*cur_scale_shell, gammaGuess(3), gammaGuess(4), gammaGuess(5), real_mu, real_z);

J_gamma = jacobian_matrix_CR(gammaGuess(1)*cur_scale_ring, gammaGuess(2)*cur_scale_shell, gammaGuess(3), gammaGuess(4), gammaGuess(5), real_mu, real_z);

J_gamma(:,1) = J_gamma(:,1)*cur_scale_ring;
J_gamma(:,2) = J_gamma(:,2)*cur_scale_shell;
Jacobian = J_gamma;

error = f -c;

error = error/norm(c);
Jacobian = Jacobian/norm(c);

%scale the values of not singularities
for i = 1:length(error)
    sc = norm(Jacobian(i,:));
    error(i) = error(i)/sc;
    Jacobian(i,:) = Jacobian(i,:)/sc;
end

end

%------------------------------------------------------------------------%  
%   c_1 = c_vector_CR(J,Lsol,R1,R2,mu,z)
%   constants in equation (solenoid side of system of equations) for case (CR) where the
%   optimal solution is the super postion of a cylindrical shell and a ring
%
%   INPUTS
%   J: J value for solenoid
%   Lsol: length of solenoid
%   R1: inner radius of solenoid
%   R2: outer radius of solenoid
%   mu: magnetic permeability 
%   z: z location in solenoid (0 if at center - origin is at center)
%
%   OUTPUTS
%   c_1 = 5x1 vector of constants based on solenoid paramters for the 5
%   equations used in case CR
%
%   REMARKS
%   equation solving for CR case [b1; b3; b5; b7; Bz];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function c_1 = c_vector_CR(J,Lsol,R1,R2,mu,z)
%C_VECTOR_CR
%    C_1 = C_VECTOR_CR(J,LSOL,R1,R2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    01-Apr-2019 08:36:36

t2 = R1.^2;
t3 = R2.^2;
t4 = Lsol.^2;
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = R1.*t2-R2.*t3;
t9 = t5.^2;
t10 = t6.^2;
t11 = Lsol.*z;
t12 = t4./4.0;
t13 = z.^2;
t14 = t2+t11+t12+t13;
t15 = sqrt(t14);
t16 = R1+t15;
t17 = log(t16);
t18 = t3+t11+t12+t13;
t19 = sqrt(t18);
t20 = R2+t19;
t21 = log(t20);
t22 = t17-t21;
t23 = t2-t11+t12+t13;
t24 = sqrt(t23);
t25 = R1+t24;
t26 = log(t25);
t27 = t3-t11+t12+t13;
t28 = sqrt(t27);
t29 = R2+t28;
t30 = log(t29);
t31 = t26-t30;
c_1 = [J.*Lsol.*t8.*pi.*(-1.0./3.0);J.*Lsol.*R1.*t5.*pi.*(-3.0./4.0e1)+J.*Lsol.*R2.*t6.*pi.*(3.0./4.0e1)+(J.*Lsol.*R1.*t2.*t4.*pi)./2.4e1-(J.*Lsol.*R2.*t3.*t4.*pi)./2.4e1;J.*Lsol.*pi.*(R1.*t2.*t5-R2.*t3.*t6).*(-1.5e1./4.48e2)-(J.*Lsol.*t7.*t8.*pi)./1.28e2+J.*Lsol.*t4.*pi.*(R1.*t5-R2.*t6).*(3.0./6.4e1);J.*Lsol.*R1.*t9.*pi.*(-1.898871527777778e-2)+J.*Lsol.*R2.*t10.*pi.*1.898871527777778e-2-J.*Lsol.*R1.*t5.*t7.*pi.*2.05078125e-2+J.*Lsol.*R2.*t6.*t7.*pi.*2.05078125e-2+J.*Lsol.*R1.*t2.*t4.*t5.*pi.*(2.5e1./5.12e2)+J.*Lsol.*R1.*t2.*t4.*t7.*pi.*1.627604166666667e-3-J.*Lsol.*R2.*t3.*t4.*t6.*pi.*(2.5e1./5.12e2)-J.*Lsol.*R2.*t3.*t4.*t7.*pi.*1.627604166666667e-3;(mu.*(J.*Lsol.*t22.*pi+J.*Lsol.*t31.*pi+J.*t22.*z.*pi.*2.0-J.*t31.*z.*pi.*2.0).*(-1.0./4.0))./pi];
end

%------------------------------------------------------------------------%  
%   f_1 = f_vector_CR(I,K,Lc,Rc,Rr,mu,z)
%   approx geometry side of system of equations with paramters for 
%   cyclinder and ring (case CR)
%
%   INPUTS
%   I: current for ring
%   K: current for cylinder
%   Lc: length of cylinder 
%   Rc: radius of cylinder 
%   Rr: radius of ring
%   mu: magnetic permeability 
%   z: z location in cylinders (0 if at center - origin is at center)
%
%   OUTPUTS
%   f1 = 5x1 vector of equations based on cylinderical shell and ring
%   paramters for the 5 equations used in case CR
%
%   REMARKS
%   equations solving for CR case [b1; b3; b5; b7; Bz];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function f_1 = f_vector_CR(I,K,Lc,Rc,Rr,mu,z)
%F_VECTOR_CR
%    F_1 = F_VECTOR_CR(I,K,LC,RC,RR,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    01-Apr-2019 08:36:36

t2 = Rr.^2;
t3 = Rc.^2;
t4 = t2.^2;
t5 = Lc.^2;
t6 = t5.^2;
t7 = t3.^2;
t8 = Lc./2.0;
t9 = t5./4.0;
t10 = z.^2;
f_1 = [I.*t2.*pi+K.*Lc.*t3.*pi;I.*t4.*pi.*(3.0./8.0)+(K.*Lc.*t3.*pi.*(t3.*3.0-t5))./8.0;I.*t2.*t4.*pi.*(1.5e1./6.4e1)+K.*Lc.*t3.*pi.*(t6+t7.*1.0e1-t3.*t5.*1.0e1).*(3.0./1.28e2);I.*t4.^2.*pi.*1.708984375e-1+K.*Lc.*t3.*pi.*(t3.*t6.*2.1e1+t3.*t7.*3.5e1-t5.*t6-t5.*t7.*7.0e1).*4.8828125e-3;(mu.*(K.*pi.*(t8+z).*1.0./sqrt(t3+t9+t10+Lc.*z).*2.0+K.*pi.*(t8-z).*1.0./sqrt(t3+t9+t10-Lc.*z).*2.0))./(pi.*4.0)+(I.*mu.*t2.*1.0./(t2+t10).^(3.0./2.0))./2.0];
end

%------------------------------------------------------------------------%  
%   Jacob1 = jacobian_matrix_CR(I,K,Lc,Rc,Rr,mu,z)
%   equations for jacobian martix for case CR
%
%   INPUTS
%   I: current for ring
%   K: current for cylinder
%   Lc: length of cylinder 
%   Rc: radius of cylinder 
%   Rr: radius of ring
%   mu: magnetic permeability 
%   z: z location in cylinders (0 if at center - origin is at center)
%
%   OUTPUTS
%   Jacob1 = 5x5 matrix of equations based on cylinderical shell and ring 
%   paramters for the 5 equations used in case CR in f_vector and 
%   derivatives with respect to the 5 parameters
%
%   REMARKS
%   case CR equations [b1; b3; b5; b7; Bz];
%   parameter order [I,K,Lc,Rc,Rr];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function Jacob1 = jacobian_matrix_CR(I,K,Lc,Rc,Rr,mu,z)
%JACOBIAN_MATRIX_CR
%    JACOB1 = JACOBIAN_MATRIX_CR(I,K,LC,RC,RR,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    01-Apr-2019 08:36:37

t2 = Rc.^2;
t3 = Rr.^2;
t4 = Lc.^2;
t5 = t2.*3.0;
t6 = t4-t5;
t7 = t3.^2;
t8 = t4.^2;
t9 = t2.^2;
t10 = t9.*1.0e1;
t12 = t2.*t4.*1.0e1;
t11 = t8+t10-t12;
t13 = t2.*t9.*3.5e1;
t14 = t2.*t8.*2.1e1;
t15 = z.^2;
t16 = Lc./2.0;
t17 = t4./4.0;
t18 = 1.0./pi;
t19 = Lc.*z;
t20 = t2+t15+t17+t19;
t21 = 1.0./sqrt(t20);
t22 = t16+z;
t23 = t16-z;
t24 = t2+t15+t17-t19;
t25 = 1.0./t20.^(3.0./2.0);
t26 = 1.0./t24.^(3.0./2.0);
t27 = t3+t15;
t28 = 1.0./t27.^(3.0./2.0);
Jacob1 = reshape([t3.*pi,t7.*pi.*(3.0./8.0),t3.*t7.*pi.*(1.5e1./6.4e1),t7.^2.*pi.*1.708984375e-1,(mu.*t3.*t28)./2.0,Lc.*t2.*pi,Lc.*t2.*t6.*pi.*(-1.0./8.0),Lc.*t2.*t11.*pi.*(3.0./1.28e2),Lc.*t2.*pi.*(t13+t14-t4.*t8-t4.*t9.*7.0e1).*4.8828125e-3,(mu.*t18.*(t23.*pi.*1.0./sqrt(t2+t15+t17-Lc.*z).*2.0+t21.*t22.*pi.*2.0))./4.0,K.*t2.*pi,K.*t2.*t4.*pi.*(-1.0./4.0)-(K.*t2.*t6.*pi)./8.0,K.*t2.*t11.*pi.*(3.0./1.28e2)-K.*Lc.*t2.*pi.*(Lc.*t2.*2.0e1-Lc.*t4.*4.0).*(3.0./1.28e2),K.*t2.*pi.*(t13+t14-t4.*t8-t4.*t9.*7.0e1).*4.8828125e-3-K.*Lc.*t2.*pi.*(Lc.*t8.*6.0+Lc.*t9.*1.4e2-Lc.*t2.*t4.*8.4e1).*4.8828125e-3,(mu.*t18.*(K.*t21.*pi+K.*1.0./sqrt(t24).*pi-K.*t22.^2.*t25.*pi-K.*t23.^2.*t26.*pi))./4.0,K.*Lc.*Rc.*pi.*2.0,K.*Lc.*Rc.*t2.*pi.*(3.0./4.0)-(K.*Lc.*Rc.*t6.*pi)./4.0,K.*Lc.*Rc.*t11.*pi.*(3.0./6.4e1)+K.*Lc.*t2.*pi.*(Rc.*t2.*4.0e1-Rc.*t4.*2.0e1).*(3.0./1.28e2),K.*Lc.*Rc.*pi.*(t13+t14-t4.*t8-t4.*t9.*7.0e1).*(5.0./5.12e2)+K.*Lc.*t2.*pi.*(Rc.*t8.*4.2e1+Rc.*t9.*2.1e2-Rc.*t2.*t4.*2.8e2).*4.8828125e-3,mu.*t18.*(K.*Rc.*t22.*t25.*pi.*2.0+K.*Rc.*t23.*t26.*pi.*2.0).*(-1.0./4.0),I.*Rr.*pi.*2.0,I.*Rr.*t3.*pi.*(3.0./2.0),I.*Rr.*t7.*pi.*(4.5e1./3.2e1),I.*Rr.*t3.*t7.*pi.*(1.75e2./1.28e2),I.*Rr.*mu.*t28-I.*Rr.*mu.*t3.*1.0./t27.^(5.0./2.0).*(3.0./2.0)],[5,5]);
end
