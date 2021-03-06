%------------------------------------------------------------------------%  
%   [K_1, K_2, Lc_1, Lc_2, Rc_1, Rc_2] = params_CC2(J_sol, L_sol, R1_sol, R2_sol)
%   solves for the cylinder parameters for 2 cylinders using the equations
%   defined by case CC
%   case CC equations [b1; b3; b5; b7; Bz; Bzd2];
%
%   INPUTS
%   J_sol: J value for solenoid
%   L_sol: length of solenoid
%   R_1: inner radius of solenoid
%   R_2: outer radius of solenoid
%
%   OUTPUTS
%   K_1: current for cylinder 1
%   K_2: current for cylinder 2
%   Lc_1: length of cylinder 1
%   Lc_2: length of cylinder 2
%   Rc_1: radius of cylinder 1
%   Rc_2: radius of cylinder 2
%
%   REMARKS
%   case CC = [b1; b3; b5; b7; Bz; Bzd2];
%   The optimal paramters of solved by using an LM algorithum
%
%   AUTHOR(S): Paige Husa
%
%                  v1.0 4/25/2019
% ----------------------------------------------------------------------- %

function [K_1, K_2, Lc_1, Lc_2, Rc_1, Rc_2, soln_found] = params_CC2(J_sol, L_sol, R1_sol, R2_sol,dispVals)
format long;
%% solenoid parameters
realS_J = J_sol; 
realS_L = L_sol; 
realS_R1 = R1_sol;
realS_R2 = R2_sol; 

delta_R = realS_R2 - realS_R1;
real_z = 0; 
real_mu = 1;%(4*pi)*1e-7;

%% solenoid side of equation
c = c_vector_CC2(realS_J, realS_L, realS_R1, realS_R2, real_mu, real_z);

%% inputs to LM function
ravg = (realS_R2+realS_R1)/2;
rdiff = (realS_R2-realS_R1)/2;
inGuess = realS_R1+rdiff/4;
outGuess = realS_R2-rdiff/4;
%initial guess the optimal paramters

cur_scale = c(1)/(pi*(realS_L*(inGuess)^2+realS_L*(outGuess)^2));

gammaGuess = [1; 1; realS_L; realS_L; inGuess; outGuess];
tol = 1e-5;

%is basically going to run forever unless you stop the algorithum or
%solution is found
maxIterations = 500;
maxTime = 60*5;
if ~exist('dispVals','var')
dispVals = 0;
end
maxChange = [];

%% LM method
[ currentEstimate, cost, soln_found ] = LM_LeastSquares( gammaGuess, @(x)costF(x,c,cur_scale), tol, maxIterations, maxTime, dispVals, maxChange );

if( any(currentEstimate <0 ))
    warning('Solution Found is non-physical. Resolving...');
    [ currentEstimate, cost, soln_found ] = LM_LeastSquares( abs(currentEstimate), @(x)costF(x,c,cur_scale), tol, maxIterations, maxTime, dispVals, maxChange );
    currentEstimate = abs(currentEstimate);
  
end
% currentEstimate
% cost

soln_found = cost < 1e-3;

%% final optimal values
K_1 = currentEstimate(1)*cur_scale;
K_2 = currentEstimate(2)*cur_scale;
Lc_1 = currentEstimate(3);
Lc_2 = currentEstimate(4);
Rc_1 = currentEstimate(5);
Rc_2 = currentEstimate(6);
end

%% cost function
function [error, Jacobian] = costF(gammaGuess,c,cur_scale)
real_z = 0; 
real_mu = 1;%(4*pi)*1e-7;

f =              f_vector_CC2(gammaGuess(1)*cur_scale, gammaGuess(2)*cur_scale, gammaGuess(3), gammaGuess(4), gammaGuess(5), gammaGuess(6), real_mu, real_z);

J_gamma = jacobian_matrix_CC2(gammaGuess(1)*cur_scale, gammaGuess(2)*cur_scale, gammaGuess(3), gammaGuess(4), gammaGuess(5), gammaGuess(6), real_mu, real_z);

J_gamma(:,1) = J_gamma(:,1)*cur_scale;
J_gamma(:,2) = J_gamma(:,2)*cur_scale;
Jacobian = J_gamma;

% error  = Jacobian*gammaGuess;
error = f -c;

error = error/norm(c);
Jacobian = Jacobian/norm(c);

%scale the jabocian appropriately so don't get a singularity
for i = 1:length(error)
   sc = norm(Jacobian(i,:));
   error(i) = error(i)/sc;
   Jacobian(i,:) = Jacobian(i,:)/sc;
end



end

%------------------------------------------------------------------------%  
%   c_1 = c_vector_CC2(J,Lsol,R1,R2,mu,z)
%   constants in equation (solenoid side of system of equations) for case (CC) where the
%   optimal solution is the super postion of 2 cylindrical shells
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
%   c_1 = 6x1 vector of constants based on solenoid paramters for the 6
%   equations used in case CC
%
%   REMARKS
%   equation solving for CC case [b1; b3; b5; b7; Bz; Bzd2;];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/25/2019
% ----------------------------------------------------------------------- %

function c_1 = c_vector_CC2(J,Lsol,R1,R2,mu,z)
%C_VECTOR_CC2
%    C_1 = C_VECTOR_CC2(J,LSOL,R1,R2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    25-Apr-2019 19:26:25

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
t32 = 1.0./pi;
t33 = z.*2.0;
t34 = Lsol+t33;
t35 = Lsol-t33;
t36 = 1.0./t16;
t37 = 1.0./sqrt(t14);
t38 = 1.0./t20;
t39 = 1.0./sqrt(t18);
t40 = t34.^2;
t41 = 1.0./t25;
t42 = 1.0./sqrt(t23);
t43 = 1.0./t29;
t44 = 1.0./sqrt(t27);
t45 = t35.^2;
t46 = t36.*t37;
t47 = 1.0./t16.^2;
t48 = 1.0./t14;
t49 = 1.0./t14.^(3.0./2.0);
t50 = 1.0./t20.^2;
t51 = 1.0./t18;
t52 = (t40.*t50.*t51)./4.0;
t53 = 1.0./t18.^(3.0./2.0);
t54 = (t38.*t40.*t53)./4.0;
t55 = t46+t52+t54-t38.*t39-(t36.*t40.*t49)./4.0-(t40.*t47.*t48)./4.0;
t56 = t41.*t42;
t57 = 1.0./t25.^2;
t58 = 1.0./t23;
t59 = 1.0./t23.^(3.0./2.0);
t60 = 1.0./t29.^2;
t61 = 1.0./t27;
t62 = (t45.*t60.*t61)./4.0;
t63 = 1.0./t27.^(3.0./2.0);
t64 = (t43.*t45.*t63)./4.0;
t65 = t56+t62+t64-t43.*t44-(t41.*t45.*t59)./4.0-(t45.*t57.*t58)./4.0;
c_1 = [J.*Lsol.*t8.*pi.*(-1.0./3.0);J.*Lsol.*R1.*t5.*pi.*(-3.0./4.0e1)+J.*Lsol.*R2.*t6.*pi.*(3.0./4.0e1)+(J.*Lsol.*R1.*t2.*t4.*pi)./2.4e1-(J.*Lsol.*R2.*t3.*t4.*pi)./2.4e1;J.*Lsol.*pi.*(R1.*t2.*t5-R2.*t3.*t6).*(-1.5e1./4.48e2)-(J.*Lsol.*t7.*t8.*pi)./1.28e2+J.*Lsol.*t4.*pi.*(R1.*t5-R2.*t6).*(3.0./6.4e1);J.*Lsol.*R1.*t9.*pi.*(-1.898871527777778e-2)+J.*Lsol.*R2.*t10.*pi.*1.898871527777778e-2-J.*Lsol.*R1.*t5.*t7.*pi.*2.05078125e-2+J.*Lsol.*R2.*t6.*t7.*pi.*2.05078125e-2+J.*Lsol.*R1.*t2.*t4.*t5.*pi.*(2.5e1./5.12e2)+J.*Lsol.*R1.*t2.*t4.*t7.*pi.*1.627604166666667e-3-J.*Lsol.*R2.*t3.*t4.*t6.*pi.*(2.5e1./5.12e2)-J.*Lsol.*R2.*t3.*t4.*t7.*pi.*1.627604166666667e-3;mu.*t32.*(J.*Lsol.*t22.*pi+J.*Lsol.*t31.*pi+J.*t22.*z.*pi.*2.0-J.*t31.*z.*pi.*2.0).*(-1.0./4.0);mu.*t32.*(J.*pi.*((t34.*t36.*t37)./2.0-(t34.*t38.*t39)./2.0).*4.0+J.*pi.*((t35.*t41.*t42)./2.0-(t35.*t43.*t44)./2.0).*4.0+J.*Lsol.*t55.*pi+J.*Lsol.*t65.*pi+J.*t55.*z.*pi.*2.0-J.*t65.*z.*pi.*2.0).*(-1.0./4.0)];
end

%------------------------------------------------------------------------%  
%   f_1 = f_vector_CC2(K_1,K_2,Lc_1,Lc_2,Rc_1,Rc_2,mu,z)
%   approx geometry side of system of equations with paramters for 2 cyclinders (case CC)
%
%   INPUTS
%   K_1: K value for cylinder 1
%   K_2: K value for cylinder 2
%   Lc_1: length of cylinder 1
%   Lc_2: length of cylinder 2
%   Rc_1: radius of cylinder 1
%   Rc_2: radius of cylinder 2
%   mu: magnetic permeability 
%   z: z location in cylinders (0 if at center - origin is at center)
%
%   OUTPUTS
%   f1 = 6x1 vector of equations based on cylinder paramters for the 6
%   equations used in case CC
%
%   REMARKS
%   equations solving for CC case [b1; b3; b5; b7; Bz; Bzd2];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/25/2019
% ----------------------------------------------------------------------- %

function f_1 = f_vector_CC2(K_1,K_2,Lc_1,Lc_2,Rc_1,Rc_2,mu,z)
%F_VECTOR_CC2
%    F_1 = F_VECTOR_CC2(K_1,K_2,LC_1,LC_2,RC_1,RC_2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    25-Apr-2019 19:26:28

t2 = Rc_1.^2;
t3 = Rc_2.^2;
t4 = Lc_1.^2;
t5 = Lc_2.^2;
t6 = t4.^2;
t7 = t2.^2;
t8 = t5.^2;
t9 = t3.^2;
t10 = Lc_1./2.0;
t11 = t4./4.0;
t12 = z.^2;
t13 = Lc_2./2.0;
t14 = t5./4.0;
t15 = t10-z;
t16 = Lc_1.*z;
t17 = t2+t11+t12+t16;
t18 = z.*2.0;
t19 = t2+t11+t12-t16;
t20 = 1.0./t19.^(3.0./2.0);
t21 = t10+z;
t22 = 1.0./t17.^(3.0./2.0);
t23 = Lc_1-t18;
t24 = Lc_1+t18;
t25 = t13-z;
t26 = Lc_2.*z;
t27 = t3+t12+t14+t26;
t28 = t3+t12+t14-t26;
t29 = 1.0./t28.^(3.0./2.0);
t30 = t13+z;
t31 = 1.0./t27.^(3.0./2.0);
t32 = Lc_2-t18;
t33 = Lc_2+t18;
f_1 = [K_1.*Lc_1.*t2.*pi+K_2.*Lc_2.*t3.*pi;(K_1.*Lc_1.*t2.*pi.*(t2.*3.0-t4))./8.0+(K_2.*Lc_2.*t3.*pi.*(t3.*3.0-t5))./8.0;K_1.*Lc_1.*t2.*pi.*(t6+t7.*1.0e1-t2.*t4.*1.0e1).*(3.0./1.28e2)+K_2.*Lc_2.*t3.*pi.*(t8+t9.*1.0e1-t3.*t5.*1.0e1).*(3.0./1.28e2);K_1.*Lc_1.*t2.*pi.*(t2.*t6.*2.1e1+t2.*t7.*3.5e1-t4.*t6-t4.*t7.*7.0e1).*4.8828125e-3+K_2.*Lc_2.*t3.*pi.*(t3.*t8.*2.1e1+t3.*t9.*3.5e1-t5.*t8-t5.*t9.*7.0e1).*4.8828125e-3;mu.*(K_1.*1.0./sqrt(t17).*t21.*pi.*2.0+K_1.*t15.*pi.*1.0./sqrt(t2+t11+t12-Lc_1.*z).*2.0).*7.957747154594767e-2+mu.*(K_2.*1.0./sqrt(t27).*t30.*pi.*2.0+K_2.*t25.*pi.*1.0./sqrt(t3+t12+t14-Lc_2.*z).*2.0).*7.957747154594767e-2;mu.*(K_1.*t15.*t20.*pi.*2.0+K_1.*t20.*t23.*pi.*2.0+K_1.*t21.*t22.*pi.*2.0+K_1.*t22.*t24.*pi.*2.0-K_1.*t15.*1.0./t19.^(5.0./2.0).*t23.^2.*pi.*(3.0./2.0)-K_1.*1.0./t17.^(5.0./2.0).*t21.*t24.^2.*pi.*(3.0./2.0)).*(-7.957747154594767e-2)-mu.*(K_2.*t25.*t29.*pi.*2.0+K_2.*t29.*t32.*pi.*2.0+K_2.*t30.*t31.*pi.*2.0+K_2.*t31.*t33.*pi.*2.0-K_2.*t25.*1.0./t28.^(5.0./2.0).*t32.^2.*pi.*(3.0./2.0)-K_2.*1.0./t27.^(5.0./2.0).*t30.*t33.^2.*pi.*(3.0./2.0)).*7.957747154594767e-2];
end

%------------------------------------------------------------------------%  
%   Jacob1 = jacobian_matrix_CC2(K_1,K_2,Lc_1,Lc_2,Rc_1,Rc_2,mu,z)
%   equations for jacobian martix for case CC
%
%   INPUTS
%   K_1: K value for cylinder 1
%   K_2: K value for cylinder 2
%   Lc_1: length of cylinder 1
%   Lc_2: length of cylinder 2
%   Rc_1: radius of cylinder 1
%   Rc_2: radius of cylinder 2
%   mu: magnetic permeability 
%   z: z location in cylinders (0 if at center - origin is at center)
%
%   OUTPUTS
%   Jacob1 = 6x6 matrix of equations based on cylinder paramters for the 6
%   equations used in case CC in f_vector and derivatives with respect to the
%   6 parameters
%
%   REMARKS
%   case CC equations [b1; b3; b5; b7; Bz; Bzd2];
%   parameter order [K_1, K_2, Lc_1, Lc_2, Rc_1, Rc_2];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/25/2019
% ----------------------------------------------------------------------- %

function Jacob1 = jacobian_matrix_CC2(K_1,K_2,Lc_1,Lc_2,Rc_1,Rc_2,mu,z)
%JACOBIAN_MATRIX_CC2
%    JACOB1 = JACOBIAN_MATRIX_CC2(K_1,K_2,LC_1,LC_2,RC_1,RC_2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    25-Apr-2019 19:26:27

t2 = Rc_1.^2;
t3 = Rc_2.^2;
t4 = Lc_1.^2;
t5 = t2.*3.0;
t6 = t4-t5;
t7 = Lc_2.^2;
t8 = t3.*3.0;
t9 = t7-t8;
t10 = t4.^2;
t11 = t2.^2;
t12 = t11.*1.0e1;
t18 = t2.*t4.*1.0e1;
t13 = t10+t12-t18;
t14 = t7.^2;
t15 = t3.^2;
t16 = t15.*1.0e1;
t19 = t3.*t7.*1.0e1;
t17 = t14+t16-t19;
t20 = t2.*t11.*3.5e1;
t21 = t2.*t10.*2.1e1;
t22 = t3.*t15.*3.5e1;
t23 = t3.*t14.*2.1e1;
t24 = Lc_1./2.0;
t25 = t4./4.0;
t26 = z.^2;
t27 = Lc_2./2.0;
t28 = t7./4.0;
t29 = Lc_1.*z;
t30 = t2+t25+t26+t29;
t31 = 1.0./sqrt(t30);
t32 = t24+z;
t33 = t24-z;
t34 = t2+t25+t26-t29;
t35 = Lc_2.*z;
t36 = t3+t26+t28+t35;
t37 = 1.0./sqrt(t36);
t38 = t27+z;
t39 = t27-z;
t40 = t3+t26+t28-t35;
t41 = 1.0./t30.^(3.0./2.0);
t42 = 1.0./t34.^(3.0./2.0);
t43 = 1.0./t36.^(3.0./2.0);
t44 = 1.0./t40.^(3.0./2.0);
t45 = z.*2.0;
t46 = Lc_1+t45;
t47 = Lc_1-t45;
t48 = Lc_2+t45;
t49 = Lc_2-t45;
t50 = t46.^2;
t51 = 1.0./t30.^(5.0./2.0);
t52 = t47.^2;
t53 = 1.0./t34.^(5.0./2.0);
t54 = t32.^2;
t55 = t33.^2;
t56 = Lc_1.*2.0;
t57 = t48.^2;
t58 = 1.0./t36.^(5.0./2.0);
t59 = t49.^2;
t60 = 1.0./t40.^(5.0./2.0);
t61 = t38.^2;
t62 = t39.^2;
t63 = z.*4.0;
t64 = Lc_2.*2.0;
t65 = 1.0./t30.^(7.0./2.0);
t66 = 1.0./t34.^(7.0./2.0);
t67 = 1.0./t36.^(7.0./2.0);
t68 = 1.0./t40.^(7.0./2.0);
Jacob1 = reshape([Lc_1.*t2.*pi,Lc_1.*t2.*t6.*pi.*(-1.0./8.0),Lc_1.*t2.*t13.*pi.*(3.0./1.28e2),Lc_1.*t2.*pi.*(t20+t21-t4.*t10-t4.*t11.*7.0e1).*4.8828125e-3,mu.*(t33.*pi.*1.0./sqrt(t2+t25+t26-Lc_1.*z).*2.0+t31.*t32.*pi.*2.0).*7.957747154594767e-2,mu.*(t32.*t41.*pi.*2.0+t33.*t42.*pi.*2.0+t41.*t46.*pi.*2.0+t42.*t47.*pi.*2.0-t32.*t50.*t51.*pi.*(3.0./2.0)-t33.*t52.*t53.*pi.*(3.0./2.0)).*(-7.957747154594767e-2),Lc_2.*t3.*pi,Lc_2.*t3.*t9.*pi.*(-1.0./8.0),Lc_2.*t3.*t17.*pi.*(3.0./1.28e2),Lc_2.*t3.*pi.*(t22+t23-t7.*t14-t7.*t15.*7.0e1).*4.8828125e-3,mu.*(t39.*pi.*1.0./sqrt(t3+t26+t28-Lc_2.*z).*2.0+t37.*t38.*pi.*2.0).*7.957747154594767e-2,mu.*(t38.*t43.*pi.*2.0+t39.*t44.*pi.*2.0+t43.*t48.*pi.*2.0+t44.*t49.*pi.*2.0-t38.*t57.*t58.*pi.*(3.0./2.0)-t39.*t59.*t60.*pi.*(3.0./2.0)).*(-7.957747154594767e-2),K_1.*t2.*pi,K_1.*t2.*t4.*pi.*(-1.0./4.0)-(K_1.*t2.*t6.*pi)./8.0,K_1.*t2.*t13.*pi.*(3.0./1.28e2)-K_1.*Lc_1.*t2.*pi.*(Lc_1.*t2.*2.0e1-Lc_1.*t4.*4.0).*(3.0./1.28e2),K_1.*t2.*pi.*(t20+t21-t4.*t10-t4.*t11.*7.0e1).*4.8828125e-3-K_1.*Lc_1.*t2.*pi.*(Lc_1.*t10.*6.0+Lc_1.*t11.*1.4e2-Lc_1.*t2.*t4.*8.4e1).*4.8828125e-3,mu.*(K_1.*t31.*pi+K_1.*1.0./sqrt(t34).*pi-K_1.*t41.*t54.*pi-K_1.*t42.*t55.*pi).*7.957747154594767e-2,mu.*(K_1.*t41.*pi.*-3.0-K_1.*t42.*pi.*3.0+K_1.*t50.*t51.*pi.*(3.0./4.0)+K_1.*t51.*t54.*pi.*3.0+K_1.*t52.*t53.*pi.*(3.0./4.0)+K_1.*t53.*t55.*pi.*3.0+K_1.*t32.*t51.*pi.*(t56+t63).*(3.0./2.0)+K_1.*t32.*t46.*t51.*pi.*3.0+K_1.*t33.*t47.*t53.*pi.*3.0-K_1.*t50.*t54.*t65.*pi.*(1.5e1./4.0)-K_1.*t52.*t55.*t66.*pi.*(1.5e1./4.0)+K_1.*t33.*t53.*pi.*(t56-z.*4.0).*(3.0./2.0)).*7.957747154594767e-2,K_2.*t3.*pi,K_2.*t3.*t7.*pi.*(-1.0./4.0)-(K_2.*t3.*t9.*pi)./8.0,K_2.*t3.*t17.*pi.*(3.0./1.28e2)-K_2.*Lc_2.*t3.*pi.*(Lc_2.*t3.*2.0e1-Lc_2.*t7.*4.0).*(3.0./1.28e2),K_2.*t3.*pi.*(t22+t23-t7.*t14-t7.*t15.*7.0e1).*4.8828125e-3-K_2.*Lc_2.*t3.*pi.*(Lc_2.*t14.*6.0+Lc_2.*t15.*1.4e2-Lc_2.*t3.*t7.*8.4e1).*4.8828125e-3,mu.*(K_2.*t37.*pi+K_2.*1.0./sqrt(t40).*pi-K_2.*t43.*t61.*pi-K_2.*t44.*t62.*pi).*7.957747154594767e-2,mu.*(K_2.*t43.*pi.*-3.0-K_2.*t44.*pi.*3.0+K_2.*t57.*t58.*pi.*(3.0./4.0)+K_2.*t58.*t61.*pi.*3.0+K_2.*t59.*t60.*pi.*(3.0./4.0)+K_2.*t60.*t62.*pi.*3.0+K_2.*t38.*t58.*pi.*(t63+t64).*(3.0./2.0)+K_2.*t38.*t48.*t58.*pi.*3.0+K_2.*t39.*t49.*t60.*pi.*3.0-K_2.*t57.*t61.*t67.*pi.*(1.5e1./4.0)-K_2.*t59.*t62.*t68.*pi.*(1.5e1./4.0)-K_2.*t39.*t60.*pi.*(t63-t64).*(3.0./2.0)).*7.957747154594767e-2,K_1.*Lc_1.*Rc_1.*pi.*2.0,K_1.*Lc_1.*Rc_1.*t2.*pi.*(3.0./4.0)-(K_1.*Lc_1.*Rc_1.*t6.*pi)./4.0,K_1.*Lc_1.*Rc_1.*t13.*pi.*(3.0./6.4e1)+K_1.*Lc_1.*t2.*pi.*(Rc_1.*t2.*4.0e1-Rc_1.*t4.*2.0e1).*(3.0./1.28e2),K_1.*Lc_1.*t2.*pi.*(Rc_1.*t10.*4.2e1+Rc_1.*t11.*2.1e2-Rc_1.*t2.*t4.*2.8e2).*4.8828125e-3+K_1.*Lc_1.*Rc_1.*pi.*(t20+t21-t4.*t10-t4.*t11.*7.0e1).*(5.0./5.12e2),mu.*(K_1.*Rc_1.*t32.*t41.*pi.*2.0+K_1.*Rc_1.*t33.*t42.*pi.*2.0).*(-7.957747154594767e-2),mu.*(K_1.*Rc_1.*t32.*t51.*pi.*6.0+K_1.*Rc_1.*t33.*t53.*pi.*6.0+K_1.*Rc_1.*t46.*t51.*pi.*6.0+K_1.*Rc_1.*t47.*t53.*pi.*6.0-K_1.*Rc_1.*t32.*t50.*t65.*pi.*(1.5e1./2.0)-K_1.*Rc_1.*t33.*t52.*t66.*pi.*(1.5e1./2.0)).*7.957747154594767e-2,K_2.*Lc_2.*Rc_2.*pi.*2.0,K_2.*Lc_2.*Rc_2.*t3.*pi.*(3.0./4.0)-(K_2.*Lc_2.*Rc_2.*t9.*pi)./4.0,K_2.*Lc_2.*Rc_2.*t17.*pi.*(3.0./6.4e1)+K_2.*Lc_2.*t3.*pi.*(Rc_2.*t3.*4.0e1-Rc_2.*t7.*2.0e1).*(3.0./1.28e2),K_2.*Lc_2.*t3.*pi.*(Rc_2.*t14.*4.2e1+Rc_2.*t15.*2.1e2-Rc_2.*t3.*t7.*2.8e2).*4.8828125e-3+K_2.*Lc_2.*Rc_2.*pi.*(t22+t23-t7.*t14-t7.*t15.*7.0e1).*(5.0./5.12e2),mu.*(K_2.*Rc_2.*t38.*t43.*pi.*2.0+K_2.*Rc_2.*t39.*t44.*pi.*2.0).*(-7.957747154594767e-2),mu.*(K_2.*Rc_2.*t38.*t58.*pi.*6.0+K_2.*Rc_2.*t39.*t60.*pi.*6.0+K_2.*Rc_2.*t48.*t58.*pi.*6.0+K_2.*Rc_2.*t49.*t60.*pi.*6.0-K_2.*Rc_2.*t38.*t57.*t67.*pi.*(1.5e1./2.0)-K_2.*Rc_2.*t39.*t59.*t68.*pi.*(1.5e1./2.0)).*7.957747154594767e-2],[6,6]);
end



