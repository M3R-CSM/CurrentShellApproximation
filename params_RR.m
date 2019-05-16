%------------------------------------------------------------------------%  
%   [Ir_1, Ir_2, Rr_1, Rr_2] = params_RR(J_sol, L_sol, R1_sol, R2_sol)
%   solves for the ring paramters for 2 rings using the equations
%   defined by case RR
%   case RR equations [b1; b3; b5; Bz];
%
%   INPUTS
%   J_sol: J value for solenoid
%   L_sol: length of solenoid
%   R_1: inner radius of solenoid
%   R_2: outer radius of solenoid
%
%   OUTPUTS
%   I_1: current for ring 1
%   I_2: current for ring 2
%   Rr_1: radius of ring 1
%   Rr_2: radius of ring 2
%
%   REMARKS
%   case RR = [b1; b3; b5; Bz];
%   The optimal paramters of solved by using an LM algorithum
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: Brandon Saunders
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function [Ir_1, Ir_2, Rr_1, Rr_2, soln_found] = params_RR(J_sol, L_sol, R1_sol, R2_sol)
format long;
%% solenoid parameters
realS_J = J_sol; 
realS_L = L_sol; 
realS_R1 = R1_sol;
realS_R2 = R2_sol; 

delta_R = realS_R2 - realS_R1;
real_z = 0; 
real_mu = 1;%(4*pi)*1e-7;

c = c_vector_RR(realS_J, realS_L, realS_R1, realS_R2, real_mu, real_z);

%% inputs to LM function
ravg = (realS_R2+realS_R1)/2;
rdiff = (realS_R2-realS_R1)/2;
% radd = realS_R2 - ravg;
inGuess = realS_R1 + rdiff/2;
outGuess = realS_R2 - rdiff/2;

% gammaGuess = [realS_J*delta_R/2; realS_J*delta_R/2; inGuess; outGuess];
cur_scale = c(1)/(pi*((inGuess)^2+(outGuess)^2));
gammaGuess = [1; 1; inGuess; outGuess];
% gamma params I_1,I_2,Rr_1,Rr_2

tol = 1e-5; 



%basically runs forever unless solution is found or manually stopped
maxIterations = 500;
maxTime = 60*5;
dispVals = 0;
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


%% final values
Ir_1 = currentEstimate(1)*cur_scale;
Ir_2 = currentEstimate(2)*cur_scale;
Rr_1 = currentEstimate(3);
Rr_2 = currentEstimate(4);

end

%% cost function
function [error, Jacobian] = costF(gammaGuess,c,scale)
real_z = 0; 
real_mu = 1;%(4*pi)*1e-7;

f =              f_vector_RR(gammaGuess(1)*scale, gammaGuess(2)*scale, gammaGuess(3), gammaGuess(4), real_mu, real_z);

J_gamma = jacobian_matrix_RR(gammaGuess(1)*scale, gammaGuess(2)*scale, gammaGuess(3), gammaGuess(4), real_mu, real_z);

J_gamma(:,1) = J_gamma(:,1)*scale;
J_gamma(:,2) = J_gamma(:,2)*scale;

Jacobian = J_gamma;

error = f -c;

error = error/norm(c);
Jacobian = Jacobian/norm(c);

%scale vlaues so no singularities
for i = 1:length(error)
   sc = norm(Jacobian(i,:));
   error(i) = error(i)/sc;
   Jacobian(i,:) = Jacobian(i,:)/sc;
end

end

%------------------------------------------------------------------------%  
%   c_1 = c_vector_RR(J,Lsol,R1,R2,mu,z)
%   constants in equation (solenoid side of system of equations) for case (RR) where the
%   optimal solution is the super postion of 2 rings
%
%   INPUTS
%   J: current for solenoid
%   Lsol: length of solenoid
%   R1: inner radius of solenoid
%   R2: outer radius of solenoid
%   mu: magnetic permeability 
%   z: z location in solenoid (0 if at center - origin is at center)
%
%   OUTPUTS
%   c_1 = 4x1 vector of constants based on solenoid paramters for the 4
%   equations used in case RR
%
%   REMARKS
%   equation solving for CR case [b1; b3; b5; Bz];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function c_1 = c_vector_RR(J,Lsol,R1,R2,mu,z)
%C_VECTOR_RR
%    C_1 = C_VECTOR_RR(J,LSOL,R1,R2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    01-Apr-2019 08:39:48

t2 = R1.^2;
t3 = R2.^2;
t4 = Lsol.^2;
t5 = t2.^2;
t6 = t3.^2;
t7 = t4.^2;
t8 = R1.*t2-R2.*t3;
t9 = Lsol.*z;
t10 = t4./4.0;
t11 = z.^2;
t12 = t2+t9+t10+t11;
t13 = sqrt(t12);
t14 = R1+t13;
t15 = log(t14);
t16 = t3+t9+t10+t11;
t17 = sqrt(t16);
t18 = R2+t17;
t19 = log(t18);
t20 = t15-t19;
t21 = t2-t9+t10+t11;
t22 = sqrt(t21);
t23 = R1+t22;
t24 = log(t23);
t25 = t3-t9+t10+t11;
t26 = sqrt(t25);
t27 = R2+t26;
t28 = log(t27);
t29 = t24-t28;
c_1 = [J.*Lsol.*t8.*pi.*(-1.0./3.0);J.*Lsol.*R1.*t5.*pi.*(-3.0./4.0e1)+J.*Lsol.*R2.*t6.*pi.*(3.0./4.0e1)+(J.*Lsol.*R1.*t2.*t4.*pi)./2.4e1-(J.*Lsol.*R2.*t3.*t4.*pi)./2.4e1;J.*Lsol.*pi.*(R1.*t2.*t5-R2.*t3.*t6).*(-1.5e1./4.48e2)-(J.*Lsol.*t7.*t8.*pi)./1.28e2+J.*Lsol.*t4.*pi.*(R1.*t5-R2.*t6).*(3.0./6.4e1);(mu.*(J.*Lsol.*t20.*pi+J.*Lsol.*t29.*pi+J.*t20.*z.*pi.*2.0-J.*t29.*z.*pi.*2.0).*(-1.0./4.0))./pi];
end

%------------------------------------------------------------------------%  
%   f_1 = f_vector_RR(I_1,I_2,Rr_1,Rr_2,mu,z)
%   approx geometry side of system of equations with paramters for 
%   2 rings (case RR)
%
%   INPUTS
%   I_1: current for ring 1
%   I_2: current for ring 2 
%   Rr_1: radius of ring 1
%   Rr_2: radius of ring 2
%   mu: magnetic permeability 
%   z: z location in cylinders (0 if at center - origin is at center)
%
%   OUTPUTS
%   f1 = 4x1 vector of equations based on the ring
%   paramters for the 4 equations used in case RR
%
%   REMARKS
%   equations solving for RR case [b1; b3; b5; Bz];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function f_1 = f_vector_RR(I_1,I_2,Rr_1,Rr_2,mu,z)
%F_VECTOR_RR
%    F_1 = F_VECTOR_RR(I_1,I_2,RR_1,RR_2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    01-Apr-2019 08:39:49

t2 = Rr_1.^2;
t3 = Rr_2.^2;
t4 = t2.^2;
t5 = t3.^2;
t6 = z.^2;
f_1 = [I_1.*t2.*pi+I_2.*t3.*pi;I_1.*t4.*pi.*(3.0./8.0)+I_2.*t5.*pi.*(3.0./8.0);I_1.*t2.*t4.*pi.*(1.5e1./6.4e1)+I_2.*t3.*t5.*pi.*(1.5e1./6.4e1);(I_1.*mu.*t2.*1.0./(t2+t6).^(3.0./2.0))./2.0+(I_2.*mu.*t3.*1.0./(t3+t6).^(3.0./2.0))./2.0];
end

%------------------------------------------------------------------------%  
%   Jacob1 = jacobian_matrix_RR(I_1,I_2,Rr_1,Rr_2,mu,z)
%   equations for jacobian martix for case RR
%
%   INPUTS
%   I_1: current for ring 1
%   I_2: current for ring 2 
%   Rr_1: radius of ring 1
%   Rr_2: radius of ring 2
%   mu: magnetic permeability 
%   z: z location in cylinders (0 if at center - origin is at center)
%
%   OUTPUTS
%   Jacob1 = 4x4 matrix of equations based on cylinderical shell and ring 
%   paramters for the 4 equations used in case RR in f_vector and 
%   derivatives with respect to the 4 parameters
%
%   REMARKS
%   case RR equations [b1; b3; b5; Bz];
%   parameter order [I_1,I_2,Rr_1,Rr_2];
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS: 
%                  v1.0 4/1/2019
% ----------------------------------------------------------------------- %

function Jacob1 = jacobian_matrix_RR(I_1,I_2,Rr_1,Rr_2,mu,z)
%JACOBIAN_MATRIX_RR
%    JACOB1 = JACOBIAN_MATRIX_RR(I_1,I_2,RR_1,RR_2,MU,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    01-Apr-2019 08:39:49

t2 = Rr_1.^2;
t3 = Rr_2.^2;
t4 = t2.^2;
t5 = t3.^2;
t6 = z.^2;
t7 = t2+t6;
t8 = 1.0./t7.^(3.0./2.0);
t9 = t3+t6;
t10 = 1.0./t9.^(3.0./2.0);
Jacob1 = reshape([t2.*pi,t4.*pi.*(3.0./8.0),t2.*t4.*pi.*(1.5e1./6.4e1),(mu.*t2.*t8)./2.0,t3.*pi,t5.*pi.*(3.0./8.0),t3.*t5.*pi.*(1.5e1./6.4e1),(mu.*t3.*t10)./2.0,I_1.*Rr_1.*pi.*2.0,I_1.*Rr_1.*t2.*pi.*(3.0./2.0),I_1.*Rr_1.*t4.*pi.*(4.5e1./3.2e1),I_1.*Rr_1.*mu.*t8-I_1.*Rr_1.*mu.*t2.*1.0./t7.^(5.0./2.0).*(3.0./2.0),I_2.*Rr_2.*pi.*2.0,I_2.*Rr_2.*t3.*pi.*(3.0./2.0),I_2.*Rr_2.*t5.*pi.*(4.5e1./3.2e1),I_2.*Rr_2.*mu.*t10-I_2.*Rr_2.*mu.*t3.*1.0./t9.^(5.0./2.0).*(3.0./2.0)],[4,4]);
end