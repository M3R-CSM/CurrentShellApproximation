function [k,e,p] = ellipkep(m,alSq,tol)
%ELLIPKE Complete elliptic integral.
%   [K,E,Pi] = ELLIPKE(M,A) returns the value of the complete elliptic
%   integrals of the first, second, and third kinds, evaluated for each
%   element of M and A.  As currently implemented, M is limited to 0 <= M <= 1.
%   
%   [K,E,Pi] = ELLIPKE(M,A,TOL) computes the complete elliptic integrals to
%   the accuracy TOL instead of the default TOL = EPS(CLASS(M)).
%
%   Some definitions of the complete elliptic integrals use the modulus
%   k instead of the parameter M.  They are related by M = k^2.
%
%   Class support for input M:
%      float: double, single
%

%   Modified to include the second kind by Bjorn Bonnevier
%   from the Alfven Laboratory, KTH, Stockholm, Sweden
%   Copyright 1984-2013 The MathWorks, Inc. 

%   ELLIPKE uses the method of the arithmetic-geometric mean
%   described in [1].

%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, 17.6.

if nargin<1
  error(message('MATLAB:ellipke:NotEnoughInputs')); 
end

classin = superiorfloat(m);

if nargin<3, tol = 100*eps(classin); end
if ~isreal(m) || ~isreal(tol) || ~isreal(alSq)
    error(message('MATLAB:ellipke:ComplexInputs'))
end
if isempty(m), k = zeros(size(m),classin); e = k; return, end
if any(m(:) < 0) || any(m(:) > 1)
  error(message('MATLAB:ellipke:MOutOfRange'));
end
if ~isscalar(tol) || tol < 0 || ~isfinite(tol)
  error(message('MATLAB:ellipke:NegativeTolerance'));
end

a0 = 1;
b0 = sqrt(1-m);
c0 = NaN;
p0 = sqrt(1-alSq);
Q0 = 1;
Qsum = 1;
s0 = m;
i1 = 0; mm = Inf;
% count = 0;
while mm > tol
%     count = count +1;
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    p1 = (p0.^2+b1.^2)./(2*p0);
    e0 = (p0.^2-b1.^2)./(p0.^2+b1.^2);
    Q1 = 1/2.*Q0.*e0;
    Qsum = Qsum+Q1;
    
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    
    %Qfac = alSq./(1-alSq).*Qsum;
    
    mm = max(max(w1(:)),max(abs(Q1(:))));
    
    % test for stagnation (may happen for TOL < machine precision)
    if isequal(c0, c1) && max(abs(Q1))<tol && mm >tol
        error(message('MATLAB:ellipke:FailedConvergence'));
    end
    
    s0 = s0 + w1;  
    a0 = a1;  b0 = b1;  c0 = c1; p0 = p1; Q0 = Q1;
end
% disp(['Counts: ' int2str(count)])
k = pi./(2*a1);
e = k.*(1-s0/2);
p = 1/2*k.*(2+alSq./(1-alSq).*Qsum);
im = find(m==1);
if ~isempty(im)
    e(im) = ones(length(im),1);
    k(im) = inf;
end
