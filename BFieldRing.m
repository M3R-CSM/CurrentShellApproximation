%------------------------------------------------------------------------%
%   [B,G,K] = BFieldRing(pos,current,radius,geometricCenter,axisDirection)
%   returns the magnetic field produced by a ring of charge using
%   the Elliptic Integral formulation.
%
%   INPUTS
%   pos: the position of the point of interest (3x1) [m]
%   current: the current density [A/m]
%   radius: the radius of the shell [m]
%   geometricCenter: the center position of the shell (3x1) [m]
%   axisDirection: the axial direction of the shell (3x1)
%
%   OUTPUTS
%   B: The magnetic field (3x1) [T]
%   G: NOT CALCULATED IN THIS FUNCTION - output is 0 (The magnetic
%   gradient (5x1) [T/m])
%   K: NOT CALCULATED IN THIS FUNCTION - output is 0 (The magnetic gradient
%   jacobian (7x1) [T/m^2])
%
%   REMARKS
%   The calculation of Br and Bz are not of the same form that A. Petruska
%   uses in the function BFieldCurrentShell.  The equations were developed
%   from a different paper.
%   Paper used: "simple analytic expressions for the magnetic field of a
%   curcilar current loop" by J. C. Simpson
%   When the Bfield is ploted, the BFieldCurrentShell function approaches
%   the BFieldRing as the length of the current shell approaches 0.
%
%   The lines of code used to calculate G and K in the BFieldCurrentShell
%   function were left in this script but commented out with the thought that
%   they could be modified to fit the ring if needed later, but now we only
%   have use for the output of B
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS:
%                  v1.0 4/2/2019
% ----------------------------------------------------------------------- %

function [B,G,K] = BFieldRing(pos,current,radius,geometricCenter,axisDirection)
% Find Displacement
disp = pos-geometricCenter;

% Normalize Axis Direction
zHat = reshape(axisDirection/norm(axisDirection),3,1);

% Find Radial Distance and z Distance
zDist = zHat'*disp; %z component center to pos
rVec = disp-zDist*zHat; %vector from center to pos in the x and y direction
rMag = norm(rVec); %magnitude of the rVec

rHat = rVec/rMag; %unit vector of rVec
if any(isnan(rHat))
    rHat = cross(rand(3,1),zHat);
    rHat = rHat/norm(rHat);
end

% rHat = rVec/rMag;

% Define Coordiante Rotation Matrix
R = [rHat cross(zHat,rHat) zHat];
Rt = R';

zeta = zDist;
hSq = 4*radius*rMag/(radius+rMag)^2;
kSq = 4*radius*rMag./((radius+rMag).^2+zeta.^2);

% if abs(hSq-1)<5e-7
%     warning('BFieldCurrentShell:Alpha Squared close to unit, ellipticPI poorly defined');
% end
% if hSq>1
%     warning('BFieldCurrentShell:Alpha Squared greator than 1, ellipticPI poorly defined');
% end

if abs(hSq-1)<1e-5
    [K,E] = ellipke(kSq);
else
    [K,E,Pi] = ellipkep(kSq,hSq);
    
end
z=zDist;
I = 1e-7*current*[1 -1];
alphaSq = radius^2 + rMag^2 + z^2 - 2*radius*rMag;
betaSq = radius^2 + rMag^2 + z^2 + 2*radius*rMag;
mu0 = 4e-7 * pi;
C = mu0 * current / pi;

if( rMag > 1e-5)
    Br = ((C * z)/(2*alphaSq*sqrt(betaSq)*rMag)) * ((radius^2 + rMag^2 + z^2)*E - alphaSq * K);
    
    Bz = (C/(2 * alphaSq*sqrt(betaSq))) * ((radius^2 - rMag^2 - z^2)*E + alphaSq * K);
    
    % Calculate Gradient
    %     dBrdK = 2*sqrt(radius/rMag)*(-(kSq-2).*E+2*(kSq-1).*K)./(kSq.*(kSq-1));
    %     dKdr = (radius.*((rMag + radius).^2 + zeta.^2).^(1./2).*(radius.^2 - rMag.^2 + zeta.^2))./((rMag.*radius).^(1./2).*(2.*rMag.*radius + zeta.^2 + rMag.^2 + radius.^2).^2);
    %     dKdz = (-(2.*rMag.*radius.*zeta)./((rMag.*radius).^(1./2).*(zeta.^2 + (rMag + radius).^2).^(3./2)));
    
    %     dBrdr = I*(dBrdK.*dKdr) -1/(2*rMag)*Br;
    
    %     dBrdz = I*(dBrdK.*dKdz);
    dBrdr = ((-C.*z).*((radius^6 + (rMag^2 + z^2)^2*(2*rMag^2 + z^2) + radius^4*(3*z^2 - 8*rMag^2) + radius^2 .* (5.*rMag^2 - 4.*rMag^2.*z^2 + 3.*z^4)).* E - alphaSq.*(radius^4 - 3*radius^2 * rMag^2 + 2*rMag^4 + (2*radius^2 + 3.*rMag^2).*z^2 + z^4).*K)) / (2*rMag^2.* alphaSq^2 .* (betaSq)^(3/2));

    a1 = (radius^6 + (rMag^2 + z^2)^2*(2*rMag^2 + z^2) + radius^4*(3*z^2 - 8*rMag^2) + radius^2 .* (5.*rMag^2 - 4.*rMag^2.*z^2 + 3.*z^4));
    a2 = (radius^4 - 3*radius^2 * rMag^2 + 2*rMag^4 + (2*radius^2 + 3.*rMag^2).*z^2 + z^4);
    a3 = (2*rMag^2.* alphaSq^2 .* (betaSq)^(3/2));
    da1dr = 
    da2dr = 
    da3dr = 
    dEdk = 
    dKdk
    dkdr = 
    
    dBrdrr = ((-C.*z).*(da1dr* E + a1* dEdk*dkdr - ddr_alphaSq.*a2.*K - alphaSq.*da2dr.*K - alphaSq.*a2.*dKdkdkdr)) / a3 -...
        ((-C.*z).*(a1* E  - alphaSq.*a2.*K))/a3^2*da3dr ;

        
    dBrdz = (C .*(((radius.^2 + rMag.^2) .* (z.^4 + (radius.^2 - rMag.^2).^2) + 2 .* z.^2 *(radius.^4 - 6.*radius.^2.*rMag + rMag.^4)).*E - alphaSq .* ((radius.^2 - rMag.^2).^2 + (radius.^2 + rMag.^2).* z.^2)*K)) / (2 .* rMag .* alphaSq^2 .* betaSq^(3/2));
    Gr = [dBrdr   0   dBrdz;  0   1/rMag*Br   0;   dBrdz  0  (-dBrdr-1/rMag*Br)];
    
    % Calculate Gradient Derivative
    %     dBrdKK = 2*sqrt(radius/rMag)*((4-7*kSq+kSq.^2).*E+(-4+9*kSq-5*kSq.^2).*K)./(kSq.^(3/2).*(kSq-1).^2);
    %     dKdrr = -(radius.^2.*((rMag + radius).^2 + zeta.^2).^(3./2).*(10.*rMag.^2.*radius.^2 + 10.*rMag.^2.*zeta.^2 + 2.*radius.^2.*zeta.^2 + 8.*rMag.*radius.^3 - 3.*rMag.^4 + radius.^4 + zeta.^4 + 8.*rMag.*radius.*zeta.^2))./(2.*(rMag.*radius).^(3./2).*(2.*rMag.*radius + zeta.^2 + rMag.^2 + radius.^2).^4);
    %     dKdrz = (rMag.*radius.^2.*zeta.*((rMag + radius).^2 + zeta.^2).^(3./2).*(4.*rMag.*radius + 5.*rMag.^2 - radius.^2 - zeta.^2))./((rMag.*radius).^(3./2).*(2.*rMag.*radius + zeta.^2 + rMag.^2 + radius.^2).^4);
    
    %     dBrdrr = I*((dBrdK.*dKdrr+dBrdKK.*dKdr.^2 -2*(1/(2*rMag)).*dBrdK.*dKdr)) +3/(4*rMag^2)*Br;
    
    %     dBrdrz = I*(dBrdKK.*dKdr.*dKdz+dBrdK.*dKdrz) - 1/(2*rMag)*dBrdz;
    
        Kr = [dBrdrr 0 dBrdrz;0   (1/rMag*dBrdr-1/rMag^2*Br)  0; dBrdrz 0 -(dBrdrr+(1/rMag*dBrdr-1/rMag^2*Br))];
    %     Kt = [0  (1/rMag*dBrdr-1/rMag^2*Br) 0; (1/rMag*dBrdr-1/rMag^2*Br) 0 1/rMag*dBrdz; 0 1/rMag*dBrdz 0];
    %     Kz = [Kr(3,1) 0 -Kr(1,1)-Kr(2,2); 0 Kt(2,3)    0; (-Kr(1,1)-Kr(2,2))  0  (-Kr(3,1)-Kt(2,3))];
    
    
else %are close to the center
    Br = ((C * z)/(2*alphaSq*sqrt(betaSq)*rMag)) * ((radius^2 + rMag^2 + z^2)*E - alphaSq * K);
    
    %     dBrdrr = 0;
    %     dBrzr3 = I*(45/4*radius^2*pi*zeta.*(3*radius^2-4*zeta.^2).*(1./(radius^2+zeta.^2).^(9/2)));
    %     dBrdrz = I*(3*radius^2*pi*zeta.*(1./(radius^2+zeta.^2)).^(5/2)) + 1/2*dBrzr3*rMag^2;
    %
    %     dBrdz = 0+dBrdrz*rMag;
    
    %     Bz = I*(2*pi*zeta.*sqrt(1./(radius^2+zeta.^2)))+1/2*dBrdrz*rMag^2;
    Bz = (mu0 * current * radius^2) / (2*(radius^2 + z^2)^(3/2));
    
%     Gr = [dBrdr   0   dBrdz;  0   dBrdr  0;   dBrdz  0  -2*dBrdr];
    
    %     Kr = [0 0 dBrdrz;0   0  0; dBrdrz 0 0];
    %     Kt = [0  0 0; 0 0 dBrdrz; 0 dBrdrz 0];
    %     Kz = [Kr(3,1) 0 -Kr(1,1)-Kr(2,2); 0 Kt(2,3)    0; (-Kr(1,1)-Kr(2,2))  0  (-Kr(3,1)-Kt(2,3))];
    
end

% Br1 = cos(angle)*Br;
% Br2 = sin(angle)*Br;
% Br = [Br1 ; Br2 ; 0];
% Bz = [0;Bz;0];
% Put field in world coordinates
B = R*[Br;0;Bz];


% G = R*Gr*Rt;
% G = [G(1:3,1);G(2:3,2);];
G= [0;0;0];

% KK1 = R*[Kr*Rt(:,1) Kt*Rt(:,1) Kz*Rt(:,1)]*Rt;
% KK2 = R*[Kr*Rt(:,2) Kt*Rt(:,2) Kz*Rt(:,2)]*Rt;

% K = [KK1(1:3,1);KK1(2:3,2);KK2(2:3,2)];
K = [0;0;0];
end