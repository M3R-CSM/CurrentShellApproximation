% [B,G,S] = BFieldCurrentShell(p,K,L,R,p_cent,z_hat)returns the magnetic
%       field produced by a surface current using the Elliptic Integral
%       formulation.
%
%   p - the position of the point of interest (3x1) [m]
%   K - the current density [A/m]
%   L - the Length of the shell [m]
%   R - the radius of the shell [m]
%   p_cent - the center position of the shell (3x1) [m]
%   z_hat  - the axial direction of the shell (3x1)
%-
%  The function returns
%   B - The magnetic field (3x1) [T]
%   G - The magnetic gradient (3x3) [T/m]
%   S - The magnetic gradient jacobian (3x3x3) [T/m^2]
%
% Author: Andrew J Petruska
% Last Mod: 2 Aprl 2017

function [B,G,S] = BFieldCurrentShell(pos,current,length,radius,geometricCenter,axisDirection)

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

% Calculate Field based on Elliptic Integral Formulation
% See: https://en.wikipedia.org/wiki/Solenoid
% Solution on Wikipedia is off by a factor of 2 and a minus sign
zeta = zDist+[length/2;-length/2]; %ends of the CS is the POI was in the center of the CS z axis?
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

I = 1e-7*current*[1 -1];

if( rMag > 1e-5)
    Br = I*(2*sqrt(radius/rMag)*((kSq-2)./sqrt(kSq).*K+2./sqrt(kSq).*E));
    
    %maybe don't need this if statement - nope you need it
     if abs(hSq-1)<1e-5
             Pi_term = (pi./(2*sqrt(1 - kSq)))+(E+(kSq-1).*K)./(2-2.*kSq).*(rMag/radius-1);
     else
         Pi_term = (radius-rMag)/(radius+rMag)*Pi;
     end
        Bz = I*(2*zeta./sqrt((radius+rMag).^2+zeta.^2).*(K+Pi_term));
%     else
%         Bz = I*(2*zeta./sqrt((radius+rMag).^2+zeta.^2).*(K));
%         
%     end
    
    % Calculate Gradient
    dBrdK = 2*sqrt(radius/rMag)*(-(kSq-2).*E+2*(kSq-1).*K)./(kSq.*(kSq-1));
    dKdr = (radius.*((rMag + radius).^2 + zeta.^2).^(1./2).*(radius.^2 - rMag.^2 + zeta.^2))./((rMag.*radius).^(1./2).*(2.*rMag.*radius + zeta.^2 + rMag.^2 + radius.^2).^2);
    dKdz = (-(2.*rMag.*radius.*zeta)./((rMag.*radius).^(1./2).*(zeta.^2 + (rMag + radius).^2).^(3./2)));
    
    dBrdr = I*(dBrdK.*dKdr) -1/(2*rMag)*Br;
    
    dBrdz = I*(dBrdK.*dKdz);
    
    Gr = [dBrdr   0   dBrdz;  0   1/rMag*Br   0;   dBrdz  0  (-dBrdr-1/rMag*Br)];
    
    % Calculate Gradient Derivative
    dBrdKK = 2*sqrt(radius/rMag)*((4-7*kSq+kSq.^2).*E+(-4+9*kSq-5*kSq.^2).*K)./(kSq.^(3/2).*(kSq-1).^2);
    dKdrr = -(radius.^2.*((rMag + radius).^2 + zeta.^2).^(3./2).*(10.*rMag.^2.*radius.^2 + 10.*rMag.^2.*zeta.^2 + 2.*radius.^2.*zeta.^2 + 8.*rMag.*radius.^3 - 3.*rMag.^4 + radius.^4 + zeta.^4 + 8.*rMag.*radius.*zeta.^2))./(2.*(rMag.*radius).^(3./2).*(2.*rMag.*radius + zeta.^2 + rMag.^2 + radius.^2).^4);
    dKdrz = (rMag.*radius.^2.*zeta.*((rMag + radius).^2 + zeta.^2).^(3./2).*(4.*rMag.*radius + 5.*rMag.^2 - radius.^2 - zeta.^2))./((rMag.*radius).^(3./2).*(2.*rMag.*radius + zeta.^2 + rMag.^2 + radius.^2).^4);
    
    dBrdrr = I*((dBrdK.*dKdrr+dBrdKK.*dKdr.^2 -2*(1/(2*rMag)).*dBrdK.*dKdr)) +3/(4*rMag^2)*Br;
    
    dBrdrz = I*(dBrdKK.*dKdr.*dKdz+dBrdK.*dKdrz) - 1/(2*rMag)*dBrdz;
    
    Kr = [dBrdrr 0 dBrdrz;0   (1/rMag*dBrdr-1/rMag^2*Br)  0; dBrdrz 0 -(dBrdrr+(1/rMag*dBrdr-1/rMag^2*Br))];
    Kt = [0  (1/rMag*dBrdr-1/rMag^2*Br) 0; (1/rMag*dBrdr-1/rMag^2*Br) 0 1/rMag*dBrdz; 0 1/rMag*dBrdz 0];
    Kz = [Kr(3,1) 0 -Kr(1,1)-Kr(2,2); 0 Kt(2,3)    0; (-Kr(1,1)-Kr(2,2))  0  (-Kr(3,1)-Kt(2,3))];
    
    
else
    %calc dBrdr and then do a linearization
    dBrdr = -I*(radius^2*pi*1./(radius^2+zeta.^2).^(3/2));
    Br = dBrdr*rMag;
    
    %dBrdrr = 0;
    dBrzr3 = I*(45/4*radius^2*pi*zeta.*(3*radius^2-4*zeta.^2).*(1./(radius^2+zeta.^2).^(9/2)));
    dBrdrz = I*(3*radius^2*pi*zeta.*(1./(radius^2+zeta.^2)).^(5/2)) + 1/2*dBrzr3*rMag^2;
    
    dBrdz = 0+dBrdrz*rMag;
    
    Bz = I*(2*pi*zeta.*sqrt(1./(radius^2+zeta.^2)))+1/2*dBrdrz*rMag^2;
    
    Gr = [dBrdr   0   dBrdz;  0   dBrdr  0;   dBrdz  0  -2*dBrdr];
    
    Kr = [0 0 dBrdrz;0   0  0; dBrdrz 0 0];
    Kt = [0  0 0; 0 0 dBrdrz; 0 dBrdrz 0];
    Kz = [Kr(3,1) 0 -Kr(1,1)-Kr(2,2); 0 Kt(2,3)    0; (-Kr(1,1)-Kr(2,2))  0  (-Kr(3,1)-Kt(2,3))];
    
end

% Put field in world coordinates
B = R*[Br;0;Bz];


G = R*Gr*Rt;
%G = [G(1:3,1);G(2:3,2);];


S1 = R*[Kr*Rt(:,1) Kt*Rt(:,1) Kz*Rt(:,1)]*Rt;
S2 = R*[Kr*Rt(:,2) Kt*Rt(:,2) Kz*Rt(:,2)]*Rt;

%K = [KK1(1:3,1);KK1(2:3,2);KK2(2:3,2)];
S = zeros(3,3,3);
S(:,:,1) = S1;
S(:,:,2) = S2;
S(:,:,3) = [ S1(3,1) S1(3,2) -(S1(1,1)+S1(2,2));
             S1(3,2) S2(2,2) -(S1(1,2)+S2(2,2));
  -(S1(1,1)+S1(2,2)) -(S1(1,2)+S2(2,2)) -(S1(3,1)+S2(2,2))];



end