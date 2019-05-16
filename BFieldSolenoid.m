function [B, G, K] = BFieldSolenoid(pos,I,R1,R2,L,geometricCenter, axisDir)
%% BFieldSolenoid(pos,I,R1,R2,L,axisDir) returns the magnetic field B and
%%  magnetic field jacobian G at the point specified by the vector pos,
%%  given the current density I (A/meter^2), Inner Radius R1, outer
%%  Radious R2, length L, solenoid's geometric center, and a winding axis
%%  vector given by axisDir.
%%   - The magetic field is valid in all space, including inside the coils.
%%   - The magetic field jacobian is only valid outside of the coils.
%%   - All parameters are required in mks units.  The field values returned
%%     are in Tesla.
%%
%% The method is given by: Conway, J.T., "Trigonometric Integrals for the
%% Magnetic Field of the Coil of Rectangular Cross Section", IEEE Trans
%% Magn, Vol 42/5, 2006
%%
%% %% Example:  A solenoid with an inner diameter of 100mm a winding thickness
%% %%  of 10mm and a length of 100mm.  We wish to know the field 10mm infront
%% %%  of the coil when 1A is applied to 16AWG wire.
%%
%%  Length = 0.100;
%%  Rinner = 0.100;
%%  Router = Rinner + 0.010;
%%  AreaWire = 1.31e-6; % http://en.wikipedia.org/wiki/American_wire_gauge
%%  CurrentDensity = 1/AreaWire;
%%  FieldPosition = [0;0;Length/2+0.010];
%%  WindingAxis = [0;0;1];
%%  SolenoidCenter = [0;0;0];
%%  Field = BFieldSolenoid( FieldPosition, CurrentDensity, Rinner, Router, Length, SolenoidCenter, WindingAxis)
%%
%% %% Example:  Same Solenoid as above.  Wish to plot the field along the axis
%% %%    of the solenoid from 10 mm behind to 10mm infrot of the face.
%%    zPos = linspace(-Length/2-.010, Length/2+.010, 500); %% 500 points along center of the line
%%    ZMagnitude = zeros( size(zPos,2), 1); %% Initialize magnitude array
%%    for index = 1:size(zPos,2)
%%       ZMagnitude(index) = [0 0 1]*BFieldSolenoid( [0; 0; zPos(index)], CurrentDensity, Rinner, Router, Length, SolenoidCenter, WindingAxis);
%%    end
%%    figure(1)
%%    plot(zPos,ZMagnitude)
%%    title('Field Magnitude Along Axis');
%%    xlabel('Position (m)');
%%    ylabel('Field Magnitude (T)');
%%
%% %% Example:  Same Solenoid as above.  Make a field vector plot.
%%    zPos = linspace(-1.5*Length, 1.5*Length, 30); %% 30 points gird
%%    xPos = zPos;
%%    XX = zeros( size(zPos,1), size(xPos,2) );
%%    ZZ = zeros( size(XX) );
%%    U = zeros( size(XX) );
%%    W = zeros( size(XX) );
%%    for indexZ = 1:size(zPos,2)
%%      for indexX = 1:size(xPos,2)
%%       pos = [xPos(indexX);0;zPos(indexZ)];
%%       XX(indexX,indexZ) = pos(1);
%%       ZZ(indexX,indexZ) = pos(3);
%%       B = BFieldSolenoid( pos, CurrentDensity, Rinner, Router, Length, SolenoidCenter, WindingAxis);
%%       U(indexX,indexZ) = B(1);
%%       W(indexX,indexZ) = B(3);
%%      end
%%    end
%%    figure(2)
%%    quiver(XX,ZZ,U,W,1);
%%    title('Field Vector Plot');
%%    xlabel('X Position (m)');
%%    ylabel('Z Position (m)');
%%    rectangle('Position', [-Router  -Length/2  Router-Rinner Length],'FaceColor',[.5;.5;.5])
%%    rectangle('Position', [ Rinner  -Length/2  Router-Rinner Length],'FaceColor',[.5;.5;.5])
%%
%% %% Example:  Same Solenoid as above.  Make a Gradient Contour Map
%% %%  This will take some time to run..... you can reduce the resolution and
%% %%  spead it up significantly
%% zPos = linspace(-2*Length, 2*Length, 100); %% 50 points gird
%% xPos = linspace(-2*Router, 2*Router, 100);
%% XX = zeros( size(zPos,1), size(xPos,2) );
%% ZZ = zeros( size(XX) );
%% maxG = zeros( size(XX) );
%% minG = zeros( size(XX) );
%% for indexZ = 1:size(zPos,2)
%%     for indexX = 1:size(xPos,2)
%%         pos = [xPos(indexX);0;zPos(indexZ)];
%%         XX(indexX,indexZ) = pos(1);
%%         ZZ(indexX,indexZ) = pos(3);
%%         [~,G] = BFieldSolenoid( pos, CurrentDensity, Rinner, Router, Length, SolenoidCenter, WindingAxis);
%%         S = svd(G);
%%         maxG(indexX,indexZ) = S(1);
%%         minG(indexX,indexZ) = S(3);
%%     end
%% end
%%
%% figure(3)
%% subplot(2,1,1)
%% contour(XX,ZZ,log(maxG));
%% title('Field Jacobian Plot - log of Max Gradient Magnitude');
%% xlabel('X Position (m)');
%% ylabel('Z Position (m)');
%% rectangle('Position', [-Router  -Length/2  Router-Rinner Length],'FaceColor',[.5;.5;.5])
%% rectangle('Position', [ Rinner  -Length/2  Router-Rinner Length],'FaceColor',[.5;.5;.5])
%% view(2)
%% subplot(2,1,2)
%% contour(XX,ZZ,log(minG));
%% title('Field Jacobian Plot - log of Min Gradient Magnitude');
%% xlabel('X Position (m)');
%% ylabel('Z Position (m)');
%% rectangle('Position', [-Router  -Length/2  Router-Rinner Length],'FaceColor',[.5;.5;.5])
%% rectangle('Position', [ Rinner  -Length/2  Router-Rinner Length],'FaceColor',[.5;.5;.5])
%% view(2)
%%
%%
%% Written By:
%%   Andrew J Petruska   25-Aug-2014
%%   Updated: 5/2/2017 To include K vector and streeline other calcs
%%
%% #eml

% Use integral of elliptic integral solution...  This is slower in matlab
% BGK = integral(@(r) integral_helper(pos,I,L,r,geometricCenter, axisDir), R1, R2, 'ArrayValued',true);
% 
% B1 = BGK(1:3);
% G1 = BGK(4:9);
% K1 = BGK(10:end);

%make sure axisDir is normalized
axisDir = reshape(axisDir,3,1)/norm(axisDir);
% make sure vectors are column vectors...
pos = reshape(pos,3,1);
geometricCenter = reshape(geometricCenter,3,1);


%First we need to construct displacement from the point to the solenoid.
disp = pos - (geometricCenter - L/2*axisDir); % the paper defines the coordinate system from the center of one face...

%Now find the Z position...
z = disp'*axisDir;
r = norm(disp-z*axisDir);

% Define directions for later...
zHat = axisDir;
if(r > 0)
    rHat = (disp-z*axisDir)/r;
else
    % r = 0 is on axis of solenoid and this direction does not matter
    rHat = cross(rand(3,1),axisDir);
    rHat = rHat/norm(rHat);
end
%th = atan2(rHat(2),rHat(1));



% Calculate Bz, Equation (56)-(58)
if( z < 0 )
    Bz =    T2(R2,R1,r,abs(z))              - T2(R2,R1,r,abs(L-z));
    %     dBzdz = sign(z)*dT2dA(R2,R1,r,abs(z))   + sign(L-z)*dT2dA(R2,R1,r,abs(L-z)) ;
    % dBzdr = dT2dR(R2,R1,r,abs(z))           - dT2dR(R2,R1,r,abs(L-z));
elseif( 0 <= z && z <= L )
    Bz = TB(R2,R1,r) - T2(R2,R1,r,abs(z))           - T2(R2,R1,r,abs(L-z));
    %     dBzdz =          -sign(z)*dT2dA(R2,R1,r,abs(z)) + sign(L-z)*dT2dA(R2,R1,r,abs(L-z)) ;
    % dBzdr = -1       - dT2dR(R2,R1,r,abs(z))        - dT2dR(R2,R1,r,abs(L-z));
else
    Bz = T2(R2,R1,r,abs(L-z)) - T2(R2,R1,r,abs(z));
    %     dBzdz = - sign(L-z)*dT2dA(R2,R1,r,abs(L-z)) - sign(z)*dT2dA(R2,R1,r,abs(z))  ;
    % dBzdr = dT2dR(R2,R1,r,abs(L-z)) - dT2dR(R2,R1,r,abs(z)) ;
end

% Calculte Br, Equation (59)
Br =        T3(R2,R1,r,abs(L-z)) -    T3(R2,R1,r,abs(z));
dBrdr  = dT3dR(R2,R1,r,abs(L-z)) - dT3dR(R2,R1,r,abs(z));
dBrdr2 = dT3dR2(R2,R1,r,abs(L-z)) - dT3dR2(R2,R1,r,abs(z));
dBrdrz = -sign(L-z)*dT3dRA(R2,R1,r,abs(L-z)) - sign(z)*dT3dRA(R2,R1,r,abs(z));
if r > 0
    dBrdz  = -sign(L-z)*dT3dA(R2,R1,r,abs(L-z)) - sign(z)*dT3dA(R2,R1,r,abs(z));
else
    dBrdz = 0;
end

% construct B in world frame
B = Bz*zHat + Br*rHat;

% Br = Br*4*pi*1e-7*I
% Bz = Bz*4*pi*1e-7*I
% Correct scaling
B = 4*pi*1e-7*I*B; % mu0*I*B

% Construct the gradient matrix

% if r > 10*eps
%     Imllt = eye(3)-zHat*zHat';
%     drdx = rHat(1);%1/2*1/sqrt(disp'*Imllt*disp)*([1 0 0]*Imllt*disp+disp'*Imllt*[1;0;0]);
%     drdy = rHat(2);%1/2*1/sqrt(disp'*Imllt*disp)*([0 1 0]*Imllt*disp+disp'*Imllt*[0;1;0]);
%     drdz = rHat(3);%1/2*1/sqrt(disp'*Imllt*disp)*([0 0 1]*Imllt*disp+disp'*Imllt*[0;0;1]);
%     dzdx = zHat(1);
%     dzdy = zHat(2);
%     %dzdz = zHat(3);
%     drHatdx = Imllt*([1;0;0]/r - drdx*disp/r^2);
%     drHatdy = Imllt*([0;1;0]/r - drdy*disp/r^2);
%     %drHatdz = Imllt*([0;0;1]/r - drdz*disp/r^2);
% else
%     drdx = 1;
%     %drdy = [0;1;0];
%     drdy = 1;
%     dzdx = zHat(1);
%     dzdy = zHat(2);
%     %dzdz = zHat(3);
%     drHatdx = 0;%not correct but it will be multiplied by zero anyway
%     drHatdy = 0;%not correct but it will be multiplied by zero anyway
%     %drHatdz = 0;%not correct but it will be multiplied by zero anyway
% end


% if r > R1 && r < R2 && z > 0 && z < L
%     %in Coil not valid
%     dBdx = zeros(3,1);
%     dBdy = zeros(3,1);
% else
%     dBzdr = dBrdz;
%     if r > 0
%         dBdx = rHat*dBrdr*drdx + zHat*dBzdr*drdx + zHat*dBzdz*dzdx + rHat*dBrdz*dzdx + Br*drHatdx;
%         dBdy = rHat*dBrdr*drdy + zHat*dBzdr*drdy + zHat*dBzdz*dzdy + rHat*dBrdz*dzdy + Br*drHatdy;
%     else
%         dBdx = [1;0;0]*dBrdr*drdx + zHat*dBzdr*drdx + zHat*dBzdz*dzdx + rHat*dBrdz*dzdx + Br*drHatdx;
%         dBdy = [0;1;0]*dBrdr*drdy + zHat*dBzdr*drdy + zHat*dBzdz*dzdy + rHat*dBrdz*dzdy + Br*drHatdy;
%     end
%     %dBdz = rHat*dBrdr*drdz + zHat*dBzdr*drdz + zHat*dBzdz*dzdz + rHat*dBrdz*dzdz + Br*drHatdz;
% end

R = [rHat cross(zHat,rHat) zHat];
Rt = R';

if r > sqrt(eps)
    Gr = [dBrdr 0 dBrdz; 0 1/r*Br 0; dBrdz 0 -(dBrdr+1/r*Br)];
else
    Gr = [dBrdr 0 0; 0 dBrdr 0; 0 0 -2*dBrdr];
end

G =4*pi*1e-7*I*R*Gr*Rt;

% pack into vec
G = [G(1:3,1);G(2:3,2);];


if r > sqrt(eps)
    KK1 = [dBrdr2      0              dBrdrz;
        0     (1/r*dBrdr-1/r^2*Br) 0;
        dBrdrz      0              -dBrdr2-(1/r*dBrdr-1/r^2*Br)];
    
    KK2 = [0        KK1(2,2)     0;
        KK1(2,2)    0        1/r*dBrdz;
        0          1/r*dBrdz 0;];
    
    KK3 = [KK1(3,1)             0       -KK1(1,1)-KK1(2,2);
        0                    KK2(2,3)    0
        -KK1(1,1)-KK1(2,2)    0      -KK1(3,1)-KK2(2,3) ];
else
    
    KK1 = [0           0           dBrdrz;
        0           0           0;
        dBrdrz      0           0];
    
    KK2 = [0        0        0;
        0        0        dBrdrz;
        0        dBrdrz   0;];
    
    KK3 = [KK1(3,1)             0       -KK1(1,1)-KK1(2,2);
        0                    KK2(2,3)    0
        -KK1(1,1)-KK1(2,2)    0      -KK1(3,1)-KK2(2,3) ];
    
end
KK1_n = R*[KK1*Rt(:,1) KK2*Rt(:,1) KK3*Rt(:,1)]*Rt;
KK2_n = R*[KK1*Rt(:,2) KK2*Rt(:,2) KK3*Rt(:,2)]*Rt;

% KK1 = KK1*4*pi*1e-7*I
% KK2 = KK2*4*pi*1e-7*I;

K = 4*pi*1e-7*I*[KK1_n(1:3,1);KK1_n(2:3,2);KK2_n(2:3,2)];

% G = 4*pi*1e-7*I*[dBdx dBdy [dBdx(3);dBdy(3);-(dBdx(1)+dBdy(2))]];
% toc
% norm(G2-G)

end

% Define sub functions
function tB = TB(R2,R1,r)
% Equation (27)
if( r < R1 )
    tB = (R2 - R1);
elseif( R1 <= r && r <= R2 )
    tB = (R2 - r);
else
    tB = 0;
end

end
function t2 = T2(R2,R1,r,a)
% Equation (39)
fun = @(phi) a.*phi.*sin(phi)./(R2-r.*cos(phi)+X(R2,r,a,phi)).*(1+R2./X(R2,r,a,phi))...
    -a.*phi.*sin(phi)./(R1-r.*cos(phi)+X(R1,r,a,phi)).*(1+R1./X(R1,r,a,phi))...
    -atan(r.*sin(phi)/a).*(R2./X(R2,r,a,phi)-R1./X(R1,r,a,phi)).*sin(phi);

t2 = 1/2*(R2-R1)-a/2.*log( (R2 + r + X(R2,r,a,pi))/(R1 + r + X(R1,r,a,pi)) );


if r*a >= 0 && r ~= 0
    if verLessThan('matlab', '8')
        t2 = t2 + r/(2*pi)*quad(fun ,0,pi,1e-4);
    else
        t2 = t2 + r/(2*pi)*integral(fun ,0,pi,'AbsTol',1e-4,'RelTol',0);
    end
end
end
function t3 = T3(R2,R1,r,a)
% Equation (35)

fun1 = @(phi) 1/(2*pi) * cos(phi).*(X(R2,r,a,phi)-X(R1,r,a,phi));
fun2 = @(phi)-r^2/(4*pi) * (phi + sin(phi).*cos(phi)).*sin(phi)./(R2-r.*cos(phi)+X(R2,r,a,phi)).*(1+R2./X(R2,r,a,phi))...
    +r^2/(4*pi) * (phi + sin(phi).*cos(phi)).*sin(phi)./(R1-r.*cos(phi)+X(R1,r,a,phi)).*(1+R1./X(R1,r,a,phi));

t3 = r/(4)*log( (R2 + r + X(R2,r,a,pi))/(R1 + r + X(R1,r,a,pi)) );
if a >= 0
    if verLessThan('matlab', '8')
        t3 = t3 + quad( fun1,0, pi);
        if r > 0
            t3 = t3 + quad(fun2, 0,pi,1e-4);
        end
    else
        t3 = t3 + integral( fun1,0, pi);
        if r > 0
            t3 = t3 + integral(fun2, 0,pi,'AbsTol',1e-4,'RelTol',0);
        end
    end
end

end
function x = X(R,r,a,phi)
% Equation (32)
x = sqrt(R.^2 + r.^2 + a.^2 - 2.*R.*r.*cos(phi));
end

function dx = dXdr(R,r,a,phi)
dx = (r-R*cos(phi))./X(R,r,a,phi);
end
function dx = dXda(R,r,a,phi)
dx = a./X(R,r,a,phi);
end

function dt3dr = dT3dR(R2,R1,r,a)

dfun = @(phi) (r.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2))) - cos(phi).*((5734161139222659.*(2.*r - 2.*R1.*cos(phi)))./(72057594037927936.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) - (5734161139222659.*(2.*r - 2.*R2.*cos(phi)))./(72057594037927936.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2))) - (r.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2))) + (r.^2.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2) - (r.^2.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2) - (R1.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R1.*cos(phi)))./(8.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R2.*cos(phi)))./(8.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2));
dt3dr = (r.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) - r.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + log((R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))./(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))./(4.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2));


if r>0
    if verLessThan('matlab','8')
        dt3dr = dt3dr+quad( dfun, 0, pi,1e-4);
    else
        dt3dr = dt3dr+integral( dfun, 0, pi,'AbsTol',1e-4,'RelTol',0);
    end
else
    dt3dr = dt3dr+(5734161139222659*pi*(R1/(R1^2 + a^2)^(1/2) - R2/(R2^2 + a^2)^(1/2)))/72057594037927936;
end


if isnan(dt3dr)
    dt3dr = (r.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) - r.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + log((R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))./(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))./(4.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2));
    if verLessThan('matlab','8')
        dt3dr = dt3dr+quad( dfun, 1e-3, pi-1e-3,1e-4);
    else
        dt3dr = dt3dr+integral( dfun, 1e-3, pi-1e-3,'AbsTol',1e-4,'RelTol',0);
    end
end

end

function dt3dr2 = dT3dR2(R2,R1,r,a)

dfun = @(phi) (sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2))) - cos(phi).*(5734161139222659./(36028797018963968.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) - 5734161139222659./(36028797018963968.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)) - ((5734161139222659.*r - 5734161139222659.*R1.*cos(phi)).*(2.*r - 2.*R1.*cos(phi)))./(72057594037927936.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + ((5734161139222659.*r - 5734161139222659.*R2.*cos(phi)).*(2.*r - 2.*R2.*cos(phi)))./(72057594037927936.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2))) - (sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2))) - (R1.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2)) + (r.^2.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2)./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^3) - (r.^2.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2)./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^3) + (r.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)))./(pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2) - (r.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2) - (r.^2.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(R1.^2./2 - (R1.^2.*cos(2.*phi))./2 + a.^2))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (r.^2.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(R2.^2./2 - (R2.^2.*cos(2.*phi))./2 + a.^2))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2)) - (R1.*r.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R1.*cos(phi)))./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*r.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R2.*cos(phi)))./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2)) + (3.*R1.*r.^2.*sin(phi).*(r - R1.*cos(phi)).^2.*(phi + sin(2.*phi)./2))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(5./2)) - (3.*R2.*r.^2.*sin(phi).*(r - R2.*cos(phi)).^2.*(phi + sin(2.*phi)./2))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(5./2)) - (R1.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R1.*cos(phi)).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R2.*cos(phi)).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2));
dt3dr2 = R2.^2./(2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) - R1.^2./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) - a.^2./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) + a.^2./(2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) - r.^2./(4.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) + r.^2./(4.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) - (3.*R1.*r)./(4.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) + (3.*R2.*r)./(4.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2));


if r>0
    if verLessThan('matlab','8')
        dt3dr2 = dt3dr2+quad( dfun, 0, pi,1e-4);
    else
        dt3dr2 = dt3dr2+integral( dfun, 0, pi,'AbsTol',1e-4,'RelTol',0);
    end
else
    dt3dr2 = dt3dr2+(5734161139222659*pi*(R1/(R1^2 + a^2)^(1/2) - R2/(R2^2 + a^2)^(1/2)))/72057594037927936;
end


if isnan(dt3dr2)
    dt3dr2 = R2.^2./(2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) - R1.^2./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) - a.^2./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) + a.^2./(2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) - r.^2./(4.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) + r.^2./(4.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) - (3.*R1.*r)./(4.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) + (3.*R2.*r)./(4.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2));
    if verLessThan('matlab','8')
        dt3dr2 = dt3dr2+quad( dfun, 1e-3, pi-1e-3,1e-4);
    else
        dt3dr2 = dt3dr2+integral( dfun, 1e-3, pi-1e-3,'AbsTol',1e-4,'RelTol',0);
    end
end

end

function dt3da = dT3dA(R2,R1,r,a)

dfun = @(phi)(a.*r.^2.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)) - (a.*r.^2.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) - cos(phi).*((5734161139222659.*a)./(36028797018963968.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) - (5734161139222659.*a)./(36028797018963968.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2))) - (R1.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2));

dt3da = (r.*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)));

if verLessThan('matlab','8')
    
    if  a*r > 0
        dt3da = dt3da + quad( dfun, 0, pi,1e-4);
    end
    if isnan(dt3da)
        dt3da = (r.*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)));
        dt3da = dt3da + quad( dfun, 1e-3, pi-1e-3,1e-4);
    end
else
    
    if  a*r > 0
        dt3da = dt3da + integral( dfun, 0, pi,'AbsTol',1e-4,'RelTol',0);
    end
    if isnan(dt3da)
        dt3da = (r.*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)));
        dt3da = dt3da + integral( dfun, 1e-3, pi-1e-3,'AbsTol',1e-4,'RelTol',0);
    end
end

end

function dt3dra = dT3dRA(R2,R1,r,a)

dfun = @(phi)cos(phi).*((5734161139222659.*a.*(2.*r - 2.*R1.*cos(phi)))./(72057594037927936.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) - (5734161139222659.*a.*(2.*r - 2.*R2.*cos(phi)))./(72057594037927936.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2))) - (R1.*a.*r.*sin(phi).*(phi + sin(2.*phi)./2))./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*a.*r.*sin(phi).*(phi + sin(2.*phi)./2))./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2)) - (a.*r.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) + (a.*r.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2))./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)) + (a.*r.^2.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R1.*cos(phi)))./(8.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) - (a.*r.^2.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R2.*cos(phi)))./(8.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2)) - (a.*r.^2.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)))./(2.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^3.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) + (a.*r.^2.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(2.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^3.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)) + (R1.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R1.*cos(phi)))./(8.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^2) + (3.*R1.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R1.*cos(phi)))./(8.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(5./2)) - (R2.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R2.*cos(phi)))./(8.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^2) - (3.*R2.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(2.*r - 2.*R2.*cos(phi)))./(8.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(5./2)) - (R1.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)))./(4.*pi.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2.*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*a.*r.^2.*sin(phi).*(phi + sin(2.*phi)./2).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(4.*pi.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2.*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2));

dt3dra = ((a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) + (r.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1).*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) - (r.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*((a.*(2.*R2 + 2.*r))./(2.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) + (a.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) + (a.*((R2 + r)./(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + 1))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)) - (a.*(2.*R1 + 2.*r).*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./(2.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) - (2.*a.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1).*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^3.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) - (r.*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*((4.*R2 + 4.*r)./(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + 4).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(16.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)).^2);

if verLessThan('matlab','8')
    
    if  a*r > 0
        dt3dra = dt3dra + quad( dfun, 0, pi,1e-4);
    else
        dt3dra = dt3dra-(5734161139222659*a*pi*(R1/(R1^2 + a^2)^(3/2) - R2/(R2^2 + a^2)^(3/2)))/72057594037927936;
    end
    if isnan(dt3da)
        dt3dra = ((a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) + (r.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1).*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) - (r.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*((a.*(2.*R2 + 2.*r))./(2.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) + (a.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) + (a.*((R2 + r)./(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + 1))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)) - (a.*(2.*R1 + 2.*r).*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./(2.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) - (2.*a.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1).*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^3.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) - (r.*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*((4.*R2 + 4.*r)./(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + 4).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(16.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)).^2);
        dt3dra = dt3dra + quad( dfun, 1e-3, pi-1e-3,1e-4);
    end
else
    
    if  a*r > 0
        dt3dra = dt3dra + integral( dfun, 0, pi,'AbsTol',1e-4,'RelTol',0);
    else
        dt3dra = dt3dra-(5734161139222659*a*pi*(R1/(R1^2 + a^2)^(3/2) - R2/(R2^2 + a^2)^(3/2)))/72057594037927936;
        
    end
    if isnan(dt3dra)
        dt3dra = ((a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) + (r.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1).*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) - (r.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*((a.*(2.*R2 + 2.*r))./(2.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(3./2)) + (a.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) + (a.*((R2 + r)./(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + 1))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)) - (a.*(2.*R1 + 2.*r).*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./(2.*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(3./2)) - (2.*a.*((R1 + r)./(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) + 1).*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^3.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))))./(4.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2))) - (r.*(a./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)) - (a.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./((R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)).^2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2))).*((4.*R2 + 4.*r)./(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2) + 4).*(R1 + r + (2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2)))./(16.*(R2 + r + (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)).^2);
        dt3dra = dt3dra + integral( dfun, 1e-3, pi-1e-3,'AbsTol',1e-4,'RelTol',0);
    end
end

end
%
% function dt2dr = dT2dR(R2,R1,r,a)
%
% % der of Equation (39)
% fun = @(phi) atan2((r.*sin(phi)),a).*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) - R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)) - r.*(atan2((r.*sin(phi)),a).*sin(phi).*((R1.*(r - R1.*cos(phi)))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2) - (R2.*(r - R2.*cos(phi)))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2)) - (a.*sin(phi).^2.*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) - R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(r.^2.*sin(phi).^2 + a.^2) + (a.*phi.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1).*(cos(phi) - (r - R1.*cos(phi))./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)))./(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).^2 - (a.*phi.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1).*(cos(phi) - (r - R2.*cos(phi))./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)))./(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).^2 - (R1.*a.*phi.*sin(phi).*(2.*r - 2.*R1.*cos(phi)))./(2.*(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)).*(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(3./2)) + (R2.*a.*phi.*sin(phi).*(2.*r - 2.*R2.*cos(phi)))./(2.*(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2)).*(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(3./2))) - (a.*phi.*sin(phi).*(R1./(R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2) + 1))./(R1 - r.*cos(phi) + (R1.^2 + a.^2 + r.^2 - 2.*R1.*r.*cos(phi)).^(1./2)) + (a.*phi.*sin(phi).*(R2./(R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2) + 1))./(R2 - r.*cos(phi) + (R2.^2 + a.^2 + r.^2 - 2.*R2.*r.*cos(phi)).^(1./2));
%
% dt2dr = -(a.*((2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) - (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2));
%
% if verLessThan('matlab','8')
%     dt2dr = dt2dr + 1/(2*pi)*quad(fun ,0,pi,1e-4);
%
%     if isnan(dt2dr)
%         dt2dr = -(a.*((2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) - (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2));
%         dt2dr = dt2dr + 1/(2*pi)*quad(fun ,1e-3,pi-1e-3,1e-4);
%     end
% else
%     dt2dr = dt2dr + 1/(2*pi)*integral(fun ,0,pi,'AbsTol',1e-4,'RelTol',0);
%
%     if isnan(dt2dr)
%         dt2dr = -(a.*((2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2) - (2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2)))./(2.*(2.*R1.*r + R1.^2 + a.^2 + r.^2).^(1./2).*(2.*R2.*r + R2.^2 + a.^2 + r.^2).^(1./2));
%         dt2dr = dt2dr + 1/(2*pi)*integral(fun ,1e-3,pi-1e-3,'AbsTol',1e-4,'RelTol',0);
%     end
% end
%
% end


% function [val] = integral_helper(pos,current,length,radius,geometricCenter, axisDirection)
% 
% [B,G,K] = BFieldCurrentShell(pos,current,length,radius,geometricCenter,axisDirection);
% 
% val = [B;G;K];
% end