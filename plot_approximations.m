%------------------------------------------------------------------------%
%   Script to solve for the parameters of an approximate geometry to represent a solenoid
%   given the geometry of a solenoid. 
%   Three different solution cases (CC, CR, RR).  
%   Then use the paramters of the approximate geometry to calculate the Bfield and plot the
%   error when compared to the real solenoid.
%
%   REMARKS
%   cc solves for paramters of 2 cylindrical shell solution
%   cr solves for the paramters of cylindrical shell and ring solution
%   rr solves for the paramters of 2 ring soltions
%
%   A C solution case and R solution case can be written using the same
%   pattern and the params_C or params_R function.
%
%   AUTHOR(S): Paige Husa
%
%   MODIFICATIONS:
%                  v1.0 5/15/2019
% ----------------------------------------------------------------------- %
function plot_approximations(J,L,R1,R2, type)
%%
% clear;
% clc;
% close all;

%% define the solenoid geometry
%J = 1; %current density
% L = 5; %length
% 
% ratioX = 0.15; %length to inner radius ratio
% ratioY = 0.6; %length to thickness ratio
% 
% %define radii from the ratios
% R1 = L/ratioX;
% R2 = (L / ratioY) + R1;

%or define radii directly
% R1 = 1;
% R2 = 2;

%% define which solution geometry you want to use to approximate the solenoid
cc = all(type == 'cc');
cr = all(type == 'cr');
rr = all(type == 'rr');

%% offest value to plot
offsetValue = [0.1, 0.5]; %10% and 50% offset will be plotted on the graph
O = length(offsetValue);

minSphere = sqrt(R2^2+L^2/4);

dist = minSphere * offsetValue(2);
newSphere = minSphere + dist;
L_50 = L/2 + dist;
R_50 = R2 + dist;
endDim = max([L_50 R_50]); %finds largest dimension the magnetic field needs to be calculated at

%% Spacing of data points
%this may need to be adjusted based on what geometry you are trying to plot
N = 50;
X = linspace(0,endDim*1.5,N);
Y = linspace(0,endDim*1.5,N);

%% CC solution
if cc == true
    %solves for parameters of CC solution
    [K_1, K_2, Lc_1, Lc_2, Rc_1, Rc_2] = params_CC2(J,L,R1,R2)
    
    %b-field
    Bcc = BFieldCurrentShell([0;0;0],K_1,Lc_1,Rc_1,[0;0;0],[0;1;0]) + BFieldCurrentShell([0;0;0],K_2,Lc_2,Rc_2,[0;0;0],[0;1;0]);
    
    %params to be inside coils
    inSolenoid = @(x,y) y<=L/2 && y>= -L/2 && (x <=R2 && x>=R1 || x <=-R1 && x>=-R2);
    
    BccX = zeros(N,N);
    BccY = zeros(N,N);
    BX = zeros(N,N);
    BY = zeros(N,N);
    Pe = zeros(N,N);
    XX = zeros(N,N);
    YY = zeros(N,N);
    
    for i=1:length(X) %for all x positions in the grid
        for j = 1:length(Y) %for all y positions in the grid
            XX(i,j) = X(i);
            YY(i,j) = Y(j);
            if ~inSolenoid(X(i),Y(j))
                %b-field of CC
                Bcc = BFieldCurrentShell([X(i);Y(j);0],K_1,Lc_1,Rc_1,[0;0;0],[0;1;0]) + BFieldCurrentShell([X(i);Y(j);0],K_2,Lc_2,Rc_2,[0;0;0],[0;1;0]);
                
                %b-field of solenoid
                B = BFieldSolenoid([X(i);Y(j);0],J,R1,R2,L,[0;0;0],[0;1;0]);
                BccX(i,j) = Bcc(1);
                BccY(i,j) = Bcc(2);
                BX(i,j) = B(1);
                BY(i,j) = B(2);
                
                %calculate error
                Pe(i,j) = (norm(B-Bcc)/norm(B))*100;
            end
        end
    end
end

%% CR
if cr == true
    %solves for parameters of CR solution
    [Ir, Kc, Lc, Rc, Rr] = params_CR(J,L,R1,R2)
    
    Bcr = BFieldCurrentShell([0;0;0],Kc,Lc,Rc,[0;0;0],[0;1;0]) + BFieldRing([0;0;0],Ir,Rc,[0;0;0],[0;1;0]);
    
    inSolenoid = @(x,y) y<=L/2 && y>= -L/2 && (x <=R2 && x>=R1 || x <=-R1 && x>=-R2);
    
    BcrX = zeros(N,N);
    BcrY = zeros(N,N);
    BX = zeros(N,N);
    BY = zeros(N,N);
    Pe = zeros(N,N);
    XX = zeros(N,N);
    YY = zeros(N,N);
    
    %Bfield at every grid point
    for i=1:length(X)
        for j = 1:length(Y)
            XX(i,j) = X(i);
            YY(i,j) = Y(j);
            if ~inSolenoid(X(i),Y(j))
                %b-field of CR
                Bcr = BFieldCurrentShell([X(i);Y(j);0],Kc,Lc,Rc,[0;0;0],[0;1;0]) + BFieldRing([X(i);Y(j);0],Ir,Rr,[0;0;0],[0;1;0]);
                
                %b-field of solenoid
                B = BFieldSolenoid([X(i);Y(j);0],J,R1,R2,L,[0;0;0],[0;1;0]);
                BcrX(i,j) = Bcr(1);
                BcrY(i,j) = Bcr(2);
                BX(i,j) = B(1);
                BY(i,j) = B(2);
                Pe(i,j) = (norm(B-Bcr)/norm(B))*100;
            end
        end
    end
end

%% RR
if rr == true
    %%solves for parameters of RR solution
    [Ir_1, Ir_2, Rr_1, Rr_2] = params_RR(J,L,R1,R2)
    
    Brr = BFieldRing([0;0;0],Ir_1,Rr_1,[0;0;0],[0;1;0]) + BFieldRing([0;0;0],Ir_2,Rr_2,[0;0;0],[0;1;0]);
    
    inSolenoid = @(x,y) y<=L/2 && y>= -L/2 && (x <=R2 && x>=R1 || x <=-R1 && x>=-R2);
    
    BrrX = zeros(N,N);
    BrrY = zeros(N,N);
    BX = zeros(N,N);
    BY = zeros(N,N);
    Pe = zeros(N,N);
    XX = zeros(N,N);
    YY = zeros(N,N);
    
    %B field at every grid point
    for i=1:length(X)
        for j = 1:length(Y)
            XX(i,j) = X(i);
            YY(i,j) = Y(j);
            if ~inSolenoid(X(i),Y(j))
                %b-field of RR
                Brr = BFieldRing([X(i);Y(j);0],Ir_1,Rr_1,[0;0;0],[0;1;0]) + BFieldRing([X(i);Y(j);0],Ir_2,Rr_2,[0;0;0],[0;1;0]);
                
                %b-field of solenoid
                B = BFieldSolenoid([X(i);Y(j);0],J,R1,R2,L,[0;0;0],[0;1;0]);
                BrrX(i,j) = Brr(1);
                BrrY(i,j) = Brr(2);
                BX(i,j) = B(1);
                BY(i,j) = B(2);
                Pe(i,j) = (norm(B-Brr)/norm(B))*100;
            end
        end
    end
end

%%
figure;

%error contours
[c,h] = contour(XX/sqrt(R2^2+L^2/4),YY/sqrt(R2^2+L^2/4),Pe,[1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, .001,.01,.1,1,5,10,25,50,100],'Linewidth',2.5,'LineColor','k');
clabel(c,h, 'FontSize',18)
hold on

%coil rectange (orange)
% coils =rectangle('Position',[R1/sqrt(R2^2+L^2/4),0,(R2-R1)/sqrt(R2^2+L^2/4),(L/2)/sqrt(R2^2+L^2/4)],'Facecolor', [0.8500 0.3250 0.0980 1]);
coils =rectangle('Position',[R1/sqrt(R2^2+L^2/4),0,(R2-R1)/sqrt(R2^2+L^2/4),(L/2)/sqrt(R2^2+L^2/4)],'Facecolor', [0.5 0.5 0.5 0.5]);
hold on

%outer block of solenoid
rectangle('Position',[0,0,(R2)/sqrt(R2^2+L^2/4),(L/2)/sqrt(R2^2+L^2/4)])
hold on

%approx geom lines
if cc == true
    plot([Rc_1 Rc_1]/minSphere, [0 Lc_1/2]/minSphere,'k-.','Linewidth',1.5);
    hold on;
    plot([Rc_2 Rc_2]/minSphere, [0 Lc_2/2]/minSphere,'k-.','Linewidth',1.5);
    title('Magnetic Field Percent Error From CC Solution');
end

if cr == true
    plot([Rc Rc]/minSphere, [0 Lc/2]/minSphere,'k-.','Linewidth',1.5);
    hold on;
    scatter(Rr/minSphere, 0,40,'k','filled');
    title('Magnetic Field Percent Error From CR Solution');
end

if rr == true
    scatter(Rr_1/minSphere, 0,40,'k','filled');
    hold on;
    scatter(Rr_2/minSphere, 0,40,'k','filled');
    title('Magnetic Field Percent Error From RR Solution');
end

%offset lines
for i=1:O
    dist = minSphere * offsetValue(i);
    newSphere = minSphere + dist;
    newL = L/2 + dist;
    newR = R2 + dist;
    plot([0, newR]/minSphere,[newL, newL]/minSphere, 'k--','Linewidth',1.5);
    hold on;
    plot([newR, newR]/minSphere,[0, newL]/minSphere, 'k--','Linewidth',1.5);
    hold on;
end
    hold off;
% text(0.03, 0.97, '10% offset', 'FontSize', 18);
% text(0.03, 1.36, '50% offset', 'FontSize', 18);

set(gca, 'FontSize', 18);
pbaspect([1 1 1])
xlabel('Minimum Bounding Sphere Radius');
ylabel('Minimum Bounding Sphere Radius');
ylim([0 2]);
xlim([0 2]);

end

