clc
clear
close all
%% LBeam example
%**************************************************************************
% Material properties
% E = Modulus of elasticity
% v = Poisson ratio
MP = struct("E",1,"v",0.3);
%**************************************************************************
% Geometry properties
% ngp = number of gaussian points at each direction
% nor = Number of rectangle
% els = Element size for each direction
% X0, Y0 and Z0 are coordiate of left-bottom-back points for each rectangle
% cube
% lx, ly and lz are lengths of each rectangle cube
% EBC is essential boundary conditions [x0,x1,y0,y1,z0,z1,Disx,Disy,Disz]
% NBC is necessary boundary conditions [x0,x1,y0,y1,z0,z1,fx,fy,fz]
% BF Specifies boundary faces  [x1 x2 y1 y2 z1 z2]
% scale and order are used for plot results. read doc inside function
GP = struct("ngp",3,...
            "nor",2,...
            "els",0.05,...
            "X0",[0,0],"Y0",[0,0.4],"Z0",[0,0],...
            "lx",[1,0.4],"ly",[0.4,0.6],"lz",[0.1,0.1],...
            "EBC",[0,0.4,1,1,0,0.1,0,0,0],...
            "NBC",[1,1,0.2,0.2,0.05,0.05,0,-1,0],...
            "BF",[-inf inf -inf inf 0 0;... % Plane XY at Z=0
                  -inf inf -inf inf 0.1 0.1;...% Plane XY at Z=0.1
                  1 1 -inf inf -inf inf;... % Plane YZ at X=1
                  0 0 -inf inf -inf inf;... % Plane YZ at X=0
                  0 1 0 0 -inf inf;... % Plane XZ at Y=0 and X>=0 and X<=1
                  -inf inf 1 1 -inf inf; % Plane XZ at Y=1
                  0.4 1 0.4 0.4 -inf inf;... % Plane XZ at Y=0.4, X>=0.4, X<=1
                  0.4 0.4 0.4 1 -inf inf... % Plane YZ at X=0.4, Y>=0.4, Y<=1
                ],'scale',5,'order',[3,1,2]);
%**************************************************************************
%% FE Analysis
SolidFEA
%% Plot results
PlotI3DVonmisesCountour(GP.nor,NodeCoord,IC,nelx,nely,nelz,5,Phi,...
    Vonmises,[3,1,2])
%% Animation
for i = 1:360
camorbit(1,0,'data',[0 0 1])
drawnow
end