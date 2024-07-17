%**************************************************************************
% MeshGeneration
[ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,...
    nely,nelz] = MeshGeneration(GP.nor,GP.lx,GP.ly,GP.lz,GP.els,GP.X0,GP.Y0,GP.Z0,GP.BF);
%**************************************************************************
% Boundary conditions implementation
[F,Dis,ukdis,kdis] = BoundaryConditionsImplementation(ndof,GP.EBC,GP.NBC,NodeCoord,GP.els);
%**************************************************************************
% Indexes which is used in stiffness assembly
iK = reshape(kron(EDofMat,ones(24,1))',24*24*noe,1);
jK = reshape(kron(EDofMat,ones(1,24))',24*24*noe,1);
%**************************************************************************
% Element stiffness calculation
[Ke,Vole,Bmat,Emat] = BrickElementStiffnessCalculation(GP.ngp,GP.els,MP.v);
%**************************************************************************
% Assembling of stiffness matrix
sK = MP.E*reshape(Ke(:)*ones(1,noe),24*24*noe,1);
K = sparse(iK,jK,sK);
K = (K+K')/2; % for assuring symmetry
%**************************************************************************
% Solving equilibrium equations
Dis(ukdis) = K(ukdis,ukdis)\(F(ukdis)-K(ukdis,kdis)*Dis(kdis));
F(kdis) = K(kdis,:)*Dis;
%**************************************************************************
% Bmat calculation at vertices
[Bmat_node] = Bmat_builder(GP.els);
%**************************************************************************
% Vonmises stress calculation at vertices
[Vonmises] = Vonmises_stress(noe,Bmat_node,Dis,EDofMat,Emat,...
    MP.E,NodeRep,ElNodes);
%**************************************************************************
%% Mesh generation function
function [ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,...
    nely,nelz] = MeshGeneration(nor,lx,ly,lz,els,X0,Y0,Z0,BF)
%==========================================================================
    % Inputs
    % nor = Number of rectangles
    % lx = array which has x_length of rectangles
    % ly = array which has y_length of rectangles
    % lz = array which has z_length of rectangles
    % els = Element size at each direction
%==========================================================================
    % Preallocation arrays for the 
    % number of elements for each rectangle at each direction
    nelx = zeros(nor,1); % nelx(i) = int64(lx(i)/els); 
    nely = zeros(nor,1); % nely(i) = int64(ly(i)/els); 
    nelz = zeros(nor,1); % nelz(i) = int64(lz(i)/els);
%==========================================================================
    % coordinates{i} = [x,y,z] where x,y and z are node coordinates for
    % i'th rectangle
    coordinates = cell(nor,1);
%==========================================================================
    % NodeCoord = [coordinates{1};...,coordinates{nor}] by removing
    % duplicate nodes
    NodeCoord = [];
%==========================================================================
    % Loop over rectangles
    for i=1:nor
        nelx(i) = int64(lx(i)/els); 
        nely(i) = int64(ly(i)/els); 
        nelz(i) = int64(lz(i)/els);
        nonx = nelx(i)+1; 
        nony = nely(i)+1; 
        nonz = nelz(i)+1;
%==========================================================================
        % X0 = X-coordinate of left-bottom-back for each rectangle
        % Y0 = Y-coordinate of left-bottom-back for each rectangle
        % Z0 = Z-coordinate of left-bottom-back for each rectangle
        % nonx, nony and nonz are number of nodes for i'th rectangle for
        % x-y-z direction, respectively which can be calculated by
        x_node = linspace(X0(i),X0(i)+lx(i),nonx);
        y_node = linspace(Y0(i),Y0(i)+ly(i),nony);
        z_node = linspace(Z0(i),Z0(i)+lz(i),nonz);
%==========================================================================
        x = repmat(x_node',nony*nonz,1);
        y = repmat(kron(y_node',ones(nonx,1)),nonz,1);
        z = kron(z_node',ones(nonx*nony,1));
        coordinates{i} = [x,y,z];
        NodeCoord = [NodeCoord;coordinates{i}];
    end
%==========================================================================
    % removing duplicate nodes
    NodeCoord = unique(NodeCoord,'stable','rows');
%==========================================================================
    % Phi specifies if material exist at nodes. size(Phi)=number of rows in
    % NodeCoord
    % Phi>=0 => material exist
    % Phi = 0=>boundary representation
    % Phi<0 => no material
    Phi= 0.1*ones(size(NodeCoord,1),1);
%==========================================================================
    % Loop over BF for boundary specification(where Phi=0)
        for j=1:size(BF,1)
            ind = NodeCoord(:,1)>=BF(j,1) & NodeCoord(:,1)<=BF(j,2)...
                & NodeCoord(:,2)>=BF(j,3) & NodeCoord(:,2)<=BF(j,4)...
                & NodeCoord(:,3)>=BF(j,5) & NodeCoord(:,3)<=BF(j,6);
            Phi(ind)=0;
        end 
%==========================================================================
    % ElNodes = global indexes for nodes of each element 
    ElNodes = [];
%==========================================================================
    % IC variable for prohibiting duplicate indexing is used, because 
    % nodes with the same coordinates in different rectangles must have 
    % unique global indexes     
    IC = cell(nor,1);
%==========================================================================
    % Loop over rectangles
    for i=1:nor
        nodes = (1:size(coordinates{i},1))';
        [~,ic,iC] = intersect(coordinates{i},NodeCoord,'stable','rows');
        nodes(ic) = iC;IC{i}=iC;
        nonx = nelx(i)+1; nony = nely(i)+1;
        a = repmat([0,1,nelx(i)+2,nelx(i)+1],nelx(i),1);
        b = repmat(a,nely(i),1)+kron((0:nely(i)-1)',ones(nelx(i),1));
        c = repmat(b,nelz(i),1)+kron((0:nelz(i)-1)',nonx*nony*ones(nelx(i)*nely(i),1));
        d = repmat((1:nelx(i)*nely(i))',nelz(i),4)+c;
        eleNode = [d,d+nonx*nony];
        if i>1
            eleNode = nodes(eleNode);
        end
        ElNodes = [ElNodes;eleNode]; %ok
    end
%==========================================================================
    % EDofMat = Degrees of freedom for each element
    EDofMat = kron(ElNodes,[3,3,3])+repmat([-2,-1,0],1,8);
%==========================================================================
    % noe = Number of elements
    % ndof = Number of degrees of freedom
    % nonodes = Number of nodes
    noe = sum(nelx.*nely.*nelz);
    nonodes = max(ElNodes,[],'all');
    ndof = 3*nonodes;
%==========================================================================
    % NodeRep = reppetition of each node in all element;
    NodeRep=groupcounts(ElNodes(:));
end
%% Boundary conditions implementation function
function [F,Dis,ukdis,kdis] = BoundaryConditionsImplementation(ndof,EBC,NBC,NodeCoord,els)
%==========================================================================
    % Dis = Displacement vector
    % F = Force vector
    Dis = nan(ndof,1);
    F = zeros(ndof,1);
%==========================================================================
    for s=1:2
        if s==1
            BC = EBC;
            ND = Dis;
        else
            BC = NBC;
            ND = F;
        end
        for i=1:size(BC,1)
            Cond = NodeCoord(:,1)>=BC(i,1) & ...
                   NodeCoord(:,1)<=BC(i,2) & ...
                   NodeCoord(:,2)>=BC(i,3) & ...
                   NodeCoord(:,2)<=BC(i,4) & ...
                   NodeCoord(:,3)>=BC(i,5) & ...
                   NodeCoord(:,3)<=BC(i,6);
            Cond = find(Cond);
            if ~isempty(Cond)
                ND(Cond*3-2) = BC(i,7);
                ND(Cond*3-1) = BC(i,8);
                ND(Cond*3) = BC(i,9);
            else
                distance = sum((NodeCoord-[BC(i,1),BC(i,3),BC(i,7)]).^2,2);
                Cond = find(distance<=els);
                ND(Cond*3-2) = BC(i,7)/numel(Cond);
                ND(Cond*3-1) = BC(i,8)/numel(Cond);
                ND(Cond*3) = BC(i,9)/numel(Cond);
            end
        end
        if s==1
            Dis = ND;
        else
            F = ND;
        end
    end
%==========================================================================
    % ukdis = Unknown displacement index
    % kdis = Known displacement index
    ukdis = isnan(Dis);
    kdis = ~isnan(Dis);
end
%% Brick element stiffness calculation function
function [Ke,Vole,Bmat,Emat] = BrickElementStiffnessCalculation(ngp,els,v)

% Notes!!!
% This function is usable for structres with elements with equal size at 
% each direction
%==========================================================================
% ngp = Number of gauss points for each direction
% els = Element size at each direction
%==========================================================================
    % Global coordinates for each element
    % because the element size dose not change in the domain of structre,
    % this global coordinate is used for all element.
    eXcor = [0,els,els,0,0,els,els,0]';
    eYcor = [0,0,els,els,0,0,els,els]';
    eZcor = [0,0,0,0,els,els,els,els]';
%==========================================================================
 % Elasticity matrix by cosidering 1 as modulus of elasticity
Emat = 1/((1+v)*(1-2*v))*...
    [1-v,v,v,0,0,0;...
    v,1-v,v,0,0,0;...
    v,v,1-v,0,0,0;
    (1-2*v)/2*[zeros(3),eye(3)]];
%==========================================================================
    % Strain-displacement matrix
    Bmat = zeros(6,24);
%==========================================================================
% Volume at gaussian points
    Volg = zeros(ngp,ngp,ngp);
%==========================================================================
    % Vole = Volume for each element
    Vole = 0;
%==========================================================================
    % Ke = stiffness matrix for element
    Ke = zeros(24);
%==========================================================================
% gauss points and their weights
    [gp,wgp]=makegaussianpoint(ngp);
%==========================================================================
   % Loop over gauss points for Stiffness calcualtion
     for i=1:ngp
        t = gp(i);
        for j=1:ngp
            s = gp(j);
            for k=1:ngp
                r = gp(k);
%==========================================================================
% 8-node shape functions for brick element
%                 N=0.125*[(1-r)*(1-s)*(1-t),(1+r)*(1-s)*(1-t),...
%                          (1+r)*(1+s)*(1-t),(1-r)*(1+s)*(1-t),...
%                         (1-r)*(1-s)*(1+t),(1+r)*(1-s)*(1+t),...
%                         (1+r)*(1+s)*(1+t),(1-r)*(1+s)*(1+t)];
%==========================================================================
                % Calculation of shape function derivative wrt local
                % coordinates
                N_r = 0.125*[-(1-s)*(1-t),(1-s)*(1-t),...
                             (1+s)*(1-t),-(1+s)*(1-t),...
                             -(1-s)*(1+t),(1-s)*(1+t),...
                            (1+s)*(1+t),-(1+s)*(1+t)];
    
                N_s = 0.125*[-(1-r)*(1-t),-(1+r)*(1-t),...
                               (1+r)*(1-t),(1-r)*(1-t),...
                             -(1-r)*(1+t),-(1+r)*(1+t),...
                              (1+r)*(1+t),(1-r)*(1+t)];
    
                N_t=0.125*[-(1-r)*(1-s),-(1+r)*(1-s),...
                           -(1+r)*(1+s),-(1-r)*(1+s),...
                            (1-r)*(1-s),(1+r)*(1-s),...
                            (1+r)*(1+s),(1-r)*(1+s)];
%==========================================================================
                % Calculation derivative of global coordinates 
                % wrt local coordinates
                X_r = N_r*eXcor; Y_r = N_r*eYcor; Z_r = N_r*eZcor;
                X_s = N_s*eXcor; Y_s = N_s*eYcor; Z_s = N_s*eZcor;
                X_t = N_t*eXcor; Y_t = N_t*eYcor; Z_t = N_t*eZcor;
%==========================================================================
                % Jacobian matrix calculation
                J = [X_r,Y_r,Z_r;X_s,Y_s,Z_s;X_t,Y_t,Z_t];
                N_X_Y_Z = J\[N_r;N_s;N_t];
%==========================================================================
                % For Normal strain at x direction
                    Bmat(1,1:3:end)=N_X_Y_Z(1,:);

                % For Normal strain at y direction
                    Bmat(2,2:3:end)=N_X_Y_Z(2,:);

                % For Normal strain at z direction
                    Bmat(3,3:3:end)=N_X_Y_Z(3,:); 

                % For Shear strain at x and y direction
                    Bmat(4,1:3:end)=N_X_Y_Z(2,:);
                    Bmat(4,2:3:end)=N_X_Y_Z(1,:);

                % For Shear strain at x and z direction
                    Bmat(5,1:3:end)=N_X_Y_Z(3,:);
                    Bmat(5,3:3:end)=N_X_Y_Z(1,:);

                % For Shear strain at y and z direction
                    Bmat(6,2:3:end)=N_X_Y_Z(3,:);
                    Bmat(6,3:3:end)=N_X_Y_Z(2,:);
%==========================================================================
                % Volume at gaussian points
                    Volg(k,j,i) = det(J)*wgp(k)*wgp(j)*wgp(i);

                % Volume for each element = sum(Arg)
                    Vole = Vole+Volg(k,j,i);
%==========================================================================
                % Stiffness for each element = sum(Stiffness at each
                % gaussian point)
                    Ke = Ke+Bmat'*Emat*Bmat*Volg(k,j,i);
%==========================================================================
            end
        end
    end
end
%% Bmat_builder function
function [Bmat_node] = Bmat_builder(els)
%==========================================================================
    % Global coordinates for each element
    % because the element size dose not change in the domain of structre,
    % this global coordinate is used for all element.
    eXcor = [0,els,els,0,0,els,els,0]';
    eYcor = [0,0,els,els,0,0,els,els]';
    eZcor = [0,0,0,0,els,els,els,els]';
%==========================================================================
    % Preallocation of Bmat at 8 nodes of each element
    Bmat_node = zeros(6,24,8);
%==========================================================================
    % Local coordinates of vertices for master element
    % Node1 = [-1,-1,-1]
    % Node2 = [1,-1,-1]
    % Node3 = [-1,1,-1]
    % Node4 = [1,1,-1]
    % Node5 = [-1,-1,1]
    % Node6 = [1,-1,1]
    % Node7 = [-1,1,1]
    % Node8 = [1,1,1]
    rn = [-1,1];
    sn = rn;
    tn = rn;
%==========================================================================
    for i=1:2
        t = tn(i);
        for j=1:2
            s = sn(j);
            for k=1:2
                r = rn(k);
%==========================================================================                
                % ind = local index of the current node in loop
                ind = k+(j-1)*2+(i-1)*4;
%==========================================================================
% 8-node shape functions for brick element
%                 N=0.125*[(1-r)*(1-s)*(1-t),(1+r)*(1-s)*(1-t),...
%                          (1+r)*(1+s)*(1-t),(1-r)*(1+s)*(1-t),...
%                         (1-r)*(1-s)*(1+t),(1+r)*(1-s)*(1+t),...
%                         (1+r)*(1+s)*(1+t),(1-r)*(1+s)*(1+t)];
%==========================================================================
                % Calculation of shape function derivative wrt local
                % coordinates
                N_r = 0.125*[-(1-s)*(1-t),(1-s)*(1-t),...
                             (1+s)*(1-t),-(1+s)*(1-t),...
                             -(1-s)*(1+t),(1-s)*(1+t),...
                            (1+s)*(1+t),-(1+s)*(1+t)];
    
                N_s = 0.125*[-(1-r)*(1-t),-(1+r)*(1-t),...
                               (1+r)*(1-t),(1-r)*(1-t),...
                             -(1-r)*(1+t),-(1+r)*(1+t),...
                              (1+r)*(1+t),(1-r)*(1+t)];
    
                N_t=0.125*[-(1-r)*(1-s),-(1+r)*(1-s),...
                           -(1+r)*(1+s),-(1-r)*(1+s),...
                            (1-r)*(1-s),(1+r)*(1-s),...
                            (1+r)*(1+s),(1-r)*(1+s)];
%==========================================================================
                % Calculation derivative of global coordinates 
                % wrt local coordinates
                X_r = N_r*eXcor; Y_r = N_r*eYcor; Z_r = N_r*eZcor;
                X_s = N_s*eXcor; Y_s = N_s*eYcor; Z_s = N_s*eZcor;
                X_t = N_t*eXcor; Y_t = N_t*eYcor; Z_t = N_t*eZcor;
%==========================================================================
                % Jacobian matrix calculation
                J = [X_r,Y_r,Z_r;X_s,Y_s,Z_s;X_t,Y_t,Z_t];
                N_X_Y_Z = J\[N_r;N_s;N_t];
%==========================================================================
                % For Normal strain at x direction
                    Bmat_node(1,1:3:end,ind)=N_X_Y_Z(1,:);

                % For Normal strain at y direction
                    Bmat_node(2,2:3:end,ind)=N_X_Y_Z(2,:);

                % For Normal strain at z direction
                    Bmat_node(3,3:3:end,ind)=N_X_Y_Z(3,:); 

                % For Shear strain at x and y direction
                    Bmat_node(4,1:3:end,ind)=N_X_Y_Z(2,:);
                    Bmat_node(4,2:3:end,ind)=N_X_Y_Z(1,:);

                % For Shear strain at x and z direction
                    Bmat_node(5,1:3:end,ind)=N_X_Y_Z(3,:);
                    Bmat_node(5,3:3:end,ind)=N_X_Y_Z(1,:);

                % For Shear strain at y and z direction
                    Bmat_node(6,2:3:end,ind)=N_X_Y_Z(3,:);
                    Bmat_node(6,3:3:end,ind)=N_X_Y_Z(2,:);
%==========================================================================
            end
        end
    end
end
%% Vonmises stress calculation function
function [Vonmises] = Vonmises_stress(noe,Bmat_node,Dis,EDofMat,...
    Emat,E,NodeRep,ElNodes)
%==========================================================================
    % Inputs
    % noe = Number of elements
    % Bmat_node = Strain-displacement matrix at vertices
    % Dis = Displacement vector
    % EDofMat = Degrees of freedom for each element
    % Emat = Elasticity matrix by cosidering 1 as modulus of elasticity
    % E = Modulus of elasticity
%==========================================================================
    % Vonmises = Vonmises stress at vertices
    Vonmises = zeros(8,noe);
%==========================================================================
    % Loop over vertices
    for i=1:8
        Stress = E*Emat*Bmat_node(:,:,i)*Dis(EDofMat');
        a = (Stress(1,:)-Stress(2,:)).^2;
        b = (Stress(2,:)-Stress(3,:)).^2;
        c = (Stress(3,:)-Stress(1,:)).^2;
        d = 6*(Stress(4,:).^2);
        e = 6*(Stress(5,:).^2);
        f = 6*(Stress(6,:).^2);
        Vonmises(i,:) = sqrt(0.5*(a+b+c+d+e+f));
    end
%==========================================================================
%     Averaging vonmises stress at common nodes of elements; 
    Vonmises =(1./NodeRep).*full(sparse(ElNodes,...
        ones(size(ElNodes)),Vonmises'));
end