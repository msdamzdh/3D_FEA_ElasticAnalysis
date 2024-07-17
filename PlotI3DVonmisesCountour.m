function [] = PlotI3DVonmisesCountour(nor,NodeCoord,IC,nelx,nely,nelz,...
    scale,Phi,Vonmises,order)
%% Description
% This function is used for plot. users can change views 
% or axis order for getting desired view
%% Inputs
% nor = Number of rectangle
% NodeCoord = Node coordinates
% IC variable is used for prohibitting duplicate indexing, because 
    % nodes with the same coordinates in different rectangle cubes must 
    % have unique global indexes
% nelx,nely and nelz are nummber of elements for each rectangle cube at 
    % each direction
% scale is a number which multiples number of nodes for plotting result
% Phi specifies if material exist at nodes. size(Phi)=number of rows in
    % NodeCoord
% Vonmises stress at vertices
%%
arguments
    nor,NodeCoord,IC,nelx,nely,nelz,scale,Phi,Vonmises
    order = [3,2,1];
end
%==========================================================================
    % InputCoord is used for saving inputted coordiantes
    % QueryCoord is (nor X 3) cell which saves query coordinates of scaled
    % coordinates
    % p1 and p2 are patches used for plot
    QueryCoord = cell(nor,3);
    InputCoord = cell(nor,3);
    p = cell(nor,1);
    camlight 'right' 
    view(3) 
    axis 'equal'
    axis off
    colorbar
    for i=1:nor
    
        y = NodeCoord(:,2);
        x = NodeCoord(:,1);
        z = NodeCoord(:,3);
        
        InputCoord{i,1} = reshape(y(IC{i}), (nelx(i)+1), (nely(i)+1),...
            (nelz(i)+1));
        InputCoord{i,2} = reshape(x(IC{i}), (nelx(i)+1), (nely(i)+1),...
            (nelz(i)+1));
        InputCoord{i,3} = reshape(z(IC{i}), (nelx(i)+1), (nely(i)+1),...
            (nelz(i)+1));
    
        y2 = linspace(min(InputCoord{i,1}(:)), max(InputCoord{i,1}(:)),...
            scale*(nely(i)+1));
        x2 = linspace(min(InputCoord{i,2}(:)), max(InputCoord{i,2}(:)),...
            scale*(nelx(i)+1));
        z2 = linspace(min(InputCoord{i,3}(:)), max(InputCoord{i,3}(:)),...
            scale*(nelz(i)+1));
    
        [QueryCoord{i,1}, QueryCoord{i,2}, QueryCoord{i,3}] = ...
            meshgrid(y2, x2, z2);
    
        p{i} = patch('EdgeColor','none','FaceColor' ,'interp');    
    end
    if all(order==[3,1,2])
        order2 = [2,3,1];
    elseif all(order==[2,3,1])
        order2 = [3,1,2];
    else
        order2 = [1,2,3];
    end
    for i=1:nor
        phi = reshape(Phi(IC{i}), (nelx(i)+1), (nely(i)+1), (nelz(i)+1));
        vm = reshape(Vonmises(IC{i}), (nelx(i)+1), (nely(i)+1), (nelz(i)+1));
        phiq = interpn(InputCoord{i,2},...
                       InputCoord{i,1},...
                       InputCoord{i,3},...
                       phi,...
                       QueryCoord{i,2},...
                       QueryCoord{i,1},...
                       QueryCoord{i,3});
        vmq = interpn(InputCoord{i,2},...
                       InputCoord{i,1},...
                       InputCoord{i,3},...
                       vm,...
                       QueryCoord{i,2},...
                       QueryCoord{i,1},...
                       QueryCoord{i,3});
        for j=1:3
            QueryCoord{i,j}=permute(QueryCoord{i,j},order);
        end

        phiq = permute(phiq,order);
        vmq = permute(vmq,order);
        [ff,vf] = isosurface(QueryCoord{i,order2(1)},...
                             QueryCoord{i,order2(2)},...
                             QueryCoord{i,order2(3)},...
                             phiq,0);
        set(p{i}, 'Faces',ff,'Vertices',vf)
        isonormals(QueryCoord{i,order2(1)},...
                             QueryCoord{i,order2(2)},...
                             QueryCoord{i,order2(3)},...
                             vmq,p{i})
        isocolors(QueryCoord{i,order2(1)},...
                             QueryCoord{i,order2(2)},...
                             QueryCoord{i,order2(3)},...
                             vmq,p{i});
        
    end
end


