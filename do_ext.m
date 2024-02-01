function [xPosExt,yPosExt,dirExt,rExt]=do_ext(xPos, yPos, dir, r, Lx, Ly, amtExt)

% copy particles to account for periodic BCs
    
    % who is near a boundary
    isNearBoundary_L = xPos<=amtExt;
    isNearBoundary_R = xPos>=Lx-amtExt;
    isNearBoundary_D = yPos<=amtExt;
    isNearBoundary_U = yPos>=Ly-amtExt;
    
    % who is near a corner
    isNearEdge_LU=isNearBoundary_L & isNearBoundary_U;
    isNearEdge_LD=isNearBoundary_L & isNearBoundary_D;
    isNearEdge_RU=isNearBoundary_R & isNearBoundary_U;
    isNearEdge_RD=isNearBoundary_R & isNearBoundary_D;
        
    % extend
    xPosExt=[xPos; xPos(isNearBoundary_L)+Lx; xPos(isNearBoundary_R)-Lx; xPos(isNearBoundary_D); xPos(isNearBoundary_U);...
            xPos(isNearEdge_LU)+Lx; xPos(isNearEdge_LD)+Lx; xPos(isNearEdge_RU)-Lx; xPos(isNearEdge_RD)-Lx];
    yPosExt=[yPos; yPos(isNearBoundary_L); yPos(isNearBoundary_R); yPos(isNearBoundary_D)+Ly; yPos(isNearBoundary_U)-Ly;...
            yPos(isNearEdge_LU)-Ly; yPos(isNearEdge_LD)+Ly; yPos(isNearEdge_RU)-Ly; yPos(isNearEdge_RD)+Ly];
    dirExt=[dir; dir(isNearBoundary_L); dir(isNearBoundary_R); dir(isNearBoundary_D); dir(isNearBoundary_U);...
            dir(isNearEdge_LU); dir(isNearEdge_LD); dir(isNearEdge_RU); dir(isNearEdge_RD)];
    rExt = [r; r(isNearBoundary_L); r(isNearBoundary_R); r(isNearBoundary_D); r(isNearBoundary_U);...
            r(isNearEdge_LU); r(isNearEdge_LD); r(isNearEdge_RU); r(isNearEdge_RD)];
    