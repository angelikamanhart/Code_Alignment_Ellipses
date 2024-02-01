function [xPos, yPos, dir, r, junctionDetails] = fullSim(p,xPos0,yPos0,dir0,r0)

% function to simulate model
nrTP=floor(p.T/p.dt);

% pre-define and initialise position
xPos = zeros(nrTP,length(xPos0)); % assign random centre points of ellipses (x)
yPos = zeros(nrTP,length(xPos0)); % assign random centre points of ellipses (y)

xPos(1,:) = xPos0;
yPos(1,:) = yPos0;

% pre-define and initialise orientation
dir = zeros(nrTP,length(xPos0));
dir(1,:) = dir0;

% pre-define and initalise r
r = zeros(nrTP,length(xPos0));

r(1,:) = r0;

% initialise cell array with extra cell-cell junction info for testing purposes
junctionDetails = cell(nrTP,3); 

thetasJoining = []; % initialise thetaJunction array (empty)
cellNumberJoining = []; % initialise cellNumberJunction array (empty)

iXpos=xPos0(:);
iYpos=yPos0(:);
iDir=dir0(:);
iR=r0(:);

% loop over time points
for i = 1:nrTP-1
    i
    
    % pre-define arrays for rotation, translation and shape change arising 
    % from potential V and cell cell junctions
    dirChangeV = zeros(p.N,1);
    xChangeV = zeros(p.N,1);
    yChangeV = zeros(p.N,1);
    rChangeV = zeros(p.N,1);
    
    dirChangeJunction = zeros(length(iXpos),1);
    xChangeJunction = zeros(length(iXpos),1);
    yChangeJunction = zeros(length(iXpos),1);
    dirChangeActin = zeros(length(iXpos),1);
    
    % define a and b from r
    iA = sqrt(iR);
    iB = 1./sqrt(iR);
    
    [thetaOverlap,thetaJunction,xJunction,yJunction,dirJunction,cellNumberJunction] = findOverlap(iXpos,iYpos,iDir,iR,p);
    
    for k = 1:p.N
        % calculate how position and direction (and aspect ratio) are effected by V
        if ~isempty(thetaOverlap{k})
            dirChangeV(k) = sum((sin(thetaOverlap{k}(2,:))).^2 - (sin(thetaOverlap{k}(1,:))).^2);
            xComponent = sum(iB(k)*(sin(thetaOverlap{k}(2,:)) - sin(thetaOverlap{k}(1,:))));
            yComponent = sum(-iA(k)*(cos(thetaOverlap{k}(2,:)) - cos(thetaOverlap{k}(1,:))));
            xChangeV(k) = cos(iDir(k))*xComponent - sin(iDir(k))*yComponent;
            yChangeV(k) = sin(iDir(k))*xComponent + cos(iDir(k))*yComponent;
            switch p.optCellShape
                case 1
                    rChangeV(k) = sum((sin(2*thetaOverlap{k}(2,:)) - sin(2*thetaOverlap{k}(1,:))));
            end 
        end     
    end
    
    % if model with cell-cell junctions is being implemented:
    switch p.optModel
        case 1
               
        if ~isempty(cellNumberJunction)
                
            max_cell = max(cellNumberJunction,[],'all');
            
            if max_cell > length(dirChangeJunction)
                dirChangeJunction = zeros(max_cell,1);
                xChangeJunction = zeros(max_cell,1);
                yChangeJunction = zeros(max_cell,1);
                dirChangeActin = zeros(max_cell,1);
            end 
       
            for j = 1:size(cellNumberJunction,2) 
                xChangeJunction(cellNumberJunction(1,j)) = xChangeJunction(cellNumberJunction(1,j)) + xJunction(2,j) - xJunction(1,j);
                yChangeJunction(cellNumberJunction(1,j)) = yChangeJunction(cellNumberJunction(1,j)) + yJunction(2,j) - yJunction(1,j);
                xChangeJunction(cellNumberJunction(2,j)) = xChangeJunction(cellNumberJunction(2,j)) + xJunction(1,j) - xJunction(2,j);
                yChangeJunction(cellNumberJunction(2,j)) = yChangeJunction(cellNumberJunction(2,j)) + yJunction(1,j) - yJunction(2,j);
                dirChangeJunction(cellNumberJunction(1,j)) = dirChangeJunction(cellNumberJunction(1,j)) + dot([xJunction(2,j) - xJunction(1,j);yJunction(2,j) - yJunction(1,j)],[-sin(dirJunction(1,j));cos(dirJunction(1,j))])*sign(sin(pi/2 - thetaJunction(1,j)));
                dirChangeJunction(cellNumberJunction(2,j)) = dirChangeJunction(cellNumberJunction(2,j)) + dot([xJunction(1,j) - xJunction(2,j);yJunction(1,j) - yJunction(2,j)],[-sin(dirJunction(2,j));cos(dirJunction(2,j))])*sign(sin(pi/2 - thetaJunction(2,j)));
                dirChangeActin(cellNumberJunction(1,j)) = dirChangeActin(cellNumberJunction(1,j)) + sin(dirJunction(1,j) - dirJunction(2,j))*sign(sin(pi/2 - thetaJunction(1,j))*sin(pi/2 - thetaJunction(2,j)));
                dirChangeActin(cellNumberJunction(2,j)) = dirChangeActin(cellNumberJunction(2,j)) + sin(dirJunction(2,j) - dirJunction(1,j))*sign(sin(pi/2 - thetaJunction(1,j))*sin(pi/2 - thetaJunction(2,j)));
            end 
            
            switch p.optTesting
                case 1
                    junctionDetails{i,1} = [xJunction(1,:);yJunction(1,:)];
                    junctionDetails{i,2} = [xJunction(2,:);yJunction(2,:)];
            end 
        end 
    end 
    
    % define rhs
    dxPos = (p.nu)*cos(iDir) - xChangeV + p.kappa*xChangeJunction(1:p.N);
    dyPos = (p.nu)*sin(iDir) - yChangeV + p.kappa*yChangeJunction(1:p.N);
    dDir =  - (2*iR)./(1 + iR.^2).*dirChangeV + 4*p.kappa.*iR.^(3/2)./(1 + iR.^2).*dirChangeJunction(1:p.N) + p.mu*sqrt(iR)./(1 + iR.^2).*dirChangeActin(1:p.N);
    dr =  - (4*iR.^2./(iR.^2+1)).*rChangeV - (16*p.gamma)./(iR.^2+1)./p.rBar.*(iR - p.rBar).*(iR.^3.*p.rBar + 1);
    
    % update
    iXpos = mod(iXpos + dxPos*p.dt,p.sq);
    iYpos = mod(iYpos + dyPos*p.dt,p.sq);
    iDir = mod(iDir+ dDir.*p.dt,2*pi);
    iR = iR+dr*p.dt;
    
    % store
    xPos(i+1,:) = iXpos;
    yPos(i+1,:) = iYpos;
    dir(i+1,:) = iDir;
    r(i+1,:) = iR;
    
end  