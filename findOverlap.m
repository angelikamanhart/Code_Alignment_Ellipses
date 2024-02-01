function [thetaOverlap,thetaJunction,xJunction,yJunction,dirJunction,cellNumberJunction] = findOverlap(iXpos,iYpos,iDir,iR,p)

% function which finds theta values corresponding to points of overlap and
% cell cell junctions

% INPUT: xPos (1xn array), yPos (1xn array), dir (1xn array), p (parameters)
% OUTPUT: thetaOverlap (in order of theta_1, theta_2, etc. with
% corresponding x and y values
% thetaJoining for theta values corresponding to links that join between
% cells
% cellNumberJoining which contains the cell number corresponding to
% thetaJoining

% NOTE: doesn't work for aspect ratio less than 2 - need to think of a
% different method here if this is ever needed


% extend domain to account for periodic boundary conditions
[xPos,yPos,dir,r]=do_ext(iXpos, iYpos, iDir, iR, p.sq, p.sq, 2*max(max(sqrt(iR)),max(1./sqrt(iR))));
    
a = sqrt(r);
b = 1./sqrt(r);

switch p.optModel
    case 0
       [cellsInNbhd,dist] = rangesearch([xPos yPos],[iXpos iYpos], 2*max(max(a),max(b)),'SortIndices',false);% find nearby ellipses to avoid loopoing through all
    case 1
       [cellsInNbhd,dist] = rangesearch([xPos yPos], [iXpos iYpos], 2*max(max(a),max(b)) + p.junction_range,'SortIndices',true);% find nearby ellipses to avoid loopoing through all
end 

% initialise arrays for junctions and corresponding cell numbers
junctionPair1 = [];
junctionPair2 = [];
cellNumberJunction = [];
xCellJunction1 = [];
xCellJunction2 = [];
yCellJunction1 = [];
yCellJunction2 = [];
junctionThetaVals1 = [];
junctionThetaVals2 = [];

% SECTION 1: find overlap theta vals 

thetaOverlap = cell(1,p.N); % pre-define cell to contain all theta values

% for saving found overlap points
overlapPts = cell(p.N,p.N);
overlapClusterSize = cell(p.N,p.N);

% define boundaries
xBoundary=(xPos+a.*cos(dir)*cos(p.theta)- b.*sin(dir)*sin(p.theta))';
yBoundary=(yPos+a.*sin(dir)*cos(p.theta)+ b.*cos(dir)*sin(p.theta))';


% loop through each cell
for i = 1:p.N
    cellsInNbhd{i} = cellsInNbhd{i}(dist{i}<2*max(a(i),b(i)));
    % if there are cells within neighbourhood
    if length(cellsInNbhd{i}) > 1
        % index of cells in neighbouhood (omitting current cell)
        idxCellsInNbhd = cellsInNbhd{i}(cellsInNbhd{i}~=i);
        % loop through cells in neighbourhood
        for k = 1:length(idxCellsInNbhd)
            
            % if this pair has not already been examined
            if idxCellsInNbhd(k)>i
                
                % find points on neighbouring ellipse edges within small distance of each
                % other
                
                idxPointsInRange = rangesearch([xBoundary(:,i),yBoundary(:,i)],[xBoundary(:,idxCellsInNbhd(k)),yBoundary(:,idxCellsInNbhd(k))],p.search_range,'SortIndices',false);
                overlap = find(~cellfun(@isempty,idxPointsInRange)); % get index of non empty cells
                
                if isempty(overlap)==0; isOL=1; else isOL=0; end
                
                switch p.optModel
                case 1                    
                        % define theta = 0 and theta = pi points on
                        % ellipses 1 and 2
                        endPoints1 = [xBoundary(1,i),yBoundary(1,i);xBoundary(p.n_t/2+1,i),yBoundary(p.n_t/2+1,i)];
                        endPoints2 = [xBoundary(1,idxCellsInNbhd(k)),yBoundary(1,idxCellsInNbhd(k));xBoundary(p.n_t/2+1,idxCellsInNbhd(k)),yBoundary(p.n_t/2+1,idxCellsInNbhd(k))];
                        thetaEndPoints1 = boundaryPosToTheta(endPoints1(:,1),endPoints1(:,2),xPos(i),yPos(i),dir(i),a(i),b(i));
                        thetaEndPoints2 = boundaryPosToTheta(endPoints2(:,1),endPoints2(:,2),xPos(idxCellsInNbhd(k)),yPos(idxCellsInNbhd(k)),dir(idxCellsInNbhd(k)),a(idxCellsInNbhd(k)),b(idxCellsInNbhd(k)));
                        idxJunction = rangesearch(endPoints1,endPoints2,p.junction_range); % find all points on ellipse within junction_range
                        % make sure that junction pairs are not repeated
                        % and consider idxJunction{1}
                        if cellsInNbhd{i}(1) < idxCellsInNbhd(k) && ~isempty(idxJunction{1})                  
%                             junctionPair1 = [junctionPair1, endPoints1(idxJunction{1},:)']; % extract neighbouring pairs from rangesearch
%                             junctionPair2 = [junctionPair2, endPoints2(1,:)']; % extract neighbouring pairs from rangesearch
                            cellNumberJunction = [cellNumberJunction, [cellsInNbhd{i}(1);idxCellsInNbhd(k)]]; % create array of corresponding cell numbers
                            xCellJunction1 = [xCellJunction1, endPoints1(idxJunction{1},1)'];
                            xCellJunction2 = [xCellJunction2, endPoints2(1,1)'];
                            yCellJunction1 = [yCellJunction1, endPoints1(idxJunction{1},2)'];
                            yCellJunction2 = [yCellJunction2, endPoints2(1,2)'];
                            junctionThetaVals1 = [junctionThetaVals1, thetaEndPoints1(idxJunction{1},:)'];
                            junctionThetaVals2 = [junctionThetaVals2, thetaEndPoints2(1,:)'];
                        end
                        % make sure that junction pairs are not repeated
                        % and consider idx_j{2}
                        if cellsInNbhd{i}(1) < idxCellsInNbhd(k) && ~isempty(idxJunction{2})                  
%                             junctionPair1 = [junctionPair1, endPoints1(idxJunction{2},:)']; % extract neighbouring pairs from rangesearch
%                             junctionPair2 = [junctionPair2, endPoints2(2,:)']; % extract neighbouring pairs from rangesearch
                            cellNumberJunction = [cellNumberJunction, [cellsInNbhd{i}(1);idxCellsInNbhd(k)]]; % create array of corresponding cell numbers
                            xCellJunction1 = [xCellJunction1, endPoints1(idxJunction{2},1)'];
                            xCellJunction2 = [xCellJunction2, endPoints2(2,1)'];
                            yCellJunction1 = [yCellJunction1, endPoints1(idxJunction{2},2)'];
                            yCellJunction2 = [yCellJunction2, endPoints2(2,2)'];
                            junctionThetaVals1 = [junctionThetaVals1, thetaEndPoints1(idxJunction{2},:)'];
                            junctionThetaVals2 = [junctionThetaVals2, thetaEndPoints2(2,:)'];
                        end 
                        
                end
                % if overlap array is not empty
                if isOL==1
                    
                    % work out difference between adjacent entries
                    overlap_difference = diff(overlap);
                    jump = [0]; % initialise jump array to zero
                    % store entries which jump more than one (corresponds to
                    % different clusters of points)
                    for w = 1:length(overlap_difference)
                        if overlap_difference(w) > 1
                            jump = [jump,w];
                        end
                    end
                    % add end point onto jump array
                    jump = [jump,length(overlap_difference)+1];
                    averageOverlapPoint = zeros(length(jump)-1,2); % initialise cluster average array
                    clusterSize = zeros(length(jump)-1,1);
                    for u = 2:length(jump)
                        clusters = [xBoundary([idxPointsInRange{overlap(jump(u-1)+1:jump(u))}],i),yBoundary([idxPointsInRange{overlap(jump(u-1)+1:jump(u))}],i);
                            xBoundary(overlap(jump(u-1)+1:jump(u)),idxCellsInNbhd(k)),yBoundary(overlap(jump(u-1)+1:jump(u)),idxCellsInNbhd(k))];
                        clusterSize(u-1) = length(clusters);
                        averageOverlapPoint(u-1,:) = mean(clusters); % find mean of all points in cluster
                    end
                    
                    % save found overlap points
                    overlapPts{i,idxCellsInNbhd(k)}=averageOverlapPoint;
                    overlapClusterSize{i,idxCellsInNbhd(k)}=clusterSize;
                    
                end
              
            % otherwise just recall value to determine thetas    
            else
                
                averageOverlapPoint=overlapPts{idxCellsInNbhd(k),i};
                clusterSize = overlapClusterSize{idxCellsInNbhd(k),i};
                if isempty(averageOverlapPoint)==0; isOL=1; else isOL=0; end
                
            end
            
            % if there is overlap
            if isOL==1
                if length(averageOverlapPoint(:,1)) > 1 && length(averageOverlapPoint(:,1)) < 6
                    Y = a(i)*(-(averageOverlapPoint(:,1) - xPos(i)).*sin(dir(i)) + ((averageOverlapPoint(:,2)-yPos(i)).*cos(dir(i))));
                    X = b(i)*((averageOverlapPoint(:,1) - xPos(i)).*cos(dir(i)) + ((averageOverlapPoint(:,2)-yPos(i)).*sin(dir(i))));
                    thetas = mod(atan2(Y,X),2*pi);
                    if size(averageOverlapPoint,1) == 5
                        [thetas,sort_idx] = sort(thetas);
                        averageOverlapPoint = averageOverlapPoint(sort_idx,:);
                        diff_thetas = diff([thetas; thetas(1)]);
                        [~,diffsort_idx] = sort(diff_thetas);
                        thetas(max(mod(diffsort_idx(1:2),5))) = [];
                        thetas = [thetas(1:2),thetas(3:4)];
%                         averageOverlapPoint(max(mod(diffsort_idx(1:2),5))) = [];
                    elseif size(averageOverlapPoint,1) == 4
                        [thetas,sort_idx] = sort(thetas);
                        averageOverlapPoint = averageOverlapPoint(sort_idx,:);
                        overlap_diff = diff(averageOverlapPoint);
                        overlap_diff_angle = atan2(overlap_diff(:,2),overlap_diff(:,1));
                        [~,min_idx] = min(abs(wrapToPi(2*overlap_diff_angle - 2*dir(i))));
                        final_idx = mod(min_idx-1:min_idx+2,4)+1;
                        thetas = thetas(final_idx);
                        thetas = [thetas(1:2),thetas(3:4)];
%                         averageOverlapPoint = averageOverlapPoint(final_idx,:);
                    elseif size(averageOverlapPoint,1) == 3
                        [~,del_idx] = max(clusterSize);
                        thetas(del_idx) = [];
                        averageOverlapPoint(del_idx,:) = [];
                        [~,min_idx] = min([mod(thetas(2) - thetas(1),2*pi),mod(thetas(1) - thetas(2),2*pi)]);
                        final_idx = mod(min_idx-1:min_idx,2)+1;
                        thetas = thetas(final_idx);
%                         averageOverlapPoint = averageOverlapPoint(final_idx,:);
                    elseif size(averageOverlapPoint,1) == 2
                        [~,min_idx] = min([mod(thetas(2) - thetas(1),2*pi),mod(thetas(1) - thetas(2),2*pi)]);
                        final_idx = mod(min_idx-1:min_idx,2)+1;
                        thetas = thetas(final_idx);
%                         averageOverlapPoint = averageOverlapPoint(final_idx,:);
                    end
                    if r(i) < 1 && size(thetas,2) == 2
                        thetas = [thetas(2),thetas(4);thetas(3),thetas(1)];
                        thetaOverlap{i} = [thetaOverlap{i},thetas];
                    else
                        thetaOverlap{i} = [thetaOverlap{i},thetas];
                    end 
                else
                    thetas = [];
                    thetaOverlap{i} = [thetaOverlap{i},thetas];
                end
            end
        end
    end
end

% SECTION 3: discard cell-cell junctions according to p.optJunction

thetaJunction = [junctionThetaVals1;junctionThetaVals2];
xJunction = [xCellJunction1; xCellJunction2];
yJunction = [yCellJunction1; yCellJunction2];
dirJunction = dir(cellNumberJunction);

switch p.optJunction
    case 2
        idx_keep = find(abs(mod(sum(thetaJunction,1),2*pi)- pi) < 0.1); %index of ccj to keep if we only want head to head or tail to tail junctions
        thetaJunction = thetaJunction(:,idx_keep);
        xJunction = xJunction(:,idx_keep);
        yJunction = yJunction(:,idx_keep);
        dirJunction = dirJunction(:,idx_keep);
        cellNumberJunction = cellNumberJunction(:,idx_keep);
    case 1
        idx_keep = find(abs(mod(sum(thetaJunction,1),2*pi)- pi) > 3); %index of ccj to keep if we only want head to tail junctions
        thetaJunction = thetaJunction(:,idx_keep);
        xJunction = xJunction(:,idx_keep);
        yJunction = yJunction(:,idx_keep);
        dirJunction = dirJunction(:,idx_keep);
        cellNumberJunction = cellNumberJunction(:,idx_keep);
end 

end
