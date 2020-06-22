function [Dmax,Dind] = max_haversine(Lp,str) % Lp stands for Location points

    % Initialization code
    numberOfCoords = size(Lp,1); %20;
%     x = Lp(:,1); %rand(1, numberOfCoords);
%     y = Lp(:,2); %rand(1, numberOfCoords);
    maxDistance = zeros(numberOfCoords, 1); %zeros(1, numberOfCoords);
    indexOfMax = zeros(numberOfCoords, 1); %zeros(1, numberOfCoords, 'int32');
    %-----------------------------------------------
    % Main engine:
    % Find the furthest away points.
    for k = 1 : numberOfCoords
        distances = haversine(Lp(:,1),Lp(:,2),Lp(k,1),Lp(k,2)); %sqrt((x-x(k)).^2 + (y-y(k)).^2);
        if(strcmp(str,'max')) [maxDistance(k), indexOfMax(k)] = max(distances);
        else [maxDistance(k), indexOfMax(k)] = min(distances);
        end
    end
    %-----------------------------------------------
    % Done!  Now show out results in the command window.
    % Display in command window.
%     for k = 1 : numberOfCoords
%     thisDistance = maxDistance(k);
%     thisIndex = indexOfMax(k);
%     fprintf('Point #%d at (%.1f, %.1f) is farthest away from point #%d at (%.1f, %.1f) and is at a distance of %f\n',...
%       thisIndex, x(thisIndex), y(thisIndex),...
%       k, x(k), y(k), thisDistance);
%     end
%     [maxDistance indexOfMax]

%     Dmax = max(maxDistance);
    Dmax = maxDistance; % Vectore of Max distance from each Lp to anyother Lp
    Dind = indexOfMax;
end