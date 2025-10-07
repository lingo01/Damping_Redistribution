function xRe = funStateTrans(x, flag)
    if size(x,2) == 1
        x = x';
    end
    
    if ~exist('flag', 'var')
       flag = 'Reduce'; 
    end
    
    % [delta omega V] --> [delta V]
    if strcmp(flag, 'Reduce')
        xRe = zeros(size(x,1), size(x,2)/3*2);
        for ii = 1:(size(x,2)/3)
            xRe(:, 2*ii-1) = x(:, 3*ii-2);
            xRe(:, 2*ii)   = x(:, 3*ii);
        end
    % [delta omega V] <-- [delta V]
    elseif strcmp(flag, 'Expand')
        xRe = zeros(size(x,1), size(x,2)/2*3);
        for ii = 1:(size(x,2)/2)
            xRe(:, 3*ii-2) = x(:, 2*ii-1);
            xRe(:, 3*ii)   = x(:, 2*ii);
        end
    end
end

