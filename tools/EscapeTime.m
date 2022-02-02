function T = EscapeTime( t,xt,threshold )
%[T,X] = EscapeTime(t,xt,threshold)
% The function returns the time when the series is for the first time below the threshold.
%   Inputs:
%       t: time vector of the size 1xN
%       xt: entries of timeseries (for instance from MyIVP) of the size AxN (single A-dimensional initial conditions) or AxBxN (A-dimensional, B-initial points)
%       threshold: value below which the system is assumed to escape
%   Outputs:
%       T: escape times of each series (size: either Ax1 or AxB) 

if size(xt,3)==1 %xt is of the size AxN
    T = nan(size(xt,1),1);
    for k=1:size(xt,1)
        temp = t(find(xt(k,:)<threshold,1));
        if ~isempty(temp)
            T(k,1) = temp;
        end
    end
else %xt is of the size AxBxN
    T = nan(size(xt,1),size(xt,2));
    for k=1:size(xt,1)
        for l=1:size(xt,2)
            temp = t(find(xt(k,l,:)<threshold,1));
            if ~isempty(temp)
                T(k,l) = temp;
            end
        end
    end
end

end
