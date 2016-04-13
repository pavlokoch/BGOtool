function time = linear_time_scale(counter, fcounter)
% convert SCDP time and fine time counters to linear time
time = zeros([numel(counter),1]);
% start from first element
time(1) = counter(1);
% 6.1.4.1 Time-Stamp Correction algorithm
for i = 2:numel(counter)
    dT = counter(i) - counter(i-1);
    time(i) = time(i-1) + dT;
    if dT <= -500000
        time(i) = time(i) + 1000000;
    elseif dT > 500000
        time(i) = time(i) - 1000000;
    end
end

% add fine time
time = time + fcounter.*027.77777e-3;