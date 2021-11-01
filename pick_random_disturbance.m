function w = pick_random_disturbance(W)
% pick disturbance form uniform distribution
verts = W.V;
b_max = max(verts)';
b_min = min(verts)';

% generate random until it will be inside of W
while true
    w0 = rand(1) .* (b_max - b_min) + b_min; % ensure that w lies in [b_min,b_max]
    w = w0(1);
    if W.contains(w0)
        break
    end
end
end