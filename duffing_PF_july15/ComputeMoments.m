function Mom = ComputeMoments(f,pw,winteg,N);
% keyboard
Mom(1) = sum(f.*pw.*winteg);

if N > 1
    for ct = 2:N
        expr = (f-Mom(1)).^ct;
        Mom(ct) = sum(expr.*pw.*winteg);
    end
end

% Mom = Mom(:);
