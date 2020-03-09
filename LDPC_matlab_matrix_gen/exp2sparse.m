function H = exp2sparse(q, factor, H_exp, H_q)
    [m, n] = size(H_exp); 
    init_matrix = sparse(eye(factor));
    zero_matrix = sparse(zeros(factor, factor));
    
    % convert nb values from degree representation
    zero_indexes = (H_q == -1);
    if q > 2
        alpha = gf(2, log2(q));
        temp = alpha.^H_q;
        H_q = double(temp.x);
    else
        H_q = ones(m, n);
    end
    H_q(zero_indexes) = 0;
    
    H = [];
    
    for ii = 1:m
        layer = [];
        for jj = 1:n
            shift = H_exp(ii, jj);
            if shift == -1
                layer = [layer zero_matrix];
            else
                temp = circshift(init_matrix, -shift);
                layer = [layer H_q(ii, jj)*temp];
            end
        end
        H = [H; layer];
    end
end

