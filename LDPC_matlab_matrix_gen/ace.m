function [H_exp, H_q] = ace(q, factor, H_base, max_depth)
    [m, n] = size(H_base);
    tests = min(100, (q-1)*factor);
    rows = [1:m]; %randperm(m);
    
    % matrix with offsets (first row can be chosen arbitrarily)
    H_exp = -1*ones(m, n);
    H_exp(rows(1), :) = H_base(rows(1), :) - 1;
    
    % matrix with non-binary (nb) values
    H_q = -1*ones(m, n);
    H_q(rows(1), :) = H_base(rows(1), :) - 1;
    col_weight = sum(H_base);
    
    for ii = 2:m
        row = rows(ii);
        disp(sprintf('row = %d (%d from %d)', row, ii, m));
        cols = find(H_base(row, :));
        cols = cols(randperm(length(cols)));
        for jj = 1:length(cols)
            col = cols(jj);
            disp(sprintf('\tcol = %d', col));
            if jj == 1
                % first element in the row => no new cycles yet
                H_exp(row, col) = 0;
                H_q(row, col) = randi(q-1)-1;
                continue;
            end
            p = randperm(factor*(q-1), tests)-1;
            offsets = p / (q-1);
            nb_values = mod(p,q-1);
            
            all_aces = get_all_ace(H_exp(rows(1:ii), :), H_q(rows(1:ii), :), factor, max_depth, full(H_base(rows(1:ii), :)), q, col_weight, col, offsets, nb_values, -1);
            if all(all(all_aces == -1))
                H_exp(row, col) = randi(factor)-1;
                H_q(row, col) = randi(q-1)-1;
                continue;
            end
            
            max_ace = all_aces(:, 1);
            max_offset = offsets(1);
            max_nb = nb_values(1);
            for i=2:tests
                if ace_less(max_ace, all_aces(:,i))
                    max_ace = all_aces(:,i);
                    max_offset = offsets(i);
                    max_nb = nb_values(i);
                end
            end

            H_exp(row, col) = max_offset;
            H_q(row, col) = max_nb;
        end
    end    
end

