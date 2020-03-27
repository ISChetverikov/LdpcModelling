function matrices_gen(l, n0, m, q, output_path, is_h_sparse, is_only_h)
    Hmask=ones(l,n0);
    [H_exp, H_q] = ace(q, m, Hmask, 10);
    H_temp = exp2sparse(q, m, H_exp, H_q);
    [H, G] = ldpc_h2g(H_temp, q);
    
    h_filename = sprintf('h_%d_%d_%d', l, n0, m);
    
    if is_h_sparse
        [colH, rowH, valH] = find(H');
        h_filename = fullfile(output_path, h_filename + ".sprs");      
        dlmwrite(h_filename, [rowH colH valH], ' ');
    else
        H = full(H);
        h_filename = fullfile(output_path, h_filename + ".csv");
        csvwrite(h_filename, H);
    end 
    
    if ~is_only_h
        g_filename = sprintf('g_%d_%d_%d', l, n0, m);
        g_filename = fullfile(output_path, g_filename + ".csv");
        csvwrite(g_filename, G); % G is dense matrix
    end
end