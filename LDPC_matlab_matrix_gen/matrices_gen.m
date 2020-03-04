function matrices_gen(l, n0, m, q, output_path)
    Hmask=ones(l,n0);
    [H_exp, H_q] = ace(q, m, Hmask, 10);
    H_temp = exp2sparse(q, m, H_exp, H_q);
    [H, G] = ldpc_h2g(H_temp, q);
    
    csvwrite(fullfile(output_path, 'H.csv'), full(H))
    csvwrite(fullfile(output_path, 'G.csv'), G)
end