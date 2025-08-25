mean_gene_Gtex = mean(gene_exp_Gtex); 

normalized_gene_exp_one = zeros(size(gene_exp));
maxi= max(max(gene_exp));

for g = 1:numel(gene_ids)
    
    ref_gene_name = gene_ids_Gtex(1, g);
    des_idx = find(strcmp(gene_ids, ref_gene_name));
%%   normalized_gene_exp_one_col = gene_exp(:, des_idx)./mean_gene_Gtex(1, g);
    
    if mean_gene_Gtex(1, g) == 0
        gene_exp_one = gene_exp(:, des_idx);
        for i =1:numel(gene_exp_one)
            aa = gene_exp_one(i);
            if aa == 0
                gene_exp_one(i) = 0;
            else
                gene_exp_one(i) = max(gene_exp_one);
            end
        end
        normalized_gene_exp_one_col = gene_exp_one;
    else
           normalized_gene_exp_one_col = gene_exp(:, des_idx)./mean_gene_Gtex(1, g);   
    end

    normalized_gene_exp_one(:, des_idx) = normalized_gene_exp_one_col;
end

a=1
