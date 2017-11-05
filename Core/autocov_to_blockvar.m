function [A_dash SIG_dash] = autocov_to_blockvar(G,blocks)

A_dash = [];
SIG_dash = [];
for blockIndx=1:length(blocks)
    
    this_block_autocov = G(blocks{blockIndx},blocks{blockIndx},:);
    [A_dash(chanIndx,chanIndx,:) SIG_dash(chanIndx,chanIndx)] = autocov_to_var(this_block_autocov);
end

end