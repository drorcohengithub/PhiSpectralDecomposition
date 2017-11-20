function [SIG_r A_r As SIGs] = autocov_to_blockvar(G,blocks)

As = {};
A_r = [];
SIGs = {};
SIG_r = [];
for blockIndx=1:length(blocks)
    this_block_autocov = G(blocks{blockIndx},blocks{blockIndx},:);
    [As{blockIndx} SIGs{blockIndx}] = autocov_to_var(this_block_autocov);
    
    try
        [G1, info] = var_to_autocov(As{blockIndx}, SIGs{blockIndx},1000);
        var_info(info,true);

    catch ME
        disp(ME)
    end
    SIG_r(blocks{blockIndx},blocks{blockIndx}) = SIGs{blockIndx};
    A_r(blocks{blockIndx},blocks{blockIndx},:) = As{blockIndx};

end

end