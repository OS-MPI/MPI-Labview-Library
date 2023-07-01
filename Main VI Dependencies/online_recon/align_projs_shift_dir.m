function [RXmag_projs] = align_projs_shift_dir(num_proj, bidirectional_shift, shift_projs, RXmag_projs)

for p = 1:num_proj
    shift_slope = shift_projs(end,p)-shift_projs(1,p);
    if bidirectional_shift == 1 && shift_slope < 0
        RXmag_projs(:,p) = flip(RXmag_projs(:,p));
    end
    
end