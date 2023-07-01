function [RXmag_projs_crop,RXphase_projs_crop,shift_projs_crop,angle_projs_crop] = reshape_projs(RXmag,RXphase,shift,angles,num_proj,bidirectional_shift, cw,cropramppts)
%reshape Summary of this function goes here
%   Detailed explanation goes here

input_size = size(RXmag);
readouts_per_image = input_size(1,2);
num_pts = readouts_per_image/num_proj;
disp(num_pts);
%RXmag_projs = zeros(num_pts,num_proj);
%RXphase_projs = zeros(num_pts,num_proj);
%shift_projs = zeros(num_pts,num_proj);
%angle_projs = zeros(num_pts,num_proj);


%for p = 1:num_proj
    %tempstart = (p-1)*num_pts+1;
    %tempstop = p*num_pts;
    
    %RXmag_projs(:,p) = RXmag(tempstart:tempstop);
    %RXphase_projs(:,p) = RXphase(tempstart:tempstop);
    %shift_projs(:,p) = shift(tempstart:tempstop);
    %angle_projs(:,p) = angles(tempstart:tempstop);
 
RXmag_projs = reshape(RXmag,num_pts,num_proj);
RXphase_projs = reshape(RXphase,num_pts,num_proj);
shift_projs = reshape(shift,num_pts,num_proj);
angle_projs = reshape(angles,num_pts,num_proj);
	

a = cropramppts;
b = cropramppts;
RXmag_projs_crop = RXmag_projs(1+a:end-b,:); 
RXphase_projs_crop = RXphase_projs(1+a:end-b,:);
shift_projs_crop = shift_projs(1+a:end-b,:);
angle_projs_crop = angle_projs(1+a:end-b,:);

%if cw == 1 && bidirectional_shift == 1
%    RXmag_projs = flipud(RXmag_projs);
%    RXphase_projs = flipud(RXphase_projs);
%end
end

