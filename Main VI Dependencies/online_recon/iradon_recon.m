function [IR_crop, proj_mag_crop_bc, theta] = iradon_recon(num_proj,proj_mag_crop_bc,angles,LPF_cut,enable_2nd_bc)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%%%% baseline correction used here uses a linear correction between mean(first 2 pts) and mean(last 2 pts) for each projection

theta = mean(angles,1);

%projections_complex_crop = projections_mag_crop.*exp(1i*projections_phase_crop);

%% baseline_correction = mean(projections_complex_crop(1:4,:));
%% proj_mag_crop_bc  = abs(projections_complex_crop - repmat(baseline_correction, size(projections_complex_crop,1) ,1)); 

%% baseline correction using the first n points of the data
%n = 2;
%bc_start = mean(projections_complex_crop(1:n,:));
%bc_end = mean(projections_complex_crop(end-n:end,:));

%for i = 1:length(bc_start)
%    bc(:,i) = linspace(bc_start(i),bc_end(i),size(projections_complex_crop,1));
%end
%proj_mag_crop_bc = abs(projections_complex_crop - bc);

% Filter baseline corrected data
%proj_mag_crop_bc = filter(filter_width, 1, proj_mag_crop_bc);
proj_mag_crop_bc = LPF(proj_mag_crop_bc,100,LPF_cut);

if enable_2nd_bc ==1
    %Second baseline corretion
    xmax = length(proj_mag_crop_bc(:,1));
    m = proj_mag_crop_bc(end,:)/xmax;
    x = linspace(0,xmax-1,xmax);
    m = repmat(m,xmax,1);
    x = repmat(x,num_proj,1);
    bc2 = (m.*x');
    proj_mag_crop_bc = (proj_mag_crop_bc - bc2);
end;

IR_crop = iradon(proj_mag_crop_bc,theta,'linear',...
    'Ram-Lak',1,length(proj_mag_crop_bc(:,1)));

end

