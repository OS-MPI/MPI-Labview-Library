function [proj_mag_crop_bc] = bc_linear(projections_mag_crop, projections_phase_crop, n)

projections_complex_crop = projections_mag_crop.*exp(1i*projections_phase_crop);

% baseline correction using the first & last n points of the data
bc_start = mean(projections_complex_crop(1:n,:));
bc_end = mean(projections_complex_crop(end-n:end,:));

for i = 1:length(bc_start)
    bc(:,i) = linspace(bc_start(i),bc_end(i),size(projections_complex_crop,1));
end
proj_mag_crop_bc = abs(projections_complex_crop - bc);