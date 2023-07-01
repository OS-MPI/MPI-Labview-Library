function [angular_pos_to_subtract, AngularHomingDetected] = angular_homing_edge(angHomingTrigFirst2Images, angularPosFirst2Images, degUpright, cw)
% for angular homing, find the location of homing within the first two
% images, that becomes the new "zero" position:
%%% in CW mode, use Falling edge
%%% in CCW mode, use Rising edge

if cw  
    inds_edge =  find(diff(angHomingTrigFirst2Images) == -1); % falling edge
else
    inds_edge =  find(diff(angHomingTrigFirst2Images) == 1); % rising edge
end

% hold on, plot(angHomingTrigFirst2Images*max(angularPosFirst2Images));
% hold on, plot(angularPosFirst2Images);

if numel(inds_edge) > 1 % if you have two falling edges, choose the first one
    inds_edge = inds_edge(1);
end

if isempty(inds_edge)
    disp('angular homing not detected within first two images. No angular correction will be used.');
    angular_pos_to_subtract = 0;
    AngularHomingDetected = 'no';
else
    angular_pos_to_subtract = angularPosFirst2Images(inds_edge)
    AngularHomingDetected = 'yes';
end

angular_pos_to_subtract = angular_pos_to_subtract + degUpright;
disp(['angular correction: ',num2str(angular_pos_to_subtract),' deg']);

end
