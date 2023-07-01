function [angular_pos_to_subtract, AngularHomingDetected] = angular_homing_falledge(angHomingTrigFirst2Images, angularPosFirst2Images, degUpright)
% for angular homing, find the location of homing within the first two
% images, that becomes the new "zero" position:

inds_falledge =  find(diff(angHomingTrigFirst2Images) == -1);
% hold on, plot(angHomingTrigFirst2Images*max(angularPosFirst2Images));
% hold on, plot(angularPosFirst2Images);

if numel(inds_falledge) > 1 % if you have two falling edges, choose the first one
    inds_falledge = inds_falledge(1);
end

if isempty(inds_falledge)
    disp('angular homing not detected within first two images. No angular correction will be used.');
    angular_pos_to_subtract = 0;
    AngularHomingDetected = 'no';
else
    angular_pos_to_subtract = angularPosFirst2Images(inds_falledge)
    AngularHomingDetected = 'yes';
end

angular_pos_to_subtract = angular_pos_to_subtract + degUpright;
disp(['angular correction: ',num2str(angular_pos_to_subtract),' deg']);

end
