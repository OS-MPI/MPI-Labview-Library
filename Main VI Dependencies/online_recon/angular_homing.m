%% Angular homing:

ndigits = 3; % Trigger line n digits output ( = n cols)
trig000 = de2bi(raw_data.image000.trigger_code,ndigits);
angHomingtrig000A = trig000(:,1); % first column is top angular homing trigger line (second column is hyper/hypo switch trig)
angHomingTrigAFirst2Images = angHomingtrig000A;
angularPosFirst2Images = raw_data.image000.angular_position;
if isfield(raw_data,'image001') % if image001 exists
    trig001 = de2bi(raw_data.image001.trigger_code,ndigits);
    angHomingtrig001A = trig001(:,1); % first column is top angular homing trigger line (second column is hyper/hypo switch trig)
    angHomingTrigAFirst2Images = [angHomingTrigAFirst2Images(:); angHomingtrig001A(:)];
    angularPosFirst2Images = [angularPosFirst2Images(:); raw_data.image001.angular_position(:)];
end
[angular_pos_to_subtract, AngularHomingDetected] = angular_homing_falledge(angHomingTrigAFirst2Images, angularPosFirst2Images,-90);

if strcmp(AngularHomingDetected,'no')
    angHomingtrig000B = trig000(:,3); % third column is bottom angular homing trigger line
    angHomingTrigBFirst2Images = angHomingtrig000B;
    if isfield(raw_data,'image001') % if image001 exists
        angHomingtrig001B = trig001(:,3);  % third column is bottom angular homing trigger line
        angHomingTrigBFirst2Images = [angHomingTrigBFirst2Images(:); angHomingtrig001B(:)];
    end
    [angular_pos_to_subtract, AngularHomingDetected] = angular_homing_falledge(angHomingTrigBFirst2Images, angularPosFirst2Images,80.8594);
end



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
