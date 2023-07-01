%%Mpi continuous mode simulation  / recon
%%3/9/2019
%%Clarissa Zimmerman Cooley
% changed variable names when switched to CC mode - EM

%%%basically iradon in gen recon (does not use current waveform or
%%%filtering)

%%% 4/9/2019
%%%option to load A matrix instead of regenerating


clear;
% close all;
displaysim = 0;
loadA = 0;
spiosize = '18';

load(['sample_',spiosize,'nm']);
% load('image_I_0,5mg_ml_CCmode_SRS_G1000_Rcs_2');
shift_current_scale_factor = 20;

% DECIDE WHICH CORRECTION TO USE:
a = 1;
b = 0;
pr = proj.lin_amp3_corr(1+a:end-b,:); 
proj_crop = pr; % OPTION A
pr = abs(proj.lin_complex3 - mean(proj.lin_complex3(:))+1);  proj_crop = -1*(pr - max(pr(:))); % OPTION B

% P = sum(pr,1);
% P = pr(end,:);
% for ii = 1:26
% proj_crop(:,ii) = pr(:,ii) - P(ii); 
% end
% proj_crop = pr;
% proj_crop = -1*(pr - max(pr(:))); % OPTION B
% 
% proj_crop = -1*(pr - max(pr(:))); % OPTION B

% figure; imagesc(proj_crop)
% figure; plot(pr); title([spiosize,'nm projections']);
figure; plot(proj_crop); title([spiosize,'nm projections']);
figure; plot(proj_crop(:,[1,13,26])); title([spiosize,'nm projections']);


theta = proj.theta;
Nproj = size(pr,2);
Nsamp = size(pr,1);


% temp measured shift curren test:
measured_shift_I_avg = repmat(proj.shift_amp_vals',1,26)/proj.shift_amp_vals(1)*2.4986*10;

% measured_shift_I_avg = shift_current_scale_factor*squeeze(mean(proj.measured_shift_VRcs,1));
measured_shift_I_avg = measured_shift_I_avg(1+a:end-b,:);
Ishift = measured_shift_I_avg(:);
Ishift_crop = reshape(Ishift, [], Nproj);


% proj_crop = repmat(pr(:,14),1,26);
figure;
FOV = 2;
xyi = linspace(-FOV/2,FOV/2,length(proj_crop));
IR = iradon(abs(proj_crop), -theta+180, 'linear','Ram-Lak',1,size(proj_crop,1));
IR = fliplr(IR);
imagesc(xyi,xyi, ((IR))); axis square; colormap gray;
title([spiosize,'nm iradon']);
xlabel('X (cm)');  ylabel('Y (cm)');

% figure; plot(proj_crop); title([spiosize,'nm projections']);


%% set up model for reconstruction


%%%IMPORTANT - this determines how far the FFL is swept in simulation by
%%%scaling the measured current in the shift coils
Ishift_scaling = 1/30;  %% assumes 30 amps results in 1 cm shift
FFLshift = flipud(Ishift_crop)*Ishift_scaling;  %% FFL shift distance during acq.
negshift = abs(min(Ishift_crop(:,1)))*Ishift_scaling; %cm
posshift = abs(max(Ishift_crop(:,1)))*Ishift_scaling; %cm

%estimated FFL starting offset from isocenter (given that that projections
%are centered)
shift_offset = posshift - negshift;  %cm in neg direction

disp(['shifted from ',num2str(posshift), 'cm to ', num2str(negshift),'cm']);
disp(['estimated FFL offset ',num2str(10*shift_offset), 'mm']);



% encoded FOV
FOV = posshift + negshift;

xyi = linspace(-negshift, posshift, Nsamp)-shift_offset/2;
% temp2 = repmat(xyi', 1, Nproj);
% FFLshiftoffset = temp2;

FFLshiftoffset = FFLshift - shift_offset/2;




res_orig = FOV/(Nsamp-1);

%% sample size decrease factor for simution
res_fact = 2;
res_sim = res_orig/res_fact;
xyisim1 = [xyi(1):res_sim:xyi(end)];


%% FOV scale factor for simulation
fov_fact = 1.5;
FOVsim = fov_fact*FOV;
temp = [xyisim1(end)+res_sim:res_sim:FOVsim/2];
xyisim2 = [-fliplr(temp), xyisim1, temp];
xyisim = xyisim2;
xyiorig_ind = [numel(temp)+1:numel(temp)+numel(xyisim1)];


Nsim = numel(xyisim);
Nrecon = Nsamp; %%matrix dimension for FFL sweep simulation


%%% relates meaured FFL shift to matrix index for simulation
for pp = 1:Nproj
    for ss = 1:Nsamp
        [val,I] = min(abs(xyisim2-FFLshiftoffset(ss,pp)));
        vals(ss,pp) = val;
        fflpos_ind(ss,pp) = I;
    end
end



%%gaussian kernal to smear FFL -- this should be be replaced with
stdcm = 0.15;  %%cm
ceil(stdcm/res_sim)/2;
H = fspecial('gaussian', [Nsim,Nsim],ceil(stdcm/res_sim)/2);


% %%check FFL
% figure; subplot(2,1,1);
% imagesc(xyisim,xyisim,H); axis square;
% title('FFL simulation');
% subplot(2,1,2); plot(xyisim,H);



if loadA == 0
    
    %% calculate model
    
    FFL_mat = zeros(Nrecon, Nrecon, Nsamp, Nproj);
    FFLnow = zeros(Nsim,Nsim);
    
    
    h = figure;
    for pp = 1:Nproj
        display(['simulating sweep ', num2str(pp),' out of ',num2str(Nproj)]);
        for ss = 1:Nsamp
            
            %%creates FFL as a single line of 1s in a zero matrix
            
            FFLnow = zeros(Nsim,Nsim);
            FFLnow(: , fflpos_ind(ss,pp) ) = ones(Nsim, 1);
            
            %apply gaussian filter to simulate imperfect FFL
            FFLnow = imfilter(FFLnow,H,'replicate');
            
            FFLnowrot = imrotate(FFLnow,theta(pp),'crop');
            FFLnowrot_crop = FFLnowrot(xyiorig_ind,xyiorig_ind);
            FFLinterp = smart_interp(FFLnowrot_crop, Nrecon);
            
            FFL_mat(:,:,ss,pp) = FFLinterp;
            
            
            if displaysim == 1
                imagesc(xyisim2,xyisim2,FFL_mat(:,:,ss,pp));
                title('selection field'); axis square; colorbar; caxis([0,0.5/10])
                %             pause(0.01);
                
                drawnow
                filename = 'modelbasedrecon_gif';
                % Capture the plot as an image
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                if pp == 1
                    imwrite(imind,cm,filename,'gif','DelayTime',0, 'Loopcount',inf);
                else
                    imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
                end
                
            end
            
        end
    end
    
    %%reshape FFL matrix for each data point to form model A
    
    
    A = reshape(FFL_mat, Nrecon^2 , []).';
    
else
    load('A_matrix.mat');
end

%% Reconstruction

%% cg reconstuction

b = proj_crop(:);  %% data

%%A is not square so solve: A'Ax = A'b
Asquare =  A'*A;
B = A'*b;
% Bsim = A'*b2noise;

figure;
% subplot(2,1,1)
for iter = 1:18
    X = pcg(Asquare,B,[],iter);
    datarecon(:,:,iter) = reshape(X, Nrecon,Nrecon);
    imagesc(xyi,xyi,fliplr(abs(datarecon(:,:,iter)))); axis square;
    title(['cg iter ', num2str(iter)]); xlabel('X (cm)');  ylabel('Y (cm)'); colormap gray;
    pause(0.5)
end
