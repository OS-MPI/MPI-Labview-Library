
clear; clc; close all;


file = '2019_04_29_18_40_14_test.h5';
RXmag = hdf5read(file, '/RXmag');
RXphase  = hdf5read(file, '/RXphase');
shift = hdf5read(file, '/shift');

cropramppts = 4;

[num_proj, RXmag_projs, RXphase_projs, shift_projs] = reshape_projs(RXmag, RXphase, shift);

[IR_crop] = recon(cropramppts,num_proj,RXmag_projs,RXphase_projs);


figure, imagesc(IR_crop)