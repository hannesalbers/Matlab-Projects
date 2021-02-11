%%
close all
clear all
clc

%addpath('./Images')
nImag = 100;
curr_dir = pwd;
sig_folder = [curr_dir,'/Images/'];                                   %\specifiy the directory of your images 
Def_name = 'VimbaImage_';
Fs =206.9898              
%Specify the sampling frequency of the camera
for i = 1:nImag
    im_raw = rgb2gray(imread([sig_folder,Def_name,num2str(i)],'bmp'));       %convert to gray scale since we only care about the intensity
   im_raw = imresize(im_raw,0.6);                                           %decrease the resolution.
   im_take = im2double(im_raw);                                             %convert the data into double.
   SM(:,i) = im_take(:);                                                    %Vectorize the image data, assign it to matrix SM (Snapshot Matrix)
end
disp('Read data')
%% The following block is an example of how to output a video from the snapshot matrix
outputVid = VideoWriter([sig_folder,['VideoOut','.avi']]);
outputVid.FrameRate = 5;                                                    %specify the fps
open(outputVid)
SM = MinMaxNormalize(SM);                                                   %necessary to normalize the data from 0 to 1!
for i = 1:nImag
     writeVideo(outputVid,reshape(SM(:,i),size(im_take,1),size(im_take,2))); %you should reshape the image back to the original image dimension
end  
close(outputVid)
%%
[POD_modes, lam,eta]=PODsvd_ExMI(SM);                                       %get the POD modes, eigenvalues, and time coefficient matrix
%keep in mind the n-th mode is the n-th column vector of the POD_modes
%e.g the first mode is POD_modes(:,1), whereas the corresponding time
%coefficient of the mode is eta(1,:)
%% plot the relative strength of each lambda write the physical meaning of each plot and add the suitable axis label
figure;
plot(lam/sum(lam))
figure;
plot(cumsum(lam)/sum(lam))
%% plot the PSD from mode 1 to mode N (up to your choice)
%plot the PSD of the time coefficient of mode 1 to mode N, shift the y-axis value of each
%mode's PSD such that the graphs do not overlap!! 

%% Mode visualisation example
%for easier visualisation you can output a video of some modes e.g
%take the first POD_modes, take the corresponding time coefficient vector
%for the first mode.

outputVid = VideoWriter([sig_folder,['Mode_1','.avi']]);
outputVid.FrameRate = 5;                                                    %specify the fps of the video 
open(outputVid)
modTmodTake = 1;                                                                %take the first mode
Reconstruct(:,i) = zeros(size(POD_modes,1),nImag);    
for i = 1:nImag
    Reconstruct(:,i) = eta(modTake,i) * POD_modes(:,modTake)...
        + Reconstruct(:,i);                                                 %create a new matrix for the chosen mode, i-th column = POD_modes*i-th element of the time coefficient vector of that mode
end 
Reconstruct = MinMaxNormalize(Reconstruct);                                 %normalize the values to 0 and 1
for i = 1:nImag
     writeVideo(outputVid,reshape(Reconstruct(:,i),size(im_take,1),size(im_take,2)));

end 
close(outputVid)

%% Reconstruct a new video composed of the two modes whose peak frequency is the same!


