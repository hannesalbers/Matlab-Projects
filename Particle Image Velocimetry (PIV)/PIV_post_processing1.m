%% Group 4:  PIV post processing

close all;
clear all;
clc;

%% Import images  

 load backup_run.mat         % time delay between pictures: 1.9s 
%load run_.mat                 % time delay between pictures: 0.8s

%% Background subtraction

N = size(images,3);             % number of images

% figure
% imshow(images(:,:,2),[])
% title('With Background')
%%
minimage = min(images,[],3);    % get minimage

for i=1:N         
    images(:,:,i) = images(:,:,i)-minimage(:,:);
end

figure
imshow(minimage,[])
title('Background')

figure
imshow(images(:,:,2),[])
title('Background substracted')

%% Evaluate with PIV_base without correlation averaging

for i=1:N-1       % loop through all image pairs
    disp(i)
    [xy_grid,uv_vecs(:,:,:,i),peaks,valid,loopdata] = PIV_base(images(:,:,i),images(:,:,i+1),0,[64,64],[10,10],[0.5,0.5],[16,16],[],0.3);

end

%% Evaluate with PIV_base with correlation averaging
%function[xy_grid,uv_vecs,peaks,valid,loopdata] = PIV_base (I1,I2,mode,wxy,sxy,oxy,gxy,loopdata,threshold)

for i=1:N-1       % loop through all image pairs
    disp(i)
    if i == 1
        [xy_grid,uv_vecs(:,:,:,i),peaks,valid,loopdata] = PIV_base(images(:,:,i),images(:,:,i+1),0,[64,64],[10,10],[0.5,0.5],[16,16],[],0.3);  
    else
        [xy_grid,uv_vecs(:,:,:,i),peaks,valid,loopdata] = PIV_base(images(:,:,i),images(:,:,i+1),0,[64,64],[10,10],[0.5,0.5],[16,16],loopdata,0.3);

    end
end


%% Convert pixel shifts to velocity

xy_grid = (xy_grid./248).*10;

uv_vecs = (uv_vecs./248)*10;

%% Mean velocity map (correlation map averaging) usinq quiver() 
uv_vecs_mean=mean(uv_vecs,4)

figure()
quiver(xy_grid(:,:,1),xy_grid(:,:,2),uv_vecs_mean(:,:,1),uv_vecs_mean(:,:,2),2)
axis ij
axis equal
ylabel ('Höhe [mm]','Fontsize',18)
xlabel ('Breite [mm]','Fontsize',18)
title('Mean velocity map with correlation averaging','Fontsize',18)

%% Mean velocity map (direct vector map averaging) usinq quiver() 

figure()
quiver(xy_grid(:,:,1), xy_grid(:,:,2),uv_vecs(:,:,1),uv_vecs(:,:,2),2)

axis ij
axis equal
set(gca,'Fontsize',18)
ylabel ('Höhe [mm]','Fontsize',18)
xlabel ('Breite [mm]','Fontsize',18)
title('Mean velocity map without correlation averaging','Fontsize',18)

%% Streamline using streamslice() OK 

figure()
set(gca,'Fontsize',18)
streamslice(xy_grid(:,:,1),xy_grid(:,:,2),uv_vecs(:,:,1,9),uv_vecs(:,:,2))
set(gca,'Fontsize',18)
title('Streamlines without correlation averaging','Fontsize',18)
axis ij
axis equal
ylabel ('Höhe [mm]','Fontsize',18)
xlabel ('Breite [mm]','Fontsize',18)
%%
figure()
set(gca,'Fontsize',18)
streamslice(xy_grid(:,:,1),xy_grid(:,:,2),uv_vecs_mean(:,:,1),uv_vecs_mean(:,:,2))
title('Streamlines with correlation averaging','Fontsize',18)
axis ij
axis equal
ylabel ('Höhe [mm]','Fontsize',18)
xlabel ('Breite [mm]','Fontsize',18)

%% Peak locking using hist() OK
figure()
set(gca,'Fontsize',18)
hist(uv_vecs(:,:,1),10)
title('Peak locking','Fontsize',18)

%% Second order quantities OK

div = divergence(xy_grid(:,:,1),xy_grid(:,:,2),uv_vecs_mean(:,:,1),uv_vecs_mean(:,:,2));
figure()
contourf(div)
set(gca,'Fontsize',18)
title('Second order quantities')
ylabel ('Höhe [mm]')
xlabel ('Breite [mm]')
axis ij
axis equal
