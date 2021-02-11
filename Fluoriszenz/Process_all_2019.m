%% analyze fluorescenc data for different pressures
%
%	3-Oct-2019
warning('off', 'all')
clear all;
close all;
warning('off', 'all')

samplename='G04_';

files= dir([samplename '*.mat']);

%% compute pressures & intensities
for m=1:numel(files)
    load (files(m).name)
    run_no(m)= sscanf (files(m).name,[samplename '%d']);
    pressure(m)= (sum(p_periodic(:))+p_pulsed)/5;
%     figure (10)
%     subplot(4,3,m)
%     imshow(img_pulsed,[]);
%     title(sprintf('p: %f',1000*(1.0-(pressure(m)-0.5)/4.0)));
%    figure (11)
%     for n=1:4
%         subplot(4,numel(files),(n-1)*numel(files)+m)
%         imshow(img_periodic(:,:,n),[])
%     end
    mean_pulsed(m)= stats_pulsed.img_mean;
    rms_pulsed(m)= stats_pulsed.img_rms;
    mean_periodic(:,m)= stats_periodic.img_mean;
    rms_periodic(:,m)= stats_periodic.img_rms;
end
pressure= 1000 * (1.0 - (pressure - 0.5)/4.0);  % in mbar

figure (1)
clf
plot(run_no,pressure,'*')
xlabel ('Run #')
ylabel('Pressure [mbar]')
title ('Recorded Pressures')

figure (2)
clf
errorbar (pressure,mean_pulsed,rms_pulsed,'*')
xlabel ('Pressure [mbar]')
ylabel ('Pulsed Intensity []')
title ('Pressure Dependent Intensity');

figure (3)
clf
hold on
for n=1:4
    errorbar (pressure,mean_periodic(n,:),rms_periodic(n,:),'*')
end
xlabel ('Pressure [mbar]')
ylabel ('Periodic Intensity []')
hold off

%% determine ROI for image analysis
load (files(5).name)
disp('please use your cursor to select the area of interest in the figure and double click')
[small,rect]= imcrop(img_periodic(:,:,1),[]);

%% find image with dark reference: minimum intensity values
[min_val,dark_index]= min(mean_pulsed);

%% find image with atmospheric pressure reference
p_max= 0;
for m=1:numel(files)
    if m ~= dark_index & pressure(m) > p_max
        p_max= pressure(m);
        ref_index= m;
    end
end

%% compute intensity ratio (requires reference image ...)
%figure (4)
load (files(dark_index).name);
dark= imcrop(img_pulsed(601:1200,:),rect);
load (files(ref_index).name);
ref= imcrop(img_pulsed(601:1200,:),rect);
n= 0;
i_ratio= [];
p_ratio= [];
for m=1:numel(files)
    if m ~= dark_index
        n= n + 1;
        load (files(m).name);
        small= imcrop(img_pulsed(601:1200,:),rect);
        i_ratio(n)= nanmedian(nanmedian((small-dark)./(ref-dark)))
        p_ratio(n)= pressure(m)
    end
end

figure (5)
clf
subplot (2,1,1)
plot (p_ratio,i_ratio,'*');
xlabel ('Pressure [mbar]');
ylabel ('Intensity Ratio []');
subplot (2,1,2)
plot (p_ratio,1./i_ratio,'*');
xlabel ('Pressure [mbar]');
ylabel ('Inverse Ratio []');

%% compute lifetime from pulse ratio (no reference image necessary !)
load (files(dark_index).name);
dark_1= imcrop(img_pulsed(1:600,:),rect);
dark_2= imcrop(img_pulsed(601:1200,:),rect);
n= 0;
p_ratio= [];
lifetime= [];
for m=1:numel(files)
    load (files(m).name);    
    if m ~= dark_index
        n= n + 1;
        small_1= imcrop(img_pulsed(1:600,:),rect);
        small_2= imcrop(img_pulsed(601:1200,:),rect);
        ratio= median(median((small_1-dark_1)./(small_2-dark_2)))
        lifetime(n)= FindTimeConstant_us(ratio,ExposureTime/1000);
        p_ratio(n)= pressure(m);
    end
end

figure (6)
clf
subplot (2,1,1)
plot (p_ratio,lifetime,'*');
xlabel ('Pressure [mbar]');
ylabel ('Lifetime [us]');
subplot (2,1,2)
plot (p_ratio,1./lifetime,'*');
xlabel ('Pressure [mbar]');
ylabel ('Inverse Lifetime [1/us]');

%% compute phase shift for periodic excitation (no reference, no dark reference !)
n= 0;
p_ratio= [];
phase= [];
for m=1:numel(files)
    load (files(m).name);
    if m ~= dark_index
        n= n + 1;
        for k=1:4
            small_avg(k)= median(median(imcrop(img_periodic(:,:,k),rect)));
        end
        small_avg
        ratio(n)= (small_avg(2)-small_avg(4))./(small_avg(1)-small_avg(3));
        p_ratio(n)= pressure(m);
    end
end

figure (7)
omega= pi/(ExposureTime/1000);
plot (p_ratio,ratio/omega,'*');
xlabel ('Pressure [mbar]');
ylabel ('Liftime [us]');

figure (8)
plot (p_ratio,omega./ratio,'*');
xlabel ('Pressure [mbar]');
ylabel ('Inverse Liftime [1/us]');
