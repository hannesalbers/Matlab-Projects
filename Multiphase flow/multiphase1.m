%% Tabula rasa
clear variables
close all
clc

%% Enter constants
freq = 10000;                                   %sampling frequency [Hz]
duration = 10;                                  %measuring duration [s]
d = 1.5e-3;                                     %distance between measuring planes [m]
nx = 32;                                        %total number of receivers [-]
ny = 16;                                        %total number of transmitters [-]
A_ch = 0.05^2;                                  %channel cross section [m2]
A_mp = A_ch/ny^2;                               %measuring point cross section [m2]
gflow = [0.3 1.5 7.5 25.0];                    %gas flow rate [l/min]
filename = {'20191107_Trial1_0p3lnminAir_50p0lminWater' '20191107_Trial2_1p5lnminAir_51p0lminWater' ...
    '20191107_Trial3_7p5lnminAir_54p0lminWater' '20191107_Trial4_25p0lnminAir_61p0lminWater'};     %filenames
n = length(filename);                           %number of measurements [-]

fid = fopen(['20191107_Trial5_0lnminAir_51p0lminWater' '.dat'],'r');      %open file
M = fread(fid,nx*ny*freq*duration,'uint16'); %data read in
fclose(fid);                                %closes opened file
M = reshape(M,[nx ny duration*freq]);       %arranges data in a 3D-matrix with dimensions [ 32_receiving_lines(sensor_plane1&2) 16_transitting_lines frames_over_time]
M = mean(M, 3);

%% Calculate void fraction distribution
for i = 1:n
    fid = fopen([filename{i} '.dat'],'r');      %open file
    A = fread(fid,nx*ny*freq*duration,'uint16'); %data read in
    fclose(fid);                                %closes opened file
    A = reshape(A,[nx ny duration*freq]);       %arranges data in a 3D-matrix with dimensions [ 32_receiving_lines(sensor_plane1&2) 16_transitting_lines frames_over_time]
    
    M = max(A,[],3);                            %maximum for each measuring point
                                                %calculates void fraction time resolved
                  
    A = 1 - A./repmat(M,[1 1 duration*freq]);
    A1 = A(17:32,:,:);                          %dumps data of sensor plane 2 (17-32)
    
    void(i,:,:) =mean(A1,3);                    %average over time for each measuring point
    
    prof1(:,i)=mean(void(i,5:12,:),2);          %average over centerlines 5 to 12 in transmitting wire direction
    prof2(:,i)=mean(void(i,:,5:12),3);          %average over centerlines 5 to 12 in receiving wire direction
end

% PLOT RESULTS
c=[0 0 0];

set(0,'DefaultAxesLineStyleOrder','-s|-o|-v|-+|-x|-*')
set(0,'defaultaxescolororder',c)

figure('position',[100 100 1000 800],'color',[1 1 1])
plot(prof1);
grid on;
set(gca,'fontsize',20)
xlabel('Measuring Point No.','fontsize',24,'fontweight','bold')
ylabel('Void Fraction[-]','fontsize',24,'fontweight','bold')
title('Void fraction profile trans')
legend('Test 1','Test 2','Test 3','Test 4')

figure('position',[100 100 1000 800],'color',[1 1 1])
plot(prof2);
grid on; 
set(gca,'fontsize',20)
xlabel('Measuring Point','fontsize',24,'fontweight','bold')
ylabel('Void Fraction[-]','fontsize',24,'fontweight','bold')
title('Void fraction profile rec')
legend('Test 1','Test 2','Test 3','Test 4')
figure('position',[100 100 1000 800],'color',[1 1 1])

%% Calculation of Bubble Velocity
clear A M;

figure('position',[100 100 350 300],'color',[1 1 1])%creates a figure

for k=1:n                                           %number for gas flow rates (that are required)
    fid = fopen([filename{k} '.dat'],'r');          %open file
    A = fread(fid,nx*ny*freq*duration,'int16');     %data read in
    fclose (fid);                                   %closes opened file
    A = reshape(A,[nx ny duration*freq]);           %arranges data in a 3D-matrix with dimensions [ 32_receiving_lines(sensor_plane1&2) 16_transitting_lines frames_over_time]
        
    M = max(A,[],3);                                %maximum for each measuring point
    
                                       %calculates void fraction time resolved
    A =1 - A./repmat(M,[1 1 duration*freq]);       %BONUS
    A(find(A<0.2))=0;                               %threshold to avoid crosscorrelation of void fractions below 0.2
   
    A1 = permute(A(1:16,:,:),[3 1 2]);              %cut data of sensorplane 1 and shift time dimension to the first matrix dimension (for later functions use)
    A2 = permute(A(17:32,:,:),[3 1 2]);             %cut data of sensorplane 2 and shift time dimension to the first matrix dimension (for later functions use)
    
    for i = 1:16                                    %loop through each measuring point individually
        for j=1:16
            
            if max(A1(:,i,j))>0.2 && max(A2(:,i,j)>0.2)     %checks if bubbles crossed the measuring point
                [c lag] = xcorr(A1(:,i,j),A2(:,i,j));       %cross-correlation through time dimension
                [val id] = max(c);                          %find maximum of the cross-correlation function
                vel(i,j) = d/(lag(id)/freq);                %velocity calculation
                if(abs(vel(i,j))> 10 || vel(i,j) <-1)       %threshold cut-off velocities
                    vel(i,j) = 0;
                end
            else
                vel(i,j)=0;                                 %blank maximum time shift, if no bubbles crossed
            end
        end
    end
    
    vj(k,:,:)=vel;                                  %saves the velocity to calculate superficial velocity

    %% Calculate superficial velocity
    Jvelmp(k,:,:)=void(k,:,:).*vj(k,:,:);
    Jvel(k)=nansum(Jvelmp(k,:).*A_mp)./A_ch;        %Superficial Velocity [m/s]
    Jvel_theor(k) = gflow(k)*(293/273)./(60*1000*A_ch);
    
    %% PLOT RESULTS
    subplot(2,2,k)
    contourf(vel); colormap jet; shading interp;
    caxis([0 1.5]);
    xlabel('Measuring Point No.')
    ylabel('Measuring Point No.')
    colorbar
    title([ 'Abb. 4 k ' 'Superficial gas velocity ' sprintf('%.3f',Jvel(k)) ' (m/s)'])
    axis equal
    xlim([1 16])
    ylim([1 16])
end

%% Superficial Velocity Theory vs. Measured
figure('position',[100 100 1000 800],'color',[1 1 1])
plot(Jvel_theor,Jvel,[0 0.4],[0 0.4]);
set(gca,'fontsize',20)
grid on;
axis([0 0.4 0 0.4]);
xlabel('Superficial Velocity [m/s]','fontsize',24,'fontweight','bold')
ylabel('Measured Superficial Velocity [m/s]','fontsize',24,'fontweight','bold')