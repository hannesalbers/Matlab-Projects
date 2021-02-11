function [xy_grid,uv_vecs,peaks,valid,loopdata] = PIV_base (I1,I2,mode,wxy,sxy,oxy,gxy,loopdata,threshold)

%   Last Update: 01-AUG-17

%   ------ List of Changes ------
%
%   inserted "padarray" call in slidesum
%   moved subpixel interpolation routines into main function
%   removed global variables
%   removed try/catch construct

%   -----------------------------

%   PIV 2-component, 2-D Particle Image Velocimetry
%
%   I1,I2:          images (no stacked slices allowed - use frame averaging to combine data)
%
%   mode:           0  = PIV with SINC subpixel interpolation
%                   1  = PIV with GAUSS subpixel interpolation
%                   -1 = PIV WITHOUT subpixel interpolation
%
%   wxy:            [wx,wy]     interrogation window size
%   sxy:            [sx,sy]     max +/- search size in x and y
%   oxy:            [ox,oy]     window offset (constant across image) 
%                   [nx,ny,2]	window offset for individual grid points
%   gxy:            [gx,gy]     window shift
%                   [nx,ny,2]	window positions in x and y
%   loopdata        structure containing two arrays nmaps, cmaps (leave empty if not needed)
%       cmaps:      composite image of correlation maps to be averaged with current data 
%       nmaps:      [nx,ny]     array containing number of correlation maps summed at every grid point
%   threshold:      threshold value for accepting correlation map (0...1)
%
%   xy_grid:        [nx,ny,2]   grid location of output results
%   uv_vecs:        [nx,ny,2]   (map averaged) u- and v- velocities
%   peaks:          [nx,ny]     map of peak heights 
%   valid:          [nx,ny]     map of flags: -1=invalid, 0=integer shifts ok, 1=subpixel shifts ok
%   loopdata:       output structure for loop processing (contents see above)
%
%   (c) T.Roesgen IFD/ETHZ
%

%----------------------- global variables -----------------------------------

CSZX= 2*sxy(1) + 1;                     % x-size of correlation map, always odd
CSZY= 2*sxy(2) + 1;                     % y-size of correlation map, always odd

MAXITER= 20;

%-----------------------------------------------------------------------------

if size(I1) ~= size(I2)
    error ('Images should be same size ...');
end

if abs(mode) > 1
    error ('Select proper subpixel processing mode: 0=Sinc, 1=Gauss, -1=None ...');
end

if wxy(1) < 1 || wxy(2) < 1
    error ('Wrong window size wxy ...');
end

if sxy(1) < 1 || sxy(2) < 1
    error ('Search range sxy too small ...');
end

if 2*sxy(1) >= wxy(1) || 2*sxy(2) >= wxy(2)
    error ('Search range sxy too large ...');
end

if gxy(1) < 1 || gxy(2) < 1
    error ('Wrong window spacing dxy ...');
end

[height,width,nslice]=size(I1);

if nslice > 1
    disp ('Only processing first slice of image stack !');
end

% make sure image data are floating point
I1=double(I1);
I2=double(I2);

% fill in arrays of interrogation grid points
if numel(gxy) ==2
    % grid spacing is prescribed
    nx= floor((width-wxy(1))/gxy(1)) + 1;
    ny= floor((height-wxy(2))/gxy(2)) + 1;
    [xgrid,ygrid]= meshgrid((0:nx-1)*gxy(1)+1,(0:ny-1)*gxy(2)+1);
    xy_grid(:,:,1)= xgrid + wxy(1)/2 - 0.5;
    xy_grid(:,:,2)= ygrid + wxy(2)/2 - 0.5;
else
    % full grid is prescribed
    [ny,nx,nz]= size(gxy);
    xy_grid= gxy;
    xgrid= round(gxy(:,:,1) - wxy(1)/2 + 0.5);
    ygrid= round(gxy(:,:,2) - wxy(2)/2 + 0.5);
end

% size of smaller interrogation window
wszx= wxy(1)-CSZX+1;
wszy= wxy(2)-CSZY+1;
n_wsz= wszx*wszy;
            
% initialize output arrays
uv_vecs= NaN(ny,nx,2);
peaks= NaN(ny,nx);
valid= -int8(ones(ny,nx));
if isempty(loopdata)
    loopdata.cmaps= zeros(ny*CSZY,nx*CSZX);
    loopdata.nmaps= zeros(ny,nx);
end

% compute edge coordinates for PIV grid
if  numel(oxy) == 2
    ox2= repmat(round(oxy(1)/2),ny,nx);
    oy2= repmat(round(oxy(2)/2),ny,nx);
elseif  numel(oxy) == numel(gxy)
    ox2= round(oxy(:,:,1)/2);
    oy2= round(oxy(:,:,2)/2);
else
    error ('Incompatible specification for oxy ...');
end

% precompute sliding means and variances (is faster than for individual
% subwindows)
sumAn= slidesum(I1,wxy(2)-CSZY+1,wxy(1)-CSZX+1);
varAn= slidesum(I1.^2,wxy(2)-CSZY+1,wxy(1)-CSZX+1) - (sumAn.^2)./n_wsz;
sumBn= slidesum(I2,wszy,wszx);
varBn= slidesum(I2.^2,wszy,wszx) - (sumBn.^2)./n_wsz;

% initialize coefficient matrices for Gaussian fit
QMAT= [ [ 1 -1 -1  1  1  1]; ...
        [ 1  0 -1  0  0  1]; ...
        [ 1  1 -1  1 -1  1]; ...
        [ 1 -1  0  1  0  0]; ...
        [ 1  0  0  0  0  0]; ...
        [ 1  1  0  1  0  0]; ...
        [ 1 -1  1  1 -1  1]; ...
        [ 1  0  1  0  0  1]; ...
        [ 1  1  1  1  1  1] ];
QQ= QMAT' * QMAT;

% initialize phase masks for FFT shifts
PHASEX= zeros(1,CSZX);
PHASEX(1:sxy(1)+1)= 2*pi*(0:sxy(1))./CSZX;
PHASEX(sxy(1)+2:CSZX)= 2*pi*(-sxy(1):-1)./CSZX;
PHASEX= repmat(PHASEX,CSZY,1);

PHASEY= zeros(1,CSZY);
PHASEY(1:sxy(2)+1)= 2*pi*(0:sxy(2))./CSZY;
PHASEY(sxy(2)+2:CSZY)= 2*pi*(-sxy(2):-1)./CSZY;
PHASEY= repmat(PHASEY',1,CSZX);

PHASEXX= PHASEX.*PHASEX;
PHASEXY= PHASEX.*PHASEY;
PHASEYY= PHASEY.*PHASEY;

% loop over all subwindows
for rcindx=1:nx*ny
    
    [rows,cols]= ind2sub([ny,ny],rcindx);
    
    % extract smaller subwindow from first image
    x1= xgrid(rows,cols)-ox2(rows,cols)+sxy(1);
    y1= ygrid(rows,cols)-oy2(rows,cols)+sxy(2);
    A= I1(y1:y1+wxy(2)-CSZY,x1:x1+wxy(1)-CSZX);
    
    % extract full subwindow from second image
    x2= xgrid(rows,cols)+ox2(rows,cols);
    y2= ygrid(rows,cols)+oy2(rows,cols);
    B= I2(y2:y2+wxy(2)-1,x2:x2+wxy(1)-1);
    
    % compute normalizing denominator
    denom= sqrt(varAn(y1,x1) .* varBn(y2:y2+CSZY-1,x2:x2+CSZX-1));
    
    % compute covariance for windows
    fft_A= fft2(A,wxy(2),wxy(1));   % padded for size of subwindow B
    fft_B= fft2(B);
    corrAB= ifft2(conj(fft_A).*fft_B,'symmetric');
    
    % compute (co-)variances
    covAB= corrAB(1:CSZY,1:CSZX) - sumAn(y1,x1).*sumBn(y2:y2+CSZY-1,x2:x2+CSZX-1)./n_wsz;
    
    % compute correlation map
    normcorr= covAB ./ denom;
    
    % insert previous map estimate if available
    rr= (rows-1)*CSZY;
    cc= (cols-1)*CSZX;
    if not(isempty(loopdata.cmaps))
        newmap= loopdata.cmaps(rr+1:rr+CSZY,cc+1:cc+CSZX);
        newcount= loopdata.nmaps(rows,cols);
    else
        newmap= zeros(CSZY,CSZX);
        newcount= 0;
    end
    
    % add current correlation map if peak threshold is exceeded
    if max(normcorr(:)) >= threshold
        newmap= newmap + normcorr;
        newcount= newcount + 1;
    end
    
    % store updated correlation map
    loopdata.cmaps(rr+1:rr+CSZY,cc+1:cc+CSZX)= newmap;
    loopdata.nmaps(rows,cols)= newcount;
    
    % peak analysis only if there are non-zero elements in cmaps
    if newcount > 0
        
        % Find value and position (xpeak, ypeak) of the correlation peak
        [peaks(rows,cols),peak_ind]= max(newmap(:)/newcount);
        [ypeak,xpeak]= ind2sub (size(newmap),peak_ind);
        
        % compute shift index
        uv_vecs(rows,cols,1)= xpeak - sxy(1) - 1 + 2*ox2(rows,cols);
        uv_vecs(rows,cols,2)= ypeak - sxy(2) - 1 + 2*oy2(rows,cols);
        valid(rows,cols)= 0;
        
        % correlation maximum is not an a edge
        if mode >= 0 && xpeak>1 && xpeak<CSZX && ypeak>1 && ypeak<CSZY
            
            % subpixel interpolation
            if mode == 0
            
                % sinc subpixel interpolation            
                % must subtract [1,1] to work with Fourier shift indices
                fcorr= fft2(newmap/newcount);
                low= [xpeak,ypeak] - 1 - 0.75;
                high= [xpeak,ypeak] - 1 + 0.75;
                
                % Newton iteration to find subpixel shift
                new_shift= [1,1];
                shift= [xpeak,ypeak] - 1;
                iter=0;
                while norm(new_shift) > 1.e-6 && iter < MAXITER
                    iter= iter + 1;

                    % compute x-, y- and total phase factors
                    totphase= PHASEY .* repmat(shift(2),CSZY,CSZX) + ...
                              PHASEX .* repmat(shift(1),CSZY,CSZX);
                    
                    % compute product of spectrum and phase factors
                    ff= fcorr .* exp(i*totphase);
                    realf= real(ff);
                    imagf= imag(ff);
                    
                    % function evaluation
                    intpeak= sum(realf(:));
                    
                    % evaluate function gradient: vector
                    grad(1)= dot(imagf(:),PHASEX(:));
                    grad(2)= dot(imagf(:),PHASEY(:));
                    
                    % evaluate Hessian matrix
                    hess(1,1)= dot(realf(:),PHASEXX(:));
                    hess(1,2)= dot(realf(:),PHASEXY(:));
                    hess(2,1)= hess(1,2);
                    hess(2,2)= dot(realf(:),PHASEYY(:));
                    
                    % compute determinant of Hessian (local curvature)
                    det= hess(1,1)*hess(2,2)-hess(1,2)*hess(2,1);
                    
                    % compute inverse Hessian (determinant inserted later)
                    ihess(1,1)= hess(2,2);
                    ihess(1,2)= -hess(1,2);
                    ihess(2,1)= -hess(1,2);
                    ihess(2,2)= hess(1,1);
                    
                    % update shift estimate (abs(det) used to prevent runaways)
                    new_shift= (ihess*(grad'))/abs(det);
                    shift= shift - new_shift';
                    
                    % limit to acceptable coordinate range
                    shift= min(high,max(low,shift));
                end
                intpeak= intpeak / (CSZY*CSZX);
                
                % correct for Fourier indexing
                shift= shift + 1;

            else
                
                % Gaussian subpixel interpolation
                intpeak= 0;
                shift= [0,0];
                
                ind= 1;
                for jj=-1:1
                    for ii=-1:1
                        rhs(ind)= log(max(newmap(ypeak+jj,xpeak+ii)/newcount,1.e-12));
                        ind= ind + 1;
                    end
                end
                if var(rhs) == 0
                    % right hand side was clipped everywhere: no good
                    iter= MAXITER;
                else
                    qr= QMAT' * rhs';
                    coeffs= QQ \ qr;
                    mmat= [ [2*coeffs(4) coeffs(5)]; [coeffs(5) 2*coeffs(6)]];
                    mrhs= [-coeffs(2); -coeffs(3)];
                    qvec= (mmat \ mrhs)';
                    if norm(qvec) > 1
                        % interpolated displacement is too large: no good
                        iter= MAXITER;
                    else
                        % interpolation is valid
                        iter= 1;
                        shift= qvec + [xpeak,ypeak];
                        % compute surface curvature criterion for Gaussian peak
                        intpeak= exp(coeffs(1) + coeffs(2)*qvec(1) + coeffs(3)*qvec(2) + ...
                            coeffs(4)*qvec(1)*qvec(1) + coeffs(5)*qvec(1)*qvec(2) + coeffs(6)*qvec(2)*qvec(2));
                    end
                end
            end
            if iter < MAXITER && intpeak/peaks(rows,cols) > 0.98 && intpeak/peaks(rows,cols) < 2.0
                % good data: final interpolated shifts superseedes integer shift estimate
                uv_vecs(rows,cols,1)= shift(1) - sxy(1) - 1 + 2*ox2(rows,cols);
                uv_vecs(rows,cols,2)= shift(2) - sxy(2) - 1 + 2*oy2(rows,cols);
                peaks(rows,cols)= intpeak;
                valid(rows,cols)= 1;
            end
        end
    end
end

return

%----------------------------------------------------------------------------------------

function ss= slidesum(a,n,m)
% sliding window summation; averaging window size is [n,m]
[na,ma]= size(a);
a1= cumsum(padarray(a,[1,1],0,'pre'),1);
a2 = a1(n+1:na+1,:)-a1(1:na+1-n,:);
a3 = cumsum(a2,2);
ss = a3(:,m+1:ma+1)-a3(:,1:ma+1-m);
return

