%OCTAVE Video Conversion VERS 7 (FINAL)
%Clear workspace
clear; close all; clc

%Load video package for Octave (remove if using Matlab)
pkg load video;

pause(0.1);

fprintf('\nImporting Video...\n');
tic
  %Imports video file into matlab
  vid = VideoReader("monte_carlo_low.mp4");

  %Collects video dimensions
  vidHeight = vid.Height;
  vidWidth = vid.Width;
  vidNFrames = vid.NumberOfFrames;
toc

%Movie Struct to collect all original frames to run movie using movie function
mov = struct('cdata', zeros(vidHeight,vidWidth,1,'uint8'), 'colormap',gray(255));

%Movie struct for reshaping and conversions
mov2 = struct('cdata', zeros(vidHeight*vidWidth,1,'uint8'), 'colormap',gray(255));

%Movie struct for postprocessing results
mov_LowRank = struct('cdata', zeros(vidHeight*vidWidth,1,'uint8'), 'colormap',gray(255));
mov_Back = struct('cdata', zeros(vidHeight,vidWidth,1,'uint8'), 'colormap',gray(255));
mov_Fore = struct('cdata', zeros(vidHeight,vidWidth,1,'uint8'), 'colormap',gray(255));
mov_all = struct('cdata', zeros(vidHeight,vidWidth,1,'uint8'), 'colormap',gray(255));

%Initializing (All In One) matrix to hold all frames into a 2D space-time matrix
DATAMovie = zeros(vidWidth*vidHeight,vid.NumberOfFrames);

%figure number tracking
nfigure = 1;

[M,N] = size(DATAMovie);

%Collects each frame into cdata
fprintf('\nConverting video to Grayscale...\n');
tic
  k = 1;
  while hasFrame(vid)
      grayImage = rgb2gray(readFrame(vid));

      mov(k).cdata = grayImage; %Stores the gray frames into the structure;
      mov(k).colormap = gray(255); %Applies Gray Colormap to every frame
      k = k+1;

  end
toc

%View original video
fprintf('\nDisplaying Gray Video...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth vidHeight]);
   movie(hf, mov, 1, vid.FrameRate);
   nfigure = nfigure + 1;
toc


fprintf('\nReshaping Gray Video into 2D Matrix...\n');
tic
  for k = 1:vidNFrames
      mov2(k).cdata = reshape(mov(k).cdata,[],1);
      mov2(k).colormap = gray(255);
      DATAMovie(:,k) = [mov2(k).cdata];
  end
toc

%Computing SVD for DATAMovie
fprintf('\nCalculating SVD...\n');
tic
  [U_DATA,S_DATA,V_DATA] = svd(DATAMovie,'econ');
toc

%initializes variables for low rank calculations with energy
s_ = 1;
EnergyTarget = 0.9999;

sum_Sigma_s = 0;
sum_Sigma_all = 0;
Energy = 0;

%Searches for a rank value that will make the Energy bigger than or equal to 99.99% to have good approximation.
fprintf('\nSearching for low rank s number achieving 99.99 percent Energy...\n');
tic
  while(Energy < 0.9999)
    s_ = s_+1;
    sum_Sigma_s = 0;
    sum_Sigma_all = 0;
    for i = 1:s_
      sum_Sigma_s = sum_Sigma_s + S_DATA(i,i)^2;
    end
    for i = 1:vid.NumberOfFrames
      sum_Sigma_all = sum_Sigma_all + S_DATA(i,i)^2;
    end

    Energy = sum_Sigma_s/sum_Sigma_all;
  end

  fprintf('\n\tLow Rank s number achieved: %d', s_);
  fprintf('\n\twith Energy level = %d\n', Energy);
toc
clear('sum_Sigma_s');
clear('sum_Sigma_all');


fprintf('\nComputing and displaying Low Rank movie...\n');
tic
  %Computes the low rank approximation s = s_
  DATAMovie_s = U_DATA(:,1:s_)*S_DATA(1:s_,1:s_)*V_DATA(:,1:s_)';

  % View the low rank approximated video
  for k = 1:vidNFrames
      mov2(k).cdata = [cast(DATAMovie_s(:,k),'uint8')];
      mov_LowRank(k).cdata = reshape(mov2(k).cdata,vidHeight,vidWidth);
      mov_LowRank(k).colormap = gray(255);
  end

  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  movie(hf, mov_LowRank, 1, vid.FrameRate);
  nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Low Rank image at 1/4 in the movie...\n');
tic
  FrameInterest = round(vidNFrames*0.25);
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  imshow(mov_LowRank(FrameInterest).cdata,mov_LowRank(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc

%------------ Starting Separation DMD into background and foreground ------------%

DATA_SxN = U_DATA(:,1:s_)'*DATAMovie_s;

fprintf('\nCreating Original and Low Rank (s) X and Y matrices...\n');
tic
  for k = 1:vidNFrames-1
      Xs(:,k) = [DATA_SxN(:,k)];
      Ys(:,k) = [DATA_SxN(:,k+1)];
  end
toc

fprintf('\nObtaining DMD Matrix A...\n');
tic
  %Calculating SVD(X) to find pseudoinv(X) (same as pinv(X))
  %[U,S,V] = svd(Xs, 'econ');
  %S_inv = inv(S);

  %Low rank DMD Matrix using Ys * pseudoinv(X)
  A = Ys*pinv(Xs);%V*S_inv*U';
toc

%------------Uncomment to use alternate DMD Matrix and comment above ----------------%
%%fprintf('\nCreating Alternate DMD Matrix...\n');
%%tic
%%  %Creating X and Y matrices out of the original Data Matrix
%%  for k = 1:vidNFrames-1
%%      X(:,k) = [DATAMovie(:,k)];
%%      Y(:,k) = [DATAMovie(:,k+1)];
%%  end
%%
%%  [U,S,V] = svd(X,'econ');
%%  A = U(:,1:s_)'*Y*V(:,1:s_)*inv(S(1:s_,1:s_));
%%  Xs = U(:,1:s_)'*X;
%%toc



dt = 1.0/vid.FrameRate;
t = (0:dt:dt*(vidNFrames-1))';

%%dt = 1.0;
%%t = (0:1:N-1)'; %Same result as using time except scaled smaller.

fprintf('\nObtaining Eigen Vectors and values and Omega of DMD Matrix A...\n');
tic
  [e_Vect, e_Val] = eig(A); %finds eigenvectors and eigenvalues of A
  e_Val_D = diag(e_Val); %Take all the diagonals of eigenvalues
  omega = log(e_Val_D)/dt;
  absOmega = abs(omega);
  absOmegaSorted = sort(absOmega);
toc

fprintf('\nPlotting Omega sorted...\n');
tic
  figure(nfigure);
  plot(absOmegaSorted, '.', 'MarkerSize', 20);
  nfigure = nfigure + 1;
toc

%Lambda
Threshold = 0.01;

%Initializing Matrices for Results
e_Val_Back = zeros([size(e_Val)]);
e_Vect_Back = zeros([size(e_Vect)]);
e_Val_Fore = zeros([size(e_Val)]);
e_Vect_Fore = zeros([size(e_Vect)]);

D_Back = zeros(M,N);
D_Fore = zeros(M,N);

Ds_Back = zeros(s_,N);
Ds_Fore = zeros(s_,N);
Ds_all = zeros(s_,N);

fprintf('\nSeparating Background and Foreground Matrices with Threshold: %d...\n',Threshold);
  %Constants of eigensolutions using t = 0 and the first column of Xs
  b = pinv(e_Vect)*Xs(:,1);
  Percentage = 1/length(A)*100;
tic
  for i = 1:length(A);

    if(absOmega(i) < Threshold)
      %Background Separation
      e_Val_Back(i,i) = e_Val(i,i);
      e_Vect_Back(:,i) = e_Vect(:,i);
      Ds_Back = Ds_Back + b(i)*e_Vect_Back(:,i).*exp(omega(i)'*t)';  %exp(1i*imag(omega(i))*t);
    else
      %Foreground Separation
      e_Val_Fore(i,i) = e_Val(i,i);
      e_Vect_Fore(:,i) = e_Vect(:,i);
      Ds_Fore = Ds_Fore + b(i)*e_Vect_Fore(:,i).*exp(omega(i)'*t)';  %exp(1i*imag(omega(i))*t);
    end

    %Full Reconstrucion
    Ds_all = Ds_all + b(i)*e_Vect(:,i).*exp(omega(i)'*t)';

    Completed = rem(i,10);
    if(Completed == 0)
      fprintf('\n Completed: %d', round(i*Percentage));
      disp('%');
    end

    %Calculate and display estimated time to complete Loop
    if(i == 1)
      FirstLoop_t = toc;
      fprintf('\nEstimated Duration = %d minutes\n',FirstLoop_t*length(A)/60);
    end

  end
toc

fprintf('\nConverting Back to MxN matrices...\n');
tic
  D_Back_og = U_DATA(:,1:s_)*Ds_Back;
  D_Fore_og = U_DATA(:,1:s_)*Ds_Fore;
  D_all_og = U_DATA(:,1:s_)*Ds_all;

  %D_all = D_Back_og + D_Fore_og;
  D_Back = real(D_Back_og);
  D_Fore = real(D_Fore_og);
  D_all = real(D_all_og);
toc


%Comment this section to view foreground separation without using
%Background substraction with residual
##D_Fore = DATAMovie_s - abs(D_Back_og);
##D_Fore_Copy = D_Fore;
##
##fprintf('\nCreating Residual Matrix...\n');
##tic
##  R = zeros([size(D_Fore)]);
##
##  D_Fore_Copy(D_Fore_Copy<0) = 0;
##  R = D_Fore - D_Fore_Copy;
##  R(R<-40) = 3*R(R<-40);%-40 for test data
##
##  D_Fore = D_Fore - R;
##
##  D_Fore(D_Fore>180) = 255;
##toc
%----------------------------------------------------------------%



%--------- View Reconstruction Videos ---------%

fprintf('\nConverting and Displaying Background Video...\n');
tic
  for k = 1:vidNFrames
      mov2(k).cdata = [cast(D_Back(:,k),'uint8')];
      mov_Back(k).cdata = reshape(mov2(k).cdata,vidHeight,vidWidth);
      mov_Back(k).colormap = gray(255);
  end

  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  movie(hf, mov_Back, 1, vid.FrameRate);
  nfigure = nfigure + 1;
toc

fprintf('\nConverting and Displaying Foreground Video...\n');
tic
  for k = 1:vidNFrames
      mov2(k).cdata = [cast(D_Fore(:,k),'uint8')];
      mov_Fore(k).cdata = reshape(mov2(k).cdata,vidHeight,vidWidth);
      mov_Fore(k).colormap = gray(255);
  end

  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  movie(hf, mov_Fore, 1, vid.FrameRate);
  nfigure = nfigure + 1;
toc

fprintf('\nConverting and Displaying All Video...\n');
tic
  for k = 1:vidNFrames
      mov2(k).cdata = [cast(D_all(:,k),'uint8')];
      mov_all(k).cdata = reshape(mov2(k).cdata,vidHeight,vidWidth);
      mov_all(k).colormap = gray(255);
  end

  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  movie(hf, mov_all, 1, vid.FrameRate);
  nfigure = nfigure + 1;
toc


%--------- View Reconstruction Images ---------%

fprintf('\nDisplaying Background image at 1/4 in the movie...\n');
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  imshow(mov_Back(FrameInterest).cdata,mov_Back(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Foreground image at 1/4 in the movie...\n');
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  imshow(mov_Fore(FrameInterest).cdata,mov_Fore(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc

fprintf('\nDisplaying all image at 1/4 in the movie...\n');
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth vidHeight]);
  imshow(mov_all(FrameInterest).cdata,mov_all(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc

%Clear Variables to relieve memory
clear('D_Back_og');
clear('D_Fore_og');
clear('D_all_og');
clear('Ds_Back');
clear('Ds_Fore');
clear('Ds_all');
clear('Xs');
clear('Ys');
clear('DATAMovie');
clear('DATAMovie_s');
clear('DATA_SxN');
clear('U_DATA');
clear('FirstLoop_t');
clear('Completed');
clear('Percentage');
clear('S_DATA');
clear('V_DATA');
clear('grayImage');

