clear all; close all; clc

%pkg load video


pause(0.1);



nfigure = 1;


fprintf('\nImporting Video...\n');
tic
  %Load data into Octave
  load cam2_1.mat;
  load cam2_2.mat;

  %Collects video dimensions
  [vidHeight1, vidWidth1, nRGB1, vidNFrames1] = size(vidFrames2_1);
  [vidHeight2, vidWidth2, nRGB2, vidNFrames2] = size(vidFrames2_2);
  FrameRate = 60;

  %Movie Struct to collect all original frames to run movie using movie function
  mov1 = struct('cdata', zeros(vidHeight1,vidWidth1,3,'uint8'), 'colormap',[]);
  mov1_Gray = struct('cdata', zeros(vidHeight1,vidWidth1,1,'uint8'), 'colormap',gray(255));
  mov2 = struct('cdata', zeros(vidHeight2,vidWidth2,3,'uint8'), 'colormap',[]);
  mov2_Gray = struct('cdata', zeros(vidHeight2,vidWidth2,1,'uint8'), 'colormap',gray(255));
  mov1_BW = struct('cdata', zeros(vidHeight1,vidWidth1,1,'uint8'), 'colormap',gray(255));
  mov2_BW = struct('cdata', zeros(vidHeight2,vidWidth2,1,'uint8'), 'colormap',gray(255));
  mov1_int = struct('cdata', zeros(vidHeight1*vidWidth1,1,'uint8'), 'colormap',gray(255));
  mov1_LowRank = struct('cdata', zeros(vidHeight1,vidWidth1,1,'uint8'), 'colormap',gray(255));


  for k = 1:vidNFrames1
    mov1(k).cdata = vidFrames2_1(:,:,:,k);
  end
  for k = 1:vidNFrames2
    mov2(k).cdata = vidFrames2_2(:,:,:,k);
  end
toc

FrameInterest = 1 ;% round(vidNFrames1*0.5);

%View original video1
fprintf('\nDisplaying Original Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1, 1, FrameRate);
   nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Video1 image at 1/2 in the movie...\n');
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1(FrameInterest).cdata,mov1(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc



%View original video
##fprintf('\nDisplaying Original Video2...\n');
##tic
##   hf = figure(nfigure);
##   set(hf, 'position', [150 150 vidWidth2 vidHeight2]);
##   movie(hf, mov2, 1, FrameRate);
##   nfigure = nfigure + 1;
##toc

%Converting to grayscale1
tic
  for k = 1:vidNFrames1

    mov1_Gray(k).cdata = rgb2gray(mov1(k).cdata);
    mov1_Gray(k).colormap = gray(255);

  end
toc
%Converting to grayscale2
##tic
##  for k = 1:vidNFrames2
##
##    mov2_Gray(k).cdata = rgb2gray(mov2(k).cdata);
##    mov2_Gray(k).colormap = gray(255);
##
##  end
##toc


%View GrayScale video1
fprintf('\nDisplaying Gray Scale Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1_Gray, 1, FrameRate);
   nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Video1 image at 1/2 in the movie...\n');
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1_Gray(FrameInterest).cdata,mov1_Gray(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc


%View GrayScale video2
##fprintf('\nDisplaying Gray Scale Video2...\n');
##tic
##   hf = figure(nfigure);
##   set(hf, 'position', [150 150 vidWidth2 vidHeight2]);
##   movie(hf, mov2_Gray, 1, FrameRate);
##   nfigure = nfigure + 1;
##toc

Data1 = zeros(vidWidth1*vidHeight1, vidNFrames1);

%---------- Tracking algorithm -------------%
fprintf('\nCreating Data1...\n');
tic
  for k = 1:vidNFrames1
    mov1_int(k).cdata = reshape(mov1_Gray(k).cdata,[],1);
    mov1_int(k).colormap = gray(255);
    Data1(:,k) = [mov1_int(k).cdata];
  end
toc
[M N] = size(Data1);


Data_Test = Data1;% - 50;
Data_Test(Data_Test<0) = 0;

%Substraction Method -- Works on some frames not all frames 50
##for i = 226:-1:1
##  if(i == 1)
##    break;
##  end
##  Data_Test(:,i) = Data_Test(:,i) - Data_Test(:,i-1);
##
##end
##Data_Test(Data_Test<0) = 0;
IntensityThreshold = 20;
##Data_Test(Data_Test>IntensityThreshold) = 255;
##Data_Test(Data_Test<=IntensityThreshold) = 0;

%Mode Background Method
fprintf('\nCreating Mode background on Data1...\n');
tic
  Mode_Background = zeros(M,N);

  for i = 1:M
    n=0;
    for j = 1:(vidNFrames1/4-2)
      Mode_Background(i,(1+4*n):(4+4*n)) = median(Data1(i,(1+4*n):(4+4*n)));
      n=n+1;
    end
  end
  Data_Test2 = abs(Data1 - Mode_Background);

  Data_Test2(Data_Test2>IntensityThreshold) = 255;
  Data_Test2(Data_Test2<=IntensityThreshold) = 0;
toc

tic
  for k = 1:vidNFrames1
    mov1_int(k).cdata = [cast(Mode_Background(:,k),'uint8')];
    mov1_BW(k).cdata = reshape(mov1_int(k).cdata,vidHeight1, vidWidth1);
    mov1_BW(k).colormap = gray(255);

  end
toc

fprintf('\nDisplaying Black and White Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1_BW, 1, FrameRate);
   nfigure = nfigure + 1;
toc




%Horizontal and Vertical Intensities
fprintf('\nCreating Horizontal and Vertical Intensities...\n');
cc_F_Interest = 1;
x_pos = [];
y_pos = [];
x_Intensity = zeros(vidNFrames1,vidWidth1);
y_Intensity = zeros(vidNFrames1,vidHeight1);
for i = 1:vidNFrames1

  tic
  Frame_Analysis = reshape(Data_Test2(:,i),vidHeight1, vidWidth1);
  Gray_Frame = reshape(Data1(:,i),vidHeight1, vidWidth1);

  if(i ~= 1)
    Y_UpperBound = round(vidHeight1-y_pos(i-1))+60;
    Y_LowerBound = round(vidHeight1-y_pos(i-1))-60;
    if(Y_UpperBound > vidHeight1)
      Y_UpperBound = vidHeight1;
    end
    if(Y_LowerBound < 1)
      Y_LowerBound = 1;
    end
    Frame_Analysis(1:Y_LowerBound,:) = 0;
    Frame_Analysis(Y_UpperBound:vidHeight1,:) = 0;

    X_UpperBound = round(x_pos(i-1))+50;
    X_LowerBound = round(x_pos(i-1))-50;
    if(X_UpperBound > vidWidth1)
      X_UpperBound = vidWidth1;
    end
    if(X_LowerBound < 1)
      X_LowerBound = 1;
    end
    Frame_Analysis(:,1:X_LowerBound) = 0;
    Frame_Analysis(:,X_UpperBound:vidWidth1) = 0;
  end
  if(i == 1)
    Frame_Analysis(:,1:260) = 0;
    Frame_Analysis(:,340:vidWidth1) = 0;
    Frame_Analysis(1:260,:) = 0;
    Frame_Analysis(370:vidHeight1,:) = 0;
  end

  for j = 1:vidWidth1
      x_Intensity(i,j) = mean(Frame_Analysis(:,j));
  end
  %x_Intensity(x_Intensity < 20) = 0;

  for j = 1:vidHeight1
    y_Intensity(i,vidHeight1+1-j) = mean(Frame_Analysis(j,:));
  end
  %y_Intensity(y_Intensity < 15) = 0;



  SumWeightPos = 0;
  SumWeight = 0;
  for j = 1:vidWidth1
    SumWeightPos = SumWeightPos + j*x_Intensity(i,j);
    SumWeight = SumWeight + x_Intensity(i,j);
  end

  x_pos = [x_pos,SumWeightPos/SumWeight];

  SumWeightPos = 0;
  SumWeight = 0;
  for j = 1:vidHeight1
    SumWeightPos = SumWeightPos + j*y_Intensity(i,j);
    SumWeight = SumWeight + y_Intensity(i,j);
  end

  y_pos = [y_pos,SumWeightPos/SumWeight];
  Frame_Analysis(:,round(x_pos(i))) = 180;
  Gray_Frame(:,round(x_pos(i))) = 0;
  Frame_Analysis(vidHeight1 - round(y_pos(i)),:) = 180;
  Gray_Frame(vidHeight1 - round(y_pos(i)),:) = 0;

  Data_Test2(:,i) = reshape(Frame_Analysis,vidHeight1*vidWidth1,1);
  Data1(:,i) = reshape(Gray_Frame,vidHeight1*vidWidth1,1);
  toc
end



figure(nfigure);
plot(x_Intensity(cc_F_Interest,:),'-o');
nfigure = nfigure + 1;

y = linspace(1,vidHeight1,vidHeight1);
figure(nfigure);
plot(y_Intensity(cc_F_Interest,:),y,'-o');
nfigure = nfigure + 1;







%Data_Test(Data_Test>50) = 255;
%Data_Test(Data_Test<=50) = 0;

%Threshhold Method

% Data1(Data1>IntensityThreshold) = 255;
% Data1(Data1<=IntensityThreshold) = 0;


tic
  for k = 1:vidNFrames1
    mov1_int(k).cdata = [cast(Data_Test2(:,k),'uint8')];
    mov1_BW(k).cdata = reshape(mov1_int(k).cdata,vidHeight1, vidWidth1);
    mov1_BW(k).colormap = gray(255);

  end
toc

fprintf('\nDisplaying Black and White Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1_BW, 1, FrameRate);
   nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Video1 image at 1/2 in the movie...\n');
tic
  FrameInterest = cc_F_Interest;
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1_BW(FrameInterest).cdata,mov1_BW(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc


tic
  for k = 1:vidNFrames1
    mov1_int(k).cdata = [cast(Data1(:,k),'uint8')];
    mov1_Gray(k).cdata = reshape(mov1_int(k).cdata,vidHeight1, vidWidth1);
    mov1_Gray(k).colormap = gray(255);

  end
toc
fprintf('\nDisplaying Reference Gray Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1_Gray, 1, FrameRate);
   nfigure = nfigure + 1;
toc





fprintf('\nDisplaying Video1 image at frame 50 in the movie...\n');
tic
  FrameInterest = cc_F_Interest;
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1_BW(50).cdata,mov1_BW(50).colormap);
  nfigure = nfigure + 1;
toc
fprintf('\nDisplaying Video1 image at frame 50 in the movie...\n');
tic
  FrameInterest = cc_F_Interest;
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1_Gray(50).cdata,mov1_Gray(50).colormap);
  nfigure = nfigure + 1;
toc

%--------------- SINDy Method Start ----------------%

X(1,:) = x_pos;
X(2,:) = y_pos;

dt = 1/FrameRate;

for i = 1:2
  for j = 1:vidNFrames1-1
    Y(i,j) = (X(i,j+1)-X(i,j))/dt;
  end
  Y(i,vidNFrames1) = Y(i,vidNFrames1-1);
end


##%Computing SVD for Data1
##fprintf('\nCalculating SVD...\n');
##tic
##  [U_DATA,S_DATA,V_DATA] = svd(Data1,'econ');
##toc
##
##%initializes variables for low rank calculations with energy
##s_ = 1;
##EnergyTarget = 0.9999;
##
##sum_Sigma_s = 0;
##sum_Sigma_all = 0;
##Energy = 0;
##
##%Searches for a rank value that will make the Energy bigger than or equal to 99.99% to have good approximation.
####fprintf('\nSearching for low rank s number achieving 99.99 percent Energy...\n');
####tic
####  while(Energy < 0.50)
####    s_ = s_+1;
####    sum_Sigma_s = 0;
####    sum_Sigma_all = 0;
####    for i = 1:s_
####      sum_Sigma_s = sum_Sigma_s + S_DATA(i,i)^2;
####    end
####    for i = 1:vidNFrames1
####      sum_Sigma_all = sum_Sigma_all + S_DATA(i,i)^2;
####    end
####
####    Energy = sum_Sigma_s/sum_Sigma_all;
####  end
####
####  fprintf('\n\tLow Rank s number achieved: %d', s_);
####  fprintf('\n\twith Energy level = %d\n', Energy);
####toc
####clear('sum_Sigma_s');
####clear('sum_Sigma_all');
##
##fprintf('\nComputing and displaying Low Rank movie...\n');
##tic
##  %Computes the low rank approximation s = s_
##  Data1_s = U_DATA(:,1:s_)*S_DATA(1:s_,1:s_)*V_DATA(:,1:s_)';
##
##  % View the low rank approximated video
##  for k = 1:vidNFrames1
##      mov1_int(k).cdata = [cast(Data1_s(:,k),'uint8')];
##      mov1_LowRank(k).cdata = reshape(mov1_int(k).cdata,vidHeight1,vidWidth1);
##      mov1_LowRank(k).colormap = gray(255);
##  end
##
##  hf = figure(nfigure);
##  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
##  movie(hf, mov1_LowRank, 1, FrameRate);
##  nfigure = nfigure + 1;
##toc

