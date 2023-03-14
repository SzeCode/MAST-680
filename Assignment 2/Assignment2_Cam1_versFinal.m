%Mass position extracting on Cam 1
%Will Sze
%40096561

clear all; close all; clc
pause(0.1);

nfigure = 1;
FrameInterest = 1;

fprintf('\nImporting Video...\n');
tic
  %Load data into Octave
  load cam1_1.mat;

  %Collects video dimensions
  [vidHeight1, vidWidth1, nRGB1, vidNFrames1] = size(vidFrames1_1);
  FrameRate = 60;

  %Movie Struct to collect all original frames to run movie using movie function
  mov1 = struct('cdata', zeros(vidHeight1,vidWidth1,3,'uint8'), 'colormap',[]);
  mov1_Gray = struct('cdata', zeros(vidHeight1,vidWidth1,1,'uint8'), 'colormap',gray(255));
  mov1_BW = struct('cdata', zeros(vidHeight1,vidWidth1,1,'uint8'), 'colormap',gray(255));
  mov1_int = struct('cdata', zeros(vidHeight1*vidWidth1,1,'uint8'), 'colormap',gray(255));

  for k = 1:vidNFrames1
    mov1(k).cdata = vidFrames1_1(:,:,:,k);
  end
toc

%View original video1
fprintf('\nDisplaying Original Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1, 1, FrameRate);
   nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Video1 image at frame %d in the movie...\n',FrameInterest);
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1(FrameInterest).cdata,mov1(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc


%Converting to grayscale1
tic
  for k = 1:vidNFrames1
    mov1_Gray(k).cdata = rgb2gray(mov1(k).cdata);
    mov1_Gray(k).colormap = gray(255);
  end
toc


%View GrayScale video1
fprintf('\nDisplaying Gray Scale Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1_Gray, 1, FrameRate);
   nfigure = nfigure + 1;
toc

fprintf('\nDisplaying Video1 image at frame %d in the movie...\n',FrameInterest);
tic
  hf = figure(nfigure);
  set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
  imshow(mov1_Gray(FrameInterest).cdata,mov1_Gray(FrameInterest).colormap);
  nfigure = nfigure + 1;
toc



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


IntensityThreshold = 50;
%Mode Background Method
fprintf('\nCreating background on Data1...\n');
tic
  Background = zeros(M,N);
  for i = 1:M
    Background(i,:) = mode(Data1(i,:));
  end
  Data_Test2 = Data1 - Background;

  Data_Test2(Data_Test2>IntensityThreshold) = 255;
  Data_Test2(Data_Test2<=IntensityThreshold) = 0;
toc

%Horizontal and Vertical Intensities
fprintf('\nLocating the mass...\n');
cc_F_Interest = 1;
x_pos = [];
y_pos = [];
x_Intensity = zeros(vidNFrames1,vidWidth1);
y_Intensity = zeros(vidNFrames1,vidHeight1);
for i = 1:vidNFrames1

  tic
  Frame_Analysis = reshape(Data_Test2(:,i),vidHeight1, vidWidth1);
  Gray_Frame = reshape(Data1(:,i),vidHeight1, vidWidth1);

  %Set Boundary on first frame
  if(i == 1)
    Frame_Analysis(:,1:300) = 0;
    Frame_Analysis(:,390:vidWidth1) = 0;
    Frame_Analysis(1:210,:) = 0;
    Frame_Analysis(320:vidHeight1,:) = 0;
  end
  %Use previous position to locate next boundary
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

  %Calculating Centroid
  for j = 1:vidWidth1
      x_Intensity(i,j) = mean(Frame_Analysis(:,j));
  end
  for j = 1:vidHeight1
    y_Intensity(i,vidHeight1+1-j) = mean(Frame_Analysis(j,:));
  end

  %X-centroid
  SumWeightPos = 0;
  SumWeight = 0;
  for j = 1:vidWidth1
    SumWeightPos = SumWeightPos + j*x_Intensity(i,j);
    SumWeight = SumWeight + x_Intensity(i,j);
  end
  x_pos = [x_pos,SumWeightPos/SumWeight];

  %Y-centroid
  SumWeightPos = 0;
  SumWeight = 0;
  for j = 1:vidHeight1
    SumWeightPos = SumWeightPos + j*y_Intensity(i,j);
    SumWeight = SumWeight + y_Intensity(i,j);
  end
  y_pos = [y_pos,SumWeightPos/SumWeight];

  %Visualize Point in video
  Frame_Analysis(:,round(x_pos(i))) = 180;
  Gray_Frame(:,round(x_pos(i))) = 0;
  Frame_Analysis(vidHeight1 - round(y_pos(i)),:) = 180;
  Gray_Frame(vidHeight1 - round(y_pos(i)),:) = 0;

  Data_Test2(:,i) = reshape(Frame_Analysis,vidHeight1*vidWidth1,1);
  Data1(:,i) = reshape(Gray_Frame,vidHeight1*vidWidth1,1);
  toc
end


fprintf('\nPlotting intensity curves at frame %d in the movie...\n',cc_F_Interest);
tic
  figure(nfigure);
  plot(x_Intensity(cc_F_Interest,:),'-o');
  nfigure = nfigure + 1;

  y = linspace(1,vidHeight1,vidHeight1);
  figure(nfigure);
  plot(y_Intensity(cc_F_Interest,:),y,'-o');
  nfigure = nfigure + 1;
toc


%View Results
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

fprintf('\nDisplaying Video1 image at frame %d in the movie...\n',FrameInterest);
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

fprintf('\nDisplaying Gray with centroid of Video1...\n');
tic
   hf = figure(nfigure);
   set(hf, 'position', [150 150 vidWidth1 vidHeight1]);
   movie(hf, mov1_Gray, 1, FrameRate);
   nfigure = nfigure + 1;
toc


