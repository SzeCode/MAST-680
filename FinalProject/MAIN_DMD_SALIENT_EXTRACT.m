%Main program for DMD Saliency Algorithm
clear all; close all; clc;

pkg load image;
Funcdir='./Functions/';
addpath(Funcdir);
Funcdir='./Functions/Sub1';
addpath(Funcdir);

pause(0.1);

tic

OriginalImage = imread("./Input/0_1_1650.jpg");

%Display original
figure();
imshow(OriginalImage);
title("RGB Image");

%Color Space Conversion and Data matrix Creation
Display = 0;
  [GrayscaleImage,N,M,X_color1,X_color2] = ColorDataf(OriginalImage,Display);
Display = 0;
  X_Lum = LumiDataf(OriginalImage,GrayscaleImage,M,N,Display);
  [M_xLum, N_xLum] = size(X_Lum);

%DMD Low and Sparse Separation
[X_DLow1, X_DSparse1] = ColorDMD2f(X_color1, GrayscaleImage, M, N);
[X_DLow2, X_DSparse2] = ColorDMD2f(X_color2, GrayscaleImage, M, N);
[X_DLow_Lum, X_DSparse_Lum] = LumiDMDf(X_Lum, M_xLum, N_xLum);

%Post-processing and final saliency map
Display = 0;
  SM_colorview = ColorFinal2f(X_DLow1, X_DLow2, X_DSparse1, X_DSparse2, M, N, Display);
  SM_LUMVIEW = LumiFinalf(X_DLow_Lum, X_DSparse_Lum, M, N);

Display = 0;
  SM_FIN = SMFinalMapf(SM_colorview, SM_LUMVIEW, M, N, Display);

%Segmentation
T = 100;
IMSeg = Segmentation(OriginalImage, SM_FIN, M, N, T);


%Display final saliency map and segmentation results
figure();
  SM_FIN = cast(SM_FIN, 'uint8');
  imshow(SM_FIN);
  title("SM FINAL");

figure();
  IMSeg = cast(IMSeg, 'uint8');
  imshow(IMSeg);
  title("Image Segmented");

toc
