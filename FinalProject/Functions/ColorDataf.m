function [GrayscaleImage,N,M,X_color1,X_color2] = ColorDataf(OriginalImage,display)
%Color Space Conversion and DMD Matrix Creation

  GrayscaleImage = rgb2gray(OriginalImage);
  [M,N] = size(GrayscaleImage);

  if(display == 1)
    figure();
    imshow(GrayscaleImage);
    title("Gray Image");
  end

  %Image stored as R, G, B
  R = OriginalImage(:,:,1);
  G = OriginalImage(:,:,2);
  B = OriginalImage(:,:,3);

  %----------------------- COLOR SPACE CONVERSION -----------------------%
  %YUL
  R_ = cast(R,'double');
  G_ = cast(G,'double');
  B_ = cast(B,'double');

  Y = 0.299*R_+0.587*G_+0.114*B_;

  U = abs((B_-Y)*0.492);
  Color_max = cast(max(max(U)),'double');
  Color_min = cast(min(min(U)),'double');
  IntCast =  cast(U,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  U = cast(IntCast,'uint8');

  V = abs((R_-Y)*0.877);
  Color_max = cast(max(max(V)),'double');
  Color_min = cast(min(min(V)),'double');
  IntCast =  cast(V,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  V = cast(IntCast,'uint8');

  if(display == 1)
    Y = cast(Y,'uint8');
    figure();
    imshow(Y);
    title("Y Image");

    U = cast(U,'uint8');
    figure();
    imshow(U);
    title("U Image");

    V = cast(V,'uint8');
    figure();
    imshow(V);
    title("V Image");
  end


  %LAB
  Lab = rgb2lab(OriginalImage);
%%  L = Lab(:,:,1);
%%  Color_max = cast(max(max(L)),'double');
%%  Color_min = cast(min(min(L)),'double');
%%  IntCast =  cast(L,'double');
%%  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
%%  L_ = cast(IntCast,'uint8');

  a = abs(Lab(:,:,2));
  Color_max = cast(max(max(a)),'double');
  Color_min = cast(min(min(a)),'double');
  IntCast =  cast(a,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  a = cast(IntCast,'uint8');

  b = abs(Lab(:,:,3));
  Color_max = cast(max(max(b)),'double');
  Color_min = cast(min(min(b)),'double');
  IntCast =  cast(b,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  b = cast(IntCast,'uint8');

  if(display == 1)
%%    figure();
%%    imshow(L_);
%%    title("L Image");

    figure();
    imshow(a);
    title("a Image");

    figure();
    imshow(b);
    title("b Image");
  end

  %YCBCR
  ycbcr = rgb2ycbcr(OriginalImage);
%%  y = ycbcr(:,:,1);
%%  Color_max = cast(max(max(y)),'double');
%%  Color_min = cast(min(min(y)),'double');
%%  IntCast =  cast(y,'double');
%%  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
%%  y_ = cast(IntCast,'uint8');


  cb = abs(ycbcr(:,:,2));
  Color_max = cast(max(max(cb)),'double');
  Color_min = cast(min(min(cb)),'double');
  IntCast =  cast(cb,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  cb_ = cast(IntCast,'uint8');

  cr = abs(ycbcr(:,:,3));
  Color_max = cast(max(max(cr)),'double');
  Color_min = cast(min(min(cr)),'double');
  IntCast =  cast(cr,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  cr_ = cast(IntCast,'uint8');

  if(display == 1)
%%    figure();
%%    imshow(y);
%%    title("y Image");

    figure();
    imshow(cb);
    title("cb Image");

    figure();
    imshow(cr);
    title("cr Image");
  end


  %----------------------- DATA MATRIX CREATION -----------------------%

  X_color1 = zeros(M*N,40);

  X_color1(:,0+1) = reshape(b,M*N,1);
  X_color1(:,0+2) = reshape(V,M*N,1);
  X_color1(:,0+3) = reshape(cr,M*N,1);
  X_color1(:,0+4) = X_color1(:,1)+X_color1(:,2)+X_color1(:,3);

  X_color1(:,4+1) = X_color1(:,0+2);
  X_color1(:,4+2) = X_color1(:,0+3);
  X_color1(:,4+3) = X_color1(:,0+1);
  X_color1(:,4+4) = X_color1(:,0+4);

  X_color1(:,8+1) = X_color1(:,0+4);
  X_color1(:,8+2) = X_color1(:,0+3);
  X_color1(:,8+3) = X_color1(:,0+2);
  X_color1(:,8+4) = X_color1(:,0+1);

  X_color1(:,12+1) = X_color1(:,0+2);
  X_color1(:,12+2) = X_color1(:,0+3);
  X_color1(:,12+3) = X_color1(:,0+4);
  X_color1(:,12+4) = X_color1(:,0+1);

  X_color1(:,16+1) = X_color1(:,0+2);
  X_color1(:,16+2) = X_color1(:,0+1);
  X_color1(:,16+3) = X_color1(:,0+4);
  X_color1(:,16+4) = X_color1(:,0+3);

  X_color1(:,21:40) = X_color1(:,1:20);

  X_color2 = zeros(M*N,40);

  X_color2(:,1) = reshape(a,M*N,1);
  X_color2(:,2) = reshape(U,M*N,1);
  X_color2(:,3) = reshape(cb,M*N,1);
  X_color2(:,4) = X_color2(:,1)+X_color2(:,2)+X_color2(:,3);

  X_color2(:,4+1) = X_color2(:,0+2);
  X_color2(:,4+2) = X_color2(:,0+3);
  X_color2(:,4+3) = X_color2(:,0+1);
  X_color2(:,4+4) = X_color2(:,0+4);

  X_color2(:,8+1) = X_color2(:,0+4);
  X_color2(:,8+2) = X_color2(:,0+3);
  X_color2(:,8+3) = X_color2(:,0+2);
  X_color2(:,8+4) = X_color2(:,0+1);

  X_color2(:,12+1) = X_color2(:,0+2);
  X_color2(:,12+2) = X_color2(:,0+3);
  X_color2(:,12+3) = X_color2(:,0+4);
  X_color2(:,12+4) = X_color2(:,0+1);

  X_color2(:,16+1) = X_color2(:,0+2);
  X_color2(:,16+2) = X_color2(:,0+1);
  X_color2(:,16+3) = X_color2(:,0+4);
  X_color2(:,16+4) = X_color2(:,0+3);

  X_color2(:,21:40) = X_color2(:,1:20);

end
