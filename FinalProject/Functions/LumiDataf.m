function X_Lum = LumiDataf(OriginalImage,GrayscaleImage,M,N,display)
%Luminosity Image Conversion and DMD Matrix Creation

  %Image stored as R, G, B
  R = OriginalImage(:,:,1);
  G = OriginalImage(:,:,2);
  B = OriginalImage(:,:,3);

  %----------------------- Luminosity SPACE CONVERSION -----------------------%
  %YUL
  R_ = cast(R,'double');
  G_ = cast(G,'double');
  B_ = cast(B,'double');

  Y = 0.299*R_+0.587*G_+0.114*B_;

  if(display == 1)
    Y = cast(Y,'uint8');
    figure();
    imshow(Y);
    title("Y Image");
  end


  %LAB
  Lab = rgb2lab(OriginalImage);
  L = Lab(:,:,1);
  Color_max = cast(max(max(L)),'double');
  Color_min = cast(min(min(L)),'double');
  IntCast =  cast(L,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  L = cast(IntCast,'uint8');

  if(display == 1)
    figure();
    imshow(L);
    title("L Image");
  end

  %YCBCR
  ycbcr = rgb2ycbcr(OriginalImage);
  y = ycbcr(:,:,1);
  Color_max = cast(max(max(y)),'double');
  Color_min = cast(min(min(y)),'double');
  IntCast =  cast(y,'double');
  IntCast = (IntCast-Color_min)/(Color_max-Color_min)*(255.0);
  y_ = cast(IntCast,'uint8');

  if(display == 1)
    figure();
    imshow(y);
    title("y Image");
  end

  %----------------------- LUMINOSITY IMAGES -----------------------%
  L = cast(L,'double');
  Y = cast(Y,'double');
  y = cast(y,'double');

  [U_L,S_L,V_L] = svd(L,'econ');
  [U_Y1,S_Y1,V_Y1] = svd(Y,'econ');
  [U_Y2,S_Y2,V_Y2] = svd(y,'econ');

  X_L = [];
  X_Y1 = [];
  X_Y2 = [];
  X_Lum = [];

  EnergyTarget = 0.95;
  StartSingular = 3;
  EndSingular = 5;
  while(EndSingular-StartSingular<20)
    EndSingular = SRankApprox(EnergyTarget, S_L);
    EnergyTarget = EnergyTarget+0.0001;
  end
  L_sum = zeros(M,N);

  if(display == 1)
    figure();
  end
  for i = StartSingular+2:EndSingular
    s_ = StartSingular:i;
    %s_ = i:i;

    L_SVD = U_L(:,s_)*S_L(s_,s_)*V_L(:,s_)';
    X_L = [X_L, reshape(L_SVD,M*N,1)];

    L_sum = L_sum + L_SVD;

    if(display == 1)
      L_SVD = cast(L_SVD, 'uint8');
      imshow(L_SVD);
      title("L Int Singular Reconst Image");
      pause(0.1)


      if(i == 2 || i == 5 || i == 7 || i == 10)
        %nfigure = nfigure + 1;
        figure();
      end
    end

  end
  %fprintf('iL = %d\n',i');
  X_Lum = [X_Lum, X_L];

  if(display == 1)
    figure();
  end
  for i = StartSingular+2:EndSingular
    s_ = StartSingular:i;
    %s_ = i:i;

    Y1_SVD = U_Y1(:,s_)*S_Y1(s_,s_)*V_Y1(:,s_)';
    X_Y1 = [X_Y1, reshape(Y1_SVD,M*N,1)];

    if(display == 1)
      Y1_SVD = cast(Y1_SVD, 'uint8');
      imshow(Y1_SVD);
      title("Y1 Int Singular Reconst Image");
      pause(0.1)
    end
  end
  %fprintf('Y1L = %d\n',i);
  X_Lum = [X_Lum, X_Y1];

  if(display == 1)
    figure();
  end
  for i = StartSingular+2:EndSingular
    s_ = StartSingular:i;
    %s_ = i:i;

    Y2_SVD = U_Y2(:,s_)*S_Y2(s_,s_)*V_Y2(:,s_)';
    X_Y2 = [X_Y2, reshape(Y2_SVD,M*N,1)];

    if(display == 1)
      Y2_SVD = cast(Y2_SVD, 'uint8');
      imshow(Y2_SVD);
      title("Y2 Int Singular Reconst Image");
      pause(0.1)
    end

  end
  %fprintf('Y2L = %d\n',i);
  X_Lum = [X_Lum, X_Y2];


  %View aggregate SVD reconstruction of single singular values
  if(display == 1)
    L_sum = cast(L_sum,'uint8');
    figure();
    imshow(L_sum);
    title("L Summed Int Singular Reconst Image");
  end

end
