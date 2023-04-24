function SM_FIN = SMFinalMapf(SM_colorview, SM_LUMVIEW, M, N, Display)

  xpos = 1:1:M;
  ypos = 1:1:N;

  %Gaussian Distribution
  xposn = xpos-M/2;
  yposn = ypos-N/2;

  Wg = zeros(M,N);
  sigma = 175;
  for i  = 1:M
    for j = 1:N
      Wg(i,j) = exp(-(xposn(i)^2+yposn(j)^2)/(2*sigma^2))/(2*pi*sigma^2);
    end
  end
  Const = 1/max(max(Wg));
  Wg = Const*Wg;


  SM_colorview = cast(SM_colorview, 'double');
  SM_colorview = Wg.*SM_colorview;
  SM_LUMVIEW = cast(SM_LUMVIEW, 'double');
  SM_LUMVIEW = Wg.*SM_LUMVIEW;

  if(Display == 1)
    figure();
      SM_colorview = cast(SM_colorview, 'uint8');
      imshow(SM_colorview);
      title("SM colorview Weighted");

    figure();
      SM_LUMVIEW = cast(SM_LUMVIEW, 'uint8');
      imshow(SM_LUMVIEW);
      title("SM LUMVIEW Weighted");
  end


  nhood = [0 1 0; 1 1 1; 0 1 0];
  SE = strel('square',3);

  %Morphological Erosion and Dilation
  N_Erode = 1;

  for i = 1:N_Erode
    SM_colorview = imdilate(SM_colorview,nhood);
    SM_LUMVIEW = imdilate(SM_LUMVIEW,nhood);
  end
  for i = 1:N_Erode*2
    SM_colorview = imerode(SM_colorview,nhood);
    SM_LUMVIEW = imerode(SM_LUMVIEW,nhood);
  end
  for i = 1:N_Erode
    SM_colorview = imdilate(SM_colorview,nhood);
    SM_LUMVIEW = imdilate(SM_LUMVIEW,nhood);
  end

  if(Display == 1)
    figure();
    imshow(SM_colorview);
    title("SM colorview Eroded");

    figure();
    imshow(SM_LUMVIEW);
    title("SM LUMVIEW Eroded");
  end

  %Sigmoid Intensity adjustment
  SM_colorview = cast(SM_colorview, 'double');
  SMa = SM_colorview - min(min(SM_colorview));
  SMb = max(max(SM_colorview)) - min(min(SM_colorview));
  SM_colorview = SMa/SMb;
  b = 30;
  C = 1+exp(-b*(SM_colorview-.3));
  SM_colorview = C.^(-1);

  SMa = SM_colorview - min(min(SM_colorview));
  SMb = max(max(SM_colorview)) - min(min(SM_colorview));
  SM_colorview = 255*SMa/SMb;

  SM_LUMVIEW = cast(SM_LUMVIEW, 'double');
  SMa = SM_LUMVIEW - min(min(SM_LUMVIEW));
  SMb = max(max(SM_LUMVIEW)) - min(min(SM_LUMVIEW));
  SM_LUMVIEW = SMa/SMb;
  b = 30;
  C = 1+exp(-b*(SM_LUMVIEW-.5));
  SM_LUMVIEW = C.^(-1);

  SMa = SM_LUMVIEW - min(min(SM_LUMVIEW));
  SMb = max(max(SM_LUMVIEW)) - min(min(SM_LUMVIEW));
  SM_LUMVIEW = 255*SMa/SMb;

  if(Display == 1)
    figure();
    SM_colorview = cast(SM_colorview, 'uint8');
    imshow(SM_colorview);
    title("SM Color Postprocessed");

    figure();
    SM_LUMVIEW = cast(SM_LUMVIEW, 'uint8');
    imshow(SM_LUMVIEW);
    title("SM Lum Postprocessed");
  end

  SM_LUMVIEW = cast(SM_LUMVIEW, 'double');
  SM_colorview = cast(SM_colorview, 'double');

  SM_FIN = max(SM_colorview,SM_LUMVIEW);

end
