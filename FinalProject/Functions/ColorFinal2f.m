function SM_colorview = ColorFinal2f(X_DLow1, X_DLow2, X_DSparse1, X_DSparse2, M, N, Display)

  %Normalisation
  X_DLow1 = real(X_DLow1);
  a = X_DLow1 - min(min(X_DLow1));
  b = max(max(X_DLow1)) - min(min(X_DLow1));
  X_DLow1_norm = a/b;

  X_DLow2 = real(X_DLow2);
  a = X_DLow2 - min(min(X_DLow2));
  b = max(max(X_DLow2)) - min(min(X_DLow2));
  X_DLow2_norm = a/b;

  X_DSparse1 = real(X_DSparse1);
  a = X_DSparse1 - min(min(X_DSparse1));
  b = max(max(X_DSparse1)) - min(min(X_DSparse1));
  X_DSparse1_norm = a/b;

  X_DSparse2 = real(X_DSparse2);
  a = X_DSparse2 - min(min(X_DSparse2));
  b = max(max(X_DSparse2)) - min(min(X_DSparse2));
  X_DSparse2_norm = a/b;

  %Sum of low and sparse components
  X_LOWOUT = X_DLow1_norm + X_DLow2_norm;
  X_SPARSEOUT = X_DSparse1_norm + X_DSparse2_norm;


  FrameNum = 1;

  x_LOWVIEW = reshape(X_LOWOUT(:,FrameNum),M,N);
  a = x_LOWVIEW - min(min(x_LOWVIEW));
  b = max(max(x_LOWVIEW)) - min(min(x_LOWVIEW));
  x_LOWVIEW = 255*a/b;


  x_SPARSEVIEW = reshape(X_SPARSEOUT(:,FrameNum),M,N);
  a = x_SPARSEVIEW - min(min(x_SPARSEVIEW));
  b = max(max(x_SPARSEVIEW)) - min(min(x_SPARSEVIEW));
  x_SPARSEVIEW = 255*a/b;

  w = 1.5;
  SM_color = (X_LOWOUT-w*X_SPARSEOUT);
  SM_color = reshape(SM_color,M,N);
  SM_color(SM_color<0) = 0;
  a = SM_color - min(min(SM_color));
  b = max(max(SM_color)) - min(min(SM_color));
  SM_colorview = 255*(a/b).^2;

  x_LOWVIEW = cast(x_LOWVIEW,'uint8');
  x_SPARSEVIEW = cast(x_SPARSEVIEW,'uint8');

  if(Display == 1)
    figure();
    imshow(x_LOWVIEW);
    title("D_A + D_B LOWVIEW Image");

    figure();
    imshow(x_SPARSEVIEW);
    title("D_A + D_B SPARSEVIEW Image");
  end

end
