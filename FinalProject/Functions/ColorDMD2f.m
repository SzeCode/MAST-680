function [X_DLow, X_DSparse] = ColorDMD2f(X_color, GrayscaleImage, M, N)

  [U_xcolor, S_xcolor, V_xcolor] = svd(X_color,'econ');
  [M_xcolor, N_xcolor] = size(X_color);

  %Define s as high as possible
  s_ = length(S_xcolor)-1;

  %DMD Low rank and sparse separation Version 2
  [Xs_DLow, Xs_DSparse, Us_xcolor] = LowSparseDMD3(X_color, s_);

  X_DLow = zeros(M_xcolor,N_xcolor);
  X_DSparse = zeros(M_xcolor,N_xcolor);

  X_DLow = Us_xcolor(:,1:s_)*Xs_DLow;%%%%%%%%%%%%%%%%%%%%%%%%
  X_DLow = sum(X_DLow')';
  %X_DSparse = Us_xcolor(:,1:s_)*Xs_DSparse;%%%%%%%%%%%%%%%%%%
  %X_DSparse = RSparseRecon(X_color, X_DSparse, 25);
  X_Gray = reshape(GrayscaleImage,M*N,1);
  X_Gray = cast(X_Gray,'double');
  X_DSparse = RSparseRecon2(X_Gray, X_DLow, 25);

end
