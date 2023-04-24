function [X_DLow_Lum, X_DSparse_Lum] = LumiDMDf(X_Lum, M_xLum, N_xLum)

  [U_xLum, S_xLum, V_xLum] = svd(X_Lum,'econ');
  [M_xLum, N_xLum] = size(X_Lum);

  %initializes variables for low rank calculations with energy
  s_ = N_xLum-1;

  X_Lum_rank = U_xLum(:,1:s_)'*X_Lum;
  [Xs_DLow_Lum, Xs_DSparse_Lum, Us_XLum] = LowSparseDMD3(X_Lum, s_);

  X_DLow_Lum = zeros(M_xLum,N_xLum);
  X_DSparse_Lum = zeros(M_xLum,N_xLum);

  X_DLow_Lum = Us_XLum(:,1:s_)*Xs_DLow_Lum;%%%%%%%%%%%%%%%%%%%%%%%%
  %X_DSparse_Lum = Us_XLum(:,1:s_)*Xs_DSparse_Lum;%%%%%%%%%%%%%%%%%%
  X_DSparse_Lum = RSparseRecon(X_Lum, X_DLow_Lum, 25);

end
