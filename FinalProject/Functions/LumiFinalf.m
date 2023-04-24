function SM_LUMVIEW = LumiFinalf(X_DLow_Lum, X_DSparse_Lum, M, N)

  X_DLow_Lum = sum(X_DLow_Lum')';
  X_DSparse_Lum = sum(X_DSparse_Lum')';

  %------------- Luminosity Normalization -------------%

  X_DLow_Lum = real(X_DLow_Lum);

  X_DSparse_Lum = real(X_DSparse_Lum);
  a = X_DSparse_Lum - min(min(X_DSparse_Lum));
  b = max(max(X_DSparse_Lum)) - min(min(X_DSparse_Lum));
  X_DSparse_Lum_norm = 255*(a/b);

  SM_LUMVIEW = reshape(X_DSparse_Lum_norm(:,1),M,N);


end
