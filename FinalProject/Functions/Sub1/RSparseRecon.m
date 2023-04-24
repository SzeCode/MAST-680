function X_DSparseR = RSparseRecon(DATA, X_DLow, thresh)

  X_DSparseR = DATA - abs(X_DLow);
  X_DSparseR_Copy = X_DSparseR;

  R = zeros([size(X_DSparseR)]);

  X_DSparseR_Copy(X_DSparseR_Copy<0) = 0;
  R = X_DSparseR - X_DSparseR_Copy;
  R(R<-thresh) = 2.5*R(R<-thresh);%-40 for test data

  X_DSparseR = X_DSparseR - R;

  %X_DSparseR(X_DSparseR>180) = 255;

end
