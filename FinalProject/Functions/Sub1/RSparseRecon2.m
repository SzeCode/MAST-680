function X_DSparseR = RSparseRecon2(DATA, X_DLow, thresh)

 X_DLow = abs(X_DLow);

    a = X_DLow - min(min(X_DLow));
    b = max(max(X_DLow)) - min(min(X_DLow));
    X_DLow = 255*(a/b);


  X_DSparseR = DATA - abs(X_DLow);



  X_DSparseR(X_DSparseR<0) = 0;

end
