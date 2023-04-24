function [LOW, SPARSE, Us] = LowSparseDMD3(X_DATA, s_)

  [M_DATA, N_DATA] = size(X_DATA);

  %Creating X and Y matrices out of the original Data Matrix
  X(:,1:N_DATA-1) = [X_DATA(:,1:N_DATA-1)];
  Y(:,1:N_DATA-1) = [X_DATA(:,2:N_DATA)];

  [U,S,V] = svd(X,'econ');
  A = U(:,1:s_)'*Y*V(:,1:s_)*inv(S(1:s_,1:s_));
  Xs = U(:,1:s_)'*X;

  dt = 1;
  t = (0:dt:dt*(N_DATA-1));

  %fprintf('\nObtaining Eigen Vectors and values and Omega of DMD Matrix A...\n');
  [e_Vect, e_Val] = eig(A); %finds eigenvectors and eigenvalues of A
  e_Val_D = diag(e_Val); %Take all the diagonals of eigenvalues
  omega = log(e_Val_D)/dt;
  absOmega = abs(omega);
  absOmegaSorted = sort(absOmega);

  Xs_DLow1 = zeros(s_,N_DATA);
  Xs_DSparse1 = zeros(s_,N_DATA);

  %Low and sparse reconstruction
  b = pinv(e_Vect)*Xs(:,1);
  for i = 1:length(A)
      if(absOmega(i) == min(absOmega))
        Xs_DLow1 = Xs_DLow1 + b(i)*e_Vect(:,i)*exp(omega(i)'*t);
      else
        Xs_DSparse1 = Xs_DSparse1 + b(i)*e_Vect(:,i)*exp(omega(i)'*t);
      end
  end

  LOW = Xs_DLow1;
  SPARSE = Xs_DSparse1;
  Us = U(:,1:s_);

end
