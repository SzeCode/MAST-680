function s_ = SRankApprox(EnergyTarget, S_MATRIX)

  s_ = 1;
  [M_S,N_S] = size(S_MATRIX);

  sum_Sigma_s = 0;
  sum_Sigma_all = 0;
  Energy = 0;
  %Searches for a rank value that will make the Energy bigger than or equal to 99.99% to have good approximation.
    while(Energy < EnergyTarget)
      s_ = s_+1;
      sum_Sigma_s = 0;
      sum_Sigma_all = 0;
      for i = 1:s_
        sum_Sigma_s = sum_Sigma_s + S_MATRIX(i,i)^2;
      end
      for i = 1:N_S
        sum_Sigma_all = sum_Sigma_all + S_MATRIX(i,i)^2;
      end

      Energy = sum_Sigma_s/sum_Sigma_all;
    end

    %fprintf('\n\tLow Rank s number achieved: %d', s_);
    %fprintf('\n\twith Energy level = %d\n', Energy);

end

