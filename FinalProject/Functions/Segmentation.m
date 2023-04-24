function IMSegment = Segmentation(OriginalImage, SalientMap, M, N, T)

IMSegment = zeros(M,N,3);
for i = 1:M
  for j = 1:N

    if(SalientMap(i,j) > T)
      IMSegment(i,j,:) = OriginalImage(i,j,:);
    end

  end
end
