function [C] = smoothing(A,n)
B = (1/n^2) * ones(n);
xLayer = [size(A,1)+size(B,1)-1, size(A,2)+size(B,2)-1];
A2 = padarray(A,xLayer,'circular');
C2 = conv2(A2, B, 'same');
C = C2(xLayer(1)+1:end-xLayer(1), xLayer(2)+1:end-xLayer(2));
end

