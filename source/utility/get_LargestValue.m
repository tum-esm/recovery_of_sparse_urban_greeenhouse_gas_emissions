function app = get_LargestValue(X,perc)
X_sor = sort(abs(X(:)),'descend');
% should thresholding the absolut value (without positive/negative sign)
n = ceil(length(X_sor)*perc);
thr = X_sor(n+1);
app = X;
app(abs(app) < thr) = 0;
end

