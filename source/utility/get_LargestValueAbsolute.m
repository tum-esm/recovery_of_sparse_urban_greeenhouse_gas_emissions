function [app] = get_LargestValueAbsolute(X,n)
X_sor = sort(abs(X(:)),'descend');
% should thresholding the absolut value (without positive/negative sign)
thr = X_sor(n);
app = X;
app(abs(app) < thr) = 0;
while sum(app > 0) > n+1
    numRs = 1;    %can be done before any looping
    while numRs < size(app, 2)
      rand_pick = app(numRs);
      if rand_pick > 0; break; 
      else
          numRs = numRs +1;
      end
    end
    app(numRs) = 0;
end
end

