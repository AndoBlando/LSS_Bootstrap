% kurtosis estimation function 
% input:
% Xobs = (n,p) matrix distributed as X = Z \Sigma^{1/2}, where Z (n,p) i.i.d entries.
% output:
% k = estimated kurtosis for each Z (thresholded to be at least 1)
function k = kappa_est(Xobs)
n = length(Xobs(:,1));
sigHat = (1/n)*(Xobs')*(Xobs);
Frobhat = trace(sigHat*sigHat) - (1/n)*(( trace(sigHat) )^2);
sigma2 = std( diag( Xobs*(Xobs') ) )^2;
sig4hat = sum( diag(sigHat).^2 );
k = max( 3 + (sigma2 - 2*Frobhat)/(sig4hat), 1 );
end
