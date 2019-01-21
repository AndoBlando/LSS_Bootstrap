
% function: bootstrap.m
% input
%	f = function for linear spectral statistic (LSS) or Non-linear spectral statistic (NlSS)
% 	code = setting (population parameters, bootstrap replicates, etc.); see set_params.m

% e.g. bootstrap(@(x) max(x), 1)
function bootstrap(f,code)
[n,p,k,distr,Sigma,lambda_pop,N,B] = set_params(code);

pop_stat = zeros(N,1);
khat_vec = zeros(N,1);
QuEST_error = zeros(N,1);
true_eigs=sort(diag(lambda_pop));
boot_stat = zeros(N,3);

% population sqrtm
% sig_pop_sqrt = sqrtm( lambda_pop );

rng shuffle;
for i=1:N
	Xobs = distr(n,p)*sqrtm(Sigma);
	[~,~,tauhat,~,~,~,~,~,~,~,~]=QuESTimate(Xobs,0); % run QuEST algorithm
	QuEST_error(i) = mean( abs( sort(tauhat) - true_eigs ) ); % average error
	pop_eig_hat = real(eig( (Xobs')*Xobs/n ));
	pop_eig_hat = sort( pop_eig_hat );
	khat_vec(i) = kappa_est(Xobs); % kurtosis estimation
	pop_stat(i) = f(pop_eig_hat)-f(true_eigs); % population LSS / NlSS

	sig_hat_sqrt = sqrtm( diag(tauhat) );
	boot_vec=zeros(B,1);

% bootstrap stats
   for b=1:B    
      W = pearsrnd_alt(0,1,0,khat_vec(i),n,p);
      % lapack eigenvalue 
      % C = lapack('dsyev', 'N', 'U', p, sig_hat_sqrt*(W')*(W)*sig_hat_sqrt/n,...
      %   p, zeros(p,1),zeros(3*p-1,1), 3*p-1, 1);
      % eig_star = real(C{6});
      eig_star = eig(sig_hat_sqrt*(W')*(W)*sig_hat_sqrt/n);
      eig_star = sort(eig_star);
      % center by QuEST version
      boot_vec(b,:) = f(eig_star) - f(sort(tauhat));
   end
   % Calculate summary stats for given
   boot_stat(i,:) = [mean(boot_vec),std(boot_vec),prctile(boot_vec,95)];
end
population_results = [mean(pop_stat),std(pop_stat),prctile(pop_stat,95)];

% population first row, bootstrap second row.
comparison_table = [population_results;mean(boot_stat)]
end

% pearson rng function, allows for Rademacher
function mat = pearsrnd_alt(m1,m2,m3,m4,n,p)
	% generate Rademacher with kappa=1
	if m4==1
		mat = (rand(n,p) < .5)*2-1;
	else
		mat = pearsrnd(m1,m2,m3,m4,n,p);
	end
end