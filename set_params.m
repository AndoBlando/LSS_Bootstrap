% function set_params.m
% input:
%   code = (integer) specify population parameters to simulate.
% output:
%   n = sample size for each experiment
%   p = dimension of each sample
%   k = population kurtosis for population
%   Sigma = population covariance matrix
%   lambda_pop = population eigenvalues (as diagonal matrix)
%   N = number of experiements to simulate
%   B = number of bootstrap replications (Monte-Carlo)

function [n,p,k,distr,Sigma,lambda_pop,N,B] = set_params(code)
% overall parameters
% beta parameters - shape, scale
alph = 6; 
bet = alph;
% t-dist parameters
nu=9; % df 
%
N=100; 
B=50; 
%
switch code
   case 1
       n=100; p=20; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);      
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       Sigma = U*lambda_pop*U';
   case 2
       n=500; p=400; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       Sigma = U*lambda_pop*U';
   case 3
       n=500; p=600; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       Sigma = U*lambda_pop*U';
   case 4
       n=500; p=200; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( (1:p).^(-.5) );
       Sigma = Q*lambda_pop*Q';
   case 5
       n=500; p=400; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( (1:p).^(-.5) );
       Sigma = Q*lambda_pop*Q';
   case 6
       n=500; p=600; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( (1:p).^(-.5) );
       Sigma = Q*lambda_pop*Q';
   case 7
       n=500; p=200; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
      [Q,R] = qr(randn(p,p));
       lambda_pop = diag( sort((1:p).^(-.5)) );
       Sigma = Q*lambda_pop*Q';
   case 8
       n=500; p=400; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
        [Q,R] = qr(randn(p,p));
       lambda_pop = diag( sort((1:p).^(-.5)) );
       Sigma = Q*lambda_pop*Q';
   case 9
       n=500; p=600; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
        [Q,R] = qr(randn(p,p));
       lambda_pop = diag( sort((1:p).^(-.5)) );
       Sigma = Q*lambda_pop*Q';
   case 10
       n=500; p=200; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       Sigma = U*lambda_pop*U';
   case 11
       n=500; p=400; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       Sigma = U*lambda_pop*U';
   case 12
       n=500; p=600; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       Sigma = U*lambda_pop*U';
    case 13
       n=500; p=200; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       Sigma = load('realsigma200.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
    case 14
       n=500; p=400; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       Sigma = load('realsigma400.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));        
    case 15
       n=500; p=600; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       Sigma = load('realsigma600.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
    case 16
       n=500; p=200; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12); 
       Sigma = load('realsigma200.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
    case 17
       n=500; p=400; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12); 
       Sigma = load('realsigma400.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));        
    case 18
       n=500; p=600; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12); 
       Sigma = load('realsigma600.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
    % Sigma = identity, Uniform 
    case 19
       n=500; p=200; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12); 
       Sigma = eye(p);
       lambda_pop = diag(eig(Sigma));
    case 20
       n=500; p=400; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12); 
       Sigma = eye(p);
       lambda_pop = diag(eig(Sigma));        
    case 21
       n=500; p=600; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12); 
       Sigma = eye(p);
       lambda_pop = diag(eig(Sigma));
    case 22
       n=500; p=200; k=3; distr=@(n,p) randn(n,p);
       Sigma = eye(p);
       lambda_pop = diag(eig(Sigma));
    case 23
       n=500; p=400; k=3; distr=@(n,p) randn(n,p); 
       Sigma = eye(p);
       lambda_pop = diag(eig(Sigma));        
    case 24
       n=500; p=600; k=3; distr=@(n,p) randn(n,p); 
       Sigma = eye(p);
       lambda_pop = diag(eig(Sigma));
   case 25
       n=5*50; p=2*50; k=3; distr=@(n,p) randn(n,p);
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( sort((1:p).^(-.5)) );
       Sigma = Q*lambda_pop*Q';
   case 26
       n=500; p=400; k=3; distr=@(n,p) randn(n,p);
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( sort((1:p).^(-.5)));
       Sigma = Q*lambda_pop*Q';
   case 27
       n=500; p=600; k=3; distr=@(n,p) randn(n,p);
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( sort((1:p).^(-.5) ));
       Sigma = Q*lambda_pop*Q';
   case 28
       n=5*100; p=2*100; k=3; distr=@(n,p) randn(n,p);
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       Sigma = U*lambda_pop*U';
   case 29
       n=500; p=400; k=3; distr=@(n,p) randn(n,p);
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lmbda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       Sigma = U*lambda_pop*U';
   case 30
       n=500; p=600; k=3; distr=@(n,p) randn(n,p);
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       Sigma = U*lambda_pop*U';
   case 31
       n=500; p=200; k=3; distr=@(n,p) randn(n,p);
       Sigma = load('realsigma200.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
   case 32
       n=500; p=400; k=3; distr=@(n,p) randn(n,p);
       Sigma = load('realsigma400.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));        
   case 33
       n=500; p=600; k=3; distr=@(n,p) randn(n,p);
       Sigma = load('realsigma600.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
   case 34
       n=500; p=200; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( (1:p).^(-.5) );
       Sigma = Q*lambda_pop*Q';
   case 35
       n=500; p=400; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( (1:p).^(-.5) );
       Sigma = Q*lambda_pop*Q';
   case 36
       n=500; p=600; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       [Q,R] = qr(randn(p,p));
       lambda_pop = diag( (1:p).^(-.5) );
       Sigma = Q*lambda_pop*Q';
   case 37
       n=500; p=200; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       Sigma = U*lambda_pop*U';
   case 38
       n=500; p=400; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       Sigma = U*lambda_pop*U';
   case 39
       n=500; p=600; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       [U,R]=qr(randn(p));
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       Sigma = U*lambda_pop*U';
   case 40
       n=500; p=200; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       Sigma = load('realsigma200.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
   case 41
       n=500; p=400; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       Sigma = load('realsigma400.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));        
   case 42
       n=500; p=600; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       Sigma = load('realsigma600.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
   case 43
       n=500; p=256; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       lambda_pop = diag( (1:p).^(-.5) );
       U = (1/sqrt(p))*hadamard(p);
       Sigma = U*lambda_pop*U';
   case 44
       n=500; p=512; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       lambda_pop = diag( (1:p).^(-.5) );
       U = (1/sqrt(p))*hadamard(p);
       Sigma = U*lambda_pop*U';
   case 45
       n=500; p=1024; k=9/5; distr=@(n,p) (-1+2*rand(n,p))/sqrt(4/12);
       lambda_pop = diag( (1:p).^(-.5) );
       U = (1/sqrt(p))*hadamard(p);
       Sigma = U*lambda_pop*U';
   case 46
       n=500; p=200; k=4; distr=@(n,p) (trnd(10,[n,p]))./sqrt(10/(10-2));
       Sigma = ar1_sigma(p,1/2);
       lambda_pop = diag(eig(Sigma));
   case 47
       n=500; p=400; k=4; distr=@(n,p) (trnd(10,[n,p]))./sqrt(10/(10-2));
       Sigma = ar1_sigma(p,1/2);
       lambda_pop = diag(eig(Sigma));
   case 48
       n=500; p=600; k=4; distr=@(n,p) (trnd(10,[n,p]))./sqrt(10/(10-2));
       Sigma = ar1_sigma(p,1/2);
       lambda_pop = diag(eig(Sigma));
   case 49
       n=500; p=200; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       Sigma = eye(p);
       lambda_pop = Sigma; 
   case 50
       n=500; p=400; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       Sigma = eye(p);
       lambda_pop = Sigma; 
   case 51
       n=500; p=600; k=4; distr=@(n,p) (trnd(10, [n,p]))./sqrt(10/(10-2));
       Sigma = eye(p);
       lambda_pop = Sigma; 
   case 52
       n=500; p=200; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );       
       lambda_pop = diag( sort((1:p).^(-.5)) );
       [U,R]=qr(randn(p));
       Sigma = U*lambda_pop*U';
   case 53
       n=500; p=400; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = diag( sort((1:p).^(-.5)) );
       [U,R]=qr(randn(p));
       Sigma = U*lambda_pop*U';
   case 54
       n=500; p=600; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = diag( sort((1:p).^(-.5)) );
       [U,R]=qr(randn(p));
       Sigma = U*lambda_pop*U';
   case 55
       n=500; p=200; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       [U,R]=qr(randn(p));
       Sigma = U*lambda_pop*U';
   case 56
       n=500; p=400; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       [U,R]=qr(randn(p));
       Sigma = U*lambda_pop*U';
   case 57
       n=500; p=600; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = eye(p);
       lambda_pop(1:10,1:10) = 3*eye(10);
       lambda_pop = diag(sort(diag( lambda_pop )));
       [U,R]=qr(randn(p));
       Sigma = U*lambda_pop*U';
   case 58
       n=500; p=200; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       Sigma = load('realsigma200.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
    case 59
       n=500; p=400; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       Sigma = load('realsigma400.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
    case 60
       n=500; p=600; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       Sigma = load('realsigma600.txt');
       Sigma = Sigma / eigs(Sigma,1);
       lambda_pop = diag(eig(Sigma));
   case 61
       n=500; p=200; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = eye(p);
       Sigma = lambda_pop;
   case 62
       n=500; p=400; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
        lambda_pop = eye(p);
       Sigma = lambda_pop;
  case 63
       n=500; p=600; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = eye(p);
       Sigma = lambda_pop;
  case 64
       n=500; p=200; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
        lambda_pop = eye(p);
       Sigma = lambda_pop;
  case 65
       n=500; p=400; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       lambda_pop = eye(p);
       Sigma = lambda_pop;
  case 66
       n=500; p=600; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       lambda_pop = eye(p);
       Sigma = lambda_pop;
  case 67
  	   n=500; p=200; k=3; distr=@(n,p) avg_proc(n,p);
  	   lambda_pop = eye(p);
  	   Sigma = lambda_pop;
  case 68
  	   n=500; p=400; k=3; distr=@(n,p) avg_proc(n,p);
  	   lambda_pop = eye(p);
  	   Sigma = lambda_pop;
  case 69
  	   n=500; p=600; k=3; distr=@(n,p) avg_proc(n,p);
  	   lambda_pop = eye(p);
  	   Sigma = lambda_pop;
  case 70
  	   n=500; p=200; k=3; distr=@(n,p) randn(n,p);
  	   lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
	   Sigma = ar1_sigma(p,Phi); 
  case 71
	   n=500; p=400; k=3; distr=@(n,p) randn(n,p);
	   lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
	   Sigma = ar1_sigma(p,Phi);
  case 72
	   n=500; p=600; k=3; distr=@(n,p) randn(n,p);
	   lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
	   Sigma = ar1_sigma(p,Phi);
  case 73
     n=500; p=200; k=3; distr=@(n,p) randn(n,p);
	   lambda_pop = eye(p);
	   Sigma = lambda_pop;
  case 74
	   n=500; p=400; k=3; distr=@(n,p) randn(n,p);
	   lambda_pop = eye(p);
	   Sigma = lambda_pop;
  case 75
	   n=500; p=600; k=3; distr=@(n,p) randn(n,p);
	   lambda_pop = eye(p);
  	 Sigma = lambda_pop;
  case 76
       n=500; p=200; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
       Sigma = lambda_pop;
  case 77
       n=500; p=400; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
        lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
       Sigma = lambda_pop;
  case 78
       n=500; p=600; k=beta_kurt(alph,bet) ; distr=@(n,p) (betarnd( repmat(alph,[n,p]), repmat(bet,[n,p]),n,p ) - repmat(alph/(alph+bet),n,p) )./sqrt( alph*bet/((alph+bet)^2*(alph+bet+1)) );
       lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
       Sigma = lambda_pop;
  case 79
       n=500; p=200; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
       Sigma = lambda_pop;
  case 80
       n=500; p=400; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
       Sigma = lambda_pop;
  case 81
       n=500; p=600; k=3+ 6/(nu-4); distr=@(n,p) (trnd(nu, [n,p]))./sqrt(nu/(nu-2));
       lambda_pop = diag( eig( ar1_sigma(p,Phi) ) );
       Sigma = lambda_pop;
end

end

% kurtosis for beta distribution based on shape/scale parameters
function kurt=beta_kurt(a,b)
extra = 6*( (a-b)^2 *(a+b+1) - a*b*(a+b+2) )/(a*b*(a+b+2)*(a+b+3));
kurt = 3 + extra;
end

function mat = ar1_sigma(m,phi)
mat = repmat( phi, [m,m]);
for i=1:m
    for j=1:m
	mat(i,j) = mat(i,j)^( abs(i-j) );
    end
end
end

function mat = VAR(n,p)
mat = zeros(2*n,p);
mat(1,:) = randn(1,p);
for i=2:(2*n)
mat(i,:) = 0.8*mat(i-1,:) + randn(1,p);
end
mat = mat((n+1):(2*n),:);
end

function mat = avg_proc(n,p)
Xt = randn(n+1,p);
mat = (Xt(2:(n+1),:) + Xt(1:n,:))./sqrt(2);
end

function vec = sim_ar(p, phi,sig)
Xt = zeros(p+200,1);
Xt(1) = randn(1,1);
for i=2:(p+200)
Xt(i) = phi*Xt(i-1) + randn(1,1)*sig;
end
vec = Xt(201:(p+200));
end

function mat = avg_proc_ar(n,p,phi,sig)
Xt = zeros(n+1,p);
for i=1:(n+1)
	Xt(i,:) = sim_ar(p,phi,sig);
end
mat = (Xt(2:(n+1),:) + Xt(1:n,:))./sqrt(2);
end

