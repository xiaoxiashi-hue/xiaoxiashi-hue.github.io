%This is to compare nonparametric iid bootstrap and the asymptotic
%distribution in the simple location model. 
%The statistic that the bootstraps apply to is the t-statistic:
%
%                    t = sqrt(n)*(X_mean-theta0)/sigma_hat
%Let X+1 be drawn from a chi-square (1) istribution.
%
%Let theta0=0
%
%The asymptotic distribution of t is N(0,1)
%
%The nonparmetric iid bootstrap draws observations from the multinomial
%distribution (X1,...,Xn)*(1/n,...,1/n)
%
%The bootstrap t should be:
%                   t = sqrt(n)*(Xbt_mean-X_mean)/sigmabt_hat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=100; %Sample Size
B=99999; %Number of Bootstrap samples drawn

X = randn(n,1).^2-1; %generating the original sample from a chi2(1)-1 distribution
mX = mean(X);
sigma_hat = sqrt(sum((X-mX).^2)/(n-1));

%% Simulated the distribution of the Bootstrap t-statistic
rs=RandStream('mt19937ar');
tbt = zeros(B,1); %bootstrap t statistic correctly recentered (recentered at mX)
%tbt2 = zeros(B,1); %bootstrap t statistic recentered at zero
for i = 1:1:B;
    indbt = ceil(rand(rs,n,1)*n);
    Xbt = X(indbt);
    mXbt = mean(Xbt);
    sigmabt_hat = sqrt(sum((Xbt-mXbt).^2)/(n-1));
    tbt(i) = sqrt(n)*(mXbt-mX)/sigmabt_hat;
    %tbt2(i) = sqrt(n)*(mXbt-mX)/sigma_hat;
end

%% Simulate the distribution of the original t-statistic
rs=RandStream('mt19937ar');
Xsim = randn(rs,n,B).^2-1;
Xmsim = mean(Xsim);
sigsim_hat = sqrt((sum((Xsim-kron(Xmsim,ones(n,1))).^2))/(n-1));
tsim = sqrt(n)*Xmsim./sigsim_hat;

%% Report the result in a graph
x = -4:0.1:4;
f = ksdensity(tbt,x);              %Draw the density of the bootstrap t-statistic
fsim = ksdensity(tsim,x);          %Draw the density of the original t-statistic
%[f2,xi2]=ksdensity(tbt2);
z = normpdf(x,0,1);                %Draw the standard normal density
plot(x,fsim,'.',x,z,'--',x,f,'-')
