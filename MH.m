% This is the auxiliary function of Metropolis-Hastings Algorithm for PS6.
function [X,acc] = MH(mu,sigma,b,se,prior)
% b is the estimate of beta (or error). prior is the indicator of a flat
% priors. se is the standard error of betas (or error).
% parameters
% lag = 1; % iterations between successive samples
nsamp = 2000; % number of samples to draw
sig = sigma; % standard deviation of Gaussian proposal
x = -1; % start point

x_min = -3;
x_max = 3;
bins = linspace(x_min, x_max, 40);

% storage
X = zeros(nsamp,1); % samples drawn from the Markov chain
acc = [0 0]; % vector to track the acceptance rate
topplot = subplot(2,1,1);
% curve = fplot(@(x) targetdist(x,mu,sigma,b,se,prior), [x_min x_max]);
curve = ezplot(@targetdist, [x_min x_max]);
hold on

bottomplot = subplot(2,1,2);
axis([x_min x_max 1 nsamp]);

hold on

% MH routine
for i = 1:nsamp  % this loop is the main algorithm
    % for j = 1:lag
    % iterate chain one time step:
    xp = normrnd(x,sig); % generate candidate from Gaussian
    accprob = targetdist(xp,mu,sigma,b,se,prior) / targetdist(x,mu,sigma,b,se,prior); % acceptance probability
    u = rand; % uniform random number
    if u <= accprob % if accepted
        xnew = xp; % new point is the candidate
        a = 1; % note the acceptance
        acctext = '-> accepted';
    else % if rejected
        xnew = x; % new point is the same as the old one
        a = 0; % note the rejection
        acctext = '-> rejected';
    end
    acc = acc + [a 1]; % track accept-reject status
    % end
    X(i,1) = xnew; % store the i-th sample
end

% update plot
figure
subplot(2,1,1)
bars = hist(X(1:i), bins);
bars = bars./100;
bbb = bar(bins, bars, 'Facecolor', [0.8 0.8 0.8] );
%     alpha(bbb, 0.5)
uistack(bbb, 'bottom')
theta = plot(x,0,'o','MarkerFaceColor','red');
theta_new = plot(xp,0,'o','MarkerEdgeColor','red');
ttt = {['Old \theta_n: ', num2str(x), ';  L(\theta_n) = ', num2str(targetdist(x,mu,sigma,b,se,prior))], ...
    ['Proposed \theta_{n+1}: ', num2str(xp), ';  L(\theta_{n+1}) = ', num2str(targetdist(xp,mu,sigma,b,se,prior))], ...
    ['Acceptance probability: ', num2str(min(1, accprob))], ...
    ['Random number drawn: ', num2str(u), acctext]};
t = text(-2.5, 2.5, ttt);

subplot(2,1,2)
thetaline = plot(X(1:i), 1:i);
thetalinedot = plot(x,i,'o','MarkerFaceColor','red');

% (un)comment the following two lines if program should run one
% sample at a time and then pause vs run all iterations at once
pause
drawnow

% remove plot items
if i < nsamp
    delete( [t, theta, theta_new, bbb, thetalinedot] )
end

x = xnew;

% L(theta) = prior * likelihood
function probX = targetdist(x,mu,sigma,b,se,prior)
if prior == 1
    probX = prior*1/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/2*sigma^2);
else
    probX = 1/sqrt(2*pi*se^2)*exp(-(x-b)^2/2*se^2)*...
            1/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2/2*sigma^2);
end