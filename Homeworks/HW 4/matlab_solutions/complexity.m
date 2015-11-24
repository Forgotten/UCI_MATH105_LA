% script to estimate the asymptotic complexity

%% testing for the built-in solver interfaced with LAPACK
nSamples = 10;
times = zeros(nSamples,1);
sizes = 2*2.^(0:nSamples-1);
for i = 1:nSamples
    m = sizes(i);
    % creating the Matrix 
    A = rand(m,m) + m*eye(m);
    % creating the rhs
    b = rand(size(A,1),1);
    % timing 
    tic(); 
    %solving the system
    x = A \ b;
    times(i) = toc();
end

%% benchmarking the Gaussian solver

times2 = zeros(nSamples,1);
for i = 1:nSamples
    m = sizes(i);
    % creating the Matrix 
    A = rand(m,m) + m*eye(m);
    % creating the rhs
    b = rand(size(A,1),1);
    % timing 
    tic(); 
    %solving the system
    x = SolveGauss(A,b);
    times2(i) = toc();
end

%% plotting everything

%making a new figure and cleaning up frame
figure(1); clf(); 
% plotting in loglog scale
loglog(sizes, times, 'r-o');
hold on; 
loglog(sizes, times2, 'b-.');
loglog(sizes, 0.000001*sizes.^3, 'k');
xlabel('Number of unknows');
ylabel('time [s]')
legend('LAPACK', 'Homebrew Gaussian Elimination', 'third order slope')
set(gca,'FontSize', 14)