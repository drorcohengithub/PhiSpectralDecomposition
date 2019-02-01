function systems = Example_systems()

%% systems
N = 2; % # components
K = 2; % order

systems = [];

%% system 1 - disconnected only inst
SIG = [1 0.4; 0.4 0.7];


% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.4 0;
    0 0.4]; % connectivity matrix at time t-1 (Aij = from j to i (column to row)

A(:,:,2) = [-0.25 0;
    0 -0.25]; % connectivity matrix at time t-2


systems(1).A = A;
systems(1).SIG = SIG;
systems(1).nm='disconnected_with_inst';

%% system 2 - uni directional no inst
SIG = [1 0; 0 0.7];


% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.2 0;
            0.4 0.2]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.25 0;
            -0.2 0.1]; % connectivity matrix at time t-2

systems(2).A = A;
systems(2).SIG = SIG;
systems(2).nm='unidir_no_inst';



%% system 3 - unidir with inst
SIG = [1 0.65; 0.65 0.7];
normalizer = 1;%sqrt(det(SIG));
SIG = SIG/normalizer;


% coeff
A = systems(2).A;


systems(3).A = A;
systems(3).SIG = SIG;
systems(3).nm='unidir_with_inst';

%% system 4 - bidir no inst
SIG = [1 0; 0 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0.2;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2


systems(4).A = A;
systems(4).SIG = SIG;
systems(4).nm='bidir_no_inst';

%% system 5 - bidir with inst

SIG = [1 0.35; 0.35 0.9];


% coeff
A = systems(3).A;
A(1,2,1) = 0.5;
A(1,2,2) = 0.15;

systems(5).A = A;
systems(5).SIG = SIG;



% SIG = [1 0.1; 0.1 0.7];
% 
% % coeff
% A = zeros(N,N,K);
% A(:,:,1) = [0.9 0.2;
%     0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i
% 
% A(:,:,2) = [-0.5 0;
%     -0.2 -0.5]; % connectivity matrix at time t-2


systems(5).A = A;
systems(5).SIG = SIG;
systems(5).nm='bidir_with_inst';

%% system 6 - 5 dimensionsal system with  bidir with inst
N = 5;
rng(1);
off_diag = rand(N);
off_diag = (off_diag + off_diag')/4;
%create the noise cov
SIG = eye(N) + off_diag ;

A = zeros(N,N,2);

A(:,:,1) = [0.3 0.2   0.1 -0.3 -0.2;
            0.2 0.1  -0.3 -0.2  0.3;
            0    0    -0.2  0.3  0.2;
            0    0.1    0.3  0.2  0.1;
            0    0    0.2  0.1 -0.3];
         
         
A(:,:,2) = [-0.3 -0.2 0.1 0.3  0.2;
            -0.2 -0.1 0.3 0.2  0.3;
             0    0   0.2 0.3  0.2;
             0    0   0.3 0.2  0.1;
             0    0   0.2 0.1 0.3]/2;

% SIG = eye(N);
% A  = var5_test;

% coeff
systems(6).A = A;
systems(6).SIG = SIG;
systems(6).nm='5dim';


