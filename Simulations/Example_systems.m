function systems = Example_systems()

%% systems
N = 2; % # components
K = 2; % order

systems = [];

%% system 1 - disconnected only inst
SIG_f = [1 0.1; 0.1 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0;
    0 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    0 -0.5]; % connectivity matrix at time t-2


systems(1).A = A;
systems(1).SIG_f = SIG_f;
systems(1).nm='disconnected_with_inst';

%% system 2 - uni directional no inst
SIG_f = [1 0; 0 0.7];
% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2

systems(2).A = A;
systems(2).SIG_f = SIG_f;
systems(2).nm='unidir_no_inst';



%% system 3 - unidir with inst
SIG_f = [1 0.1; 0.1 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2


systems(3).A = A;
systems(3).SIG_f = SIG_f;
systems(3).nm='unidir_with_inst';

%% system 4 - bidir no inst
SIG_f = [1 0; 0 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0.2;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2


systems(4).A = A;
systems(4).SIG_f = SIG_f;
systems(4).nm='bidir_no_inst';

%% system 5 - bidir with inst
SIG_f = [1 0.1; 0.1 0.7];

% coeff
A = zeros(N,N,K);
A(:,:,1) = [0.9 0.2;
    0.16 0.8]; % connectivity matrix at time t-1 (Aij = from j to i

A(:,:,2) = [-0.5 0;
    -0.2 -0.5]; % connectivity matrix at time t-2


systems(5).A = A;
systems(5).SIG_f = SIG_f;
systems(5).nm='bidir_with_inst';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

