function [Cov_E_p, A_p, masked_Delta] = adjmat_argminGV_mask(Cov_X,Cov_XY,K,split_mask_A,split_mask_E,initial_A,iter_max,gamma,min_error)
 
%some things for the mex file
%The classes of the variables 
assert(isa(Cov_X, 'double'))
assert(isa(Cov_XY, 'double'))
assert(isa(K, 'double'))
assert(isa(split_mask_A, 'double'))
assert(isa(split_mask_E, 'double'))
assert(isa(initial_A, 'double'))
assert(isa(iter_max, 'double'))
assert(isa(gamma, 'double'))
assert(isa(min_error, 'double'))

% the maximal size of Cov_X is 50by50
assert(all(size(Cov_X) <= [50 50]));
% the maximal size of Cov_XY is 50by50by3000
assert(all(size(Cov_XY) <= [50 50 3000]));
% the maximal size of mask is 50by50
assert(all(size(split_mask_A) <= [50 50]));
% the maximal size of mask is 50by50
assert(all(size(split_mask_E) <= [50 50]));
% the maximal size of A is 50by50by3000
assert(all(size(initial_A) <= [50 50 3000]));

% these are all scalars
assert(all(size(iter_max) == [1 1]));
assert(all(size(gamma) == [1 1]));
assert(all(size(min_error) == [1 1]));


% set initial values of the connectivity matrix in disconnected model
N   = size(Cov_X,1);
A_p = initial_A; 

%make sure this executes at least once
if iter_max < 1
    error('iter_max needs to be set to at least 1. Much larger values are likely neccessary for meaningful estimates')   
    
end

if K > 1
    Cov_B = zeros(N*K,N*K);
    for i=1: K
        for j=i: K
            if i== j
                Cov_B((i-1)*N+1:i*N,(j-1)*N+1:j*N) = Cov_X;
            else
                Cov_B((i-1)*N+1:i*N,(j-1)*N+1:j*N) = Cov_XY(:,:,(j-i));
                Cov_B((j-1)*N+1:j*N,(i-1)*N+1:i*N) = Cov_XY(:,:,(j-i))';
            end
            
        end
    end
    
else
    Cov_B = Cov_X;
end

A_B = zeros(N,N*K);
Cov_XY_B = zeros(N,N*K);
for j=1: K
    A_B(:,(j-1)*N+1:j*N) = A_p(:,:,j);
    Cov_XY_B(:,(j-1)*N+1:j*N) = Cov_XY(:,:,j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make sure these are set to somthing to keep matlab/codegen happy
Cov_E_p =  nan; 
masked_Delta = nan;
%
Delta2=zeros(1,iter_max);        
fprintf('Optimzing...\n');  
for iter=1: iter_max

    if mod(iter,500)==0
        fprintf('Iteration number: %d, Delta: %d \n', int32(iter),Delta2(iter-1));
    end
    
   Cov_E_p =  Cov_X - Cov_XY_B*A_B' - A_B*Cov_XY_B' + A_B*Cov_B*A_B'; % covariance matrix of residuals in disconnected model
   Cov_E_p = split_mask_E.*Cov_E_p;
    
   Delta = zeros(N,N,K); % derivative of determinant of Cov_E_p with respect to the diagonal elements of A_p
    for k=1: K
        Delta(:,:,k) = - Cov_XY(:,:,k);
        for l=1: K
            if k > l
                Delta(:,:,k) =  Delta(:,:,k) + A_p(:,:,l)*Cov_XY(:,:,k-l);
            elseif k < l
                Delta(:,:,k) =  Delta(:,:,k) + A_p(:,:,l)*Cov_XY(:,:,l-k)';
            else
                Delta(:,:,k) =  Delta(:,:,k) + A_p(:,:,l)*Cov_X;
            end
        end
        Delta(:,:,k) = Cov_E_p\Delta(:,:,k);
    end
    
    % update A_p
    big_mask = repmat(split_mask_A,1,1,K);
    masked_Delta = big_mask.*Delta;
    A_p_tmp = A_p - gamma*masked_Delta;
    A_p = A_p_tmp;
%     
%     for k=1: K
%         for i=1: N
%             for j=1: N
%                 A_p(i,j,k) = A_p(i,j,k) - gamma*Delta(i,j,k);
%             end
%         end
%     end
    
%     for j=1: K
%         A_B(:,(j-1)*N+1:j*N) = A_p(:,:,j);
%     end
    A_B = reshape(A_p,size(A_p,2),size(A_p,2)*size(A_p,3));
    
% check convergence
%     Delta2 = 0;
%     for k=1: K
%         for i=1: N
%             for j=1: N
%                 Delta2 = Delta2 + Delta(i,j,k)^2;
%             end
%         end
%     end
    %Dynamically adjust 
    Delta2(iter) = sqrt(sum(masked_Delta(:).^2));
    if mod(iter,500)==0
        if any(Delta2(iter-100:iter-1)<=Delta2(iter))
%             keyboard()
            gamma = gamma/5;
            disp(['Changing gamma to: ' num2str(gamma)]);
        end
    end
% Delta2 = sqrt(Delta2)
% disp(iter);
    
%     fprintf('iter=%d logdet(Cov_E_p)=%f A_p(1,1)=%f A_p(2,2)=%f\n',iter,logdet(Cov_E_p),A_p(1,1),A_p(2,2));
%     fprintf('Delta2=%f Delta(1,1)=%f Delta(2,2)=%f\n',Delta2,Delta(1,1),Delta(2,2));
    
    if Delta2(iter) < min_error
       fprintf('Convergence reached at iter=%d\n', int32(iter));

%         Disp(['Convergence reached at iter ' num2str(iter)])
        break;
    end
    
end

if iter == iter_max
    fprintf('Max iter reached at iter=%d\n', int32(iter));
end
% Disp(['Max iter reached (' num2str(iter) ')'])

end