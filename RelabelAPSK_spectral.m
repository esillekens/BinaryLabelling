function C = RelabelAPSK_spectral(C)

[M,D] = size(C);

% bits
m = log2(M);

mA = floor(m/2);
mP = ceil(m/2);

% ensure more phase than amplitude
if mA==mP
    mA = mA-1;
    mP = mP+1;
end

% shift balance if single quadrant
if sum(all(C>=0))==2
    mA = mA+1;
    mP = mP-1;
end

% order
MP = 2^mP;
MA = 2^mA;

% phase and amplitude
phi = atan2(C(:,2),C(:,1));
amp = sqrt(sum(C.^2,2));

% prepare mapping
mapP = MA*graymap(MP);
mapA = graymap(MA);
mapOut = zeros(M,1);

% use Fiedler vector for amplitude clustering
idx_cluster = spectral_clustering(Polar2Cart(amp,mod(phi, pi/2/MP)), MA);

for i = 1:MA
    idx = find(idx_cluster==i);
    [~,argsortA] = sort(phi(idx)); % argsort the slice according to phase
    mapOut(idx(argsortA)) = mapP + mapA(i);
end

C(mapOut+1,:) = C;
end

function mapping = graymap(order)
    j = int32((0:order-1)');
    mapping = cast(bitxor(j,bitshift(j,-1)),'like',order);
end

function idx = spectral_clustering(X, k)

M = size(X,1);

EucD = pdist2(X,X,'squaredeuclidean');

SNR = 12;
gamma = 10.^(SNR/10);

A = exp(-gamma*EucD) - eye(M);
D = sum(A,1);

L = diag(D) - A; %lagrangian
L = 1./sqrt(D).' .* L .* 1./sqrt(D); %normalisation

[x, ~] = eig(L);
x = x(:,1:k);
vectors = 1./sqrt(D).' .* x;
[~,iii] = sort(vectors(:,2));

idx = zeros(M,1);
idx(iii) = kron((1:k).',ones(M/k,1));

end

function X = Polar2Cart(amp,phi)

X = amp.*[cos(phi),sin(phi)];

end