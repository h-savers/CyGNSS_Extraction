%   index = detect_coherence(DDM, SNR) returns the indices corresponding to
%   the coherent entries of the Delay Doppler Map DDM under the given SNR.
%
%   DDM is a matrix of N x 187 for CyGNSS Delay Doppler Map data.
%   SNR is a vector with N length.
%
%   The method is based on https://doi.org/10.1109/TGRS.2020.3009784

function index = detect_coherence(DDM,SNR)

%Pre-processing
if size(DDM,2)~=187
    error('The input matrix should be reshaped to N x 187.');
end
[N]=size(DDM,1);
DDM=reshape(DDM',[11,17,N]);

%Noise-exclusion threshold in Fig.2 of the reference paper
x=0:5:30;
y=[70 40 15 7.5 5 0.25 0];
ex_noise=polyval(polyfit(x,y,3),SNR)/100;

%Find the indices of the sharp peak in DDMs
[ind_dop,ind_tau,cin,cout]=deal(zeros(N,1));
for a=1:N
    [value,tmp]=max(DDM(:,:,a));
    [~,ind_tau(a)]=max(value);
    ind_dop(a)=tmp(ind_tau(a));
end

%Reserve edge to prevent from exception as the method use +/-2 Doppler
%and +/- 1 Delay for calculation
ind_tau(ind_tau==17)=16;ind_tau(ind_tau==1)=2;
ind_dop(ind_dop<3)=3;ind_dop(ind_dop>9)=9;

%Reduce thermal noise and calculate PR values according to (1) of the
%reference papar
for a=1:N
    tmp=DDM(:,:,a);tmp(tmp<ex_noise(a)*DDM(ind_dop(a),ind_tau(a),a))=0;
    cin(a)=sum(sum(DDM(ind_dop(a)-2:ind_dop(a)+2,ind_tau(a)-1:ind_tau(a)+1,a)));
    cout(a)=sum(sum(tmp))-cin(a);
end
pr=cin./cout;

%Return indices
% index=find(pr>=2);

if (pr>=2)
    index=1;
else
    index=0;
end
