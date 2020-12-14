function [Obs]=sleepOBS(tt,n,T)

%NumN=length(unique(n));
Obs=zeros(1,T);

Obs(tt)=n;

Obs=Obs+1;

return