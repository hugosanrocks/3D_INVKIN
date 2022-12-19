function det0=subdet(i,j,k,l,M)

det0=M(i,k,:).*M(j,l,:)-M(i,l,:).*M(j,k,:);