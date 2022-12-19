function det0=determinant_vec_de_mat(A)
%A representa matrizes puesto en vector sobre la tercera dimension A(:,:,i)
%%%%%%%%%%%%%%%%%%%%%%%  ATTENTION OJO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on suppose ici que toutes les matrices ont les zeros aux memes endroit
%on decompose avec les sous matrices
%on recherche en premier les lignes ou les colonnes ayant le moins de
%termes


[n1,n2,m]=size(A);
if m==1
    det0=det(A);
else
    if n1~=n2
        disp('err: no square matrix');
        return
    elseif n1==1
        det0=squeeze(A);
    elseif n1==2
        det0=squeeze(A(1,1,:).*A(2,2,:)-A(1,2,:).*A(2,1,:));
    elseif n1==3
        det0=squeeze(A(1,1,:).*A(2,2,:).*A(3,3,:)+A(2,1,:).*A(3,2,:).*A(1,3,:)+A(1,2,:).*A(2,3,:).*A(3,1,:)-...
            (A(3,1,:).*A(2,2,:).*A(1,3,:)+A(1,1,:).*A(3,2,:).*A(2,3,:)+A(1,2,:).*A(2,1,:).*A(3,3,:)));
    else
        B=logical(abs(A(:,:,1)));
        [mini,imin]=min(sum(B,2));
        [minj,jmin]=min(sum(B,1));
        ind =1:n1;
        %on suppose que toutes les matrices ont les zeros aux memes endroit
        if mini<=minj
            if mini==2
                ind0=ind(B(imin,:));
                A(:,ind0(2),:)=A(:,ind0(2),:)-A(:,ind0(1),:).*reshape(repmat(squeeze(A(imin,ind0(2),:)./A(imin,ind0(1),:)).',n1,1),n1,1,m);
                i=imin;
                j=ind(ind0(1));
                Aij=A([(1:i-1),(i+1:n1)],[(1:j-1),(j+1:n1)],:);
                det0=(-1)^(i+j)*squeeze(A(i,j,:)).*determinant_vec_de_mat(Aij);
            else
                i=imin;
                det0=zeros(m,1);
                for j=ind(B(imin,:))
                    Aij=A([(1:i-1),(i+1:n1)],[(1:j-1),(j+1:n1)],:);
                    det0=det0+(-1)^(i+j)*squeeze(A(i,j,:)).*determinant_vec_de_mat(Aij);
                end
            end
        else
            if minj==2
                ind0=ind(B(:,jmin));
                A(ind0(2),:,:)=A(ind0(2),:,:)-A(ind0(1),:,:).*reshape(repmat(squeeze(A(ind0(2),jmin,:)./A(ind0(1),jmin,:)).',n1,1),1,n1,m);
                i=imin;
                j=ind(ind0(1));
                Aij=A([(1:i-1),(i+1:n1)],[(1:j-1),(j+1:n1)],:);
                det0=(-1)^(i+j)*squeeze(A(i,j,:)).*determinant_vec_de_mat(Aij);
            else
                det0=zeros(m,1);
                j=jmin;
                for i=ind(B(:,jmin))
                    Aij=A([(1:i-1),(i+1:n1)],[(1:j-1),(j+1:n1)],:);
                    det0=det0+(-1)^(i+j)*squeeze(A(i,j,:)).*determinant_vec_de_mat(Aij);
                end
            end
        end
    end
end