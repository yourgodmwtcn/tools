function maskout=ZPP_process_mask(maskin)
%
maskout=maskin;
%
[M,L]=size(maskout);
Mm=M-1;
Lm=L-1;
Mmm=Mm-1;
Lmm=Lm-1;
%
neibmask=0.*maskout;
neibmask(2:Mm,2:Lm)=maskout(1:Mmm,2:Lm)+maskout(3:M,2:Lm)+...
                    maskout(2:Mm,1:Lmm)+maskout(2:Mm,3:L);
%
while sum(sum((neibmask(2:Mm,2:Lm)>=3 & maskout(2:Mm,2:Lm)==0)|...
              (neibmask(2:Mm,2:Lm)<=1 & maskout(2:Mm,2:Lm)==1)))>0
%
  maskout(neibmask>=3 & maskout==0)=1;
  maskout(neibmask<=1 & maskout==1)=0;
%
  maskout(1,2:Lm)=maskout(2,2:Lm);
  maskout(M,2:Lm)=maskout(Mm,2:Lm);
  maskout(2:Mm,1)=maskout(2:Mm,2);
  maskout(2:Mm,L)=maskout(2:Mm,Lm);
%
  maskout(1,1)=min([maskout(1,2) maskout(2,1)]);
  maskout(M,1)=min([maskout(M,2) maskout(Mm,1)]);
  maskout(1,L)=min([maskout(1,Lm) maskout(2,L)]);
  maskout(M,L)=min([maskout(M,Lm) maskout(Mm,L)]);
%
  neibmask(2:Mm,2:Lm)=maskout(1:Mmm,2:Lm)+maskout(3:M,2:Lm)+...
                      maskout(2:Mm,1:Lmm)+maskout(2:Mm,3:L);
%		    
end
%
% Be sure that there is no problem close to the boundaries
%
maskout(:,1)=maskout(:,2);
maskout(:,L)=maskout(:,Lm);
maskout(1,:)=maskout(2,:);
maskout(M,:)=maskout(Mm,:);
%
