function [ diceVal ] = Dice( X,Y )
% Y: original clusters label; X: identified clusters label
N=length(Y);
K=length(unique(Y));

diceVal=0;DiceRes=zeros(1,K);
for k1=1:K
    k1Ind=find(Y==k1);lenk1=length(k1Ind);
    dice_vec=zeros(1,K);
    for k2=1:K
        k2Ind=find(X==k2);lenk2=length(k2Ind);
        dice_vec(k2)=2*length(intersect(k1Ind,k2Ind))/(lenk1+lenk2);
    end
    DiceRes(k1)=max(dice_vec);
end
diceVal=sum(DiceRes)/K;


