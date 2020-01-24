function dydt = cosmofuncVariableRelease(t,x)

global gL gA dL dA cA cL rA KL KA nL nA

dydt = zeros(size(x));

NL=x(1);%pop density of L-A+
NA=x(2);%pop density of A-L+
SL=x(3);%conc of Lys
SA=x(4);%conc of Hyp

%%variable release rate
%%Col 1: hyp conc; Col 2: release rate of L.
R=[0.00	0.52
0.58	0.83
0.65	0.78
0.73	0.53
0.80	0.38
1.07	0.08];
Hyp=R(:,1);
rLM=R(:,2);% Measured rL
rLi=interp1(Hyp,rLM,SA); %interpolated rL.

if rLi>0
    rL=rLi;
else rL=0;
end

dydt(1)=gL*SL^nL/(SL^nL+KL^nL)*NL-dL*NL;
dydt(2)=gA*SA^nA/(SA^nA+KA^nA)*NA-dA*NA;
dydt(3)=-cL*gL*SL^nL/(SL^nL+KL^nL)*NL+rL*NA;
dydt(4)=-cA*gA*SA^nA/(SA^nA+KA^nA)*NA+rA*NL;

end


%dydt(1)=gL*SL*NL-dL*NL;
%dydt(3)=-cL*gL*SL*NL+rL*NA*dA;



