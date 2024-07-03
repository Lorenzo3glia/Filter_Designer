function [gap C1]=calcolagap(Cg,h,w,er)
%CO deve essere espresso in picofarad ed w in metri di conseguenza anche h
%Si utilizza un epsilon_r = 9.6 in base alle formule 

wh=w/h;
mo=wh*(0.619*log10(wh)-0.3853);
ko=4.26-1.453*log10(wh);
C1=0;

for n=1:5
    %C2=Cgap
    CO=2*Cg+C1;
    %g/w
    gw=(CO/w*exp(-ko))^(1/mo);
    
    %calcolo di me - ke 
    %if gw>= 0.1 && gw<=0.3
    if(gw<=0.3)
        me=0.8675;
        ke=2.043*wh^0.12;
   
    elseif gw>= 0.3 && gw<=1
        me=(1.565/(wh)^0.16 -1);
        ke=1.97-0.03/wh;
    else
        warning('g/w fuori range')
        gw=NaN;
        break
    end
    %CE calcolata in pico farad
    CE = w*gw^me * exp(ke);
    C1=CE/2;
    
    C1=C1*(er/9.6)^0.8;
end
%risultato in metri
gap=gw;
