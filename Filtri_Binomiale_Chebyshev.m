clc
close all
clear

disp('Benvenuto nel software di calcolo: Filtri')
disp('Realizzato da Lorenzo Treglia')
%% Fornisco i dati obbligatori
strz0 = sprintf('Inserire impedenza Z0: ');
z0=input(strz0);

% %Nel caso fornito dal professore é fissa a 120 ohm
% strzc = sprintf('Inserire impedenza ZC: ');
% z1=input(strzc);
z1=z0;


strfc= sprintf('Inserire frequenza di centrobanda: ');
f0=input(strfc);

strbw= sprintf('Inserire larghezza di banda bilatera: ');
df=input(strbw);

strzN= sprintf('Inserire ordine del filtro: ');
N=input(strzN);

% strzripp= sprintf('Inserire il ripple desiderato [dB]: ');
% ripple=input(strzripp);

%%
z2=z1;
y0=1/z0;
c0=3*1e8;
%Implementiamo le formule del passa basso 
%% Prototipo passa basso
k=1:N;  
strztpe= sprintf('Inserire il tipo di filtro desiderato( B per Binomiale C per Chebyshev): ');
tipo=input(strztpe,'s');

if tipo=='B'
    %% Formule per il massimamente piatto (8.100b)
    g=2*sin((2*k-1)/(2*N)*pi);
    g(N+1)=1;       %Resistenza del generatore
%% Formule per il chebyshev formule(8.101a-8.101c)
elseif tipo=='C'
        strzripp= sprintf('Inserire il ripple desiderato [dB]: ');
        ripple=input(strzripp);
        k2=10^(ripple/10)-1;        %impongo il k^2
        b=log((sqrt(1+k2)+1)/(sqrt(1+k2)-1));
        ak=sin((2*k-1)/(2*N)*pi);
        bk=sinh(b/(2*N))^2+sin(k*pi/N).^2;
        g(1)=2*ak(1)/sinh(b/(2*N));
        for k=2:N           %Siccome abbiamo già calcolato il primo elemento dobbiamo fare un for che parte dal secondo elemento
            g(k)=4*ak(k-1)*ak(k)/bk(k-1)/g(k-1);
        end
        if mod(N,2) %Assegnazione del valore dell'impedenza del generatore
            g(N+1)=1;
        else
            g(N+1)=1/(2*k2+1-2*sqrt(k2)*sqrt(1+k2));
        end

else 
    error('Tipologia di filtro errata')
end
g(N+1)=1/g(N+1);

%% Scalatura impedenza (Passaggio numero 2 nell'ordine)
%Viene effettuata una scalatura per adattarci al carico di riferimento
g
g(1:2:N)=g(1:2:N)/z2;
g(2:2:N)=g(2:2:N)*z2;
g(N+1)=g(N+1)*z2;

%% Passa Banda -LAYOUT A (terzo passaggio)
%Il layout differisce solamente l'elemento in testa al filtro se una
%capacità oppure un'induttanza
%f1 , f2 sono rispettivamente frequenza superiore ed inferiore
f1=f0-df/2;
f2=f0+df/2;
w0=2*pi*sqrt(f1*f2); %Posso scriverlo in funzione della frequenza facendo uscire fuori dalla radice 2*pi
w1= 2*pi * f1;
w2=2*pi*f2;
%Consideriamo la tipologia A come modello di riferimento le G dispari sono
%le capcita mentre le G pari solo le induttanze;


%% Risuonatore Parallelo 
C(1:2:N)=g(1:2:N)/(w2-w1);  %Le capacità dispari sono le g dispari
L(1:2:N)=1./(w0^2*C(1:2:N)); %Le L dispari corrispondenti saranno

%% Risuonatore Serie
L(2:2:N)=g(2:2:N)/(w2-w1);
C(2:2:N)=1./(w0^2*L(2:2:N));

%Le induttanze sono nano henry e le capacità pico Farad
disp(['L         C'])
disp([L'*1e9 C'*1e12])
disp(['Zgen: '])
g(N+1)

%% Passaggio in parallelo (serie parlallelo) INVERTER
%Quello che si fa in questo passaggio è andare a trasformare tutti i
%risuonatori in risuonatori parallelo, in questo passaggio è però anche
%necessario modificare i valori dei risuonatori invariati per uniformarli
%tutti alle stesse dimensioni.
rapporto=((pi*y0)/2)^2;
prodotto=1/(w0^2);
C0=sqrt(prodotto*rapporto);
L0=sqrt(prodotto/rapporto);
disp(['L0         C0'])
disp([L0*1e9 C0*1e12])


K(1)=z2;
K(2:2:N)= sqrt(sqrt(L(2:2:N)./(C(2:2:N)))*sqrt(L0/C0));
K(3:2:N+1)=sqrt(sqrt(L(2:2:N)./(C(2:2:N)))*sqrt(L0/C0));
if mod(N,2)==1
    K(N+1)=g(N+1);
end

L(2:2:N)=L0;
C(2:2:N)=C0;

disp(['L         C'])
disp([L'*1e9 C'*1e12])
disp('...........K(prima modifica)..........')
disp([K'])

% in questo passaggio avviene la conversione dei valori dei risuonatori non
% modificati per portarli a valori di L0 e C0
 K(1:2:N)=K(1:2:N).*sqrt(sqrt(L0/C0)*sqrt(C(1:2:N)./L(1:2:N)));
 K(2:2:N+1)=K(2:2:N+1).*sqrt(sqrt(L0/C0)*sqrt(C(1:2:N)./L(1:2:N)));
 L(1:2:N)=L0;
 C(1:2:N)=C0;
 K(N+1)=K(N+1)*sqrt(z1/g(N+1));
 disp('.........K(seconda modifica)............')
 disp('Da questo momento tutte le L e le C valgono L0 C0')
 disp([K'])

 %% Calcolo capacità di GAP
 %Una volta ottenuti tutti i valori di L0 e C0 bisogna implementare anche
 %gli inverter, quest'ultimi si implementano mediante dei tratti di linea
 %quest'ultimi caratterizzati da alcune capacità di gap
 J = 1./K;
 theta = atan(J/y0);
 Cg=0.5*tan(2*theta)*y0/w0
 disp('.........Cg............')
 disp([Cg*1e12])
for n=1:N
    lung(n)=pi-theta(n)-theta(n+1);
end

lung;
strer = sprintf('Inserire permittività mezzo er[Qucs]: ');
er=input(strer);
strer_eff = sprintf('Inserire permittività efficace er_eff[Qucs]: ');
er_eff=input(strer_eff);

 disp('.........Lunghezza fisica in millimetri............')
lfisica=lung*c0/w0/sqrt(er_eff) %Sono in metri,per convertire moltiplicare per mille
lfisica*1000

%% Voglio calcolare le capacità di gap a 9.6

disp('.........Cg(9.6)............')
Cg96=Cg*(9.6/er)^0.8*1e12;
disp([Cg*(9.6/er)^0.8*1e12])


%w=3.95e-3;
%% Fare attenzione ad inserire i valori di w e h
%per ottenere il valore di w inserire sul qucs impedenza della linea e h
%con er dopo di che premere su analizza per ottenre il risultato in
%lunghezza della linea
w=1.44389e-3;
h=1e-3;

for n=1:length(Cg96)
    [temp_gap, temp_C1] = calcolagap(Cg96(n), h, w,er);
    gap(n) = temp_gap * w * 1e3;
    C1(n) = temp_C1* w * 1e3;
end

disp('Gap in millimetri')
gap
disp('Accorciamenti di shunt')
deltal=C1*1e-12*z0*3e8/sqrt(er_eff)

for n=1:N
    length_finale(n)=lfisica(n)-deltal(n)-deltal(n+1);
end

disp('Lunghezze finali')
length_finale

