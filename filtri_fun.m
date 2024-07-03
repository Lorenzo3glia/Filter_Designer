function Cg=filtri_fun(df,ripple,f0,N,z0)


%% Filtri (pdf di riferimento Estratto Collin filtri 1)
  
z1=z0;
z2=z1;
y0=1/z0;

c0=3*1e8;

%Implementiamo le formule del passa basso 
%% Prototipo passa basso
k=1:N;              %vettore di interi da 1 ad N
tipo='C';
if tipo=='B'
    %% Formule per il massimamente piatto (8.100b)
    g=2*sin((2*k-1)/(2*N)*pi);
    g(N+1)=1;       %Resistenza del generatore
%% Formule per il chebyshev formule(8.101a-8.101c)
elseif tipo=='C'
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

%% Scalatura impedenza
g
g(1:2:N)=g(1:2:N)/z2;
g(2:2:N)=g(2:2:N)*z2;
g(N+1)=g(N+1)*z2;

%% Passa Banda -LAYOUT A
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

disp(['L         C'])
disp([L'*1e9 C'*1e12])
disp(['Zgen: '])
g(N+1)

%% Passaggio in parallelo (serie parlallelo) INVERTER
rapporto=((pi*y0)/2)^2;
prodotto=1/(w0^2);
C0=sqrt(prodotto*rapporto);
L0=sqrt(prodotto/rapporto);
disp(['L0         C0'])
disp([L0*1e9 C0*1e12])

K(1)=z2
K(2:2:N)= sqrt(sqrt(L(2:2:N)./(C(2:2:N)))*sqrt(L0/C0));
K(3:2:N+1)=sqrt(sqrt(L(2:2:N)./(C(2:2:N)))*sqrt(L0/C0));
if mod(N,2)==1
    K(N+1)=g(N+1);
end

L(2:2:N)=L0;
C(2:2:N)=C0;

disp(['L         C'])
disp([L'*1e9 C'*1e12])
disp('.....................')
disp([K'])

 K(1:2:N)=K(1:2:N).*sqrt(sqrt(L0/C0)*sqrt(C(1:2:N)./L(1:2:N)));
 K(2:2:N+1)=K(2:2:N+1).*sqrt(sqrt(L0/C0)*sqrt(C(1:2:N)./L(1:2:N)));
 L(1:2:N)=L0;
 C(1:2:N)=C0;
 K(N+1)=K(N+1)*sqrt(z1/g(N+1));
 disp('.........K............')
 disp([K'])

 %% Calcolo capacità di GAP
 J = 1./K;
 theta = atan(J/y0);
 Cg=0.5*tan(2*theta)*y0/w0
 disp('.........Cg............')
 disp([Cg*1e12])
for n=1:N
    lung(n)=pi-theta(n)-theta(n+1);
end

lung
lfisica=lung*c0/w0 %Sono in metri,per convertire moltiplicare per mille
 




