dflin=linspace(50,300,101)*1e6;
ripplelin=linspace(0.01,1,101);
f0=1e9;
N=4;
z0=50;

for nf=1:length(dflin)
    for nr=1:length(ripplelin)
        df=dflin(nf);
        ripple=ripplelin(nf);
        Cg=filtri_fun(df,ripple,f0,N,z0)*1e12;

        Cg1(nf,nr)=Cg(1);
        Cg2(nf,nr)=Cg(2);
        Cg3(nf,nr)=Cg(3);

    end
    h=waitbar(nf/length(dflin),'Progress')
end
