$offSymXRef
$OnMultiR
$offSymList

option limrow = 0;
option limcol = 0;
option solprint=on;
option sysout=off;

$TITLE SDP

scalar a;
scalar v;
set puntos /1*5/;
parameter m(puntos) /1 EPS,2 EPS,3 EPS,4 EPS,5 EPS/;
parameter b(puntos) /1 EPS,2 EPS,3 EPS,4 EPS,5 EPS/;

parameter caudal

Variables

u,s,g1,g2,delta,F;

POSITIVE VARIABLES

g1,g2,delta,u,s ;

Equations

funobj, eq1,eq2,eq3,eq4,eq5,eq6;

funobj.. F=e= 10*g1+20*g2+100*delta + max(m("1")*(v-u-s+a)+b("1"), m("2")*(v-u-s+a)+b("2"),m("3")*(v-u-s+a)+b("3"),m("4")*(v-u-s+a)+b("4"));

eq1.. 0 =l= v-u-s+a ;
eq2.. v-u-s+a =l= 100 ;
eq3.. u =l= 50 ;
eq4.. g1+g2+0.9*u+delta =e= 45 ;
eq5.. g1 =l= 20 ;
eq6.. g2 =l= 25 ; 

Model problema1 /funobj,eq1,eq2,eq3,eq4,eq5,eq6/ ;

parameter almacen(puntos) /1 EPS,2 25,3 50,4 75,5 100/;
parameter alpha(puntos) /1 EPS,2 EPS,3 EPS,4 EPS,5 EPS/;

set esc /1*3/;

parameter alphaesc(esc) /1 EPS,2 EPS,3 EPS/;

set t /1*4/;

Table caudales(t,esc)
    1   2   3
1   14  12  16
2   17  13  15
3   25  18  19
4   15  20  17;

* Prom de caudales
alias(puntos,puntosp);
loop(t,
    loop(puntos,
        v = almacen(puntos);
        loop(esc,
            a = caudales(t,esc);
            solve problema1 using DNLP minimizing F;
            alphaesc(esc) = F.l + EPS;
        alpha(puntos) = (alphaesc('1')+alphaesc('2')+alphaesc('3'))/3;)
    );
    display alpha;
    
    loop(puntos, 
        m(puntos) = (sum(puntosp$[ord(puntosp)=ord(puntos)+1],alpha(puntosp))-alpha(puntos))/25;
        b(puntos) = alpha(puntos)-m(puntos)*almacen(puntos);
    );
    display m,b;
);

Execute_unload 'Resultados';
