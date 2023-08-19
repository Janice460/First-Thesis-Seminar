$offSymXRef
$OnMultiR
$offSymList

option limrow = 0;
option limcol = 0;
option solprint=on;
option sysout=off;

$TITLE DP

scalar a;
scalar v;
set puntos /1*5/;
parameter m(puntos) /1 EPS,2 EPS,3 EPS,4 EPS,5 EPS/ ;
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

set esc /1*1/;
set t /1*4/;

Table caudales(t,esc)
    1
1   14
2   17
3   18
4   15;

alias(puntos,puntosp);
loop(t,
    loop(puntos,
        v = almacen(puntos);
        a = caudales(t,'1');
        solve problema1 using DNLP minimizing F;
        alpha(puntos) = F.l + EPS;
    );
    display alpha;
    
    loop(puntos, 
        m(puntos) = (sum(puntosp$[ord(puntosp)=ord(puntos)+1],alpha(puntosp))-alpha(puntos))/25;
        b(puntos) = alpha(puntos)-m(puntos)*almacen(puntos);
    );
    display m, b;
);

Execute_unload 'Resultados';
