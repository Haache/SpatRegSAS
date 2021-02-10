%macro MoranTest(Data= ,Var= , Poids=)/Des="Tests pour détecter une
autocorrélation spatiale";
options nonotes;
ods exclude all;
proc iml;
Title "Tests d'autocorrélation spatiale de MORAN et de GEARY";
/******************************************* Lecture des données dans une
matrice******************************************************************/
use &Data; read all var {&Var} into Y;
/*****************************************Matrice de poids dans le même ordre
que les Id**********************************************************/
use &Poids; read all into W;
/*Standardiser les lignes*/
W_row=W[,+];
do j=1 to nrow(W_row);
if W_row[j,1]=0 then W_row[j,1]=1;
end;
W=W/W_row;
N=nrow(Y);
Wcol=W[+,];
Wrow=W[,+];
K=(sum((Y-mean(Y))##4)/N)/((sum((Y-mean(Y))##2))/N)**2;
/******************************************* Calcul de la
statistique******************************************************************
*/
MORAN=(t(Y-mean(Y)) * W * (Y-mean(Y)))/(t((Y-mean(Y))) * (Y-mean(Y)));
EMORAN=-1/(N-1); /*Espérance de moran*/
S0=sum(W);
*S1=0.5*sum((W+t(W))##2);
S1=sum(W#W+W#t(W));
S2=sum((Wcol+t(Wrow))##2);
/*VMORAN=(N**2*S1-N*S2+3*S0**2)/((N-1)*(N+1)*S0**2)-(1/(N-1))**2; /*Variance
de moran*/
VMORAN=(((N*((N**2-3*N+3)*S1-3*S2+3*S0**2))-(K*((N**2-3)*S1-
2*N*S2+6*S0**2))))/((N-1)*(N-2)*(N-3)*S0**2)-EMORAN**2;
SMORAN=(VMORAN)**0.5;
RESULTAT=MORAN||EMORAN||VMORAN;
NAME={"Statistique","Espérance","Variance"};
OUT= nrow(Y);
OUTCN={"Nombre d'observations"};
OUT1={&VAR};
options notes;
ods select all;
print OUT1[label=" " colname={"Variable"}] OUT[label=" " colname=OUTCN] ;
print RESULTAT[label="Indice de Moran" colname=NAME rowname=" "];
options nonotes;
ods exclude all;
MORANST=(MORAN-EMORAN)/SMORAN;
PROB1=CDF("normal",MORANST,0,1); /*unilatérale à gauche*/
PROB2=1-CDF("normal",MORANST,0,1); /*unilatérale à droite*/
PROB3=2-2*CDF("normal",abs(MORANST),0,1); /*bilatérale*/
RESULTAT=MORANST || PROB1 || PROB2 ||PROB3;
NAME={"Statistique standardisée", "Prob unilatérale à gauche H1: I<0", "Prob
unilatérale à droite H1: I>0", "Prob bilatérale: Alternative générale"};
options notes;
ods select all;print RESULTAT[label="Test de Moran" colname=NAME rowname=" " format=10.5];
options nonotes;
ods exclude all;
/********************************************************************GEARY
INDEX*****************************************************/
Y2=Y##2;
Wcol=W[+,];
GEARY=(N-1)*(sum(Y2)-2*t(Y)*W*Y+Wcol*Y2)/(2*S0*sum((Y-mean(Y))##2));
EGEARY=1;
/*VGEARY=( (2*S1+S2) * (N-1) -4*S2 )/( 2*(N+1)*S0 );*/
/*VGEARY=S0/(2*S1+S2);*/
VGEARY=(((N-1)*S1*(N**2-3*N+3-(N-1)*K))-(0.25*((N-1)*S2*(N**2+3*N-6-(N**2-
N+2)*K)))+(S0**2*(N**2-3-(N-1)**2*K)))/(N*(N-2)*(N-3)*S0**2);
SGEARY=VGEARY**0.5;
RESULTAT=GEARY||EGEARY||VGEARY;
NAME={"Statistique","Espérance","Variance"};
options notes;
ods select all;
print RESULTAT[label="Indice de Geary" colname=NAME rowname=""];
options nonotes;
ods exclude all;
GEARYST=(GEARY-EGEARY)/SGEARY;
PROB1=CDF("normal",GEARYST,0,1); /*unilatérale à gauche*/
PROB2=1-CDF("normal",GEARYST,0,1); /*unilatérale à droite*/
PROB3=2-2*CDF("normal",abs(GEARYST),0,1); /*bilatérale*/
RESULTAT=GEARYST || PROB1 || PROB2 ||PROB3;
NAME={"Statistique standardisée", "Prob unilatérale à gauche H1: G<1", "Prob
unilatérale à droite H1: G>1", "Prob bilatérale: Alternative générale"};
options notes;
ods select all;
print RESULTAT[label="Test de Geary" colname=NAME rowname=" " format=10.5];
quit;
title " ";
%put
*****************************************************************************
**********;
%put * MERCI D AVOIR UTILISE LES MACROS "REGRESSIONS SPATIALES"
*;
%put * DE Elysée Aristide HOUNDETOUNGAN
*;
%put * Courriel : ariel92and@gmail.com
*;
%put * N hésitez pas à envoyer un mail en cas de besoin
*;
%put
*****************************************************************************
**********;
%mend;
