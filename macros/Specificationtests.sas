%macro SpatialTest(Data= , Modele= , Poids=)/Des="Tests pour détecter une mauvaise spécification du modèle linéaire standard";
options nonotes;
ods exclude all;
/*Créer deux variables macros intégrant la variable endogène et les variables explicatives*/
%let ENDOG=%Scan(&Modele,1,"~"); %let EXOG=%Scan(&Modele,2,"~");  %let LEXOG= %scan(&EXOG,1,"+") ; %let q=2;
%if %length(%scan(&EXOG,2,"+")) gt 0 %then %do;
%do %until(%length(%scan(&EXOG,&q,"+"))=0);
	%let LEXOG=&LEXOG, %scan(&EXOG,&q,"+") ;
	%let q=%eval(&q+1);
%end;
%end;

/*************************** Commençons les bonnes choses maintenant***********************************/
proc iml;
/******************************************* Lecture des données dans une matrice******************************************************************/
Title "Test de mauvaise spécification du modèle linéaire standard : (SARMA, MORAN, LMERR et LMLAG)";
/**** *******************************************Variables explicatives******************************************************************************/
use &Data; read all var {&LEXOG} into X;
/**************************************************Variable endogène******************************************************************************/ 
read all var {&ENDOG} into Y;
/*****************************************Matrice de poids dans le même ordre que les Id******************************************************/
use &Poids;read all into W;
/*Standardiser les lignes*/
W_row=W[,+];
do j=1 to nrow(W_row);
if W_row[j,1]=0 then W_row[j,1]=1;
end;
W=W/W_row;
/*************************************** Calcul des estimateurs sous HO *******************************************************************/
K=J(nrow(X),1,1); 
X=K||X;
BETA=inv(t(X)*X)*t(X)*Y;
M=diag(K)-X*inv(t(X)*X)*t(X); 
ERROR=M*Y;
Sigma2=(t(ERROR)*ERROR)/(nrow(X));

/***************************************** Scores sous H0*********************************************/
Sr=(t(ERROR)*W*y)/Sigma2;
Sl=(t(ERROR)*W*ERROR)/Sigma2;
/**********************************Autres matrices du test de SARMA**********************/
T=trace(W*W+t(W)*W);
D=t(W*X*BETA)*M*(W*X*BETA)+T*Sigma2;
E=(D/Sigma2)*T-T**2;
SARMA=(E**(-1)) * ( (Sl**2) * (D/Sigma2) + (Sr**2) * T - 2 * Sl* Sr* T);
PROB=1-PROBCHI(SARMA,2);
RESULTAT=SARMA || 2 || PROB;
NAME={"Statistique", "DDL", "Probabilité"};
OUT1={"&ENDOG" "Maximum de vraisemblance"};
OUT1CN={"Variable dépendante" "Méthode d'esimation"};
OUT11= nrow(X) || ncol(X)-1;
OUT11CN={"Nombre d'observations" "Nombre de régresseurs"};
options notes;
ods select all;
print OUT1[label="" colname=OUT1CN];
print OUT11[label="" colname=OUT11CN];
print RESULTAT[label="Test SARMA" colname=NAME  rowname=" " ];
options nonotes;
ods exclude all;
/************************************** TEST de MORAN********************************************************/
MORAN=(t(ERROR)*W*ERROR)/(t(ERROR)*ERROR); /*Statistique de moran*/
EMORAN=(trace(M*W))/(nrow(X)-ncol(X));                   /*Espérance de moran*/
VMORAN=(trace(M*W*M*t(W))+trace((M*W)**2)+(trace(M*W))**2)/((nrow(X)-ncol(X))*(nrow(X)-ncol(X)+2)); /*Variance de moran*/
SMORAN=(VMORAN)**0.5;
RESULTAT=MORAN||EMORAN||VMORAN;
NAME={"Statistique","Espérance","Variance"};
options notes;
ods select all;
print RESULTAT[label="Indice de Moran" colname=NAME  rowname=" " ];
options nonotes;
ods exclude all;
MORANST=(MORAN-EMORAN)/SMORAN;
PROB1=CDF("normal",MORANST,0,1);   /*unilatérale à gauche*/
PROB2=1-CDF("normal",MORANST,0,1);   /*unilatérale à droite*/
PROB3=2-2*CDF("normal",abs(MORANST),0,1);   /*bilatérale*/
RESULTAT=MORANST || PROB1 || PROB2 ||PROB3;
NAME={"Statistique standardisée", "Prob unilatérale à gauche H1: I<0", "Prob unilatérale à droite H1: I>0", "Prob bilatérale: Alternative générale"};
options notes;
ods select all;
print RESULTAT[label="Test de Moran" colname=NAME rowname=" " ];
options nonotes;
ods exclude all;
/*******************************************Test LMERR alternative autocoréélation résidu seulement*************************************************/
LMERR=(((t(ERROR)*W*ERROR)/Sigma2)**2)/(T);
PROB=1-PROBCHI(LMERR,1);
RESULTAT=LMERR || 1 || PROB;
/*Version robuste LMERR*/
RLMERR=((Sl-T*Sigma2*(D**(-1))*Sr)**2)/(T-(T**2)*Sigma2*D**(-1));
PROB=1-PROBCHI(RLMERR,1);
RESULTAT1=RLMERR || 1 || PROB;
RESULTAT=RESULTAT //  RESULTAT1;
COLNAME={"Statistique", "DDL", "Probabilité"};
ROWNAME={"Test LMERR Simple", "Test LMERR Robuste"};
options notes;
ods select all;
print RESULTAT[label="Test LMERR :  Autocorrélation spatiale des résidus (SEM) " colname=COLNAME  rowname=ROWNAME ];
options nonotes;
ods exclude all;
/**************************************Test LMMAG alternative variable endogène spatiale décalée**************************************************/
T1=D/Sigma2;
LMLAG=(((t(ERROR)*W*Y)/Sigma2)**2)/(T1);
PROB=1-PROBCHI(LMLAG,1);
RESULTAT=LMLAG || 1 || PROB;
/*Version robuste LMLAG*/
RLMLAG=((Sr-Sl)**2)/(Sigma2**(-1)*D-T);
PROB=1-PROBCHI(RLMLAG,1);
RESULTAT1=RLMLAG || 1 || PROB;
COLNAME={"Statistique", "DDL", "Probabilité"};
ROWNAME={"Test LMLAG Simple", "Test LMLAG Robuste"};
RESULTAT=RESULTAT // RESULTAT1;
options notes;
ods select all;
print RESULTAT[label="Test LMLAG :  Modèle Autorégressif (SAR)" colname=COLNAME  rowname=ROWNAME];
options notes;
ods select all;
quit;
title " ";
%put ***************************************************************************************;
%put *              MERCI D AVOIR UTILISE LES MACROS "REGRESSIONS SPATIALES"               *;                                                      
%put *                          DE Elysée Aristide HOUNDETOUNGAN                           *;  
%put *                          Courriel :  ariel92and@gmail.com                           *; 
%put *                  N hésitez pas à envoyer un mail en cas de besoin                   *; 
%put ***************************************************************************************;

%mend;
