%macro GSM(Data= , Modele= , Poids=, CORR=, Rho=0.9999, Lambda=0.9999)/Des="General Spatial Model";
/*options nonotes;
ods exclude all;*/
/*Créer deux variables macros intégrant la variable endogène et les variables explicatives*/
%let ENDOG=%Scan(&Modele,1,"~"); %let EXOG=%Scan(&Modele,2,"~");  %let LEXOG= %scan(&EXOG,1,"+") ;
%let LEXOG1= %scan(&EXOG,1,"+") ; %let LEXOG_= "&LEXOG1" ;%let q=2;
%if %length(%scan(&EXOG,2,"+")) gt 0 %then %do;
%do %until(%length(%scan(&EXOG,&q,"+"))=0);
	%let LEXOG=&LEXOG, %scan(&EXOG,&q,"+") ;
	%let LEXOG1=%scan(&EXOG,&q,"+") ; %let LEXOG_=&LEXOG_ "&LEXOG1" ;
	%let q=%eval(&q+1);
%end;
%end;


proc iml;
/******************************************* Lecture des données dans une matrice**************************************************************************/
/**** *******************************************Variables explicatives**************************************************************************************/
use &Data; read all var {&LEXOG} into X;
/**************************************************Variable endogène**************************************************************************************/ 
read all var {&ENDOG} into Y;
/*****************************************Matrice de poids dans le même ordre que les Id********************************************************************/
use &Poids;read all into W;
/*Standardiser les lignes*/
W_row=W[,+];
do j=1 to nrow(W_row);
if W_row[j,1]=0 then W_row[j,1]=1;
end;
W=W/W_row;

K=J(nrow(X),1,1); 
X=K||X;
I=diag(K);
M=I-X*inv(t(X)*X)*t(X); 
Nrow=nrow(Y);
/****************************************************************************************************************************************************************************/
/*************************** Commençons les bonnes choses maintenant***********************************/
%if &CORR="SAR" %then %do;
/*Estimation des modèles*/
Title "Modèle SAR : Variable endogène spatiale décalée";
/*SAR*/
/*****************************************Ecriture de la vraisemblance**********************************************/
start Logv(Argum) global(I,X, Y, W, Nrow);

Rho_=Argum;
A_=I-Rho_*W;
BETA_=inv(t(X)*X)*t(X)*A_*(Y);
Sigma2_=(1/nrow(X))*t(A_*Y-X*Beta_)*(A_*Y-X*Beta_);
JA_=log(det(A_));
F_=-0.5*Nrow*log(2*arcos(-1))-0.5*Nrow*log(Sigma2_)+JA_- Nrow/2;
return (F_);
finish;

/*Fournir le gradient pour faciliter le calcul des estimateurs*/
start Grad(Argum) global(I,X, Y, W, Nrow);
Rho_=Argum;
A_=I-Rho_*W;
BETA_=inv(t(X)*X)*t(X)*A_*Y;
Sigma2_=(1/Nrow)*t(A_*Y-X*Beta_)*(A_*Y-X*Beta_);
G=-trace((inv(A_))*W)+(1/Sigma2_)*t(A_*Y-X*BETA_)*W*Y;
return (G);
finish;

/*Espace des paramètres d'autocorrélation*/
EIG=(EIGVAL(W))[,1];
MAXE=MAX(EIG); MINE=MIN(EIG);

/*Trouver les paramètres*/
Cons=(MINE) //  (MAXE );
Argum=&Rho;

/*Argum[1]=(sum(t(M*Y)*M*Y))/nrow(X);*/
optn={1 0};

call nlptr(rc,Best,"Logv",Argum ,optn, Cons) grd="Grad"; /*Lancement de la procédure*/ 

/*Construire la matrice Jacobienne*/
Rho=Best; A=I-Rho*W; INVA=inv(A); BETA=inv(t(X)*X)*t(X)*A*(Y); 
Sigma2=(1/nrow(X))*t(A*Y-X*Beta)*(A*Y-X*Beta); 

JACOB=J(nrow(Beta)+2, nrow(Beta)+2, 0);
JACOB[1,1]=(nrow(X))/(2*Sigma2**2);
JACOB[1,2]=(1/Sigma2)*trace(W*INVA);
JACOB[2,1]=(1/Sigma2)*trace(W*INVA);
JACOB[2,2]=trace((W*INVA)**2)+trace(t(W*INVA)*W*INVA) + (1/Sigma2)*trace(t(W*INVA*X*Beta)*W*INVA*X*Beta);
JACOB[2,3:nrow(Beta)+2]=t((1/Sigma2)*t(X)*W*(INVA*X*Beta));
JACOB[3:nrow(Beta)+2,2]=(1/Sigma2)*t(X)*W*(INVA*X*Beta);
JACOB[3:nrow(Beta)+2,3:nrow(Beta)+2]=(1/Sigma2)*t(X)*X;

if  Det(JACOB) =0 then do;
	Print "Le Jacobien est singulier. La variance de Rho est donc biaisée. Ceci est du à l'ordre de grandeur élevé (resp. faible) des variables. Le réduire (resp. augmenter)  en divisant (resp. multipliant) les variables par un même coefficient et réestimer le modèle.";
	options notes;
	ods select all;
end;
else do;
		VAR=inv(JACOB);
		STD=(Vecdiag(VAR))##0.5;

/*Vraisemblance  du modèle contraint pour les coef Rho et Beta*/
	Beta1=inv(t(X)*X)*t(X)*Y; Sigma21=(1/Nrow)*(-t(Y)*X*Beta1+t(Y)*Y);    /*Paramètre du modèle contraint rho=0*/
	LRA=J(2,1,0); 
	LRA[1,1]=Logv(Best);
	LRA[2,1]=-0.5*Nrow*log(2*arcos(-1))-0.5*Nrow*log(Sigma21)-Nrow/2;

/*********************** Statistique des tests ************************************/
	/*Wald pour le coefficient Rho*/
	SWRHO=t(Rho)*inv(Var[2,2])*Rho;   PSWRHO=1-PROBCHI(SWRho,1);/*Test de Wald pour la nullité du paramètre de corrélation*/
	SLR=2*(LRA[1,1]-LRA[2,1]);   PSLR=1-PROBCHI(SLR,1); /*Test de rapport de vraissemblance pour le paramètre de corrélation*/
	SSTU=Rho/(STD[2,1]);   PSSTU=PSWRHO ;/*Test de significativité simple du paramètre de corrélation*/
	SWBETA=t(BETA)*inv(Var[3:ncol(Var),3:ncol(Var)])*BETA;  PSWBETA=1-PROBCHI(SWBETA,nrow(Beta));/*Test de nullité de wald pour la nullité de tous les coefficients des variables explicatives*/
/******************** Critère d'information*************************************/
	AIC=-2*LRA[1,1]+2*(nrow(Beta)+2);
	BIC=-2*LRA[1,1]+(log(nrow(X)))*(nrow(Beta)+2);
/*Test LMERR**/
	ERROR=A*Y - X*BETA;
	Sl=(t(ERROR)*W*ERROR)/Sigma2;
	LMERR2=(Sl**2)/(trace(W*W+t(W)*W)-((trace((W**2)*(inv(A))+t(W)*W*inv(A)))**2)*VAR[2,2]);
	LMPROB=1-PROBCHI(LMERR2,1);

/*Paramétrage des sorties*/
OUT1={"&ENDOG" "Maximum de vraisemblance"};
OUT1CN={"Variable dépendante" "Méthode d'esimation"};
OUT11= nrow(X) || ncol(X)-1;
OUT11CN={"Nombre d'observations" "Nombre de régresseurs"};
OUT2=(LRA[1,1]) // AIC // BIC // Sigma2 ;
OUT2RN={"Vraisemblance du modèle" "Critère AIC" "Critère BIC" "Variance des erreurs"};
OUT3=Rho // Var[2,2]**0.5;
OUT3RN={"Paramètre Rho" "Variance"};
OUT4=(SLR || PSLR) // (SSTU || PSSTU) // (SWRHO || PSWRHO);
OUT4RN={"Rapport de vraisemblance" "Student" "Wald"}; OUT4CN={"Statistique" "Probabilité"};
BETAPROB=J(nrow(BETA),2,0);
do j=1 to nrow(BETA);
BETAPROB[j,1]=BETA[j,1]/STD[j+2,1];
BETAPROB[j,2]=1-PROBCHI(BETAPROB[j,1]**2,1);
end;
OUT5= BETA || STD[3:nrow(STD),1]|| BETAPROB;
OUT5RN={"Constante" &LEXOG_};
OUT5CN={ "Coefficient" "Ecart-type" "Statistique" "Probabilité"};
OUT6=SWBETA || PSWBETA;
OUT6CN={"Statistique" "Probabilité"};
OUT7=LMERR2|| LMPROB ;
OUT7CN={"Statistique" "Probabilité"}; 
options notes;
ods select all;
print OUT1[label="" colname=OUT1CN];
print OUT11[label="" colname=OUT11CN];
print OUT2[label="Informations sur le modèle" rowname=OUT2RN];
print OUT3[label="Estimation du paramètre Rho" rowname=OUT3RN];
print OUT4[label="Significativité de Rho" rowname=OUT4RN colname=OUT4CN format=10.5];
print OUT5[label="Estimation des paramètres du modèle" rowname=OUT5RN colname=OUT5CN format=10.5];
print OUT6[label="Test de Wald de significativité globale des paramètres" colname=OUT6CN format=10.5];
print OUT7[label="Test LMERR* d'autocorrélation des résidus" colname=OUT7CN format=10.5];
end;
quit;
title " ";
%put ***************************************************************************************;
%put *              MERCI D AVOIR UTILISE LES MACROS "REGRESSIONS SPATIALES"               *;                                                      
%put *                          DE Elysée Aristide HOUNDETOUNGAN                           *;  
%put *                          Courriel :  ariel92and@gmail.com                           *; 
%put *                  N hésitez pas à envoyer un mail en cas de besoin                   *; 
%put ***************************************************************************************;
%end;
/*************************************************************************************************************************************************************************************************/

%if &CORR="SEM" %then %do;
/*Estimation des modèles SEM*/
Title "Modèle SEM :  Erreurs spatialement autocorrélées";
/*****************************************Ecriture de la vraisemblance******************************************************************************/
start Logv(Argum) global(I,X, Y, W, Nrow);
Lambda_=Argum;
B_=I-Lambda_*W;
BETA_=inv(t(B_*X)*B_*X)*t(B_*X)*B_*Y;
Sigma2_=(1/Nrow)*t(Y-X*BETA_)*t(B_)*B_*(Y-X*BETA_);
JB_=log(det(B_));

F=-0.5*Nrow*log(2*arcos(-1))-0.5*Nrow*log(Sigma2_)+JB_-Nrow/2;
return (F);
finish;

/*Fournir le gradient pour faciliter le calcul des estimateurs*/
start Grad(Argum) global(I,X, Y, W,Nrow);
Lambda_=Argum;
B_=I-Lambda_*W;
BETA_=inv(t(B_*X)*B_*X)*t(B_*X)*B_*Y;
Sigma2_=(1/Nrow)*t(Y-X*BETA_)*t(B_)*B_*(Y-X*BETA_);

G=-trace(inv(B_)*(W))+(1/Sigma2_)*(t(Y-X*BETA_))*(t(W))*B_*(Y-X*BETA_);
return (G);
finish;

/*Espace des paramètres d'autocorrélation*/
EIG=(EIGVAL(W))[,1];
MAXE=MAX(EIG); MINE=MIN(EIG);

/*Trouver les paramètres*/
Cons=(MINE) //  (MAXE );
Argum=&Lambda;

/*Argum[1]=(sum(t(M*Y)*M*Y))/nrow(X);*/
optn={1 0};
call nlptr(rc,Best,"Logv",Argum ,optn, Cons) grd="Grad"; /*Lancement de la procédure*/ 

/*Estimateur et variance*/
Lambda=Best; B=I-Lambda*W; InvOmega=t(B)*B; BETA=inv(t(X)*InvOmega*X)*t(X)*InvOmega*Y; InvB=inv(B);
Sigma2=(1/Nrow)*t(Y-X*BETA)*InvOmega*(Y-X*BETA);
VARB=inv(t(X)*InvOmega*X)*Sigma2;
STDB=(vecdiag(VARB))##0.5;

/*Construire la matrice Jacobienne, la matrice variance covariance et les écart types des estimateurs sauf Beta*/
JACOB=J(2,2, 0);
JACOB[1,1]=Nrow/(2*Sigma2**2);
JACOB[1,2]=0.5*(1/Sigma2)*trace((W*t(B)+B*t(W))*t(InvB)*InvB;
JACOB[2,1]=0.5*(1/Sigma2)*trace((W*t(B)+B*t(W))*t(InvB)*InvB;
JACOB[2,2]=trace((W*t(B)+B*t(W))*t(InvB)*InvB*W*InvB);
STDL=J(2,1,0);
	if  Det(JACOB) ^=0 then do;
		VARL=inv(JACOB);
		STDL=(vecdiag(VARL))##0.5;
	end;
	else do;
		STD=(1/vecdiag(JACOB))##0.5;
	end;

/*Vraisemblance  du modèle pour les coef Lamda et Beta*/
	Beta1=inv(t(X)*X)*t(X)*Y; Sigma21=(1/Nrow)*(-t(Y)*X*Beta1+t(Y)*Y);    /*Paramètre du modèle contraint rho=0*/
	LRA=J(2,1,0); 
	LRA[1,1]=Logv(Best);
	LRA[2,1]=-0.5*Nrow*log(2*arcos(-1))-0.5*Nrow*log(Sigma21)-Nrow/2;
	
/*********************** Statistique des tests ************************************/
	SWLAMBDA=t(LAMBDA)*inv(STDL[2,1]**2)*LAMBDA;   PSWLAMBDA=1-PROBCHI(SWLAMBDA,1);/*Test de Wald pour la nullité du paramètre de corrélation*/
	SLR=2*(LRA[1,1]-LRA[2,1]);   PSLR=1-PROBCHI(SLR,1); /*Test de rapport de vraissemblance pour le paramètre de corrélation*/
	SSTU=LAMBDA/(STDL[2,1]);   PSSTU=PSWLAMBDA ;/*Test de significativité simple du paramètre de corrélation*/
	SWBETA=t(BETA)*inv(VarB)*BETA;  PSWBETA=1-PROBCHI(SWBETA,nrow(Beta));/*Test de nullité de wald pour la nullité de tous les coefficients des variables explicatives*/
/******************** Critère d'information*************************************/
AIC=-2*LRA[1,1]+2*(ncol(X)+2);
BIC=-2*LRA[1,1]+(log(nrow(X)))*(ncol(X)+2);
OUT1={"&ENDOG" "Maximum de vraisemblance"};
OUT1CN={"Variable dépendante" "Méthode d'esimation"};
OUT11= nrow(X) || ncol(X)-1;
OUT11CN={"Nombre d'observations" "Nombre de régresseurs"};
OUT2=(LRA[1,1]) // AIC // BIC // Sigma2 ;
OUT2RN={"Vraisemblance du modèle" "Critère AIC" "Critère BIC" "Variance des erreurs"};
OUT3=Lambda // STDL[2,1];
OUT3RN={"Paramètre Lambda" "Variance"};
OUT4=(SLR || PSLR) // (SSTU || PSSTU) // (SWLAMBDA || PSWLAMBDA);
OUT4RN={"Rapport de vraisemblance" "Student" "Wald"}; OUT4CN={"Statistique" "Probabilité"};
BETAPROB=J(nrow(BETA),2,0);
do j=1 to nrow(BETA);
BETAPROB[j,1]=BETA[j,1]/STDB[j,1];
BETAPROB[j,2]=1-PROBCHI(BETAPROB[j,1]**2,1);
end;
OUT5=BETA || STDB || BETAPROB ;/*J(nrow(BETA),1,2)-2*(PROBT(abs(BETA/STDB),nrow(X)-ncol(X)));*/
OUT5RN={"Constante" &LEXOG_};
OUT5CN={"Coefficient" "Ecart-type" "Statistique" "Probabilité"};
OUT6=SWBETA || PSWBETA;
OUT6CN={"Statistique" "Probabilité"}; 
options notes;
ods select all;
print OUT1[label="" colname=OUT1CN];
print OUT11[label="" colname=OUT11CN];
print OUT2[label="Informations sur le modèle" rowname=OUT2RN];
print OUT3[label="Estimation du paramètre Lambda" rowname=OUT3RN];
if  Det(JACOB)=0 then do;
	Print "Le Jacobien est singulier. La variance de Lambda est donc biaisée. Ceci est du à l'ordre de grandeur élevé (resp. faible) des variables. Le réduire (resp. augmenter)  en divisant (resp. multipliant) les variables par un même coefficient et réestimer le modèle.";
end;
print OUT4[label="Significativité de Lambda" rowname=OUT4RN colname=OUT4CN format=10.5];
print OUT5[label="Estimation des paramètres du modèle" rowname=OUT5RN colname=OUT5CN format=10.5];
print OUT6[label="Test de Wald de significativité globale des paramètres" colname=OUT6CN format=10.5];
quit;
%put ***************************************************************************************;
%put *              MERCI D AVOIR UTILISE LES MACROS "REGRESSIONS SPATIALES"               *;                                                      
%put *                          DE Elysée Aristide HOUNDETOUNGAN                           *;  
%put *                          Courriel :  ariel92and@gmail.com                           *; 
%put *                  N hésitez pas à envoyer un mail en cas de besoin                   *; 
%put ***************************************************************************************;
%end;
title "";
options notes;
ods select all;
%mend;
