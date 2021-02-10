/*Je travaille ici sur des données simulées. J'ai deux séries de variables sur la production agricole. La première variable devrait répondre à un modèle SAR. La seconde à un modèle SEM
J'ai simulé les variables donc je savais a priori les modèles. J'ai ensuite 3 variables explicatives : Précipitation moyenne annuelle, Température  et Taxe. L'objectif est de 
modéliser la production agricole en fonction de la Précipitaiton, Température moyenne annuelle. Je représente d'abord les quatre variables sur la carte d'Afrique afin de mesurer
les différentes tendances*/
libname TP "D:\Sauvegardes\Résultats Statistiques\SAS\SAS Printemps 2017 Québec";
/*Une bibrairie pour les macros*/
options sasmstore=Mmacro mstored;
/*A cette étape il faut compiler les macros*/
goptions reset=all;
ods pdf style=ocean;
proc gmap data=TP.Base map=Maps.Africa;
format ProducAgric1 ProducAgric2 Precipitation Temperature Taxe comma10.1;
id ID;
choro ProducAgric1/coutline=red ; 
choro ProducAgric2/coutline=red ; 
choro Precipitation/coutline=red ; 
choro Temperature/coutline=red ; 
choro Taxe/coutline=red; 
run;
quit;

/*On peut noter une autocorrélation spatiale surtout pour la première variable ProducAgric1. Pour les autres variables, il n'y a pas de tendances générales. 
On va se référer au test de Moran pour confirmer. J'ai utilisé la Macro, Morantest*/
/*Dans Data on met la Base, dans Var la Variable et dans Poids la Matrice de poids*/
%Morantest(Data=TP.Base, Var=ProducAgric1, Poids=TP.Matrice);
%Morantest(Data=TP.Base, Var=ProducAgric2, Poids=TP.Matrice);
%Morantest(Data=TP.Base, Var=Precipitation, Poids=TP.Matrice);
%Morantest(Data=TP.Base, Var=Temperature, Poids=TP.Matrice);
%Morantest(Data=TP.Base, Var=Taxe, Poids=TP.Matrice);


/************************************************************************************Modèle de ProducAgric1*******************************************************************/
/*Nous allons nous comporter comme si nous ne savons pas qu'on va estimer un modèle SAR. Nous allons nous mettre dans un cadre réel. Nous allons estimer 
un modèle linéaire simple et tester s'il y a une mauvaise spécification*/
proc reg data=TP.Base;
model ProducAgric1= Precipitation Temperature Taxe;
run;

/*Tests de mauvaise spécification avec la macro Spatialtest. Dans modèle on met le modèle estimé*/
%Spatialtest(Data=TP.Base, Modele= ProducAgric1 ~ Precipitation + Temperature + Taxe, Poids=TP.Matrice);
/*Les tests ont montré une mauvaise spécification du modèle MCO. Nous allons estimer le modèle SAR avec la macro GSM. Dans la macro variable CORR on indique SAR*/
%GSM(Data=TP.Base, Modele=  ProducAgric1 ~ Precipitation + Temperature + Taxe, Poids=TP.Matrice, CORR="SAR");

/************************************************************************************Modèle de ProducAgric2*******************************************************************/
/*Reprenons le même principe sur la seconde variable*/
/*modèle linéaire simple*/
proc reg data=TP.Base;
model ProducAgric2= Precipitation Temperature Taxe;
run;

/*Tests de mauvaise spécification*/
%Spatialtest(Data=TP.Base, Modele= ProducAgric2 ~ Precipitation + Temperature + Taxe, Poids=TP.Matrice);
/*Les tests indiques un modèle SEM*/
%GSM(Data=TP.Base, Modele=  ProducAgric2 ~ Precipitation + Temperature + Taxe, Poids=TP.Matrice, CORR="SEM");
ods pdf close;
