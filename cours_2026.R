#On regarde ce qu'il y a en mémoire et on le supprime
ls()
rm(list=ls())

#répértoire courant
setwd("/home/melen/Dropbox/cours_agro")

#Packages R 
library(multcomp)
library(multcompView)
library(emmeans)

#on charge les données
data<-read.csv("lesion_asco.csv",stringsAsFactors = T)

#on regarde le jeu de données et on vérifie qu'il est bien chargé et que les variables sont dans le bon format
summary(data)

# temps : temps post-inoculation
# genotype : genotype de pois
# espece : espèce inoculée
# sfeuille : surface de la feuille en cm2
# sain : surface saine
# lesion : surface nécrosée


# nous avons la surface de la lésion au cours du temps. Si la lésion est un cercle on a le rayon=sqrt(surface/pi)
data$rayon<-sqrt(data$lesion/pi)

### sous sélectionner un data.frame

# Une façon de sélectionner un sous partie d'un tableau de données
# rappel pour un data.frame d, on peut choisir les lignes et les colonnes à sélectionner
# d[lignes,colonnes]
#ex
data[1,] # première ligne
data[3,4] # première ligne 4ème colonne
data[1:10,c(3,5)] # lignes de 1 à 10 x colonnes 3 et 5
data[c(1,10),3:5] #
##

# que les lignes où l'espèce est dp
data[data$espece=="pm",]
#que les lignes avec pm et sensible
data[data$espece=="pm" & data$genotype=="sensible",]
#que les lignes avec pm et sensible et le temps > 3
data[data$espece=="pm" & data$genotype=="sensible" & data$temps>3,]

nouveau<-data[data$espece=="dp" & data$genotype=="resistant" & data$temps<7,]
nouveau
######################################################################
######### Analyse descriptive/graphique  ########################################
####################################################################

#On visualise les données 
#regarder ?plot pour les arguments de la fonction plot

par(mfrow=c(1,2))
plot(data$lesion~data$temps,pch=20,xlab="Temps (jour)",ylab="Surface de la lésion(cm²)",xlim=c(0,10))
plot(data$rayon~data$temps,pch=20,xlab="Temps (jour)",ylab="Rayon de la lésion (cm)",xlim=c(0,10))


#coplot surface
coplot(data$lesion~data$temps|data$espece)
coplot(data$lesion~data$temps|data$genotype)
plot(data$lesion~data$temps,col=as.factor(data$espece),pch=20)
plot(data$lesion~data$temps,col=as.factor(data$genotype),pch=20)

#coplot rayon
coplot(data$rayon~data$temps|data$espece)
coplot(data$rayon~data$temps|data$genotype)
plot(data$rayon~data$temps,col=as.factor(data$espece),pch=20)
plot(data$rayon~data$temps,col=as.factor(data$genotype),pch=20)

# Si la vitesse de progression est constante v=dr/dt=cste alors on a une relation linéaire entre le rayon et le temps
# on analyse donc l'évolution du rayon de la lésion en fonction du temps

#comparons les espèces
with(data[data$espece=="dp",],plot(rayon~temps,pch=20,xlim=c(0,10)))
with(data[data$espece=="pm",],plot(rayon~temps,pch=20,xlim=c(0,10)))

#ou 
plot(rayon~temps,pch=20,xlim=c(0,10),data=data[data$espece=="dp",],col="blue")
points(rayon~temps,pch=20,xlim=c(0,10),data=data[data$espece=="pm",],col="red")
legend("topleft",legend=c("Pinodes","Phoma"),pch=20,bty='n',col=c("blue","red"))

#comparons les variétés (= génotypes)
with(data[data$genotype=="resistant",],plot(rayon~temps,pch=20,xlim=c(0,10)))
with(data[data$genotype=="sensible",],plot(rayon~temps,pch=20,xlim=c(0,10)))

plot(rayon~temps,pch=20,xlim=c(0,10),data[data$genotype=="sensible",],col="blue")
points(rayon~temps,pch=20,xlim=c(0,10),data[data$genotype=="resistant",],col="red")
legend("topleft",legend=c("sensible","résistant"),pch=20,bty='n',col=c("blue","red"))

#Pour chaque cas genotye*espece
par(mfrow=c(2,2))
with(data[data$espece=="dp" & data$genotype=="resistant",],plot(rayon~temps,pch=20,xlim=c(0,10),ylim=c(0,1.4),main="dp résistant"))
with(data[data$espece=="pm" & data$genotype=="resistant",],plot(rayon~temps,pch=20,xlim=c(0,10),ylim=c(0,1.4),main="pm résistant"))
with(data[data$espece=="dp" & data$genotype=="sensible",],plot(rayon~temps,pch=20,xlim=c(0,10),ylim=c(0,1.4),main="dp sensible"))
with(data[data$espece=="pm" & data$genotype=="sensible",],plot(rayon~temps,pch=20,xlim=c(0,10),ylim=c(0,1.4),main="pm sensible"))
dev.off()#pour réinitialiser les paramètres graphiques

#######################################################################################################################################
######################################################################################################################################
# On veut ajuster un modèle linéaire qui décrit le développement moyen du rayon de la lésion en fonction du temps
# regarder ?lm pour les arguments de la fonction lm
#On a un modèle linéaire du type rayon= intercept + a*temps
#lm nous permet (entre autres) d'estimer les paramètres intercept et a à partir des données

#En considérant toutes les données

mod<-lm(rayon~temps,data=data)
mod #coefficients estimés --> rayon= -0.1033 + 0.1165 * jours
anova(mod) #analyse de la variance, test des effets globaux
summary(mod) #coef estimés, R² etc
qqnorm(residuals(mod)) # des résidus normalement distribués?
plot(residuals(mod)~data$temps,pch=20) # écart modèle/donné (=résidus)
abline(h=0)

plot(predict(mod,data.frame(temps=seq(0:11)))~seq(0:11),type='l',lwd=2,col='red',ylim=c(0,1.5),xlab="Temps (jours)",ylab="Rayon de la lésion (cm)")
points(data$rayon~data$temps,pch=20)

# ou 
plot(data$rayon~data$temps,pch=20,ylim=c(0,1.5),xlab="Temps (jours)",ylab="Rayon de la lésion (cm)")
abline(a=-0.1033,b=0.1165,col='blue') # pour ajouter une droite ?abline


#Pour chaque espece
#dp

mod_dp<-lm(rayon~temps,data=data[data$espece=="dp",])
mod_dp
summary(mod_dp)
qqnorm(residuals(mod_dp))
plot(residuals(mod_dp)~data[data$espece=="dp",]$temps,pch=20)
abline(h=0)

plot(predict(mod_dp,data.frame(temps=seq(0:11)))~seq(0:11),type='l',lwd=2,col='red',ylim=c(0,1.5),
     xlab="Temps (jours)",ylab="Rayon de la lésion (cm)")
with(data[data$espece=="dp",],points(rayon~temps,pch=20))

# ou
plot(rayon~temps,pch=20,ylim=c(0,1.5),xlab="Temps (jours)",ylab="Rayon de la lésion (cm)",data=data[data$espece=='dp',])
abline(a=coef(mod_dp)[1],b=coef(mod_dp)[2],lwd=2,col='red')


#pm

mod_pm<-lm(rayon~temps,data=data[data$espece=="pm",])
mod_pm
summary(mod_pm)
qqnorm(residuals(mod_pm))

plot(predict(mod_pm,data.frame(temps=seq(0:11)))~seq(0:11),type='l',lwd=2,col='red',ylim=c(0,1.5),
     xlab="Temps (jours)",ylab="Rayon de la lésion (cm)")
with(data[data$espece=="pm",],points(rayon~temps,pch=20))

#comparons dp et pm

plot(predict(mod_dp,data.frame(temps=seq(0:11)))~seq(0:11),type='l',lwd=2,col='red',
     ylim=c(0,1.5),xlab="Temps (jours)",ylab="Rayon de la lésion (cm)")
lines(predict(mod_pm,data.frame(temps=seq(0:11)))~seq(0:11),type='l',lwd=2,col='orange') #lines() pour ajouter une courbe au graphique
legend("topleft",legend=c("dp","pm"),col=c("red","orange"),lwd=2)#ajout de la légende --> ?legend 

# ou, autre façon...

plot(rayon~temps,pch=5,ylim=c(0,1.5),xlab="Temps (jours)",ylab="Rayon de la lésion (cm)",data=data[data$espece=='dp',])
abline(a=coef(mod_dp)[1],b=coef(mod_dp)[2],lwd=2,col='red')
points(rayon~temps,pch=6,ylim=c(0,1.5),xlab="Temps (jours)",ylab="Rayon de la lésion (cm)",col='gray',data=data[data$espece=='pm',])
abline(a=coef(mod_pm)[1],b=coef(mod_pm)[2],lwd=2,col='orange')

########################################################################################################
# on observe une différence entre les espèces, y a t-il un effet significatif?
# Pour tester cela on réalise une Analyse de covariance ANCOVA 

# Quand on regarde "l'effet" d'une variable sur une droite (i.e. on compare deux droites) quelles peuvent être les différences ?

# des pentes différentes
plot(NULL, type="n", xlab="x", ylab="y", xlim=c(0, 10), ylim=c(0, 10))
abline(a=2,b=0.5,col='red')
abline(a=2,b=1,col='orange')

# des ordonnées à l'origine différentes
plot(NULL, type="n", xlab="x", ylab="y", xlim=c(0, 10), ylim=c(0, 10))
abline(a=2,b=0.5,col='red')
abline(a=0.5,b=0.5,col='orange')

# des ordonnées à l'origine et des pentes différentes
plot(NULL, type="n", xlab="x", ylab="y", xlim=c(0, 10), ylim=c(0, 10))
abline(a=2,b=0.5,col='red')
abline(a=0.5,b=1,col='orange')


# Revenons au problème initial.
# Pour cela on ajoute l'effet d'une variable explicative qualitative (effet espece à deux modalités: dp pm) au modèle de régression
# On passe de rayon=intercept + a*time à rayon =intercept*effet souche + a*effet souche*time 
# on peut avoir des effets espèces significatifs sur la pente et/ou sur la constante (=ordonnée à l'origine)

modele<-lm(rayon~temps+espece+temps:espece,data=data)# modèle pour ancova (rayon~temps+espece+temps:espece = rayon~temps*espece)
dummy.coef(modele) #affiche tous les coefs
anova(modele) #analyse de variance qui nous donne une p-value
summary(modele) 
qqnorm(residuals(modele))

#variance expliquée

anova(modele)$Sum/(sum(anova(modele)$Sum))
cbind(variable=row.names(anova(modele)),influence=anova(modele)$Sum/(sum(anova(modele)$Sum)))

# tracer les droites pour dp et pm à partir des coefs ??
plot(NULL, type="n", xlab="x", ylab="y", xlim=c(0, 10), ylim=c(0, 2))
#dp
abline(a=-0.119+0,b=0.1337+0,col='red')
#pm
abline(a=-0.119+0.0153,b=0.1337-0.035268,col='orange')

# comparer les pentes pour chaque modalité avec emtrends et cld 

#cld(emmeans(modele,"espece"))

emtrends(modele,"espece",var="temps") # 
cld(emtrends(modele,"espece",var="temps")) # Tukey différences de pentes significatives pour les deux espèces
plot(cld(emtrends(modele,"espece",var="temps")))


#################################################################"
###### Effet genotype?

m<-lm(rayon~temps*genotype,data=data)
dummy.coef(m)
anova(m)
summary(m)

cld(emtrends(m,"genotype",var="temps")) 
plot(cld(emtrends(m,"genotype",var="temps")))

#################################################################"
###### Interaction genotype espece?

m<-lm(rayon~temps*genotype*espece,data=data)
dummy.coef(m)
anova(m)
summary(m)
dummy.coef(m)

# 
cld(emtrends(m,~genotype|espece,var="temps")) 
plot(cld(emtrends(m,~genotype|espece,var="temps")))

cld(emtrends(m,~espece|genotype,var="temps")) 
plot(cld(emtrends(m,~espece|genotype,var="temps")))

