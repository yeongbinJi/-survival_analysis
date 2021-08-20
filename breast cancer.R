#데이터
library(dynpred) 
breast=read.csv("nki.csv",header = T)
head(breast)
str(breast); 
y.tmp = rep(1,nrow(breast))
y.tmp[which((breast[,3]=="y"))] = 1
y.tmp[which((breast[,3]=="n"))] = 0
breast[,3] = y.tmp
breast$Posnodes=as.integer(breast$Posnodes)
breast2<-breast[order(breast$SampleID),] ; head(breast2)
data(nki);help(nki);head(nki)
breast3<-cbind(breast2[,c(3,4,5,6,9)],nki[,c(4,6,8,9,11)])
y.tmp2 = rep(1,nrow(breast3))
y.tmp2[which((breast3[,8]=="Yes"))] = 1
y.tmp2[which((breast3[,8]=="No"))] = 0
breast3[,8] = y.tmp2
y.tmp3 = rep(1,nrow(breast3))
y.tmp3[which((breast3[,9]=="Yes"))] = 1
y.tmp3[which((breast3[,9]=="No"))] = 0
breast3[,9] = y.tmp3
breast3$chemotherapy=as.integer(breast3$chemotherapy)
breast3$hormonaltherapy=as.integer(breast3$hormonaltherapy)
head(breast3);str(breast3)
attach(breast3)

#기초통계량
summary(breast3)
str(breast3)
#사망여부
table(EVENTdeath)
prop.table(table(EVENTdeath))
#원격 전이 유무 
table(EVENTmeta)
prop.table(table(EVENTmeta))
#림프절 상태 
table(Posnodes)
prop.table(table(Posnodes))
#에스트로겐 수용체유무
table(ESR1)
prop.table(table(ESR1))
#조직 등급
table(histolgrade)
prop.table(table(histolgrade))
#화학요법 
table(chemotherapy)
prop.table(table(chemotherapy))
#호르몬요법 
table(hormonaltherapy)
prop.table(table(hormonaltherapy))



#카플란 마이어 추정량: 생존함수의 추정량으로 가장 널리 사용되는 방법
#개개의 생존시간들은 서로 독립적이며 중도 절단은 생존시간과 독립이라고 가정  
kmfit<-survfit(Surv(TIMEsurvival, EVENTdeath)~1)
summary(kmfit)
print(kmfit,print.rmean = T)

#생존함수 신뢰구간
plot(kmfit,main="KM Curve", xlab="year", ylab="Surv Prob")



#누적위험함수 신뢰구간
plot(kmfit$time, -log(kmfit$surv),
     xlab = "year", ylab = "hazard rate", type="s",main="Cumulative Hazard")
lines(kmfit$time, -log(kmfit$upper), lty=3)
lines(kmfit$time, -log(kmfit$lower), lty=3)



#Nelson-Aalen 누적 위험률함수와 생존함수
par(mfrow=c(2,3))
#원격전이(EVENTmeta)에 따른 
metano<-breast3[EVENTmeta==0,]
metayes<-breast3[EVENTmeta==1,]
metano.na<-survfit(Surv(metano$TIMEsurvival,metano$EVENTdeath)~1,type="fleming")
metayes.na<-survfit(Surv(metayes$TIMEsurvival,metayes$EVENTdeath)~1,type="fleming")
plot(metano.na$time,-log(metano.na$surv),type="s",col=2, main="N-A EVENTmeta",xlab="year",ylab="hazard rate")
lines(metayes.na$time,-log(metayes.na$surv),type="s",lty=2,col=4)
legend("bottomright",c('원격전이 No : 0','원격전이 Yes : 1'),lty=1:2,col=c(2,4))
#림프구상태(Posnodes)에 따른 
posneg<-breast3[Posnodes==0,]
pospos<-breast3[Posnodes==1,]
posneg.na<-survfit(Surv(posneg$TIMEsurvival,posneg$EVENTdeath)~1,type="fleming")
pospos.na<-survfit(Surv(pospos$TIMEsurvival,pospos$EVENTdeath)~1,type="fleming")
plot(posneg.na$time,-log(posneg.na$surv),type="s",col=2, main="N-A Posnodes",xlab="year",ylab="hazard rate")
lines(pospos.na$time,-log(pospos.na$surv),type="s",lty=2,col=4)
legend("bottomright",c('림프절 negative : 0','림프절 positive : 1'),lty=1:2,col=c(2,4))
#에스트로겐 수용체에 따른 
estno<-breast3[ESR1==0,]
estyes<-breast3[ESR1==1,]
estno.na<-survfit(Surv(estno$TIMEsurvival,estno$EVENTdeath)~1,type="fleming")
estyes.na<-survfit(Surv(estyes$TIMEsurvival,estyes$EVENTdeath)~1,type="fleming")
plot(estno.na$time,-log(estno.na$surv),type="s",col=2, main="N-A ESR1",xlab="year",ylab="hazard rate")
lines(estyes.na$time,-log(estyes.na$surv),type="s",lty=2,col=4)
legend("bottomright",c('에스트로겐수용체 No : 0','에스트로겐수용체 Yes : 1'),lty=1:2,col=c(2,4))
#화학요법 유무에 따른 
cheno<-breast3[chemotherapy==0,]
cheyes<-breast3[chemotherapy==1,]
cheno.na<-survfit(Surv(cheno$TIMEsurvival,cheno$EVENTdeath)~1,type="fleming")
cheyes.na<-survfit(Surv(cheyes$TIMEsurvival,cheyes$EVENTdeath)~1,type="fleming")
plot(cheno.na$time,-log(cheno.na$surv),type="s",col=2, main="N-A chemotherapy",xlab="year",ylab="hazard rate")
lines(cheyes.na$time,-log(cheyes.na$surv),type="s",lty=2,col=4)
legend("bottomright",c('화학요법 No : 0','화학요법 Yes : 1'),lty=1:2,col=c(2,4))
#호르몬치료 유무에 따른 
horno<-breast3[hormonaltherapy==0,]
horyes<-breast3[hormonaltherapy==1,]
horno.na<-survfit(Surv(horno$TIMEsurvival,horno$EVENTdeath)~1,type="fleming")
horyes.na<-survfit(Surv(horyes$TIMEsurvival,horyes$EVENTdeath)~1,type="fleming")
plot(horno.na$time,-log(horno.na$surv),type="s",col=2, main="N-A hormonaltherapy",xlab="year",ylab="hazard rate")
lines(horyes.na$time,-log(horyes.na$surv),type="s",lty=2,col=4)
legend("bottomright",c('호르몬치료 No : 0','호르몬치료 Yes : 1'),lty=1:2,col=c(2,4))
#조직 등급에 따른 
hisinter<-breast3[histolgrade=="Intermediate",]
hispoor<-breast3[histolgrade=="Poorly diff",]
hiswell<-breast3[histolgrade=="Well diff",]
hisinter.na<-survfit(Surv(hisinter$TIMEsurvival,hisinter$EVENTdeath)~1,type="fleming")
hispoor.na<-survfit(Surv(hispoor$TIMEsurvival,hispoor$EVENTdeath)~1,type="fleming")
hiswell.na<-survfit(Surv(hiswell$TIMEsurvival,hiswell$EVENTdeath)~1,type="fleming")
plot(hisinter.na$time,-log(hisinter.na$surv),type="s",col=2, main="N-A histolgrade",xlab="year",ylab="hazard rate")
lines(hispoor.na$time,-log(hispoor.na$surv),type="s",lty=2,col=4)
lines(hiswell.na$time,-log(hiswell.na$surv),type="s",lty=3,col=3)
legend("bottomright",c('Intermediate','Poorly',"Well"),lty=1:3,col=c(2,4,3))







#로그순위검정
par(mfrow=c(2,3))
#원격전이 생존함수 동일성 검정
metafit<- survdiff(Surv(TIMEsurvival, EVENTdeath)~EVENTmeta, data=breast3)
metafit
plot(survfit(Surv(TIMEsurvival, EVENTdeath)~EVENTmeta),col=c(2,4),lty=1:2,
     main="EVENTmeta", xlab="year",ylab="Surv Prob" )
legend("bottomleft",c('원격전이 No : 0','원격전이 Yes : 1'),lty=1:2,col=c(2,4))
#림프구상태 생존함수 동일성 검정
posfit<- survdiff(Surv(TIMEsurvival, EVENTdeath)~Posnodes, data=breast3)
posfit
plot(survfit(Surv(TIMEsurvival, EVENTdeath)~Posnodes),col=c(2,4),lty=1:2,
     main="Posnodes", xlab="year",ylab="Surv Prob" )
legend("bottomleft",c('음성질환 : n','양성질환 : y'),lty=1:2,col=c(2,4))
#에스트로겐 수용체 생존함수 동일성 검정
esfit<-survdiff(Surv(TIMEsurvival, EVENTdeath)~ESR1)
esfit
plot(survfit(Surv(TIMEsurvival, EVENTdeath)~ESR1),col=c(2,4),lty=1:2,
     main="ESR1", xlab="year",ylab="Surv Prob" )
legend("bottomleft",c('에스트로겐수용체 No : 0','에스트로겐수용체 Yes : 1'),lty=1:2,col=c(2,4))
#화학요법
chefit<-survdiff(Surv(TIMEsurvival, EVENTdeath)~chemotherapy)
chefit
plot(survfit(Surv(TIMEsurvival, EVENTdeath)~chemotherapy),col=c(2,4),lty=1:2,
     main="chemotherapy", xlab="year",ylab="Surv Prob" )
legend("bottomleft",c('화학요법 No : 0','화학요법 Yes : 1'),lty=1:2,col=c(2,4))
#호르몬치료
horfit<-survdiff(Surv(TIMEsurvival, EVENTdeath)~hormonaltherapy)
horfit
plot(survfit(Surv(TIMEsurvival, EVENTdeath)~hormonaltherapy),col=c(2,4),lty=1:2,
     main="hormonaltherapy", xlab="year",ylab="Surv Prob" )
legend("bottomleft",c('호르몬치료 No : 0','호르몬치료 Yes : 1'),lty=1:2,col=c(2,4))
#조직 등급 생존함수 동일성 검정 
gfit<- survdiff(Surv(TIMEsurvival, EVENTdeath)~histolgrade)
gfit
plot(survfit(Surv(TIMEsurvival, EVENTdeath)~histolgrade),col=c(2,4,3),lty=1:3,
     main="histolgrade", xlab="year",ylab="Surv Prob" )
legend("bottomleft",c('Intermediate','Poorly',"Well"),lty=1:3,col=c(2,4,3))







#cox 비례위험모형
fitcox<-coxph(Surv(TIMEsurvival, EVENTdeath)~EVENTmeta+Posnodes+ESR1+age+diameter
              +chemotherapy+hormonaltherapy+histolgrade, data=breast3)
fitcox
summary(fitcox)
cox.zph(fitcox)
extractAIC(fitcox)
#cox 비례위험모형-단계적 모형 선택
library(MASS)
coxfit=stepAIC(fitcox,trace = F)
coxfit
summary(coxfit)
cox.zph(coxfit)
extractAIC(coxfit)


#종양의 지름 층화_최종모형
fitcox2<-coxph(Surv(TIMEsurvival, EVENTdeath)~EVENTmeta+Posnodes+ESR1+age+strata(diameter)
               +chemotherapy+hormonaltherapy+histolgrade, data=breast3)
fitcox2
summary(fitcox2)
cox.zph(fitcox2)
extractAIC(fitcox2)
#cox 비례위험모형-단계적 모형 선택
library(MASS)
coxfit2=stepAIC(fitcox2,trace = F)
coxfit2
summary(coxfit2)
cox.zph(coxfit2)
extractAIC(coxfit2)


#조직 등급 층화 
fitcox3<-coxph(Surv(TIMEsurvival, EVENTdeath)~EVENTmeta+Posnodes+ESR1+age+diameter
                 +chemotherapy+hormonaltherapy+strata(histolgrade), data=breast3)
fitcox3
summary(fitcox3)
cox.zph(fitcox3)
extractAIC(fitcox3)
#cox 비례위험모형-단계적 모형 선택
library(MASS)
coxfit3=stepAIC(fitcox3,trace = F)
coxfit3
summary(coxfit3)
cox.zph(coxfit3)
extractAIC(coxfit3)

#교호작용 
fitcox4<-coxph(Surv(TIMEsurvival, EVENTdeath)~EVENTmeta+Posnodes+ESR1+age+diameter
               +chemotherapy*hormonaltherapy+histolgrade, data=breast3)
fitcox4
summary(fitcox4)
cox.zph(fitcox4)
extractAIC(fitcox4)
#cox 비례위험모형-단계적 모형 선택
library(MASS)
coxfit4=stepAIC(fitcox4,trace = F)
coxfit4
summary(coxfit4)
cox.zph(coxfit4)
extractAIC(coxfit4)




#잔차모형진단
cox.zph(coxfit2)
par(mfrow=c(1,3))
#Schoenfeld 잔차 산점도
plot(cox.zph(coxfit2))

#공변량 없는 모형
par(mfrow=c(2,2))

fit0<-coxph(Surv(TIMEsurvival,EVENTdeath)~1)
res1<-resid(fit0,type='martingale')

plot(EVENTmeta,res1,xlab="EVENTmeta",ylab="martingale")
lines(lowess(EVENTmeta,res1,iter=0),lty=2)

plot(ESR1,res1,xlab="ESR1",ylab="martingale")
lines(lowess(ESR1,res1,iter=0),lty=2)

plot(chemotherapy,res1,xlab="chemotherapy",ylab="martingale")
lines(lowess(chemotherapy,res1,iter=0),lty=2)

histolgrade=as.integer(histolgrade)
plot(histolgrade,res1,xlab="histolgrade",ylab="martingale")
lines(lowess(histolgrade,res1,iter=0),lty=2)


#공변량 있는 모형
par(mfrow=c(2,2))

res2<-resid(coxph(Surv(TIMEsurvival,EVENTdeath)~EVENTmeta+ESR1+chemotherapy+histolgrade+strata(diameter)),type='martingale')
plot(EVENTmeta,res2,xlab="EVENTmeta",ylab="martingale")
lines(lowess(EVENTmeta,res2,iter=0),lty=2)

plot(ESR1,res2,xlab="ESR1",ylab="martingale")
lines(lowess(ESR1,res2,iter=0),lty=2)

plot(chemotherapy,res2,xlab="age",ylab="martingale")
lines(lowess(chemotherapy,res2,iter=0),lty=2)

histolgrade=as.integer(histolgrade)
plot(histolgrade,res2,xlab="age",ylab="martingale")
lines(lowess(histolgrade,res2,iter=0),lty=2)


