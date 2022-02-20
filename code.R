###### SET UP #####
#Directory set up
setwd("C:/Users/adria/Dropbox/UIEM/LEAD/Proyectos/ENSANUT")

#Package load
pacman::p_load(dplyr,tidyr,ggstatsplot,readxl,tableone,car,easystats,
               patchwork,MASS,see,qqplotr,bootStepAIC,performance,
               rpart,rpart.plot,gtools,broom,lmtest,visdat,parameters,
               ggcharts,conflicted,RNHANES,rattle,mlogit,MLmetrics,transplantr,
               survey,haven,survey)

# Solving duplicate functions conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter", "dplyr")

###### DATA WRANGNLING #####
#Data upload
bioq <- read_excel("detBioCronicosAdultos.xlsx")
bioq$UPM <- as.character(bioq$UPM)
bioq$VIV_SEL <- as.character(bioq$VIV_SEL)
bioq$HOGAR <- as.character(bioq$HOGAR)
bioq$NUMREN <- as.character(bioq$NUMREN)

antro <- read_excel("CN_ANTROPOMETRIA.xlsx")
antro$UPM <- as.character(antro$UPM)
antro$VIV_SEL <- as.character(antro$VIV_SEL)
antro$HOGAR <- as.character(antro$HOGAR)
antro$NUMREN <- as.character(antro$NUMREN)

#Dataset meriging
antro <- antro %>% select(UPM,VIV_SEL,HOGAR,NUMREN,PESO1_2,TALLA4_2,P27_2_1,
                          CIRCUNFERENCIA8_2,CINTURA21_2,P27_2_2) %>% 
                  mutate(id=interaction(UPM,VIV_SEL,HOGAR,NUMREN)) %>%               
                  mutate(imc=(PESO1_2/(TALLA4_2/100)^2)) %>% 
                  mutate(imcq  = quantcut(imc,q=5)) %>% 
                  mutate(ica=CINTURA21_2/TALLA4_2) %>% 
                  mutate(icc=CINTURA21_2/CIRCUNFERENCIA8_2) %>% 
                  rename(psis=P27_2_1,pdias=P27_2_2,peso=PESO1_2,talla=TALLA4_2,
                         circcintura=CINTURA21_2,circcadera=CIRCUNFERENCIA8_2)

bioq <- bioq %>% select(UPM,VIV_SEL,HOGAR,NUMREN,EDAD,SEXO,VALOR_GLU_SUERO,
                        VALOR_HB1AC,
                        VALOR_ALBUM,VALOR_COL_HDL,VALOR_COL_LDL,VALOR_COLEST,
                        VALOR_CREAT,VALOR_INSULINA,VALOR_TRIG) %>% 
  mutate(id=interaction(UPM,VIV_SEL,HOGAR,NUMREN)) %>%
  mutate(tgtohdl=VALOR_TRIG/VALOR_COL_HDL) %>% 
  mutate(nonhdl=VALOR_COLEST-VALOR_COL_HDL) %>% 
  mutate(ldlsmpsn= (VALOR_COLEST/0.948)-(VALOR_COL_HDL/0.971)-((VALOR_TRIG/8.56)+
        ((VALOR_TRIG*nonhdl)/2140)-(VALOR_TRIG/16100))-9.44) %>% 
  mutate(rmnnt=VALOR_COLEST-ldlsmpsn-VALOR_COL_HDL) %>% 
  mutate(rmnntq  = quantcut(rmnnt,q=5)) %>% 
  mutate(homair=VALOR_INSULINA*(VALOR_GLU_SUERO/405)) %>% 
  mutate(etnia="non-black") %>% 
  rename(edad=EDAD,sexo=SEXO,glu=VALOR_GLU_SUERO,
         hba1c=VALOR_HB1AC,alb=VALOR_ALBUM,hdl=VALOR_COL_HDL,
         ldl=VALOR_COL_LDL,colest=VALOR_COLEST,cr=VALOR_CREAT,
         ins=VALOR_INSULINA,trig=VALOR_TRIG) %>% 
  mutate(tfg = ckd_epi(creat = cr, age = edad, sex = sexo,
                        eth = etnia, units = "US")) %>% 
  mutate(across(where(is.numeric), round, 2)) %>%
  filter(between(glu,40,1000),between(alb,1,6),between(hdl,14,101),
         between(cr,0.2,10),between(ins,0,50),between(ldlsmpsn,20,1000),
         between(homair,0,15)) %>% 
  select(!ldl)%>% 
  na.omit()
  
#Merging third (INSP) dataset
data <- merge(bioq,antro,by="id")
data2 <- read_dta("DATA.dta")
data2$UPM <- as.character(data2$code_upm)
data2$VIV_SEL <- as.character(data2$viv_sel)
data2$HOGAR <- as.character(data2$hogar)
data2$NUMREN <- as.character(data2$numren)
data2 <- data2 %>% mutate(id=interaction(UPM,VIV_SEL,HOGAR,NUMREN)) 
data2 <- data2 %>% select(UPM,VIV_SEL,HOGAR,NUMREN,id,dx_infarto,estrato,
                          sexo,edad,entidad,nse3,ecv,alc2,nhba1c,
                          nhdl,nldl,ncolest,ninsulina,nse3,imc_valid,fuma3,
                          encuesta,alc2,escolaridad,tas,tad,fuma3,ntrig,
                          padres_dia,padres_hta,h_ayuno,ngluc,talla5) %>% 
                    mutate(dx_infarto=ifelse(dx_infarto==1,"1","0"))
data <- merge(data,data2,by="id")

#FINAL DATASET
write.csv(data,"merged_and_filtered_data.csv")

##### DEFINING DISEASES IN THE DATASET #####
data <- data %>% select(UPM,VIV_SEL,HOGAR,NUMREN,id,edad.x,sexo.x,glu,hba1c,alb,
                        hdl,colest,cr,ldlsmpsn,ins,trig,tgtohdl,homair,tfg,psis,
                        pdias,imc,rmnnt,imcq,rmnntq,dx_infarto,estrato,
                        entidad,ecv,nse3,alc2,fuma3,nhba1c,
                        nhdl,nldl,ncolest,ninsulina,imc_valid,fuma3,encuesta,
                        alc2,ntrig) %>% 
                  mutate(dm=hba1c>=6.5) %>% 
                  mutate(erc=tfg<=90) %>% 
                  mutate(dlp=ldlsmpsn>=190|trig>=150|hdl<=40) %>% 
                  mutate(has=psis>=140|pdias>90) %>% 
                  mutate(obes=imc>=30) %>% 
                  mutate(sexo = ifelse(sexo.x=="M","Mujer","Hombre")) %>% 
                  mutate(across(where(is.numeric), round, 2)) %>% 
                  mutate(dx_infarto=as.factor(dx_infarto))
write.csv(data,"disease_definitions.csv")

##### ENCODING ENSANUT'S STRATIFIED COMPLEX DESIGN #####
design <- svydesign(id =  ~ UPM, nest = TRUE, survey.lonely.psu = "adjust", 
                    data = data)

###### DESCRIPTIVE STATS #####
tab1 <- CreateCatTable(data=data,
               vars = c("sexo","dm","dlp","has","obes","nse3","alc2","fuma3",
                        "dm","erc","dlp","has","obes","dx_infarto"))

tab2 <- CreateContTable(data=data,
                       vars = c("edad.x","glu","hba1c","alb","hdl","colest","cr",
                                "ldlsmpsn","ins","trig","tgtohdl","homair",
                                "tfg","psis","pdias","imc","rmnnt"))

### Punctual estimates accounting for ENSANUT'S complex design
### Cuantitative variables
edad <- svymean(~edad.x, design, na=TRUE, ci=TRUE)
edad_ci <- confint(svymean(~edad.x, design))
edad
edad_ci

glu <- svymean(~glu, design, na=TRUE, ci=TRUE)
glu_ci <- confint(svymean(~glu, design))
glu
glu_ci

hba1c <- svymean(~hba1c, design, na=TRUE, ci=TRUE)
hba1c_ci <- confint(svymean(~hba1c, design))
hba1c
hba1c_ci

hdl <- svymean(~hdl, design, na=TRUE, ci=TRUE)
hdl_ci <- confint(svymean(~hdl, design))
hdl
hdl_ci

colest <- svymean(~colest, design, na=TRUE, ci=TRUE)
colest_ci <- confint(svymean(~colest, design))
colest
colest_ci

ldl <- svymean(~ldlsmpsn, design, na=TRUE, ci=TRUE)
ldl_ci <- confint(svymean(~ldlsmpsn, design))
ldl
ldl_ci

trig <- svymean(~trig, design, na=TRUE, ci=TRUE)
trig_ci <- confint(svymean(~trig, design))
trig
trig_ci

imc <- svymean(~imc, design, na=TRUE, ci=TRUE)
imc_ci <- confint(svymean(~imc, design))
imc
imc_ci

tfg <- svymean(~tfg, design, na=TRUE, ci=TRUE)
tfg_ci <- confint(svymean(~tfg, design))
tfg
tfg_ci

tgtohdl <- svymean(~tgtohdl, design, na=TRUE, ci=TRUE)
tgtohdl_ci <- confint(svymean(~tgtohdl, design))
tgtohdl
tgtohdl_ci

rmnnt <- svymean(~rmnnt, design, na=TRUE, ci=TRUE)
rmnnt_ci <- confint(svymean(~rmnnt, design))
rmnnt
rmnnt_ci

homair <- svymean(~homair, design, na=TRUE, ci=TRUE)
homair_ci <- confint(svymean(~homair, design))
homair
homair_ci


### Categorical variables
sexo <- svyciprop(~I(sexo=="Mujer"), design, method="me",na.rm=TRUE)
sexo
dm <- svyciprop(~I(dm=="TRUE"), design, method="me",na.rm=TRUE)
dm
has <- svyciprop(~I(has=="TRUE"), design, method="me",na.rm=TRUE)
has
evc <- svyciprop(~I(ecv=="1"), design, method="me",na.rm=TRUE)
evc
iam <- svyciprop(~I(dx_infarto=="1"), design, method="me",na.rm=TRUE)
iam
erc <- svyciprop(~I(erc=="TRUE"), design, method="me",na.rm=TRUE)
erc
obes <- svyciprop(~I(obes=="TRUE"), design, method="me",na.rm=TRUE)
obes
dlp <- svyciprop(~I(dlp=="TRUE"), diamesign, method="me",na.rm=TRUE)
dlp
nse1 <- svyciprop(~I(nse3=="1"), design, method="me",na.rm=TRUE)
nse1
nse2 <- svyciprop(~I(nse3=="2"), design, method="me",na.rm=TRUE)
nse2
nse3 <- svyciprop(~I(nse3=="3"), design, method="me",na.rm=TRUE)
nse3
alchol <- svyciprop(~I(alc2=="1"), design, method="me",na.rm=TRUE)
alchol
tabaco <- svyciprop(~I(fuma3=="1"), design, method="me",na.rm=TRUE)
tabaco
edad <- svyquantile(~edad.x, design, c(.25,.5,.75),ci=TRUE)
edad

###### "VALIDATING" AGAINST NHANES #####
# Comparing LMHRs vs NHANES
# Getting NHANES' data
## LDL (in mg/dL) is encoded as "LBLDL" and Triglycerides (in mg/dL) as "LBXTR"
nhanes1 <-  nhanes_load_data("TRIGLY_G", "2011-2012")
## HDL (in mg/dL) is encoded as "LBDHDD"
nhanes2 <- nhanes_load_data("HDL_G","2011-2012")
## BMI (kg/m2) is encoded as "BMXBMI"
nhanes3 <- nhanes_load_data("BMX_G","2011-2012")
## AGE (in years) is encoded as "RIDAGEYR"
nhanes4 <- nhanes_load_data("DEMO_G","2011-2012")
## TCHOL (in years) is encoded as "RIDAGEYR"
nhanes5 <- nhanes_load_data("TCHOL_G","2011-2012")

### Data years and encoding can be corroborated at: 
### https://wwwn.cdc.gov/nchs/nhanes/search/default.aspx

# Merging NHANES' data sets and excluding children
nhanes_lipids <- merge(nhanes1,nhanes2,by="SEQN")
nhanes_lipids <- merge(nhanes_lipids,nhanes5,by="SEQN")
nhanes_demo <- merge(nhanes3,nhanes4,by="SEQN") %>% filter(RIDAGEYR>17)
nhanes_data <- merge(nhanes_lipids,nhanes_demo,by="SEQN") %>% 
#Calculating remnants in NHANES' IV
            mutate(rmnnt=LBXTC-LBDLDL-LBDHDD)

### Visual Comparisons

## BMI
nhanes_imc <- ggplot(data=data,mapping = aes(x=imc))+
  geom_density(inherit.aes = TRUE,color="#006341", fill="#006341",alpha=0.3)+
  geom_density(inherit.aes = FALSE, data=nhanes_data,mapping = (aes(x=BMXBMI))
               ,color="blue3",fill="blue3",alpha=0.3)+
  ggtitle("BMI in ENSANUT-18 vs NHANES-IV")+
  scale_y_continuous(name="Density",breaks=waiver())+
  scale_x_continuous(name="kg/m2",breaks=waiver(),limits = c(10,60))+
  theme(plot.title = element_text(size = 20))+
  theme_classic()
ggsave("nhanes_imc.tiff", units="cm", width=20, height=20, dpi=500, 
       compression = 'lzw')
nhanes_imc

## Cholesterol
nhanes_chol <- ggplot(data=data,mapping = aes(x=colest))+
  geom_density(inherit.aes = TRUE,color="#006341", fill="#006341",alpha=0.3)+
  geom_density(inherit.aes = FALSE, data=nhanes_data,mapping = (aes(x=LBXTC))
               ,color="blue3",fill="blue3",alpha=0.3)+
  ggtitle("Total Cholesterol in ENSANUT-18 vs NHANES-IV")+
  scale_y_continuous(name="Density",breaks=waiver())+
  scale_x_continuous(name="mg/dL",breaks=waiver(),limits = c(20,400))+
  theme(plot.title = element_text(size = 20))+
  theme_classic()
ggsave("nhanes_colest.tiff", units="cm", width=20, height=20, dpi=500, 
       compression = 'lzw')
nhanes_chol

## LDL
nhanes_ldl <- ggplot(data=data,mapping = aes(x=ldlsmpsn))+
  geom_density(inherit.aes = TRUE,color="#006341", fill="#006341",alpha=0.3)+
  geom_density(inherit.aes = FALSE, data=nhanes_data,mapping = (aes(x=LBDLDL))
               ,color="blue3",fill="blue3",alpha=0.3)+
  ggtitle("LDL in ENSANUT-18 vs NHANES-IV")+
  scale_y_continuous(name="Density",breaks=waiver())+
  scale_x_continuous(name="mg/dL",breaks=waiver(),limits = c(20,300))+
  theme(plot.title = element_text(size = 20))+
  theme_classic()
ggsave("nhanes_ldl.tiff", units="cm", width=20, height=20, dpi=500, 
       compression = 'lzw')
nhanes_ldl

## Remnants
nhanes_rmnnt <- ggplot(data=data,mapping = aes(x=rmnnt))+
  geom_density(inherit.aes = TRUE,color="#006341", fill="#006341",alpha=0.3)+
  geom_density(inherit.aes = FALSE, data=nhanes_data,mapping = (aes(x=rmnnt))
               ,color="blue3",fill="blue3",alpha=0.3)+
  ggtitle("RLP in ENSANUT-18 vs NHANES-IV")+
  scale_y_continuous(name="Density",breaks=waiver())+
  scale_x_continuous(name="mg/dL",breaks=waiver(),limits = c(20,300))+
  theme(plot.title = element_text(size = 20))+
  theme_classic()
ggsave("nhanes_rmnnts.tiff", units="cm", width=20, height=20, dpi=500, 
       compression = 'lzw')
nhanes_rmnnt

nhanes_imc+nhanes_chol+nhanes_ldl+nhanes_rmnnt

###### LIPID MARKERS HEAD-TO-HEAD COMPARISONS #####
### IAM
iam_lip <- glm(dx_infarto ~ hdl+ldlsmpsn+trig+rmnnt,
         data= data, family = "binomial")
summary(iam_lip)
confint(iam_lip)

iam_boot <- boot.stepAIC(iam_lip,data,B=100, )
iam_boot

### Stroke
evc_lip <- glm(ecv ~ hdl+ldlsmpsn+trig+rmnnt,
               data= data, family = "binomial")
summary(evc_lip)
confint(evc_lip)

evc_boot <- boot.stepAIC(evc_lip,data,B=100, )
evc_boot

### Diabetes
dm_lip <- glm(dm ~ hdl+ldlsmpsn+trig+rmnnt,
               data= data, family = "binomial")
summary(dm_lip)
confint(dm_lip)

dm_boot <- boot.stepAIC(dm_lip,data,B=100, )
dm_boot

### Hypertension
has_lip <- glm(has ~ hdl+ldlsmpsn+trig+rmnnt,
              data= data, family = "binomial")
summary(has_lip)
confint(has_lip)

has_boot <- boot.stepAIC(has_lip,data,B=100, )
has_boot

### Chronic Kidney Disease
erc_lip <- glm(erc ~ hdl+ldlsmpsn+trig+rmnnt,
               data= data, family = "binomial")
summary(erc_lip)
confint(erc_lip)

erc_boot <- boot.stepAIC(erc_lip,data,B=100, )
erc_boot

### Obesity
obes_lip <- glm(obes ~ hdl+ldlsmpsn+trig+rmnnt,
               data= data, family = "binomial")
summary(obes_lip)
confint(obes_lip)

obes_boot <- boot.stepAIC(obes_lip,data,B=100, )
obes_boot

###### LOGISTIC REGRESSIONS #####
iam_lr <- glm(dx_infarto~rmnntq+sexo+edad.x,data=data,family = "binomial")
summary(iam_lr)
confint(iam_lr)
iam_betas <- ggcoefstats(iam_lr,output = "plot",exclude.intercept = T,
                         effsize = "omega", 
                         title = "Myocardial infarction")
iam_betas

evc_lr <- glm(ecv~rmnntq+sexo+edad.x,data=data,family = "binomial")
summary(evc_lr)
confint(evc_lr)
evc_betas <- ggcoefstats(evc_lr,output = "plot",exclude.intercept = T,
                         effsize = "omega", 
                         title = "Stroke")
evc_betas

dm_lr <- glm(dm~rmnntq+sexo+edad.x,data=data,family = "binomial")
summary(dm_lr)
confint(dm_lr)
dm_betas <- ggcoefstats(dm_lr,output = "plot",exclude.intercept = T,
                        effsize = "omega", 
                        title = "Diabetes")
dm_betas

has_lr <- glm(has~rmnntq+sexo+edad.x,data=data,family = "binomial")
summary(has_lr)
confint(has_lr)
has_betas <- ggcoefstats(has_lr,output = "plot",exclude.intercept = T,
                         effsize = "omega", 
                         title = "Hypertension")
has_betas

erc_lr <- glm(erc~rmnntq+sexo+edad.x,data=data,family = "binomial")
summary(erc_lr)
confint(erc_lr)
erc_betas <- ggcoefstats(erc_lr,output = "plot",exclude.intercept = T,
                         effsize = "omega", 
                         title = "Chronic Kidney Disease")
erc_betas

obes_lr <- glm(obes~rmnntq+sexo+edad.x,data=data,family = "binomial")
summary(obes_lr)
confint(obes_lr)
obes_betas <- ggcoefstats(obes_lr,output = "plot",exclude.intercept = T,
                          effsize = "omega", 
                          title = "Obesity")
obes_betas
