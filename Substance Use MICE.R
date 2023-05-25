###Multiple Imputation for Weekly IPV, Substance Use and Childhood Maltreatment Logs###
install.packages("devtools")
library(devtools)
devtools::install_github("kkleinke/countimp", force = TRUE)
install.packages("VIM")
install.packages("sjmisc")
install.packages("AER")
install.packages('TMB', type = 'source')

#load libraries
library(countimp)
library(foreign)
library(mice)
library(VIM)
library(data.table)
library(sjmisc)
library(AER)
        
              
## Substance Use Data ##
weeklySU <- as.data.frame(read.spss("Z:\\R21 Clinical Trial\\data\\Weekly Log Data\\SU\\SU Week Sums.sav"))

#removing ppts with only baseline data
weeklySUclean <- weeklySU[-c(2,10,12,16,27,28,29,30,35,36,39,41,43,48,49,54,60,61,70),]

#removing the totals column
weeklySUclean <- weeklySU[-21]

# Remove row based on condition
summary(weeklySUclean) #53 rows with at least 2 weekly log responses

#reshape to long
weeklySUclean_long <- melt(setDT(weeklySUclean), id.vars = c("INT","PID"), variable.name = "week")

#pattern of missingness
md.pattern(weeklySUclean)

aggr_plot <- aggr(weeklySUclean, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#linear regression just to see the fit of the data compared to imputed
fitols <- lm(value ~ week + INT, data=weeklySUclean_long)
summary(fitols)

#checking for overdispersion
rd <- glm(value ~ week + INT, data = weeklySUclean_long, family = poisson)
dispersiontest(rd,trafo=1)
#seems to be evidence of overdispersion as c = 2.83 and p < .001, so zero-inflated negative binomial probabily better than poisson

## poisson imputation
imp1 <- countimp(weeklySUclean_long, meth = c(" ", " "," ", "pois"))
summary(imp1) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp1$imp$value) 

weeklySUimps <- mice::complete(imp1, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklySUmids <-as.mids(weeklySUimps)

#lm with pooled data
fitimp <- with(weeklySUmids,
               lm(value ~ week + INT))

summary(pool(fitimp)) #the p-values look very similar so imputation is successful

#merging imputations
pois_imp <- merge_imputations(
  weeklySUclean_long,
  weeklySUmids,
  ori = weeklySUclean_long,
  summary = c("den"),
  filter = NULL
)

#imputation fit looks to be okay/not the best

#zero inflated poisson
imp2 <- countimp(weeklySUclean_long, meth = c(" ", " "," ", "zip")) #the spaces tell us to skip the first 3 columns, we only want the 4th column imputed
summary(imp2) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp2$imp$value) 

weeklySUimps_zip <- mice::complete(imp2, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklySUmids_zip <-as.mids(weeklySUimps_zip)

#lm with pooled data
fitimp_zip <- with(weeklySUmids_zip,
               lm(value ~ week + INT))

summary(pool(fitimp_zip)) #the p-values look very similar so imputation is successful

#merging imputations
zip_imp <- merge_imputations(
  weeklySUclean_long,
  weeklySUmids_zip,
  ori = weeklySUclean_long,
  summary = c("hist"),
  filter = NULL
)


##zero inflated negative binomial
imp3 <- countimp(weeklySUclean_long, meth = c(" ", " "," ", "zinb"))
summary(imp3) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp3$imp$value) 

weeklySUimps_zinb <- mice::complete(imp3, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklySUmids_zinb <-as.mids(weeklySUimps_zinb)

#lm with pooled data
fitimp_zinb <- with(weeklySUmids_zinb,
                   lm(value ~ week + INT))

summary(pool(fitimp_zinb)) #the p-values look very similar so imputation is successful

#merging imputations
zinb_imp <- merge_imputations(
  weeklySUclean_long,
  weeklySUmids_zinb,
  ori = weeklySUclean_long,
  summary = c("den"),
  filter = NULL
)

View(zinb_imp$data) #contains the original df with imputed data

#based on theory alone, zinb is probably better, so I'd say its the most jutified to use ZINB data

#converting back to wide format
zinb_imp_wide <- reshape(zinb_imp$data, idvar = c("PID","INT"), timevar = "week", v.names = c("value", "value_imp"), direction = "wide")

#calculating sum totals
zinb_imp_wide$SU_total <- rowSums(zinb_imp_wide[,c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)])
zinb_imp_wide$SU_total_imp <- rowSums(zinb_imp_wide[,c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42)])

summary(zinb_imp_wide$SU_total) #mean is 18.09
summary(zinb_imp_wide$SU_total_imp) #mean is 14.54

#exporting as csv
write.csv(zinb_imp_wide, file = "C:\\Users\\Ellen Martin\\OneDrive\\Desktop\\Yale Stover Lab\\R21 Clinical Trial Data Analysis LOCAL\\weeklySUimputed.csv")



########## IPV Data ##########

weeklyIPV <- as.data.frame(read.spss("Z:\\R21 Clinical Trial\\data\\Weekly Log Data\\IPV\\IPV Week Sums_1.sav"))


weeklyIPVclean <- weeklyIPV[-c(2,10,12,16,27,28,29,30,35,36,39,41,43,48,49,54,60,61,70),]
#removing the totals column

weeklyIPVclean <- weeklyIPVclean[,-c(24,25,26,3)]

# Remove row based on condition
summary(weeklyIPVclean) #53 rows with at least 2 weekly log responses

#reshape to long

weeklyIPVclean_long <- melt(setDT(weeklyIPVclean), id.vars = c("INT","PID"), variable.name = "week")

#pattern of missingness
md.pattern(weeklyIPVclean)

aggr_plot <- aggr(weeklyIPVclean, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#linear regression just to see the fit of the data compared to imputed
fitols <- lm(value ~ week + INT, data=weeklyIPVclean_long)
summary(fitols)

#checking for overdispersion
rd <- glm(value ~ week + INT, data = weeklyIPVclean_long, family = poisson)
dispersiontest(rd,trafo=1)
#seems to be evidence of overdispersion as c = 1.83 and p = .001, so zero-inflated negative binomial probabily better than poisson

## poisson imputation
imp1_IPV <- countimp(weeklyIPVclean_long, meth = c(" "," "," ", "pois"))
summary(imp1_IPV) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp1_IPV$imp$value) 

weeklyIPVimps <- mice::complete(imp1_IPV, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklyIPVmids <-as.mids(weeklyIPVimps)

#lm with pooled data
fitimp_IPV <- with(weeklyIPVmids,
               lm(value ~ week + INT))

summary(pool(fitimp_IPV)) #the p-values look very similar so imputation is successful

#merging imputations
pois_imp_IPV <- merge_imputations(
  weeklyIPVclean_long,
  weeklyIPVmids,
  ori = weeklyIPVclean_long,
  summary = c("den"),
  filter = NULL
)

#actually looks like it fits the data super well

#zero inflated poisson - doesn't work
imp2_IPV <- countimp(weeklyIPVclean_long, meth = c(" ", " "," ", "zip"))
summary(imp2_IPV) #there are 5 datasets of imputed data which can be found in imp1$imp$value
#didn't work
head(imp2_IPV$imp$value) 

weeklyIPVimps_zip <- mice::complete(imp2_IPV, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklySUmids_zip <-as.mids(weeklySUimps_zip)

#lm with pooled data
fitimp_zip <- with(weeklySUmids_zip,
                   lm(value ~ week + INT))

summary(pool(fitimp_zip)) #the p-values look very similar so imputation is successful

#merging imputations
zip_imp <- merge_imputations(
  weeklySUclean_long,
  weeklySUmids_zip,
  ori = weeklySUclean_long,
  summary = c("hist"),
  filter = NULL
)



##zero inflated negative binomial
imp3_IPV <- countimp(weeklyIPVclean_long, meth = c(" "," "," ","zinb"))
summary(imp3) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp3$imp$value) 

weeklySUimps_zinb <- mice::complete(imp3, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklySUmids_zinb <-as.mids(weeklySUimps_zinb)

#lm with pooled data
fitimp_zinb <- with(weeklySUmids_zinb,
                    lm(value ~ week + INT))

summary(pool(fitimp_zinb)) #the p-values look very similar so imputation is successful

#merging imputations
zinb_imp <- merge_imputations(
  weeklySUclean_long,
  weeklySUmids_zinb,
  ori = weeklySUclean_long,
  summary = c("den"),
  filter = NULL
)

View(zinb_imp$data) #contains the original df with imputed data

#we will use the pois_imp_IPV becauase its the only one that converged

#back to wide format
pois_imp_IPV_wide <- reshape(pois_imp_IPV$data, idvar = c("PID","INT"), timevar = "week", v.names = c("value", "value_imp"), direction = "wide")

#calculating sum totals
pois_imp_IPV_wide$IPV_total <- rowSums(pois_imp_IPV_wide[,c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)])
pois_imp_IPV_wide$IPV_total_imp <- rowSums(pois_imp_IPV_wide[,c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42)])

summary(pois_imp_IPV_wide$IPV_total) #mean is 3.273
summary(pois_imp_IPV_wide$IPV_total_imp) #mean is 2.64

#changing names
names(pois_imp_IPV_wide) <- sub('value.', '', names(pois_imp_IPV_wide))
#exporting as csv
write.csv(pois_imp_IPV_wide, file = "C:\\Users\\Ellen Martin\\OneDrive\\Desktop\\Yale Stover Lab\\R21 Clinical Trial Data Analysis LOCAL\\weeklyIPVimputed.csv")


## Child Maltreatment Logs ##

weeklyCM <- as.data.frame(read.spss("Z:\\R21 Clinical Trial\\data\\Weekly Log Data\\CM\\CM Logs Weekly Sums.sav"))

#removing ppts with only baseline data
weeklyCMclean <- weeklyCM[-c(2,10,12,16,27,28,29,30,35,36,39,41,43,48,49,54,60,61,70),]

#removing the totals column
weeklyCMclean <- weeklyCMclean[,-c(23:44)]


# Remove row based on condition
summary(weeklyCMclean) #53 rows with at least 2 weekly log responses

#reshape to long

weeklyCMclean_long <- melt(setDT(weeklyCMclean), id.vars = c("INT","PID"), variable.name = "week")

#pattern of missingness
md.pattern(weeklyCMclean)

aggr_plot <- aggr(weeklyCMclean, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#linear regression just to see the fit of the data compared to imputed
fitols <- lm(value ~ week + INT, data=weeklyCMclean_long)
summary(fitols)

#checking for overdispersion
rd <- glm(value ~ week + INT, data = weeklyCMclean_long, family = poisson)
dispersiontest(rd,trafo=1)
#seems to be evidence of overdispersion as c = 1.18 and p < .001, so zero-inflated negative binomial probabily better than poisson

## poisson imputation
imp1 <- countimp(weeklyCMclean_long, meth = c(" ", " "," ", "pois"))
summary(imp1) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp1$imp$value) 

weeklyCMimps <- mice::complete(imp1, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklyCMmids <-as.mids(weeklyCMimps)

#lm with pooled data
fitimp <- with(weeklyCMmids,
               lm(value ~ week + INT))

summary(pool(fitimp)) #the p-values look slightly different

#merging imputations
pois_imp <- merge_imputations(
  weeklyCMclean_long,
  weeklyCMmids,
  ori = weeklyCMclean_long,
  summary = c("den"),
  filter = NULL
)


#zero inflated poisson - didn't work
imp2 <- countimp(weeklyCMclean_long, meth = c(" ", " "," ", "zip"))
summary(imp2) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp2$imp$value) 

weeklyCMimps_zip <- mice::complete(imp2, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklyCMmids_zip <-as.mids(weeklyCMimps_zip)

#lm with pooled data
fitimp_zip <- with(weeklyCMmids_zip,
                   lm(value ~ week + INT))

summary(pool(fitimp_zip)) #the p-values look very similar so imputation is successful

#merging imputations
zip_imp <- merge_imputations(
  weeklySUclean_long,
  weeklySUmids_zip,
  ori = weeklySUclean_long,
  summary = c("hist"),
  filter = NULL
)


##zero inflated negative binomial - didn't work
imp3 <- countimp(weeklyCMclean_long, meth = c(" ", " "," ", "zinb"))
summary(imp3) #there are 5 datasets of imputed data which can be found in imp1$imp$value

head(imp3$imp$value) 

weeklyCMimps_zinb <- mice::complete(imp3, action="long", include = TRUE)


# Convert back to mids type - mice can work with this type
weeklyCMmids_zinb <-as.mids(weeklyCMimps_zinb)

#lm with pooled data
fitimp_zinb <- with(weeklyCMmids_zinb,
                    lm(value ~ week + INT))

summary(pool(fitimp_zinb)) #the p-values look very similar so imputation is successful

#merging imputations
zinb_imp <- merge_imputations(
  weeklyCMclean_long,
  weeklyCMmids_zinb,
  ori = weeklyCMclean_long,
  summary = c("den"),
  filter = NULL
)


#none of the methods for CM worked except poisson so we'll stick with it

#back to wide format
pois_imp_wide <- reshape(pois_imp$data, idvar = c("PID","INT"), timevar = "week", v.names = c("value", "value_imp"), direction = "wide")

#calculating sum totals
pois_imp_wide$CM_total <- rowSums(pois_imp_wide[,c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)])
pois_imp_wide$CM_total_imp <- rowSums(pois_imp_wide[,c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42)])

summary(pois_imp_wide$CM_total) #mean is 2.54
summary(pois_imp_wide$CM_total_imp) #mean is 2.22

#exporting as csv
write.csv(pois_imp_wide, file = "C:\\Users\\Ellen Martin\\OneDrive\\Desktop\\Yale Stover Lab\\R21 Clinical Trial Data Analysis LOCAL\\weeklyCMimputed.csv")




