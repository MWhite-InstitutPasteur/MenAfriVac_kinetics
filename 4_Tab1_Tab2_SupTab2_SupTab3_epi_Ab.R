
###################################
###################################
##                               ##
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
###################################


AB_2 = read.csv("C://U//Ab_dynamics//MenAfriVac//Data//Data_002//Data_002_process.csv")

N_part_2 = nrow(AB_2)



AB_3 = read.csv("C://U//Ab_dynamics//MenAfriVac//Data//Data_003//Data_003_process.csv")

N_part_3 = nrow(AB_3)


###################################
###################################
##                               ## 
##  ######  ####  #####     ##   ##
##    ##   ##  ## ##  ##   ###   ##
##    ##   ###### #####     ##   ##
##    ##   ##  ## ##  ##    ##   ##
##    ##   ##  ## #####    ####  ##
##                               ##
###################################
###################################

TAB1 = matrix(NA, nrow=9, ncol=28)

colnames(TAB1) = c("country", "vaccine", "N", "age", "age_low", "age_high", "gender(female)", "T_follow", "T_follow_low", "T_follow_high",
                   "SBA_base", "SBA_base_low", "SBA_base_high", "SBA_peak", "SBA_peak_low", "SBA_peak_high", "SBA_boost", "SBA_boost_low", "SBA_boost_high",
                   "ELISA_base", "ELISA_base_low", "ELISA_base_high", "ELISA_peak", "ELISA_peak_low", "ELISA_peak_high", "ELISA_boost", "ELISA_boost_low", "ELISA_boost_high")
TAB1 = as.data.frame(TAB1)

TAB1[,1] = c("Gambia", "Gambia", "Gambia", "Mali", "Mali", "Mali", "Gambia",  "Mali", "Senegal")
TAB1[,2] = c("Hib-TT/PsA-TT", "PsA-TT/Hib-TT", "PsA-TT/PsA-TT", "Hib-TT/PsA-TT", "PsA-TT/Hib-TT", "PsA-TT/PsA-TT", "PsA-TT/ -",  "PsA-TT/ -",  "PsA-TT/ -") 


index_1 = intersect( which(AB_2$country=="Gambia"), intersect(which(AB_2$prim_vac=="Hib-TT"), which(AB_2$boost_vac=="PsA-TT")))

TAB1[1,3]     = length(index_1)
TAB1[1,4:6]   = quantile( AB_2$age[index_1], prob=c(0.5, 0.025, 0.975))/12  
TAB1[1,7]     = length(which(AB_2$sex[index_1]=="female"))
TAB1[1,8:10]  = quantile( (AB_2$t_v9 - AB_2$t_v1)[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[1,11:13] = quantile( AB_2$MenA_v1[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[1,14:16] = quantile( AB_2$MenA_v3[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[1,17:19] = quantile( AB_2$MenA_v6[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[1,20:22] = quantile( AB_2$elisa_v1[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[1,23:25] = quantile( AB_2$elisa_v3[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[1,26:28] = quantile( AB_2$elisa_v6[index_1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[1,11] = exp(mean( log(AB_2$MenA_v1[index_1]),na.rm=TRUE))
TAB1[1,14] = exp(mean( log(AB_2$MenA_v3[index_1]),na.rm=TRUE))
TAB1[1,17] = exp(mean( log(AB_2$MenA_v6[index_1]),na.rm=TRUE))
TAB1[1,20] = exp(mean( log(AB_2$elisa_v1[index_1]),na.rm=TRUE))
TAB1[1,23] = exp(mean( log(AB_2$elisa_v3[index_1]),na.rm=TRUE))
TAB1[1,26] = exp(mean( log(AB_2$elisa_v6[index_1]),na.rm=TRUE))




index_2 = intersect( which(AB_2$country=="Gambia"), intersect(which(AB_2$prim_vac=="PsA-TT"), which(AB_2$boost_vac=="Hib-TT")))

TAB1[2,3]     = length(index_2)
TAB1[2,4:6]   = quantile( AB_2$age[index_2], prob=c(0.5, 0.025, 0.975))/12  
TAB1[2,7]     = length(which(AB_2$sex[index_2]=="female"))
TAB1[2,8:10]  = quantile( (AB_2$t_v9 - AB_2$t_v1)[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[2,11:13] = quantile( AB_2$MenA_v1[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[2,14:16] = quantile( AB_2$MenA_v3[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[2,17:19] = quantile( AB_2$MenA_v6[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[2,20:22] = quantile( AB_2$elisa_v1[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[2,23:25] = quantile( AB_2$elisa_v3[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[2,26:28] = quantile( AB_2$elisa_v6[index_2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[2,11] = exp(mean( log(AB_2$MenA_v1[index_2]),na.rm=TRUE))
TAB1[2,14] = exp(mean( log(AB_2$MenA_v3[index_2]),na.rm=TRUE))
TAB1[2,17] = exp(mean( log(AB_2$MenA_v6[index_2]),na.rm=TRUE))
TAB1[2,20] = exp(mean( log(AB_2$elisa_v1[index_2]),na.rm=TRUE))
TAB1[2,23] = exp(mean( log(AB_2$elisa_v3[index_2]),na.rm=TRUE))
TAB1[2,26] = exp(mean( log(AB_2$elisa_v6[index_2]),na.rm=TRUE))




index_3 = intersect( which(AB_2$country=="Gambia"), intersect(which(AB_2$prim_vac=="PsA-TT"), which(AB_2$boost_vac=="PsA-TT")))

TAB1[3,3]     = length(index_3)
TAB1[3,4:6]   = quantile( AB_2$age[index_3], prob=c(0.5, 0.025, 0.975))/12  
TAB1[3,7]     = length(which(AB_2$sex[index_3]=="female"))
TAB1[3,8:10]  = quantile( (AB_2$t_v9 - AB_2$t_v1)[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[3,11:13] = quantile( AB_2$MenA_v1[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[3,14:16] = quantile( AB_2$MenA_v3[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[3,17:19] = quantile( AB_2$MenA_v6[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[3,20:22] = quantile( AB_2$elisa_v1[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[3,23:25] = quantile( AB_2$elisa_v3[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[3,26:28] = quantile( AB_2$elisa_v6[index_3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[3,11] = exp(mean( log(AB_2$MenA_v1[index_3]),na.rm=TRUE))
TAB1[3,14] = exp(mean( log(AB_2$MenA_v3[index_3]),na.rm=TRUE))
TAB1[3,17] = exp(mean( log(AB_2$MenA_v6[index_3]),na.rm=TRUE))
TAB1[3,20] = exp(mean( log(AB_2$elisa_v1[index_3]),na.rm=TRUE))
TAB1[3,23] = exp(mean( log(AB_2$elisa_v3[index_3]),na.rm=TRUE))
TAB1[3,26] = exp(mean( log(AB_2$elisa_v6[index_3]),na.rm=TRUE))



index_4 = intersect( which(AB_2$country=="Mali"), intersect(which(AB_2$prim_vac=="Hib-TT"), which(AB_2$boost_vac=="PsA-TT")))

TAB1[4,3]     = length(index_4)
TAB1[4,4:6]   = quantile( AB_2$age[index_4], prob=c(0.5, 0.025, 0.975))/12  
TAB1[4,7]     = length(which(AB_2$sex[index_4]=="female"))
TAB1[4,8:10]  = quantile( (AB_2$t_v9 - AB_2$t_v1)[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[4,11:13] = quantile( AB_2$MenA_v1[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[4,14:16] = quantile( AB_2$MenA_v3[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[4,17:19] = quantile( AB_2$MenA_v6[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[4,20:22] = quantile( AB_2$elisa_v1[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[4,23:25] = quantile( AB_2$elisa_v3[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[4,26:28] = quantile( AB_2$elisa_v6[index_4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[4,11] = exp(mean( log(AB_2$MenA_v1[index_4]),na.rm=TRUE))
TAB1[4,14] = exp(mean( log(AB_2$MenA_v3[index_4]),na.rm=TRUE))
TAB1[4,17] = exp(mean( log(AB_2$MenA_v6[index_4]),na.rm=TRUE))
TAB1[4,20] = exp(mean( log(AB_2$elisa_v1[index_4]),na.rm=TRUE))
TAB1[4,23] = exp(mean( log(AB_2$elisa_v3[index_4]),na.rm=TRUE))
TAB1[4,26] = exp(mean( log(AB_2$elisa_v6[index_4]),na.rm=TRUE))


index_5 = intersect( which(AB_2$country=="Mali"), intersect(which(AB_2$prim_vac=="PsA-TT"), which(AB_2$boost_vac=="Hib-TT")))

TAB1[5,3]     = length(index_5)
TAB1[5,4:6]   = quantile( AB_2$age[index_5], prob=c(0.5, 0.025, 0.975))/12  
TAB1[5,7]     = length(which(AB_2$sex[index_5]=="female"))
TAB1[5,8:10]  = quantile( (AB_2$t_v9 - AB_2$t_v1)[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[5,11:13] = quantile( AB_2$MenA_v1[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[5,14:16] = quantile( AB_2$MenA_v3[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[5,17:19] = quantile( AB_2$MenA_v6[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[5,20:22] = quantile( AB_2$elisa_v1[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[5,23:25] = quantile( AB_2$elisa_v3[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[5,26:28] = quantile( AB_2$elisa_v6[index_5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[5,11] = exp(mean( log(AB_2$MenA_v1[index_5]),na.rm=TRUE))
TAB1[5,14] = exp(mean( log(AB_2$MenA_v3[index_5]),na.rm=TRUE))
TAB1[5,17] = exp(mean( log(AB_2$MenA_v6[index_5]),na.rm=TRUE))
TAB1[5,20] = exp(mean( log(AB_2$elisa_v1[index_5]),na.rm=TRUE))
TAB1[5,23] = exp(mean( log(AB_2$elisa_v3[index_5]),na.rm=TRUE))
TAB1[5,26] = exp(mean( log(AB_2$elisa_v6[index_5]),na.rm=TRUE))



index_6 = intersect( which(AB_2$country=="Mali"), intersect(which(AB_2$prim_vac=="PsA-TT"), which(AB_2$boost_vac=="PsA-TT")))

TAB1[6,3]     = length(index_6)
TAB1[6,4:6]   = quantile( AB_2$age[index_6], prob=c(0.5, 0.025, 0.975))/12  
TAB1[6,7]     = length(which(AB_2$sex[index_6]=="female"))
TAB1[6,8:10]  = quantile( (AB_2$t_v9 - AB_2$t_v1)[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[6,11:13] = quantile( AB_2$MenA_v1[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[6,14:16] = quantile( AB_2$MenA_v3[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[6,17:19] = quantile( AB_2$MenA_v6[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[6,20:22] = quantile( AB_2$elisa_v1[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[6,23:25] = quantile( AB_2$elisa_v3[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[6,26:28] = quantile( AB_2$elisa_v6[index_6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[6,11] = exp(mean( log(AB_2$MenA_v1[index_6]),na.rm=TRUE))
TAB1[6,14] = exp(mean( log(AB_2$MenA_v3[index_6]),na.rm=TRUE))
TAB1[6,17] = exp(mean( log(AB_2$MenA_v6[index_6]),na.rm=TRUE))
TAB1[6,20] = exp(mean( log(AB_2$elisa_v1[index_6]),na.rm=TRUE))
TAB1[6,23] = exp(mean( log(AB_2$elisa_v3[index_6]),na.rm=TRUE))
TAB1[6,26] = exp(mean( log(AB_2$elisa_v6[index_6]),na.rm=TRUE))




index_7 = intersect( which(AB_3$country=="Gambia"), which(AB_3$prim_vac=="PsA-TT") )

TAB1[7,3]     = length(index_7)
TAB1[7,4:6]   = quantile( AB_3$age[index_7], prob=c(0.5, 0.025, 0.975)) 
TAB1[7,7]     = length(which(AB_3$sex[index_7]=="female"))
TAB1[7,8:10]  = quantile( (AB_3$t_v5 - AB_3$t_v1)[index_7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[7,11:13] = quantile( AB_3$MenA_v1[index_7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[7,14:16] = quantile( AB_3$MenA_v3[index_7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[7,20:22] = quantile( AB_3$elisa_v1[index_7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[7,23:25] = quantile( AB_3$elisa_v3[index_7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)

TAB1[7,11] = exp(mean( log(AB_3$MenA_v1[index_7]),na.rm=TRUE))
TAB1[7,14] = exp(mean( log(AB_3$MenA_v3[index_7]),na.rm=TRUE))
TAB1[7,20] = exp(mean( log(AB_3$elisa_v1[index_7]),na.rm=TRUE))
TAB1[7,23] = exp(mean( log(AB_3$elisa_v3[index_7]),na.rm=TRUE))



index_8 = intersect( which(AB_3$country=="Mali"), which(AB_3$prim_vac=="PsA-TT") )

TAB1[8,3]     = length(index_8)
TAB1[8,4:6]   = quantile( AB_3$age[index_8], prob=c(0.5, 0.025, 0.975)) 
TAB1[8,7]     = length(which(AB_3$sex[index_8]=="female"))
TAB1[8,8:10]  = quantile( (AB_3$t_v5 - AB_3$t_v1)[index_8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[8,11:13] = quantile( AB_3$MenA_v1[index_8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[8,14:16] = quantile( AB_3$MenA_v3[index_8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[8,20:22] = quantile( AB_3$elisa_v1[index_8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[8,23:25] = quantile( AB_3$elisa_v3[index_8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)


TAB1[8,11] = exp(mean( log(AB_3$MenA_v1[index_8]),na.rm=TRUE))
TAB1[8,14] = exp(mean( log(AB_3$MenA_v3[index_8]),na.rm=TRUE))
TAB1[8,20] = exp(mean( log(AB_3$elisa_v1[index_8]),na.rm=TRUE))
TAB1[8,23] = exp(mean( log(AB_3$elisa_v3[index_8]),na.rm=TRUE))




index_9 = intersect( which(AB_3$country=="Senegal"), which(AB_3$prim_vac=="PsA-TT") )

TAB1[9,3]     = length(index_9)
TAB1[9,4:6]   = quantile( AB_3$age[index_9], prob=c(0.5, 0.025, 0.975))
TAB1[9,7]     = length(which(AB_3$sex[index_9]=="female"))
TAB1[9,8:10]  = quantile( (AB_3$t_v5 - AB_3$t_v1)[index_9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[9,11:13] = quantile( AB_3$MenA_v1[index_9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[9,14:16] = quantile( AB_3$MenA_v3[index_9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[9,20:22] = quantile( AB_3$elisa_v1[index_9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)
TAB1[9,23:25] = quantile( AB_3$elisa_v3[index_9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE)


TAB1[9,11] = exp(mean( log(AB_3$MenA_v1[index_9]),na.rm=TRUE))
TAB1[9,14] = exp(mean( log(AB_3$MenA_v3[index_9]),na.rm=TRUE))
TAB1[9,20] = exp(mean( log(AB_3$elisa_v1[index_9]),na.rm=TRUE))
TAB1[9,23] = exp(mean( log(AB_3$elisa_v3[index_9]),na.rm=TRUE))



############################################################
############################################################
##                                                        ## 
##   ####  ##  ## #####    ######  ####  #####     ####   ##
##  ##     ##  ## ##  ##     ##   ##  ## ##  ##   ##  ##  ##
##   ####  ##  ## #####      ##   ###### #####       ##   ##
##      ## ##  ## ##         ##   ##  ## ##  ##     ##    ##
##   ####   ####  ##         ##   ##  ## #####     #####  ##
##                                                        ##
############################################################
############################################################
##
## A nice example of interpreting interactions can be found here:
##
## https://www.r-bloggers.com/interpreting-interaction-coefficient-in-r-part1-lm/

library(nlme)


AB_2_sub = AB_2[,c(10, 11, 2, 4, 5, 6, 7, 
                   30, 32, 33, 35, 
                   21, 23, 24, 26)]

AB_2_sub = cbind( rep("two", N_part_2), AB_2_sub )

colnames(AB_2_sub) = c("trial", "prim_vac", "boost_vac", "country", "age", "sex", "height", "weight", 
                       "SBA_base", "SBA_peak", "SBA_pre_boost", "SBA_boost", 
                       "elisa_base", "elisa_peak", "elisa_pre_boost", "elisa_boost")


AB_2_sub$age = AB_2_sub$age/12

AB_3_sub = AB_3[,c(9, 2, 4, 5, 7, 8,
                   20, 22, 15, 17)]

AB_3_sub = cbind( rep("three", N_part_3), AB_3_sub[,1], rep(NA, N_part_3),
                  AB_3_sub[,2:8],  matrix(NA, nrow=N_part_3, ncol=2),
                  AB_3_sub[,9:10], matrix(NA, nrow=N_part_3, ncol=2) )

colnames(AB_3_sub) = colnames(AB_2_sub)



AB_combi = rbind( AB_2_sub, AB_3_sub)





#############################################
## Relevel factors

AB_combi$trial <- relevel( AB_combi$trial, ref="two")

##AB_combi$prim_vac <- relevel( AB_combi$prim_vac, ref="PsA-TT")

##AB_combi$boost_vac <- relevel( AB_combi$boost_vac, ref="PsA-TT")

AB_combi$prim_vac <- relevel( AB_combi$prim_vac, ref="Hib-TT")

AB_combi$boost_vac <- relevel( AB_combi$boost_vac, ref="Hib-TT")


AB_combi$country <- relevel( AB_combi$country, ref="Gambia")

AB_combi$sex <- relevel( AB_combi$sex, ref="female")




###################
## Linear model SBA_peak

AB_mod_SBA_peak1 <- lm( log10(SBA_peak) ~ prim_vac + country + trial*age + sex + height + weight + log10(SBA_base),  
                         data=AB_combi)

summary(AB_mod_SBA_peak1)



cbind(
	summary(AB_mod_SBA_peak1)$coef[,1], 
	summary(AB_mod_SBA_peak1)$coef[,1] - 1.96*summary(AB_mod_SBA_peak1)$coef[,2],   
	summary(AB_mod_SBA_peak1)$coef[,1] + 1.96*summary(AB_mod_SBA_peak1)$coef[,2]
)



###################
## Linear model SBA_boost

AB_mod_SBA_boost1 <- lm( log10(SBA_boost) ~  boost_vac + country + age + sex + height + weight + log10(SBA_peak),  
                         data=AB_combi)

summary(AB_mod_SBA_boost1)



cbind(
	summary(AB_mod_SBA_boost1)$coef[,1], 
	summary(AB_mod_SBA_boost1)$coef[,1] - 1.96*summary(AB_mod_SBA_boost1)$coef[,2],   
	summary(AB_mod_SBA_boost1)$coef[,1] + 1.96*summary(AB_mod_SBA_boost1)$coef[,2]
)


	


###################
## Linear model elisa_peak

AB_mod_elisa_peak1 <- lm( log10(elisa_peak) ~ trial*age + prim_vac + country + sex + height + weight + log10(elisa_base),  
                         data=AB_combi)

summary(AB_mod_elisa_peak1)



cbind(
	summary(AB_mod_elisa_peak1)$coef[,1], 
	summary(AB_mod_elisa_peak1)$coef[,1] - 1.96*summary(AB_mod_elisa_peak1)$coef[,2],   
	summary(AB_mod_elisa_peak1)$coef[,1] + 1.96*summary(AB_mod_elisa_peak1)$coef[,2]
)



###################
## Linear model elisa_boost

AB_mod_elisa_boost1 <- lm( log10(elisa_boost) ~  boost_vac + country + age + sex + height + weight + log10(elisa_peak),  
                         data=AB_combi)

summary(AB_mod_elisa_boost1)



cbind(
	summary(AB_mod_elisa_boost1)$coef[,1], 
	summary(AB_mod_elisa_boost1)$coef[,1] - 1.96*summary(AB_mod_elisa_boost1)$coef[,2],   
	summary(AB_mod_elisa_boost1)$coef[,1] + 1.96*summary(AB_mod_elisa_boost1)$coef[,2]
)




############################################################
############################################################
##                                                        ## 
##   ####  ##  ## #####    ######  ####  #####     ####   ##
##  ##     ##  ## ##  ##     ##   ##  ## ##  ##   ##  ##  ##
##   ####  ##  ## #####      ##   ###### #####       ##   ##
##      ## ##  ## ##         ##   ##  ## ##  ##   ##  ##  ##
##   ####   ####  ##         ##   ##  ## #####     ####   ##
##                                                        ##
############################################################
############################################################



############################################################
## SBA_peak

AB_SBA_peak_vaccine =  lm( log10(SBA_peak) ~ prim_vac,  
                         data=AB_combi)

summary(AB_SBA_peak_vaccine)

cbind(
	summary(AB_SBA_peak_vaccine)$coef[,1], 
	summary(AB_SBA_peak_vaccine)$coef[,1] - 1.96*summary(AB_SBA_peak_vaccine)$coef[,2],   
	summary(AB_SBA_peak_vaccine)$coef[,1] + 1.96*summary(AB_SBA_peak_vaccine)$coef[,2]
)





AB_SBA_peak_study =  lm( log10(SBA_peak) ~ trial,  
                         data=AB_combi)

summary(AB_SBA_peak_study)

cbind(
	summary(AB_SBA_peak_study)$coef[,1], 
	summary(AB_SBA_peak_study)$coef[,1] - 1.96*summary(AB_SBA_peak_study)$coef[,2],   
	summary(AB_SBA_peak_study)$coef[,1] + 1.96*summary(AB_SBA_peak_study)$coef[,2]
)




AB_SBA_peak_country =  lm( log10(SBA_peak) ~ country,  
                         data=AB_combi)

summary(AB_SBA_peak_country)

cbind(
	summary(AB_SBA_peak_country)$coef[,1], 
	summary(AB_SBA_peak_country)$coef[,1] - 1.96*summary(AB_SBA_peak_country)$coef[,2],   
	summary(AB_SBA_peak_country)$coef[,1] + 1.96*summary(AB_SBA_peak_country)$coef[,2]
)






AB_SBA_peak_age_two =  lm( log10(SBA_peak) ~ age,  
                         data=AB_combi[which(AB_combi$trial=="two"),])

summary(AB_SBA_peak_age_two)

cbind(
	summary(AB_SBA_peak_age_two)$coef[,1], 
	summary(AB_SBA_peak_age_two)$coef[,1] - 1.96*summary(AB_SBA_peak_age_two)$coef[,2],   
	summary(AB_SBA_peak_age_two)$coef[,1] + 1.96*summary(AB_SBA_peak_age_two)$coef[,2]
)






AB_SBA_peak_age_three =  lm( log10(SBA_peak) ~ age,  
                         data=AB_combi[which(AB_combi$trial=="three"),])

summary(AB_SBA_peak_age_three)

cbind(
	summary(AB_SBA_peak_age_three)$coef[,1], 
	summary(AB_SBA_peak_age_three)$coef[,1] - 1.96*summary(AB_SBA_peak_age_three)$coef[,2],   
	summary(AB_SBA_peak_age_three)$coef[,1] + 1.96*summary(AB_SBA_peak_age_three)$coef[,2]
)






AB_SBA_peak_sex =  lm( log10(SBA_peak) ~ sex,  
                         data=AB_combi)

summary(AB_SBA_peak_sex)

cbind(
	summary(AB_SBA_peak_sex)$coef[,1], 
	summary(AB_SBA_peak_sex)$coef[,1] - 1.96*summary(AB_SBA_peak_sex)$coef[,2],   
	summary(AB_SBA_peak_sex)$coef[,1] + 1.96*summary(AB_SBA_peak_sex)$coef[,2]
)





AB_SBA_peak_height =  lm( log10(SBA_peak) ~ height,  
                         data=AB_combi)

summary(AB_SBA_peak_height)

cbind(
	summary(AB_SBA_peak_height)$coef[,1], 
	summary(AB_SBA_peak_height)$coef[,1] - 1.96*summary(AB_SBA_peak_height)$coef[,2],   
	summary(AB_SBA_peak_height)$coef[,1] + 1.96*summary(AB_SBA_peak_height)$coef[,2]
)




AB_SBA_peak_weight =  lm( log10(SBA_peak) ~ weight,  
                         data=AB_combi)

summary(AB_SBA_peak_weight)

cbind(
	summary(AB_SBA_peak_weight)$coef[,1], 
	summary(AB_SBA_peak_weight)$coef[,1] - 1.96*summary(AB_SBA_peak_weight)$coef[,2],   
	summary(AB_SBA_peak_weight)$coef[,1] + 1.96*summary(AB_SBA_peak_weight)$coef[,2]
)







AB_SBA_peak_base =  lm( log10(SBA_peak) ~ SBA_base,  
                         data=AB_combi)

summary(AB_SBA_peak_base)

cbind(
	summary(AB_SBA_peak_base)$coef[,1], 
	summary(AB_SBA_peak_base)$coef[,1] - 1.96*summary(AB_SBA_peak_base)$coef[,2],   
	summary(AB_SBA_peak_base)$coef[,1] + 1.96*summary(AB_SBA_peak_base)$coef[,2]
)





############################################################
## elisa_peak

AB_elisa_peak_vaccine =  lm( log10(elisa_peak) ~ prim_vac,  
                         data=AB_combi)

summary(AB_elisa_peak_vaccine)

cbind(
	summary(AB_elisa_peak_vaccine)$coef[,1], 
	summary(AB_elisa_peak_vaccine)$coef[,1] - 1.96*summary(AB_elisa_peak_vaccine)$coef[,2],   
	summary(AB_elisa_peak_vaccine)$coef[,1] + 1.96*summary(AB_elisa_peak_vaccine)$coef[,2]
)





AB_elisa_peak_study =  lm( log10(elisa_peak) ~ trial,  
                         data=AB_combi)

summary(AB_elisa_peak_study)

cbind(
	summary(AB_elisa_peak_study)$coef[,1], 
	summary(AB_elisa_peak_study)$coef[,1] - 1.96*summary(AB_elisa_peak_study)$coef[,2],   
	summary(AB_elisa_peak_study)$coef[,1] + 1.96*summary(AB_elisa_peak_study)$coef[,2]
)




AB_elisa_peak_country =  lm( log10(elisa_peak) ~ country,  
                         data=AB_combi)

summary(AB_elisa_peak_country)

cbind(
	summary(AB_elisa_peak_country)$coef[,1], 
	summary(AB_elisa_peak_country)$coef[,1] - 1.96*summary(AB_elisa_peak_country)$coef[,2],   
	summary(AB_elisa_peak_country)$coef[,1] + 1.96*summary(AB_elisa_peak_country)$coef[,2]
)






AB_elisa_peak_age_two =  lm( log10(elisa_peak) ~ age,  
                         data=AB_combi[which(AB_combi$trial=="two"),])

summary(AB_elisa_peak_age_two)

cbind(
	summary(AB_elisa_peak_age_two)$coef[,1], 
	summary(AB_elisa_peak_age_two)$coef[,1] - 1.96*summary(AB_elisa_peak_age_two)$coef[,2],   
	summary(AB_elisa_peak_age_two)$coef[,1] + 1.96*summary(AB_elisa_peak_age_two)$coef[,2]
)






AB_elisa_peak_age_three =  lm( log10(elisa_peak) ~ age,  
                         data=AB_combi[which(AB_combi$trial=="three"),])

summary(AB_elisa_peak_age_three)

cbind(
	summary(AB_elisa_peak_age_three)$coef[,1], 
	summary(AB_elisa_peak_age_three)$coef[,1] - 1.96*summary(AB_elisa_peak_age_three)$coef[,2],   
	summary(AB_elisa_peak_age_three)$coef[,1] + 1.96*summary(AB_elisa_peak_age_three)$coef[,2]
)






AB_elisa_peak_sex =  lm( log10(elisa_peak) ~ sex,  
                         data=AB_combi)

summary(AB_elisa_peak_sex)

cbind(
	summary(AB_elisa_peak_sex)$coef[,1], 
	summary(AB_elisa_peak_sex)$coef[,1] - 1.96*summary(AB_elisa_peak_sex)$coef[,2],   
	summary(AB_elisa_peak_sex)$coef[,1] + 1.96*summary(AB_elisa_peak_sex)$coef[,2]
)





AB_elisa_peak_height =  lm( log10(elisa_peak) ~ height,  
                         data=AB_combi)

summary(AB_elisa_peak_height)

cbind(
	summary(AB_elisa_peak_height)$coef[,1], 
	summary(AB_elisa_peak_height)$coef[,1] - 1.96*summary(AB_elisa_peak_height)$coef[,2],   
	summary(AB_elisa_peak_height)$coef[,1] + 1.96*summary(AB_elisa_peak_height)$coef[,2]
)




AB_elisa_peak_weight =  lm( log10(elisa_peak) ~ weight,  
                         data=AB_combi)

summary(AB_elisa_peak_weight)

cbind(
	summary(AB_elisa_peak_weight)$coef[,1], 
	summary(AB_elisa_peak_weight)$coef[,1] - 1.96*summary(AB_elisa_peak_weight)$coef[,2],   
	summary(AB_elisa_peak_weight)$coef[,1] + 1.96*summary(AB_elisa_peak_weight)$coef[,2]
)







AB_elisa_peak_base =  lm( log10(elisa_peak) ~ SBA_base,  
                         data=AB_combi)

summary(AB_elisa_peak_base)

cbind(
	summary(AB_elisa_peak_base)$coef[,1], 
	summary(AB_elisa_peak_base)$coef[,1] - 1.96*summary(AB_elisa_peak_base)$coef[,2],   
	summary(AB_elisa_peak_base)$coef[,1] + 1.96*summary(AB_elisa_peak_base)$coef[,2]
)










############################################################
## SBA_boost

AB_SBA_boost_vaccine =  lm( log10(SBA_boost) ~ boost_vac,  
                         data=AB_combi)

summary(AB_SBA_boost_vaccine)

cbind(
	summary(AB_SBA_boost_vaccine)$coef[,1], 
	summary(AB_SBA_boost_vaccine)$coef[,1] - 1.96*summary(AB_SBA_boost_vaccine)$coef[,2],   
	summary(AB_SBA_boost_vaccine)$coef[,1] + 1.96*summary(AB_SBA_boost_vaccine)$coef[,2]
)







AB_SBA_boost_country =  lm( log10(SBA_boost) ~ country,  
                         data=AB_combi)

summary(AB_SBA_boost_country)

cbind(
	summary(AB_SBA_boost_country)$coef[,1], 
	summary(AB_SBA_boost_country)$coef[,1] - 1.96*summary(AB_SBA_boost_country)$coef[,2],   
	summary(AB_SBA_boost_country)$coef[,1] + 1.96*summary(AB_SBA_boost_country)$coef[,2]
)






AB_SBA_boost_age_two =  lm( log10(SBA_boost) ~ age,  
                         data=AB_combi[which(AB_combi$trial=="two"),])

summary(AB_SBA_boost_age_two)

cbind(
	summary(AB_SBA_boost_age_two)$coef[,1], 
	summary(AB_SBA_boost_age_two)$coef[,1] - 1.96*summary(AB_SBA_boost_age_two)$coef[,2],   
	summary(AB_SBA_boost_age_two)$coef[,1] + 1.96*summary(AB_SBA_boost_age_two)$coef[,2]
)








AB_SBA_boost_sex =  lm( log10(SBA_boost) ~ sex,  
                         data=AB_combi)

summary(AB_SBA_boost_sex)

cbind(
	summary(AB_SBA_boost_sex)$coef[,1], 
	summary(AB_SBA_boost_sex)$coef[,1] - 1.96*summary(AB_SBA_boost_sex)$coef[,2],   
	summary(AB_SBA_boost_sex)$coef[,1] + 1.96*summary(AB_SBA_boost_sex)$coef[,2]
)





AB_SBA_boost_height =  lm( log10(SBA_boost) ~ height,  
                         data=AB_combi)

summary(AB_SBA_boost_height)

cbind(
	summary(AB_SBA_boost_height)$coef[,1], 
	summary(AB_SBA_boost_height)$coef[,1] - 1.96*summary(AB_SBA_boost_height)$coef[,2],   
	summary(AB_SBA_boost_height)$coef[,1] + 1.96*summary(AB_SBA_boost_height)$coef[,2]
)




AB_SBA_boost_weight =  lm( log10(SBA_boost) ~ weight,  
                         data=AB_combi)

summary(AB_SBA_boost_weight)

cbind(
	summary(AB_SBA_boost_weight)$coef[,1], 
	summary(AB_SBA_boost_weight)$coef[,1] - 1.96*summary(AB_SBA_boost_weight)$coef[,2],   
	summary(AB_SBA_boost_weight)$coef[,1] + 1.96*summary(AB_SBA_boost_weight)$coef[,2]
)







AB_SBA_boost_base =  lm( log10(SBA_boost) ~ SBA_peak,  
                         data=AB_combi)

summary(AB_SBA_boost_base)

cbind(
	summary(AB_SBA_boost_base)$coef[,1], 
	summary(AB_SBA_boost_base)$coef[,1] - 1.96*summary(AB_SBA_boost_base)$coef[,2],   
	summary(AB_SBA_boost_base)$coef[,1] + 1.96*summary(AB_SBA_boost_base)$coef[,2]
)





############################################################
## elisa_boost

AB_elisa_boost_vaccine =  lm( log10(elisa_boost) ~ boost_vac,  
                         data=AB_combi)

summary(AB_elisa_boost_vaccine)

cbind(
	summary(AB_elisa_boost_vaccine)$coef[,1], 
	summary(AB_elisa_boost_vaccine)$coef[,1] - 1.96*summary(AB_elisa_boost_vaccine)$coef[,2],   
	summary(AB_elisa_boost_vaccine)$coef[,1] + 1.96*summary(AB_elisa_boost_vaccine)$coef[,2]
)







AB_elisa_boost_country =  lm( log10(elisa_boost) ~ country,  
                         data=AB_combi)

summary(AB_elisa_boost_country)

cbind(
	summary(AB_elisa_boost_country)$coef[,1], 
	summary(AB_elisa_boost_country)$coef[,1] - 1.96*summary(AB_elisa_boost_country)$coef[,2],   
	summary(AB_elisa_boost_country)$coef[,1] + 1.96*summary(AB_elisa_boost_country)$coef[,2]
)






AB_elisa_boost_age_two =  lm( log10(elisa_boost) ~ age,  
                         data=AB_combi[which(AB_combi$trial=="two"),])

summary(AB_elisa_boost_age_two)

cbind(
	summary(AB_elisa_boost_age_two)$coef[,1], 
	summary(AB_elisa_boost_age_two)$coef[,1] - 1.96*summary(AB_elisa_boost_age_two)$coef[,2],   
	summary(AB_elisa_boost_age_two)$coef[,1] + 1.96*summary(AB_elisa_boost_age_two)$coef[,2]
)








AB_elisa_boost_sex =  lm( log10(elisa_boost) ~ sex,  
                         data=AB_combi)

summary(AB_elisa_boost_sex)

cbind(
	summary(AB_elisa_boost_sex)$coef[,1], 
	summary(AB_elisa_boost_sex)$coef[,1] - 1.96*summary(AB_elisa_boost_sex)$coef[,2],   
	summary(AB_elisa_boost_sex)$coef[,1] + 1.96*summary(AB_elisa_boost_sex)$coef[,2]
)





AB_elisa_boost_height =  lm( log10(elisa_boost) ~ height,  
                         data=AB_combi)

summary(AB_elisa_boost_height)

cbind(
	summary(AB_elisa_boost_height)$coef[,1], 
	summary(AB_elisa_boost_height)$coef[,1] - 1.96*summary(AB_elisa_boost_height)$coef[,2],   
	summary(AB_elisa_boost_height)$coef[,1] + 1.96*summary(AB_elisa_boost_height)$coef[,2]
)




AB_elisa_boost_weight =  lm( log10(elisa_boost) ~ weight,  
                         data=AB_combi)

summary(AB_elisa_boost_weight)

cbind(
	summary(AB_elisa_boost_weight)$coef[,1], 
	summary(AB_elisa_boost_weight)$coef[,1] - 1.96*summary(AB_elisa_boost_weight)$coef[,2],   
	summary(AB_elisa_boost_weight)$coef[,1] + 1.96*summary(AB_elisa_boost_weight)$coef[,2]
)







AB_elisa_boost_base =  lm( log10(elisa_boost) ~ elisa_peak,  
                         data=AB_combi)

summary(AB_elisa_boost_base)

cbind(
	summary(AB_elisa_boost_base)$coef[,1], 
	summary(AB_elisa_boost_base)$coef[,1] - 1.96*summary(AB_elisa_boost_base)$coef[,2],   
	summary(AB_elisa_boost_base)$coef[,1] + 1.96*summary(AB_elisa_boost_base)$coef[,2]
)







