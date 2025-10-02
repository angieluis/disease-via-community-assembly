#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(vegan) 
library(ape) 
library(ggbiplot)
library(tidyverse) 
library(gridExtra)
library(gclus) 
library(scales) 
library(betapart)
library(lmtest)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
webspeciesdata <- read.csv("WebSpeciesdata.csv")
envdat <- read.csv("EnvironmentalData.csv", row.names = 1)
trait.data <- read.csv("SpeciesTraitData.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate number of infected deer mice as density * prevalence
webspeciesdata$pm.inf <- webspeciesdata$pm * webspeciesdata$prevalence 

webspeciesdata <- webspeciesdata %>%
  mutate(site = str_split_i(web, "[.]", 1))

spe <- as.matrix(webspeciesdata[,3:36]) # untransformed mean abundance estimates
rownames(spe) <- webspeciesdata$web
spel <- log(spe+1) # log + 1 transform



# standardize environmental data so mean=0 and unit variance
envdat.sd <- decostand(envdat, "standardize", MARGIN=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dilution Effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


dilution.prev.plot <- ggplot(webspeciesdata,aes(x=invSimpson, y=prevalence, group=in2018pub)) +
  geom_point(size = 2.5, aes(col=in2018pub, shape = in2018pub)) +
  coord_cartesian(ylim=c(0,0.31)) +
  geom_smooth(method="lm",aes(col=in2018pub)) +
  xlab("Rodent diversity") +
  ylab("SNV prevalence") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("N" = "firebrick", "Y"="steelblue"))
dilution.prev.plot

dilmod <- lm(prevalence ~ invSimpson, data = webspeciesdata)
summary(dilmod) 
summary(lm(prevalence ~ invSimpson*in2018pub, data = webspeciesdata)) # interaction is significant



# log density of infected mice + 1 as the risk metric
dilution.loginf.plot <- ggplot(webspeciesdata,aes(x=invSimpson, y=log(pm.inf+1), group=in2018pub)) +
  geom_point(size = 2.5, aes(col=in2018pub, shape = in2018pub)) +
  coord_cartesian(ylim=c(0,max(log(webspeciesdata$pm.inf+1)))) +
  geom_smooth(method="lm",aes(col=in2018pub)) +
  xlab("Rodent diversity") +
  ylab("Infected host density") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("N" = "firebrick", "Y"="steelblue"))
dilution.loginf.plot

dilmod2 <- lm(log(pm.inf+1) ~ invSimpson, data = webspeciesdata)
summary(dilmod2) 
summary(lm(log(pm.inf+1) ~ invSimpson*in2018pub, data = webspeciesdata)) # interaction is significant

library(patchwork)
dilution.prev.plot / dilution.loginf.plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta diversity and Nestedness
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

spe.i <- replace(spe, spe>0,1)

# beta diversity by Sorensen indices
beta.multi(spe.i)

# nestedness analyses
out <- nestedtemp(spe.i)
plot(out, kind="incid", names=TRUE)
csim <- oecosimu(spe.i,nestedchecker,"r00",statistic="C.score", nsimul = 999) #alternative = "greater", 
ntemp <- oecosimu(spe.i,nestedtemp,"quasiswap", nsimul = 999)
nnodf <- oecosimu(spe.i, nestednodf, "quasiswap",  nsimul = 999) #alternative = "greater",

# make a table of results
results_df <- data.frame(
  Method = c("Matrix Temperature", "C-score", "NODF"),
  Statistic = unname(c(ntemp$statistic$statistic,
                       csim$statistic$C.score,
                       nnodf$statistic$statistic[3])),
  Mean_Null = c(mean(ntemp$oecosimu$simulated),
                mean(csim$oecosimu$simulated),
                mean(nnodf$oecosimu$simulated)),
  SES = unname(c((ntemp$oecosimu$statistic - mean(ntemp$oecosimu$simulated)) / sd(ntemp$oecosimu$simulated),  # ntemp doesn't return SES directly
                 (csim$statistic$C.score - mean(csim$oecosimu$simulated)) / sd(csim$oecosimu$simulated),
                 (nnodf$statistic$statistic[3] - mean(nnodf$oecosimu$simulated)) / sd(nnodf$oecosimu$simulated))),
  P_value = c(ntemp$oecosimu$pval,
              csim$oecosimu$pval,
              nnodf$oecosimu$pval[3])
)

results_df[, 2:4] <- round(results_df[, 2:4], digits = 2)
results_df[, 5] <- round(results_df[, 5], digits = 3)
# print(xtable(results_df), include.rownames = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Redundancy Analysis (RDA)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Hellinger-transform the species dataset
spe.hel <- decostand(spe, "hellinger")

## RDA of the Hellinger-transformed species data constrained by all the environmental variables
spe.rda <- rda(spe.hel ~ ., as.data.frame(envdat.sd))

summary(spe.rda)$cont	# Scaling 2 (default)

# Canonical coefficients from the rda object
coef(spe.rda)
# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(spe.rda)$r.squared)
# 0.8435653

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)
# 0.7601334

# Global test of the RDA result
anova(spe.rda, permutations = how(nperm = 9999))

anova(spe.rda, by = "axis", permutations = how(nperm = 999))
# first 5 canonical axes are significant.

anova(spe.rda, by = "terms", permutations = how(nperm = 999))
# all environmental variables are significant


# RDA plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

spe.rda.sum <- vegan::scores(spe.rda,tidy=TRUE)

df1 <- spe.rda.sum %>%
  filter(score=="species")
df2 <- spe.rda.sum %>%
  filter(score=="biplot")

rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2)) +
  geom_text(aes(label=label),size=4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #coord_fixed() + # makes sure x and y axes are the same lengths
  # xlim(-0.09, 0.09) +
  # ylim(-0.09,0.09) +
  theme_classic()
#rda.plot
rda.plot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=RDA1*.75, y=0, yend=RDA2*.75), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=RDA1*.75,y=RDA2*.75,label=label,
                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="red", size=4)
rda.plot


## zoomed in

rda.plot2 <- ggplot(df1, aes(x=RDA1, y=RDA2)) +
  geom_text(aes(label=label),size=4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  #coord_fixed() + # makes sure x and y axes are the same lengths
  xlim(-0.095, 0.095) +
  ylim(-0.095,0.095) +
  theme_classic() +
  geom_segment(data=df2, aes(x=0, xend=RDA1*.08, y=0, yend=RDA2*.08), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=RDA1*.08,y=RDA2*.08,label=label,
                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="red", size=4)
rda.plot2


spe.rda.pca.sum <- vegan::scores(spe.rda,choices= 9:10, display = c("sp"), tidy=TRUE)
df3 <- spe.rda.pca.sum %>%
  filter(score=="species") 

rda.pcplot <- ggplot(df3, aes(x=PC1, y=PC2)) + 
  geom_text(aes(label=label),size=4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_classic()

rda.pcplot

grid.arrange(rda.plot, rda.plot2, rda.pcplot, nrow = 1)


#### 3D plot

library(rgl)
knitr::knit_hooks$set(webgl = hook_webgl)
x<- summary(spe.rda)$species
x.env <- summary(spe.rda)$biplot
plot3d(x[,1], x[,2], x[,3], xlab="RDA1", ylab="RDA2", zlab="RDA3", type="n")
text3d(x[,1], x[,2], x[,3], texts = rownames(x))
plot3d(c(0,x.env[1,1]), c(0, x.env[1,2] ), c(0,x.env[1,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[1,1], x.env[1,2], x.env[1,3], texts = "tmin", col="blue")
plot3d(c(0,x.env[2,1]), c(0, x.env[2,2] ), c(0,x.env[2,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[2,1], x.env[2,2], x.env[2,3], texts = "tmax", col="blue")
plot3d(c(0,x.env[3,1]), c(0, x.env[3,2] ), c(0,x.env[3,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[3,1], x.env[3,2], x.env[3,3], texts = "prcp", col="blue")
plot3d(c(0,x.env[4,1]), c(0, x.env[4,2] ), c(0,x.env[4,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[4,1], x.env[4,2], x.env[4,3], texts = "snow", col="blue")
plot3d(c(0,x.env[5,1]), c(0, x.env[5,2] ), c(0,x.env[5,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[5,1], x.env[5,2], x.env[5,3], texts = "elev", col="blue")
plot3d(c(0,x.env[6,1]), c(0, x.env[6,2] ), c(0,x.env[6,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[6,1], x.env[6,2], x.env[6,3], texts = "biomass", col="blue")
plot3d(c(0,x.env[7,1]), c(0, x.env[7,2] ), c(0,x.env[7,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[7,1], x.env[7,2], x.env[7,3], texts = "tree cover", col="blue")
plot3d(c(0,x.env[8,1]), c(0, x.env[8,2] ), c(0,x.env[8,3]) , add=TRUE, type="l", col="blue")
text3d(x.env[8,1], x.env[8,2], x.env[8,3], texts = "bare ground", col="blue")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Abiotic Regressions on deer mouse abundance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# How well do environmental variables predict deer mouse abundance?

dat <- webspeciesdata %>%
  select(site, web, n_months, pm) %>%
  mutate(pm = log(webspeciesdata$pm + 1))


# PCA on correlated environmental variables ~~~~~~~~~~~~~~~~~~~~~
pca_env <- prcomp(envdat.sd)
summary(pca_env)
biplot(pca_env, col=c("gray","red"))
biplot(pca_env, choices = 3:4, col=c("gray","red"))

dat$PC1 <- pca_env$x[,1]
dat$PC2 <- pca_env$x[,2]
dat$PC3 <- pca_env$x[,3]
dat$PC4 <- pca_env$x[,4]

# web as the replicate
envintmod.wt.max <- lm(pm ~ PC1 + PC2 + PC3 + PC4, data=dat, weights = n_months)
step(envintmod.wt.max)
envintmod.wt <- lm(pm ~ PC1 + PC2 + PC3, data=dat, weights = n_months)
summary(envintmod.wt)
envintmod.wt.min <- lm(pm ~ PC1, data=dat, weights = n_months)
summary(envintmod.wt.min)

# same thing but now site as the replicate
sdat <- dat %>%
  group_by(site) %>%
  summarise(across(c(2:7), mean)) # take the means within sites
summary(lm(pm ~ PC1, data=sdat, weights = n_months))
# same results



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Biotic Regressions on deer mouse abundance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Do phylogenetic relationships or functional trait overlap predict the effect 
# of other species on deer mice?  
  
# Prepare Phylogenetic Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get phylogenetic tree 
tree <- read.tree("PhyloTree.txt")
plot.phylo(tree,no.margin=TRUE,cex=0.9)

# calculate distances
rodent.phylodistmat <- cophenetic.phylo(tree)
# re-order alphabetically
rodent.phylodistmat <- rodent.phylodistmat[sort(rownames(rodent.phylodistmat)),
                                           sort(colnames(rodent.phylodistmat))]

# phylogenetic distance from pm & reorder to match spe (order of 2 letter codes)
phylo.dist <- rodent.phylodistmat[which(rownames(rodent.phylodistmat)=="Peromyscus_maniculatus"),
                                  match(trait.data$species,colnames(rodent.phylodistmat))]



# Prepare Trait Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pm.ind <- which(trait.data$sp.code=="pm")
### How different is each species diet from P. maniculatus? 0 is total overlap, 200 is no overlap
diet.dist <- apply(trait.data[,5:14],1,function(x){sum(abs(trait.data[pm.ind,5:14]-x))})

# difference in mass between each species and deer mice
mass.dist <- abs(trait.data$mass[pm.ind] - trait.data$mass)

# 0 if also nocturnal, 1 if not
activity.dist <- 1 - trait.data$Nocturnal

distance.data <- data.frame(trait.data[,1:2], phylo.dist,
                            diet.dist, mass.dist, activity.dist)

rownames(distance.data) <- NULL
distance.data


# Convert distances in phylogeny and traits to similarity ~~~~~~~~~~~~~~~~~~~~
# ranges from 0 to 1, where 1 is complete overlap
similarity.data <- distance.data[,1:2]
similarity.data$phylo <- 1 - distance.data$phylo.dist/max(distance.data$phylo.dist)
similarity.data$diet <- 1 - distance.data$diet.dist/200
similarity.data$mass <- 1 - distance.data$mass.dist/max(distance.data$mass.dist)
similarity.data$activity <- 1 - distance.data$activity.dist
similarity.data$size <- trait.data$mass - trait.data$mass[pm.ind]
similarity.data$size <- log(similarity.data$size - min(similarity.data$size)+1) # log it
similarity.data$size <- similarity.data$size/max(similarity.data$size) #divide by highest so 0 to 1 like rest of them

rownames(similarity.data) <- NULL
similarity.data.inclpm <- similarity.data 
# remove deer mice for analysis
similarity.data <- similarity.data[-which(similarity.data$species=="Peromyscus_maniculatus"),]




# Calculate the different competition coefficients --------------------------- #

# use similarity of each species to pm as competition coefficient, alpha
# then multiply it by the abundance of each species 

pm.ind <- which(colnames(spel)=="pm")
pm <- spel[,pm.ind]
spel.nopm <- spel[,-pm.ind]
webspeciesdata$sum.phylocomp <- spel.nopm %*% similarity.data$phylo # summed phylogenetic competition
webspeciesdata$sum.dietcomp <- spel.nopm %*% similarity.data$diet # summed diet competition
webspeciesdata$sum.masscomp <- spel.nopm %*% similarity.data$mass # competition based on similarity in mass
webspeciesdata$sum.actcomp <- spel.nopm %*% similarity.data$activity # competition with other nocturnal rodents
webspeciesdata$sum.sizecomp <- spel.nopm %*% similarity.data$size # competition with larger rodents


# merge with environmental data ---------------------------------------------- #

webspeciesdata <- left_join(webspeciesdata, 
                                     mutate(envdat.sd, web = rownames(envdat.sd)))



# Average over webs within sites so can analyze with site as replicate ------- #

mean.sitedata <- webspeciesdata %>%
  group_by(site) %>%
  summarise(across(c(2:38,40:52), mean))


# create tables to hold model results
web.model.table <- data.frame(model = character(), AIC = numeric(), 
                              F_stat = numeric(), p_val = numeric(), R2 = numeric())
site.model.table <- data.frame(model = character(), AIC = numeric(), 
                               F_stat = numeric(), p_val = numeric(), R2 = numeric())

# Simpson's diversity model for comparison ----------------------------------- #

# using web as replicate
divmod <- lm(log(pm+1) ~ invSimpson, data = webspeciesdata, weights = n_months)
web.model.table <- data.frame(model = "Simpson's diversity", 
                              AIC = AIC(divmod), 
                              F_stat = summary(divmod)$fstatistic[1],
                              p_val = anova(divmod)$`Pr(>F)`[1], 
                              R2 = summary(divmod)$adj.r.squared)

# using site as replicate
sdivmod <- lm(log(pm+1) ~ invSimpson, data = mean.sitedata, weights = n_months)
site.model.table <- data.frame(model = "Simpson's diversity", 
                              AIC = AIC(sdivmod), 
                              F_stat = summary(sdivmod)$fstatistic[1],
                              p_val = anova(sdivmod)$`Pr(>F)`[1], 
                              R2 = summary(sdivmod)$adj.r.squared)


# Regression using phylogenetic similarity ----------------------------------- #

# using web as replicate
phylomod <- lm(log(pm+1) ~ sum.phylocomp, data = webspeciesdata, weights = n_months)
web.model.table <- rbind(web.model.table,
                         data.frame(model = "Phylogenetic Competition", 
                              AIC = AIC(phylomod), 
                              F_stat = summary(phylomod)$fstatistic[1],
                              p_val = anova(phylomod)$`Pr(>F)`[1], 
                              R2 = summary(phylomod)$adj.r.squared))

# using site as replicate
sphylomod <- lm(log(pm+1) ~ sum.phylocomp, data = mean.sitedata, weights = n_months)
site.model.table <- rbind(site.model.table,
                         data.frame(model = "Phylogenetic Competition", 
                                    AIC = AIC(sphylomod), 
                                    F_stat = summary(sphylomod)$fstatistic[1],
                                    p_val = anova(sphylomod)$`Pr(>F)`[1], 
                                    R2 = summary(sphylomod)$adj.r.squared))


# Regression using Diet similarity ------------------------------------------- #

# web as replicate
dietmod <- lm(log(pm+1) ~ sum.dietcomp,data = webspeciesdata, weights = n_months)
web.model.table <- rbind(web.model.table,
                         data.frame(model = "Diet Competition", 
                                    AIC = AIC(dietmod), 
                                    F_stat = summary(dietmod)$fstatistic[1],
                                    p_val = anova(dietmod)$`Pr(>F)`[1], 
                                    R2 = summary(dietmod)$adj.r.squared))

# site as replicate
sdietmod <- lm(log(pm+1) ~ sum.dietcomp, data = mean.sitedata, weights = n_months)
site.model.table <- rbind(site.model.table,
                         data.frame(model = "Diet Competition", 
                                    AIC = AIC(sdietmod),  
                                    F_stat = summary(sdietmod)$fstatistic[1],
                                    p_val = anova(sdietmod)$`Pr(>F)`[1], 
                                    R2 = summary(sdietmod)$adj.r.squared))


# Regression using Mass similarity ------------------------------------------- #

# web as replicate
massmod <- lm(log(pm+1) ~ sum.masscomp, data = webspeciesdata, weights = n_months)
web.model.table <- rbind(web.model.table,
                          data.frame(model = "Mass Competition", 
                                     AIC = AIC(massmod),  
                                     F_stat = summary(massmod)$fstatistic[1],
                                     p_val = anova(massmod)$`Pr(>F)`[1], 
                                     R2 = summary(massmod)$adj.r.squared))

# site as replicate
smassmod <- lm(log(pm+1) ~ sum.masscomp, data = mean.sitedata, weights = n_months)
site.model.table <- rbind(site.model.table,
                          data.frame(model = "Mass Competition", 
                                     AIC = AIC(smassmod),  
                                     F_stat = summary(smassmod)$fstatistic[1],
                                     p_val = anova(smassmod)$`Pr(>F)`[1], 
                                     R2 = summary(smassmod)$adj.r.squared))



# Regression using activity similarity (nocturnal or not) -------------------- #

# web as replicate
actmod <- lm(log(pm+1) ~ sum.actcomp, data = webspeciesdata, weights = n_months)
web.model.table <- rbind(web.model.table,
                          data.frame(model = "Nocturnal Competition", 
                                     AIC = AIC(actmod),  
                                     F_stat = summary(actmod)$fstatistic[1],
                                     p_val = anova(actmod)$`Pr(>F)`[1], 
                                     R2 = summary(actmod)$adj.r.squared))

# site as replicate
sactmod <- lm(log(pm+1) ~ sum.actcomp, data = mean.sitedata, weights = n_months)
site.model.table <- rbind(site.model.table,
                          data.frame(model = "Nocturnal Competition", 
                                     AIC = AIC(sactmod),  
                                     F_stat = summary(sactmod)$fstatistic[1],
                                     p_val = anova(sactmod)$`Pr(>F)`[1], 
                                     R2 = summary(sactmod)$adj.r.squared))


# Regression using larger size ----------------------------------------------- #
# i.e., is there more of an effect of individuals of larger species?

# web as replicate
sizemod <- lm(log(pm+1) ~ sum.sizecomp, data = webspeciesdata, weights = n_months)
web.model.table <- rbind(web.model.table,
                          data.frame(model = "Larger Competition", 
                                     AIC = AIC(sizemod),  
                                     F_stat = summary(sizemod)$fstatistic[1],
                                     p_val = anova(sizemod)$`Pr(>F)`[1], 
                                     R2 = summary(sizemod)$adj.r.squared))

# site as replicate
ssizemod <- lm(log(pm+1) ~ sum.sizecomp, data = mean.sitedata, weights = n_months)
site.model.table <- rbind(site.model.table,
                          data.frame(model = "Larger Competition", 
                                     AIC = AIC(ssizemod), 
                                     F_stat = summary(ssizemod)$fstatistic[1],
                                     p_val = anova(ssizemod)$`Pr(>F)`[1], 
                                     R2 = summary(ssizemod)$adj.r.squared))


# Adding mitigating effects of biomass --------------------------------------- #
# Can sites with higher carrying capacity (herbaceous biomass) mitigate effects
# of diet competition and improve model fit?

# web as replicate
dietbiomassmod <- lm(log(pm+1)  ~ sum.dietcomp + biomass, 
                      data = webspeciesdata, weights = n_months)
web.model.table <- rbind(web.model.table,
                          data.frame(model = "Diet Competition + Biomass", 
                                     AIC = AIC(dietbiomassmod), 
                                     F_stat = summary(dietbiomassmod)$fstatistic[1],
                                     p_val = anova(dietbiomassmod)$`Pr(>F)`[1], 
                                     R2 = summary(dietbiomassmod)$adj.r.squared))

# site as replicate
sdietbiomassmod <- lm(log(pm+1)  ~ sum.dietcomp + biomass, 
                     data = mean.sitedata, weights = n_months)
site.model.table <- rbind(site.model.table,
                          data.frame(model = "Diet Competition + Biomass", 
                                     AIC = AIC(sdietbiomassmod), 
                                     F_stat = summary(sdietbiomassmod)$fstatistic[1],
                                     p_val = anova(sdietbiomassmod)$`Pr(>F)`[1], 
                                     R2 = summary(sdietbiomassmod)$adj.r.squared))

# this is the best model for either web or site as replicate

# order tables of model results by AIC
web.model.table <- web.model.table[order(web.model.table$AIC),]
web.model.table
site.model.table <- site.model.table[order(site.model.table$AIC),]
site.model.table


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Structural Equation Models                                                              
# Environment and species interactions both important.
# How do they fit together?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(lavaan)
library(semPlot)

# do a separate PCA on abiotic variables -- removing biomass because want it as own variable
pca_abiotic <- prcomp(envdat.sd[,-6])

# scale all data so mean of 0 and sd of 1
data <- webspeciesdata %>%
  select(tmin,tmax,prcp,swe,elevation,biomass,pm,sum.dietcomp, prevalence, pm.inf, invSimpson) %>%
  mutate(PC1 = scale(pca_abiotic$x[,1]), 
         PC2 = scale(pca_abiotic$x[,2]),
         PC3 = scale(pca_abiotic$x[,3]),
         PC4 = scale(pca_abiotic$x[,4]),
         pm = scale(log(pm +1)),
         prevalence = scale(prevalence),
         l.pm.inf = scale(log(pm.inf+1)),
         pm.inf = scale(pm.inf),
         sum.dietcomp = scale(sum.dietcomp),
         D = scale(invSimpson))
data$abiotic <- with(data, PC1 + PC2 + PC3)


############# hypothesized model:
sem.model <- '
  # regressions
  sum.dietcomp ~ abiotic + biomass
  biomass ~ abiotic
  pm ~ sum.dietcomp +  biomass
  prevalence ~ pm 
  
  '
pathfit <- sem(sem.model, data=data)
summary(pathfit, fit.measures = TRUE, rsquare = TRUE)
semPaths(pathfit,'std', layout='tree2', nCharNodes = 9, sizeMan = 8, sizeLat = 10)


### log(infected mouse density + 1) instead of prevalence
sem.model.lihd <- '
  # regressions
  sum.dietcomp ~ abiotic + biomass
  biomass ~ abiotic
  pm ~ sum.dietcomp +  biomass
  l.pm.inf ~ pm 
  
  '
pathfit.lihd <- sem(sem.model.lihd, data=data)
summary(pathfit.lihd, fit.measures = TRUE, rsquare = TRUE)
semPaths(pathfit.lihd,'std', layout='tree2', nCharNodes = 9, sizeMan = 8, sizeLat = 10)
# nearly identical to prevalence model


############# does abiotic also directly affect pm?

sem.model2 <- '
  # regressions
  sum.dietcomp ~ abiotic + biomass
  biomass ~ abiotic
  pm ~ sum.dietcomp +  biomass + abiotic
  prevalence ~ pm 
  
  '
pathfit2 <- sem(sem.model2, data=data)
summary(pathfit2, fit.measures = TRUE, rsquare = TRUE)
semPaths(pathfit2,'std', layout='tree2', nCharNodes = 9, sizeMan = 8, sizeLat = 10)
### no the path is not significant and the overall model is worse and no longer significant



############# does diversity affect pm?
sem.model3 <- '
  # regressions
  sum.dietcomp ~ abiotic + biomass
  biomass ~ abiotic
  pm ~ sum.dietcomp +  biomass + D
  prevalence ~ pm
  
  '
pathfit3 <- sem(sem.model3, data=data)
summary(pathfit3, fit.measures = TRUE, rsquare = TRUE)
semPaths(pathfit3,'std', layout='tree2', nCharNodes = 9, sizeMan = 8, sizeLat = 10)
### no the path is not significant and the overall model is worse and no longer significant



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Community Clustering                                                              
# Which sites are most similar based on their environmental conditions? 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Euclidean distance matrix of the standardized environmental data 
env.de <- dist(envdat.sd)


#### Compare clustering algorithms  ------------------------------------------ #

# Compute single linkage agglomerative clustering using chord distances
env.de.single <- hclust(env.de, method = "single")

# Compute complete-linkage agglomerative clustering
env.de.complete <- hclust(env.de, method = "complete")

# Compute UPGMA agglomerative clustering
env.de.UPGMA <- hclust(env.de, method = "average")

# Compute centroid clustering
env.de.centroid <- hclust(env.de, method = "centroid")


### To test which is the 'best', we can look at 
# Cophenetic correlations ---------------------------------------------------- #
# This calculates R between the original dissimilarity matrix and the cophenetic matrix.
# The cophenetic distance between two objects in a dendrogram is the distance 
# where the two objects become members of the same group. Locate any two objects, 
# start from one, and “climb up the tree” to the first node leading down to the 
# second object: the level of that node along the distance scale is the 
# cophenetic distance between the two objects. A cophenetic matrix is a matrix 
# representing the cophenetic distances among all pairs of objects. 

# Single linkage clustering
env.de.single.coph <- cophenetic(env.de.single)
cor(env.de, env.de.single.coph) 
# Complete linkage clustering
env.de.comp.coph <- cophenetic(env.de.complete)
cor(env.de,env.de.comp.coph) 
# Average clustering
env.de.UPGMA.coph <- cophenetic(env.de.UPGMA)
cor(env.de, env.de.UPGMA.coph) 
# Centroid clustering
env.de.centroid.coph <- cophenetic(env.de.centroid)
cor(env.de,env.de.centroid.coph)  

# UPGMA is best, so use that

k <- 5  

env.UPGMA.g <- cutree(env.de.UPGMA, k = k)


# Reorder clusters
env.deo <- reorder.hclust(env.de.UPGMA,env.de)

# reordered dendrogram by group color
#hcoplot(env.de.UPGMA, env.de, lab = rownames(envdat.sd), k = k)

# Convert the "hclust" object into a "dendrogram" object
dend <- as.dendrogram(env.deo)

cols <- hue_pal()(5)

heatmap(
  t(envdat.sd),
  Rowv = NA,
  Colv = dend,
  col = gray.colors(10,start=1,end=0.1),
  scale = "none",
  margin = c(8, 8),
  ylab = "Standardized Environmental Variables",
  xlab = "Sites",
  ColSideColors = cols[env.UPGMA.g]
)

# add density, diversity, SNV for illustration
rods <- data.frame(pm, D = webspeciesdata$invSimpson,
                   SNV = webspeciesdata$prevalence)
rownames(rods) <- rownames(envdat.sd)

heatmap(
  t(cbind( decostand(rods,"standardize", MARGIN=2), envdat.sd)),
  cexRow=2,cexCol=1.5,
  Rowv = NA,
  Colv = dend,
  col = gray.colors(10,start=1,end=0.1),
  scale = "none",
  margin = c(12, 8),
  #ylab = "Standardized Environmental Variables",
  xlab = "Sites",
  ColSideColors = cols[env.UPGMA.g],
  RowSideColors = c(rep("darkred",dim(rods)[2]),rep("white",dim(envdat.sd)[2]))
)

# Now with infected host density instead of prevalence
rods2 <- data.frame(pm, D = webspeciesdata$invSimpson,
                   IHD = webspeciesdata$pm.inf)
rownames(rods2) <- rownames(envdat.sd)

heatmap(
  t(cbind( decostand(rods2,"standardize", MARGIN=2), envdat.sd)),
  cexRow=2,cexCol=1.5,
  Rowv = NA,
  Colv = dend,
  col = gray.colors(10,start=1,end=0.1),
  scale = "none",
  margin = c(12, 8),
  #ylab = "Standardized Environmental Variables",
  xlab = "Sites",
  ColSideColors = cols[env.UPGMA.g],
  RowSideColors = c(rep("darkred",dim(rods)[2]),rep("white",dim(envdat.sd)[2]))
)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How do the community clusters affect disease processes?                                                             
# Is the dilution effect only present in some clusters? 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add cluster to the data
webspeciesdata$cluster <- factor(env.UPGMA.g)


ggplot(webspeciesdata, aes(invSimpson, prevalence, col=cluster)) + 
  ggtitle("Dilution Effect by Environmental Clusters") +
  labs(x="Rodent Diversity", y="SNV prevalence") + 
  coord_cartesian(ylim=c(0,0.31)) +
  geom_point(fill=as.numeric(webspeciesdata$cluster)+1) +
  theme_classic(base_size = 18)

# There is a dilution effect in clusters 1, 2, and 4. 
# Now group these together and see if they vary in biological processes.

webspeciesdata <- mutate(webspeciesdata, cluster124 = factor(ifelse(cluster==1|cluster==2|cluster==4,0,1)))

grid.arrange(
  ggplot(webspeciesdata, aes(invSimpson, prevalence)) + 
    #ggtitle("Dilution Effect by Environmental Clusters") +
    labs(x="Rodent Diversity", y="SNV prevalence") + 
    coord_cartesian(ylim=c(0,0.31)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray", aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  ,
  ggplot(webspeciesdata, aes(log(pm+1), prevalence)) + 
    #ggtitle("Density Dependent prevalence by Environmental Clusters") +
    labs(x="Deer mouse density", y="SNV prevalence") + 
    coord_cartesian(ylim=c(0,0.31)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray", aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  ,
  ggplot(webspeciesdata, aes(invSimpson, log(pm+1))) + 
    #ggtitle("Effect of diversity on density") +
    labs(x="Rodent Diversity", y="Deer mouse density") + 
    coord_cartesian(ylim=c(0,4.5)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray", aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  ,
  ggplot(webspeciesdata, 
         aes(invSimpson, sum.dietcomp)) + 
    #ggtitle("How does diversity relate to competition?") +
    labs(x="Rodent Diversity", y="Diet Competition") + 
    #coord_cartesian(ylim=c(0,0.31)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray",aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  , nrow=2)

# Statistics
clustmod1 <- lm(prevalence~invSimpson*cluster124 -1, data=webspeciesdata) # xtable(tidy(clustmod1))
clustmod2 <- lm(prevalence~log(pm+1)*cluster124 -1, data=webspeciesdata)
clustmod3 <- lm(log(pm+1)~invSimpson*cluster124 -1, data=webspeciesdata)
clustmod4 <- lm(sum.dietcomp~invSimpson*cluster124 -1, data=webspeciesdata)


#### with infected host density as the disease metric
grid.arrange(
  ggplot(webspeciesdata, aes(invSimpson, log(pm.inf+1))) + 
    labs(x="Rodent Diversity", y="Infected host density") + 
    coord_cartesian(ylim=c(0,3)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray", aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  ,
  ggplot(webspeciesdata, aes(log(pm+1), log(pm.inf+1))) + 
    labs(x="Deer mouse density", y="Infected host density") + 
    coord_cartesian(ylim=c(0,3)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray", aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  ,
  ggplot(webspeciesdata, aes(invSimpson, log(pm+1))) + 
    labs(x="Rodent Diversity", y="Deer mouse density") + 
    coord_cartesian(ylim=c(0,4.5)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray", aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  ,
  ggplot(webspeciesdata, aes(invSimpson, sum.dietcomp)) + 
    labs(x="Rodent Diversity", y="Diet Competition") + 
    #coord_cartesian(ylim=c(0,0.31)) +
    geom_point(size = 2.5, aes(col=cluster)) +
    geom_smooth(method='lm', color="darkgray",aes(linetype=cluster124)) +
    theme_classic(base_size = 18)
  , nrow=2)

# Statistics
clustmod1a <- lm(log(pm.inf+1)~invSimpson*cluster124, data=webspeciesdata) # xtable(tidy(clustmod1))
clustmod2a <- lm(log(pm.inf+1)~log(pm+1)*cluster124, data=webspeciesdata)
