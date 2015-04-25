# Charles Zheng
# Code for calculating simplex probabilities

empprob <- function(a,b,res0=10000) {
  m1 <- length(a)
  m2 <- length(b)
  mat1 <- matrix(rexp(m1*res0),res0,m1)
  mat2 <- matrix(rexp(m2*res0),res0,m2)
  sums <- mat1 %*% a - mat2 %*% b
  sum(sums > 0)/res0
}

calcCs <- function(a) {
  m <- length(a)
  ans <- numeric(m)
  for (i in 1:m) {
    ans[i] = a[i]^(m-1)/prod(a[i]-a[-i])
  }
  ans
}

centmat <- function(a,b) {
  mat <- matrix(0,length(a),length(b))
  ans <- matrix(a[row(mat)]/(a[row(mat)] + b[col(mat)]),length(a),length(b))
  ans[a[row(mat)]==0] = 0
  return (ans)
}

theprob <- function(a,b) {
  calcCs(a) %*% centmat(a,b) %*% calcCs(b)
}

calcCs2 <- function(a){
    if (length(unique(a)) == length(a)) {
        return (calcCs(a))
    } else {
        o <- order(a)
        a2 <- a[order(a)]
        ua <- unique(sort(a))
        k <-  length(ua)
        for (jj in 1:k) {
            inds <- which(a2==ua[jj])
            m <- length(inds)
            delta <- 1e-3 * ua[k]
            if (jj > 1) delta <- min(c(delta, 1e-3 * (ua[jj] - ua[jj-1])))
            if (jj < k) delta <- min(c(delta, 1e-3 * (ua[jj+1] - ua[jj])))
            a2[inds] <- a2[inds] + delta * ((1:m)/(m+1) - .5)
        }
        a3 <- numeric()
        a3[o] <- a2
        return (calcCs(a3))
    }
}

# theoretical probabilty in case some values a[i] = a[j]
theprob2 <-  function(a,b) {
    calcCs2(a) %*% centmat(a,b) %*% calcCs2(b)
}

domprob <- function(v1,v2) {
    v <- v1 - v2
    if (sum(v != 0) == 0) return(0.5)
    if (sum(v > 0) == 0) return(0)
    if (sum(v < 0) == 0) return(1)
    a <- v[v > 0]
    b <- -v[v < 0]
    return (theprob2(a,b))
}

# In this example I compute the probability that 1*X1+2*X2+3*X3-4*X4-5*X5-6*X6

a <- c(1,2,3)
b <- c(4,5,6)
empprob(a,b) # empirical probability from 10000 samples
theprob(a,b) # theoretical probability

# Calculation of HDI indices

hdi <- read.csv("HDI_formatted.csv",stringsAsFactors = FALSE)
names(hdi)
health_index <- (hdi$life_expectancy - 20)/65
health_index[health_index > 1] <- 1
mean_schooling_index <- (hdi$mean_schooling)/15
mean_schooling_index[mean_schooling_index > 1] <- 1
expected_schooling_index <- (hdi$expected_schooling)/18
expected_schooling_index[expected_schooling_index > 1] <- 1
education_index <- (mean_schooling_index + expected_schooling_index)/2
income_index <- (log(hdi$gni_per_capita) - log(100))/(log(75000)-log(100))
income_index[income_index > 1] <- 1
hdi_calc <-(health_index * education_index * income_index)^(1/3)
hdi_aug <- data.frame(hdi, health_index, education_index, income_index, hdi_calc)
head(hdi_calc)
#plot(hdi_calc, hdi$hdi2014)
bad_inds <- which(abs(hdi_calc - hdi$hdi2014) > 1e-3)
ncountries <- dim(hdi_aug)[1]
rawmat <- cbind(log(health_index), log(education_index), log(income_index))

                                        # compute all pairwise fuzzy probs
nset <- ncountries
fprobs <- matrix(0, nset, nset)
hdidiff <- matrix(0, nset, nset)
for (ii in 1:nset) {
    for (jj in 1:nset) {
        v1 <- rawmat[ii, ]
        v2 <- rawmat[jj, ]
        fprobs[ii,jj] = domprob(v1,v2)
        hdidiff[ii,jj] = hdi$hdi2014[ii] - hdi$hdi2014[jj]
    }
}
names(hdi)
cnames <- paste(hdi$rank, hdi$country, sep=". ")
rownames(fprobs) <- cnames[1:nset]
colnames(fprobs) <- cnames[1:nset]
#write.csv(fprobs, "top10__hdi.csv")

pdf("HDI_plot.pdf")
xx <- as.vector(hdidiff)
yy <- as.vector(fprobs)
filt <- xx > 0
xx <- xx[filt]
yy <- yy[filt]
plot(xx,yy,pch="o",
     col=rgb(0.5,0.5,0.5),cex=.5,xlab="Difference in HDI", ylab="Truth Value")
points(xx,yy,pch="o",col=rgb(0,0,0,0.05),cex=.5)
dev.off()

plot(as.vector(hdidiff[1:50,1:50]),as.vector(fprobs[1:50,1:50]),pch="o",main="Top 50 pairwise dominance",col=rgb(0,0,0,0.5))

list.files()
epi <- read.csv("2014_epi_framework_indicator_scores_friendly.csv",sep=",")
options(width=160)
head(epi)

newnames <- c('Rank','Country','EPI','Ten.Yr.Change','Environmental.Health','Ecosystem.Vitality','EH.Health.Impacts','EH.Air.Quality','EH.Water.Sanitation','EV.Water.Resources','EV.Agriculture','EV.Forests','EV.Fisheries','EV.Biodiversity.Habitat','EV.Climate.Energy','Child.Mortality','Household.Air.Quality','Air.Exposure.PM2.5','Air.PM2.5.Exceedance','Access.to.Sanitation','Access.to.Drinking.Water','Wastewater.Treatment','Agricultural.Subsidies','Pesticide.Regulation','Change.Forest.Cover','Fish.Stocks','Fishing.Pressure','TPA.National.Biome','TPA.Global.Biome','Marine.Protected.Areas','Critical.Habitat.Protection','Trend.Carbon','Change.Trend.Carbon','Trend.CO2.Emissions','Access.to.Electricity')
cbind(names(epi),newnames)
oldnames <- names(epi)
names(epi) <- newnames

head(epi)

names(epi)

prox.epi <- (epi$Environmental.Health * .4 + epi$Ecosystem.Vitality * .6)
plot(epi$EPI, prox.epi)
prox.eh <- (epi$EH.Health.Impacts + epi$EH.Air.Quality + epi$EH.Water.Sanitation)/3
#cbind(epi$Environmental.Health, prox.eh)
nepi <- epi; nepi[is.na(epi)] <- 0
prox.ev <- .25 * epi$EV.Water.Resources + 0.05 * epi$EV.Agriculture + .1 * epi$EV.Forests + .1 * epi$EV.Fisheries + .25 * epi$EV.Biodiversity.Habitat + .25*epi$EV.Climate.Energy
wtot <- 1-(.25 * is.na(epi$EV.Water.Resources) + 0.05 * is.na(epi$EV.Agriculture) + .1 * is.na(epi$EV.Forests) + .1 * is.na(epi$EV.Fisheries) + .25 * is.na(epi$EV.Biodiversity.Habitat) + .25*is.na(epi$EV.Climate.Energy))
prox.ev <- (.25 * nepi$EV.Water.Resources + 0.05 * nepi$EV.Agriculture + .1 * nepi$EV.Forests + .1 * nepi$EV.Fisheries + .25 * nepi$EV.Biodiversity.Habitat + .25*nepi$EV.Climate.Energy)/wtot


# weight-adjusted epi components

wadj <- cbind(.1333*epi$EH.Health.Impacts, .1333*epi$EH.Air.Quality, .1333*epi$EH.Water.Sanitation,
    .15*nepi$EV.Water.Resources/wtot, .03 * nepi$EV.Agriculture/wtot, .06 * nepi$EV.Forests / wtot, .06 * nepi$EV.Fisheries / wtot,
    .15 * nepi$EV.Biodiversity.Habitat / wtot, .15 * nepi$EV.Climate.Energy / wtot)
prox.epi <- apply(wadj,1,sum)
cbind(epi$EPI, prox.epi)


                                        # compute all pairwise EPI
nset <- 10
fprobs <- matrix(0, nset, nset)
epidiff <- matrix(0, nset, nset)
for (ii in 1:nset) {
    for (jj in 1:nset) {
        v1 <- wadj[ii, ]
        v2 <- wadj[jj, ]
        fprobs[ii,jj] = domprob(v1,v2)
        epidiff[ii,jj] = epi$EPI[ii] - epi$EPI[jj]
    }
}
cnames <- paste(epi$Rank, epi$Country, sep=". ")
rownames(fprobs) <- cnames[1:nset]
colnames(fprobs) <- cnames[1:nset]
write.csv(fprobs, "top10_epi.csv")
