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
nset <- 10
fprobs <- matrix(0,nset,nset)
for (ii in 1:nset) {
    for (jj in 1:nset) {
        v1 <- rawmat[ii, ]
        v2 <- rawmat[jj, ]
        fprobs[ii,jj] = domprob(v1,v2)
    }
}
names(hdi)
cnames <- paste(hdi$rank, hdi$country, sep=". ")
rownames(fprobs) <- cnames[1:nset]
colnames(fprobs) <- cnames[1:nset]
write.csv(fprobs, "top10__hdi.csv")
