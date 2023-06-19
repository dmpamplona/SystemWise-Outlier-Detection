
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
####  Iterative Model-Free Approach in Detectin System Wise-Outliers  ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#### PART 1: SIMULATION STUDY --------------------------------------------

##   PART 1.1: Data Generating Process
# No Outliers
dat.gen <- function(n, tp, var1, var2, seed) {
  y <- vector()
  
  for(i in 1:n) {
    set.seed(seed); lambda <- rnorm(2*tp, 0, var1)
    set.seed((seed*2)*i); a <- rnorm(2*tp, 0, var2)
    y_i <- vector()
    y_i[1] <- a[1]
    
    for(t in 2:(2*tp)) {
      y_i[t] <- 8 + 0.75*y_i[t-1] + lambda[t] + a[t]
    }
    
    len <- length(y_i)
    y_i <- y_i[c((0.5*len+1):len)]
    y <- c(y, y_i)
  }
  
  time <- c(1:tp)
  data <- data.frame(Time = rep(time,n), Series = rep(c(1:n), each=tp), Value = y)
  
  return(data)
}


# Outlier: Variance
dat.gen.out.var <- function(n, tp, var, seed, n.ref) {
  y <- vector()
  
  for(i in 1:n) {
    set.seed(seed*i)
    a <- rnorm(2*tp, 0, var)
    y_i <- vector()
    y_i[1] <- a[1]
    
    for(t in 2:(2*tp)) {
      y_i[t] <- 8 + 0.75*y_i[t-1] + a[t]
    }
    
    len <- length(y_i)
    y_i <- y_i[c((0.5*len+1):len)]
    y <- c(y, y_i)
  }
  
  time <- c(1:tp)
  data <- data.frame(Time = rep(time,n), Series = rep(c((n.ref+1):(n.ref+n)), each=tp), Value = y)
  
  return(data)
}


# Outlier: Trend
dat.gen.out.trend <- function(n, tp, var, seed, n.ref) {
y <- vector()

for(i in 1:n) {
  set.seed(seed*i)
  a <- rnorm(2*tp, 0, var)
  y_i <- vector()
  y_i[1] <- a[1] + 35
  
  for(t in 2:tp) {
    y_i[t] <- y_i[t-1] - 0.9 + a[t]
  }
  y <- c(y, y_i)
}

time <- c(1:tp)
data <- data.frame(Time = rep(time,n), Series = rep(c((n.ref+1):(n.ref+n)), each=tp), Value = y)

return(data)
}


##  PART 1.2: D Measure

D.fun <- function(data, tp) {
  sdi <- vector()
  madi <- vector()
  for(i in 1:tp) {
    val <- subset(data, data$Time == data$Time[i])$Value
    sdi[i] <- sd(val, na.rm=TRUE)  
    madi[i] <- median(abs(val-mean(val, na.rm=TRUE)))
  }
  return(c(mean(sdi), mean(madi)))
}


##  PART 1.3: Monte Carlo Simulation
library(doParallel)
library(foreach)
library(dplyr)
numCores <- detectCores()
registerDoParallel(numCores)


# Set Simulation Scenario
n.sim = 25        # No. of obs   (10, 25, 50)
tp.sim = 25       # No. of time series   (10, 25, 50)
var1.sim = 1      # System variance  (1)   
var2.sim = 3      # Series variance (1, 3, 10)
M = 400           # No. of Monte Carlo runs (400)
s = c(1.5, 1.75, 2, 1.25, 2.5, 2.75, 3)  # IQR Multipliers


## Simulation (No Outliers)
foreach(k=1:length(s)) %dopar% {
  library(dplyr)
  wrong <- vector()
  for(m in 1:M) {
    
    # Part 1. Simulate the Needed Data Setting
    seed.m <- seeds*m*k
    data <- dat.gen(n=n.sim, tp=tp.sim, var1=var1.sim, var2=var2.sim, seed=seed.m)  #### Change every scenario
    data$out <- rep(0, length(data[,1]))
    
    # Part 2. Prepare The Initial Vectors and Data Frames
    
    data.run <- data %>%
      group_by(Series) %>%
      summarize(out.ob = mean(out))
    data.run <- as.data.frame(data.run)
    
    out.pred <- rep(0,length(data.run$Series))
    D.sd <- vector()
    
    # Part 3. Run the Algo per Data
    repeat
    {
      for(i in 1:dim(table(data$Series))) {
        data.i <- subset(data, data[,2] != i) # Leave out Series = i
        result <- D.fun(data=data.i, dim(table(data[,1])))
        D.sd[i] <- result[1]
      }
      
      dat <- data.frame(Series = unique(data$Series), Stat = D.sd)
      dat <- dat[order(dat$Stat, decreasing=TRUE),]
      
      Q1 <- quantile(dat$Stat, probs=0.25)
      lower <- Q1 - s[k]*(IQR(dat$Stat))
      min <- dat[nrow(dat),]$Stat
      
      if(min < lower)
      {
        out.pred[dat[nrow(dat),]$Series] <- 1
        data <- subset(data, data$Series != dat[nrow(dat),]$Series)
        D.sd <- vector()
      } else {
        break
      }
    }
    data.run$out.pred <- out.pred
    wrong[m] <- sum(data.run$out.ob != data.run$out.pred)
  }
  
  100*(sum(wrong)/(M*n.sim))
}


## Simulation (With Outlier)

# Set Simulation Scenario
n.sim.o = 1              # No. of outlier series (1,2,3)
var.sim.o = 3*var2.sim   # Variance of outlier series 

foreach(k=1:5) %dopar% {
  library(dplyr)
  correct <- vector()
  wrong <- vector()
  
  for(m in 1:M) {
    
    # Part 1. Simulate the Needed Data Setting
    seed.m <- seeds*m*k
    y <- dat.gen(n=n.sim, tp=tp.sim, var1=var1.sim, var2=var2.sim, seed=seed.m)
    #y.out <- dat.gen.out.var(n=n.sim.o, tp=tp.sim, var=var.sim.o, seed=(seed.m*2), n.ref=n.sim)  # For High Variance Outlier
    #y.out <- dat.gen.out.trend(n=n.sim.o, tp=tp.sim, var=var.sim.o, seed=(seed.m*2), n.ref=n.sim)  # For Trend Outlier
    data <- data.frame(rbind(y, y.out), Behavior = c(rep("Normal", n.sim*tp.sim), rep("Outlier", n.sim.o*tp.sim)))
    out <- as.numeric(data$Behavior == "Outlier")  # For simulation only
    data$out <- out
    
    
    # Part 2. Prepare The Initial Vectors and Data Frames
    data.run <- data %>%
      group_by(Series) %>%
      summarize(out.ob = mean(out))
    data.run <- as.data.frame(data.run)
    out.pred <- rep(0,length(data.run$Series))
    D.sd <- vector()
    
    
    # Part 3. Run the LOOD Algorithm per Simulated Data
    repeat
    {
      for(i in 1:dim(table(data$Series))) {
        data.i <- subset(data, data[,2] != i) # Leave out Series = i
        result <- D.fun(data=data.i, dim(table(data[,1])))
        D.sd[i] <- result[1]
      }
      
      dat <- data.frame(Series = unique(data$Series), Stat = D.sd)
      dat <- dat[order(dat$Stat, decreasing=TRUE),]
      
      Q1 <- quantile(dat$Stat, probs=0.25)
      lower <- Q1 - s[k]*(IQR(dat$Stat))
      min <- dat[nrow(dat),]$Stat
      
      if(min < lower)
      {
        out.pred[dat[nrow(dat),]$Series] <- 1
        data <- subset(data, data$Series != dat[nrow(dat),]$Series)
        D.sd <- vector()
      } else {
        break
      }
    }
    data.run$out.pred <- out.pred
    wrong[m] <- sum((data.run[data.run$out.ob == 0,])$out.ob != (data.run[data.run$out.ob == 0,])$out.pred)
    correct[m] <- sum((data.run[data.run$out.ob == 1,])$out.ob == (data.run[data.run$out.ob == 1,])$out.pred)
  }
  
  c(100*(sum(correct))/(M*n.sim.o), 100*(sum(wrong)/(M*n.sim)))
  
}



#### PART 2: APPLICATION TO FISHERIES DATA -------------------------------------
library(reshape2)

# Data Preparation
dat <- read.csv(".\\tiger-prawn.csv")    ## Add data here
y <- melt(dat, id.vars="Time")
colnames(y) <- c("Time", "Series", "Value")
Region <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8",
            "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16", "R17")
Time <- c(1:72)

# Standardize the values
std.val <- vector()
for(i in Region) {
  dat <- y[y$Series == i,]
  dat$Value <- (dat$Value - mean(dat$Value, na.rm=TRUE))/sd(dat$Value, na.rm=TRUE)
  std.val <- c(std.val, dat$Value)
}
y$Value <- std.val


# Plot the multiple time series
g <- ggplot(y, aes(x=Time, y=Value, group=Series)) +
  geom_line(color='#207C7E') + theme_minimal() +
  ylab("Species") + xlab("Year")
g


# Apply the Algorithm
data <- y
D.sd <- vector()

repeat
{
  for(i in unique(data$Series)) {
    data.i <- subset(data, data[,2] != i) # Leave out Series = i
    result <- D.fun(data=data.i, dim(table(data[,1])))
    D.sd[i] <- result[1]
  }
  
  dat <- data.frame(Series = unique(data$Series), Stat = D.sd)
  dat <- dat[order(dat$Stat, decreasing=TRUE),]
  
  Q1 <- quantile(dat$Stat, probs=0.25)
  lower <- Q1 - 2*(IQR(dat$Stat))
  min <- dat[nrow(dat),]$Stat
  
  if(min < lower)
  {
    print(dat[nrow(dat),]$Series)
    data <- subset(data, data$Series != dat[nrow(dat),]$Series)
    D.sd <- vector()
  } else {
    break
  }
}


