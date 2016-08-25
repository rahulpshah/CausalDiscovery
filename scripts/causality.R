# Load the libraries 
library(vars)
library(urca)
library(stats)
library(pcalg)
# Read the input data 
data_path <- "~/Documents/05.Topic-IV.Project-5.Manufacturer-Retailer-Price.CausalRelationDiscovery.distribution/Input Data/data.csv"
data <- read.csv(data_path)
# Build a VAR model 
#Choosing the best lag value under Schwartz C criterion[3]
lag <- vars::VARselect(data, lag.max = 10, type="const")$selection[3]
var.model <- vars::VAR(data, p = lag, type = "const", ic = "SC")
# Extract the residuals from the VAR model 
qty.residuals <- var.model$varresult$Move$residuals
mprice.residuals <- var.model$varresult$MPRICE$residuals	
rprice.residuals <- var.model$varresult$RPRICE$residuals

# Check for stationarity using the Augmented Dickey-Fuller test 
best_lag <- 6   #Turns out same for all 3
qty.df <- ur.df(qty.residuals, type = "none", lags = best_lag)
mprice.df <- ur.df(mprice.residuals, type = "none", lags = best_lag)
rprice.df <- ur.df(rprice.residuals, type = "none", lags = best_lag)

#Reconfirming the test using adf.test
adf.test(qty.residuals)
adf.test(mprice.residuals)
adf.test(rprice.residuals)

#Since all the test - statistic are negative therefore the there is no unit root of the characteristic polynomial..thus it is stationary series

#Check whether the variables follow a Gaussian distribution  

ks.test(qty.residuals, "pnorm", mean = mean(qty.residuals), sd = sd(qty.residuals))
ks.test(mprice.residuals, "pnorm", mean = mean(mprice.residuals), sd = sd(mprice.residuals))
ks.test(rprice.residuals, "pnorm", mean = mean(rprice.residuals), sd = sd(rprice.residuals))

#Since the p value is so less, the residuals follow non-gaussian distribution

#Write the residuals to a csv file to build causal graphs using Tetrad software

x <- cbind(qty.residuals, mprice.residuals, rprice.residuals)
colnames(x) <- c("qty","mprice","rprice")
output_path <- "~/Documents/05.Topic-IV.Project-5.Manufacturer-Retailer-Price.CausalRelationDiscovery.distribution/tetrad.csv"
write.csv(x, file = output_path, row.names = FALSE)

#PC Algorithm
suffStat <- list(C=cor(x),n=nrow(x))
pc_fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.1, labels = colnames(x), skel.method = "original", verbose = TRUE)
plot(pc_fit, main = "PC Output")
# LiNGAM algorithm
lingam_fit <- LINGAM(x, verbose = TRUE)
show(lingam_fit)
