# Clear workspace and close all plots
rm(list = ls())
graphics.off()

library(xlsx)

# Read scenario
cat('Select scenario Excel file\n')
filename <- 'hoosfield_scenario1.xlsx'
tab <- read.table(filename, header = FALSE, stringsAsFactors = FALSE)

outname <- 'hoosfield_scenario1'
dataname <- 'DATA_hoosfield\\monthly_weather_and_input_hoosfield_scenario1.xlsx'
parname <- 'DATA_hoosfield\\hoosfield_parameters.xlsx'

# Create output folder
if (!dir.exists(outname)) {
  dir.create(paste('OUTPUT_', outname, sep = ''))
}

# Read files
datatab <- read.xlsx(dataname, sheetIndex=1, header=TRUE)
partab <- read.xlsx(parname,sheetIndex=1, header = FALSE)

# Global parameters
k <- partab$X3[1:4]
r <- partab$X3[5]
eta <- partab$X3[6]
clay <- partab$X3[7]
d <- partab$X3[8]
meantemp <- partab$X3[9]

# Initial SOC content
C0 <- partab$X3[10:13]
IOM <- partab$X3[14]

# Read data
years <- datatab[1]
months <- datatab[2]
temp <- datatab[3]
rain <- datatab[4]
pet <- datatab[5]
g <- datatab[6]
f <- datatab[7]
cover <- datatab[8]

# Settings
gamma <- r/(r + 1)
avecg <- c(gamma, 1 - gamma, 0, 0)
avecf <- c(eta, eta, 0, 1 - 2 * eta)

b <- g * avecg + f * avecf
b <- t(b)

x <- 1.67 * (1.85 + 1.6 * exp(-0.0786 * clay))
biohumfac <- 1 / (x + 1)
alpha <- 0.46 * biohumfac
beta <- 0.54 * biohumfac

# Time parameters
t0 <- 1
T <- length(years)
tspan <- c(t0, T)

# Rate modifying factors
source("rhofun.R")  # Assuming rhofun.R contains the rhofun function

#ka <- kb <- kc <- acc <- rep(0, T)  # initialize variables
  result <- rhofun(as.vector(temp), as.vector(rain), as.vector(pet), as.vector(cover), clay, d, meantemp)
  ka <- result$a
  kb <- result$b
  kc <- result$c
  acc <- result$accTSMDD


rho <- ka * kb * kc

# Save rate modifying factors table
rho_table <- data.frame(years, months, ka, kb, kc, acc_TSMD = acc)
write.table(rho_table, file = sprintf('OUTPUT_%s/RHO_%s.csv', outname, outname), sep = ",", row.names = FALSE)

# Run
gamma_matrix <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, alpha, alpha, alpha, alpha, beta, beta, beta, beta), nrow = 4)
A <- -diag(k) %*% (diag(4) - gamma_matrix)
M <- diag(4) - gamma_matrix
Minv <- solve(M)
At <- -M %*% diag(k) %*% Minv
Atinv <- -M %*% diag(1 / k) %*% Minv

Cout <- matrix(0, nrow = length(C0), ncol = T)
Cout[, 1] <- C0

for (n in seq_along(tout)[-1]) {
  if (rho[n] == 0) {
    C0 <- C0 + b[, n]
  } else {
    F <- Atinv %*% (expm(rho[n] * At) - diag(4)) / rho[n]
    C0 <- C0 + F %*% (rho[n] * A %*% C0 + b[, n])
  }
  Cout[, n] <- C0
}

tout <- t0:T
Cout <- t(Cout)

soc <- rowSums(Cout) + IOM

# Save SOC table
out_table <- data.frame(years, months, Cout, IOM = rep(IOM, T), SOC = soc)
colnames(out_table) <- c('year', 'month', 'DPM', 'RPM', 'BIO', 'HUM', 'IOM', 'SOC')
write.table(out_table, file = sprintf('OUTPUT_%s/SOC_%s.csv', outname, outname), sep = ",", row.names = FALSE)

# Plot SOC
plot(years, soc, type = 'l', col = 'blue', lwd = 2, xlab = 'Year', ylab = 'SOC', main = outname)
dev.copy(png, filename = sprintf('OUTPUT_%s/SOC_%s.png', outname, outname))
dev.off()

cat(sprintf('Simulation ended with success. Results have been saved in the folder %s\n', outname))
