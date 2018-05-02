# For whatever model you're working with, load parameter values from your saved optimization outputs
# DAHe1He2 predicted
He1m = lk.DAHe1He2$par['He1m']
He2m = lk.DAHe1He2$par['He2m']
Ac1 = lk.DAHe1He2$par['Ac']
Ae1 = lk.DAHe1He2$par['Ae']
wmm = (master.g1p_11)
wmo = (master.g1p_10)
wom = (master.g1p_01)
woo = (master.g1p_00)
wnn = master.naive
uc = Uc.master
ua = Ua.master
ue = Ue.master

# DAH predicted

Hm = lk.DAH$par['Hm']
Ac = lk.DAH$par['Ac']
Ae = lk.DAH$par['Ae']

# Calculate model predictions
# Use the equation from the relevant likeklihood function for pp
DAHe1He2_predicted = p0.master*(He1m*He2m*wmm + He1m*wmo + He2m*wom + woo + wnn)*(Ac*uc+ua+Ae*ue)

DAH_predicted = p0.master*(Hm*(weights.master.1+weights.master.2)+weights.master.3+weights.master.naiive)*(Ac*Uc.master+Ua.master+Ae*Ue.master)

DHe1He2_predicted = p0.master*(He1m*He2m*wmm + He1m*wmo + He2m*wom + woo + wnn)

# Plot observed data
xx = barplot(colSums(data.master), ylab = "cases", xlab = 'birth year', main = "observed vs. predicted", col = 'gray', border = 'gray') # This gives you a barplot of cases per birth year (aggregate over all rows)

# Overlay model predictions
## DAH_predicted/rowSums(DAH_predicted) This normalizes the predicted values to give the fraction of cases predicted in each birth year
## Finally we have to multiply by the total number of cases observed in the row to get the model prediction.
lines(xx, colSums(DAH_predicted/rowSums(DAH_predicted)*rowSums(data.master)), col = 'red')
lines(xx, colSums(DAHe1He2_predicted/rowSums(DAHe1He2_predicted)*rowSums(data.master)), col = 'blue')

lines(xx, colSums(DHe1He2_predicted/rowSums(DHe1He2_predicted)*rowSums(data.master)), col = 'green')

################ Here we attempt to fix our likelihood
wmn = weights.master.1 + weights.master.2 - wmm - wmo

lk.DAHe1He2 = optim(par = c(He1m = .2, He2m = .2,  Ac = 2, Ae = 2), fn = lk.D.A.He1.He2, wmm = (master.g1p_11), wmo = (master.g1p_10), wom = (master.g1p_01), woo = (master.g1p_00), wnn = master.naive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 1, 1), upper = c(1, 1, 100, 100)); lk.DAHe1He2
