par(mfrow = c(2,1), mar = c(2.5, 2, 1, 1))
plot.new()  ## Set up a 2-panel plot and initialize the top plot
xx = 0:97
## This part initializes the top plot and plots the timeline using polygons
plot.window(xlim = c(0, 97), ylim = c(0, 1))
axis(1, at = xx[which(1918:2015 %in% rev(c(2015, 2009, 1977, 1968, 1957, 1918)))], labels = FALSE, adj = 0, line = -.19)
mtext(rev(c(2009, 1968, 1957, 1918)), side = 1, line = 1, at = xx[which(1918:2015 %in% rev(c(2009, 1968, 1957, 1918)))], cex = 1.15, adj = .75)
mtext(c(1977,2017), side = 1, line = 1, at = xx[which(1918:2017 %in% c(2015, 1977))], cex = 1.15, adj = .25)
pandemics = xx[which(1918:2017 %in% rev(c(2015, 2009, 1977, 1968, 1957, 1918)))]
cols = c('dodgerblue', 'lightskyblue', 'tomato', 'dodgerblue')
#1918:1957 H1N1
polygon(x = c(pandemics[1:2], pandemics[2:1]), y = rep(c(0, 1), each = 2), col = 'dodgerblue', border = 'dodgerblue')
text(xx[which(1918:2017 == 1940)]+.5, .5, 'H1N1', cex= 2, col = 'black')
#1957:1968 H2N2
polygon(x = c(pandemics[2:3], pandemics[3:2]), y = rep(c(0, 1), each = 2), col = 'lightskyblue', border = 'lightskyblue')
text(xx[which(1918:2017 == 1962)]+.5, .5, 'H2N2', cex= 1, col = 'black')
#1968:1977 H3N2
polygon(x = c(pandemics[3:4], pandemics[4:3]), y = rep(c(0, 1), each = 2), col = 'tomato', border = 'tomato')
#1977:2015 H3N2
polygon(x = c(pandemics[c(4,6)], pandemics[c(6,4)]), y = rep(c(.5, 1), each = 2), col = 'tomato', border = 'tomato')
text(xx[which(1918:2017 == 1990)]+.5, .75, 'H3N2', cex= 1.6, col = 'black')
#1977:2009 H1N1
polygon(x = c(pandemics[c(4,6)], pandemics[c(6,4)]), y = rep(c(0, .5), each = 2), col = 'dodgerblue', border = 'dodgerblue')
text(xx[which(1918:2017 == 1997)]+.5, .25, 'H1N1', cex= 1.6, col = 'black')
#2009:2015 H1N1
# polygon(x = c(pandemics[c(5, 6)], pandemics[c(6,5)]), y = rep(c(0, .5), each = 2), col = 'dodgerblue', border = 'dodgerblue')
# text(xx[which(1918:2015 == 2012)], .25, 'H1N1\np09', cex= .75, col = 'black')
# A label
mtext(text = 'a', side = 3, at = -2, cex = 1.5, font = 2, line = -.3)
#segments(x0 = xx[which(1918:2015 == 1968)], y0 = 0, y1 = 1, lty = 2, col = 'black')
################################################################

## This part extracts some arbitrary row from weights.master.WHATEVER, stacks it up with all the rows from the other matrices and plots it as a stacked barplot. You could do this with all 9 or 12 or whatever of your cases!
par(mar = c(2.5, 2, 2, 1))
## B. china 2015 exposure weights
xx = barplot(rbind((master.g1p_11[13,as.character(1918:2015)]), (master.g1p_1n[13,as.character(1918:2015)]), (master.g1p_01[13,as.character(1918:2015)]), (master.g1p_10[13,as.character(1918:2015)]), (master.g1p_00[13,as.character(1918:2015)]), (master.g1p_0n[13,as.character(1918:2015)]), (master.p_nn[13,as.character(1918:2015)])), col = c('dodgerblue', 'lightskyblue', 'forestgreen', 'darkolivegreen2','tomato', 'lightpink', 'grey'), xlab = 'birth year', yaxt = 'n', xaxt = 'n', border = NA, main = 'Initial IAV exposures by birth year, China 2017', ylab = 'Population fraction')
axis(2, at = c(0, .5, 1), line = -1.2, cex = 1.15)
axis(1, at = xx[which(1918:2015 %in% rev(c(2015, 2009, 1977, 1968, 1957, 1918)))], labels = FALSE, adj = 0)
mtext(rev(c(2009, 1968, 1957, 1918)), side = 1, line = 1, at = xx[which(1918:2017 %in% rev(c(2009, 1968, 1957, 1918)))], cex = 1.15, adj = .75)
mtext(c(1977,2017), side = 1, line = 1, at = xx[which(1918:2015 %in% c(2015, 1977))], cex = 1.15, adj = .25)
#legend(2, .98, c('H1N1', 'H2N2', 'H3N2', 'naive'), fill =  c('dodgerblue','lightskyblue', 'tomato',  'grey'), bg = 'black', cex = 1.1)
# B label
#mtext(text = 'b', side = 3, at = -2, cex = 1.5, font = 2, line = .3)
#par(xpd=TRUE)
#legend('right', legend = c('Group 1 mm', 'Group 1 mn', 'Group 1 om', 'Group 1 mo', 'Group 1 oo', 'Group 1 om'), pch = c(NA, NA, NA, NA, NA, NA, NA), lty = c(1, 1, 1, 1, 1, 1, 1), col = c('dodgerblue', 'lightskyblue', 'forestgreen', 'darkolivegreen2','tomato', 'lightpink', 'grey'))


#####################################################################
## Let's say the optimum parameter estimates for model DAH are stored in the variable lk.DAH:
# DAH predicted
nested_wmm = master.g1p_11
nested_wmo = master.g1p_10 + master.g1p_1n
nested_wom = master.g1p_01
nested_woo = master.g1p_00 + master.g1p_0n
nested_wnn = master.p_nn

DAH_prediction = p0.master*(lk.DAH$par['Hm']*(weights.master.1+weights.master.2)+weights.master.3+weights.master.naiive)*(lk.DAH$par['Ac']*Uc.master+Ua.master+lk.DAH$par['Ae']*Ue.master)

#DAHe1He2_prediction = p0.master*(lk.DAHe1He2$par['He1m']*lk.DAHe1He2$par['He2m']*nested_wmm + lk.DAHe1He2$par['He1m']*nested_wmo + lk.DAHe1He2$par['He2m']*nested_wom + nested_woo + nested_wnn)*(lk.DAHe1He2$par['Ac']*Uc.master+Ua.master+lk.DAHe1He2$par['Ae']*Ue.master)

DAHep_prediction = p0.master*(lk.DAHmep$par['Hm']*lk.DAHmep$par['e1']*lk.DAHmep$par['e2']*nested_wmm + lk.DAHmep$par['Hm']*lk.DAHmep$par['e1']*nested_wmo + lk.DAHmep$par['Hm']*nested_wom + nested_woo + nested_wnn)*(lk.DAHmep$par['Ac']*Uc.master+Ua.master+lk.DAHmep$par['Ae']*Ue.master)

# You should be able to use this to figure out how to make predictions for other models on your own!
#test_prediction = p0.master*(.5*lk.DAHe1He2$par['He2m']*nested_wmm + .5*nested_wmo + lk.DAHe1He2$par['He2m']*nested_wom + nested_woo + nested_wnn)*(lk.DAHe1He2$par['Ac']*Uc.master+Ua.master+lk.DAHe1He2$par['Ae']*Ue.master)



## save as a .pdf
pdf('model_predictions.pdf')
# Plot the total case counts per birth year in the observed data using a bar plot
xx = barplot(colSums(data.master), ylab = "cases", xlab = 'birth year', main = "H5N1 observed vs. predicted", col = 'gray', border = 'gray')
# Then overlay the model prediction. Dividing by row sums normalizes to give the fraction of cases predicted in each birth year for each country-year. The multiply by the row sums of data master, which gives the total number of observed cases in each country-year. This scales the prediction so that the total number of predicted cases matches the total number of observed cases.
lines(xx, colSums(DAH_prediction/rowSums(DAH_prediction)*rowSums(data.master)), col = 'darkred', lwd  = 1.5)
#lines(xx, colSums(DAHe1He2_prediction/rowSums(DAHe1He2_prediction)*rowSums(data.master)), col = 'blue', lwd  = 1.5)
lines(xx, colSums(DAHep_prediction/rowSums(DAHep_prediction)*rowSums(data.master)), col = 'blue', lwd  = 1.5)
# You can basically duplicate the above line, but plot predictions for other models that you made above to add several predictions to the same plot.
#lines(xx, colSums(test_prediction/rowSums(test_prediction)*rowSums(data.master)), col = 'orange', lwd  = 1.5)
# Sample legend
legend('topright', legend = c('data', 'DAHep', 'Old Best DAH'), pch = c(15, NA, NA), lty = c(NA, 1, 1), col = c('gray', 'blue', 'darkred'))
## 
dev.off()

