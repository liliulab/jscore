########################################################
########## simulations used in the manuscript ##########
########################################################

library(jScore)
library(aricode)
library(clevr)

## Given a "total" number of data points, create "k" clusters of different sizes. 
random.size = function(total, k) {
	pp = rbeta(k, shape1=2, shape2=2)
	ss = round(total * pp/sum(pp))
	dd = total - sum(ss)
	ss[length(ss)] = ss[length(ss)] + dd
	cat('total ', sum(ss), ' - ', ss, '\n');
	idx.0 = which(ss <= 0)
	for(idx in idx.0) {
		idx.more = which(ss >= 2)
		idx.more = sample(idx.more, abs(ss[idx]) + 1)
		ss[idx.more] = ss[idx.more] - 1
		ss[idx] = 1
	}
	cat('   total ', sum(ss), ' - ', ss, '\n'); flush.console();
	return(ss)
}
	
## Given clusters of different "sizes", create class labels.
random.simulate.truth = function(sizes) {
	y = c()
	label = c()
	for(i in 1:length(sizes)) {
		this.size = sizes[i]
		y = c(y, rnorm(this.size, mean=i, sd=0.05))
		label = c(label, rep(i, this.size))
	}
	data.df = data.frame(idx=1:length(y), y=y, label=label)
	return(data.df)
}
	
## compute J-score
my.jscore = function(truth, pred) {
	j = NA;
	if(length(truth)==length(pred)){
		est.mat <- c()
		for(i in unique(pred)){
			est.mat<- rbind(est.mat, pred==i)
		}
		est.num <- apply(est.mat,1, sum)

		truth.mat<-c()
		for(i in unique(truth)){
			truth.mat <- rbind(truth.mat, truth==i)
		}
		truth.num<- apply(truth.mat,1, sum)

		est.mat.acc <- matrix(nrow=dim(est.mat)[1], ncol=dim(truth.mat)[1])
		for(i in 1:dim(est.mat)[1]) {
			for(j in 1:dim(truth.mat)[1]) {
				est.mat.acc[i,j] <- sum(est.mat[i,] & truth.mat[j,]) / sum(est.mat[i,] | truth.mat[j,])
			}
		}

		M1 <- sum(apply(est.mat.acc,1,max)*est.num)/length(pred)
		M1.1 <- sum(apply(est.mat.acc,2,max)*truth.num)/length(pred)
		M2 <- 2*M1*M1.1/(M1+M1.1)
		j = M2		
	} else {
		cat('Truth and Pred have different lengths.\n')
	}
	return(j);
}

## compute H-score
h_score = function(truth, pred, reverse=F) {
	h = NA;
	if(length(truth)==length(pred)){
		crosstab = as.matrix.data.frame(table(truth, pred))
		mm = 0
		for(i in 1:nrow(crosstab)) {
			this.pred = crosstab[i, ]
			this.pred.mm = this.pred[-which.max(this.pred)]
			mm = mm + sum(this.pred.mm)
		}
		h = mm / length(truth)
	} else {
		cat('Truth and Pred have different lengths.\n')
	}
	if(reverse) h = 1- h;
	return(h)
}

## compute F-score
f_score = function(truth, pred) {
	ff = NA;
	if(length(truth)==length(pred)){
		crosstab = as.matrix.data.frame(table(truth, pred))
		best.f = c()
		for(i in 1:nrow(crosstab)) {
			f = c()
			for(j in 1:ncol(crosstab)) {
				m = crosstab[i, j]
				precision = m / sum(crosstab[, j])
				recall = m / sum(crosstab[i, ])
				this.f = 2 * precision * recall / (precision + recall)
				f = c(f, this.f)
			}
			best.f = c(best.f, max(f, na.rm=T))
		}
		wt = rowSums(crosstab) / sum(crosstab)
		ff = sum(wt * best.f)
	} else {
		cat('Truth and Pred have different lengths.\n')
	}
	return(ff)
}

## compute H-score, F-score, and J-score
get.h.f.j = function(truth, pred) c(h=h_score(truth, pred), f=f_score(truth, pred), j=my.jscore(truth, pred))

## compute eight different accuracy measures
all_measures = function(truth, pred) {
	jj1 = jscore(truth, pred)
	hh1 = h_score(truth, pred, T)
	aa1 = ARI(truth, pred)
	rr1 = RI(truth, pred)
	nn1 = NMI(truth, pred)
	vv1 = v_measure(truth, pred)
	nv1 = 1 - NVI(truth, pred)
	ff1 = f_score(truth, pred)

	jj2 = jscore(pred, truth)
	hh2 = h_score(pred, truth, T)
	aa2 = ARI(pred, truth)
	rr2 = RI(pred, truth)
	nn2 = NMI(pred, truth)
	vv2 = v_measure(pred, truth)
	nv2 = 1 - NVI(pred, truth)
	ff2 = f_score(pred, truth)
	
	ss = cbind(c(jj1, hh1, aa1, rr1, nn1, vv1, nv1, ff1), c(jj2, hh2, aa2, rr2, nn2, vv2, nv2, ff2))
	return(ss)
}

## plot average and standard deviation of accuracy measures
plot.by.set = function(score.wk, legend.x=5, legend.y=0.95) {
	ccs = c(ARI='orange', F='lightgoldenrod1', H='turquoise1', J='black', NMI='gray60', NVI='purple', RI='limegreen', V='blue')
	score.wk$color = ccs[score.wk$method]
	score.wk$set = paste(score.wk$method, score.wk$k, sep=':')
	score.wk.by.set = split(score.wk, score.wk$set)
	
	plot(0, 0, xlim=c(1, 50), ylim=c(0, 1), xlab='', ylab='Score', type='n', las=1, cex.lab=1.6, cex.axis=1.5); mtext('Number of Clusters', side=1, line=2.5)
	xxx = lapply(score.wk.by.set, function(x) {
		if(x$method[1] == 'NMI') x$score.1 = x$score.1 + 0.02
		avg = mean(x$score.1)
		err = sd(x$score.1)
		points(x$k[1], avg, col=x$color[1], pch=16, cex=0.8)
		arrows(x$k[1], avg-err, x$k[1], avg+err, length=0.02, angle=90, code=3, col=x$color[1])
	})	
	if(!is.na(legend.x) & ! is.na(legend.y))
		legend(legend.x, legend.y, legend=names(ccs), col=ccs, lty=1, bty='n')
}

## Simulations used in Figure 1
sizes = c(10, 30, 60)
data.df = random.simulate.truth(sizes)
plot(y ~ idx, data=data.df, ylim=c(3, 1), type='p', las=1, col=label)

pred = c(rep(1, 10), rep(2, 30), rep(3, 40), rep(4, 20))
get.h.f.j(data.df$label, pred)  ## 0.2   0.88   0.77
pred = c(rep(1, 10), rep(2, 30), rep(3, 40), rep(4, 10), rep(5, 10))
get.h.f.j(data.df$label, pred)  ## 0.2   0.88   0.75
pred = c(rep(1, 7), rep(2, 3), rep(1, 21), rep(2, 9), rep(1, 42), rep(2, 18))
get.h.f.j(data.df$label, pred)  ## 0.3   0.53   0.39
pred = c(rep(1, 7), rep(2, 3), rep(1, 21), rep(2, 9), rep(1, 42), rep(2, 13), rep(3, 5))
get.h.f.j(data.df$label, pred)  ## 0.3   0.53   0.38
		
## Simulations used in Figure 2
k.max = 10
sizes = random.size(total=1000, k=k.max)
data.df.orig = random.simulate(sizes)

ss.all = c()
for(new.k in 1:50) {
	for(iter in 1:200) {
		data.df = data.df.orig
		new.sizes = random.size(total=1000, k=new.k)
		data.df$label.new = unlist(lapply(1:new.k, function(x) rep(x, new.sizes[x])))
		table(data.df$label, data.df$label.new); all_measures(data.df$label, data.df$label.new)
		ss = all_measures(data.df$label, data.df$label.new)
		ss.all = rbind(ss.all, cbind(iter, new.k, ss))	
	}
}
ss.all = data.frame(method=rep(c('J', 'H', 'ARI', 'RI', 'NMI', 'V', 'NVI', 'F'), nrow(ss.all)/8), iter=ss.all[, 1], k=ss.all[, 2], score.1=ss.all[, 3], score.2=ss.all[, 4], stringsAsFactors=F)
plot.by.set(ss.all, legend.x=10, legend.y=0.3)
abline(v=10, lty=2)
	
## Simulations used in Figure 4
total = 1000
k.max = 50
ss.all = c()
for(k in 2:k.max) {
	for(iter in 1:50) {
		sizes = random.size(total=total, k=k)
		truth = unlist(lapply(1:k, function(x) rep(x, sizes[x])))
		table(truth)
		clusters.1 = sample(truth, length(truth))
		table(clusters.1)
		clusters.2 = sample(truth, length(truth))
		table(clusters.2)
		ss = all_measures(clusters.1, clusters.2)
		ss.all = rbind(ss.all, cbind(iter, k, ss))
	}
}
ss.all = data.frame(method=rep(c('J', 'H', 'ARI', 'RI', 'NMI', 'V', 'NVI', 'F'), nrow(ss.all)/8), iter=ss.all[, 1], k=ss.all[, 2], score.1=ss.all[, 3], score.2=ss.all[, 4], stringsAsFactors=F)
plot.by.set(ss.all, legend.x=30, legend.y=0.8)
	
