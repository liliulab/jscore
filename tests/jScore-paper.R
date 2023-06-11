########################################################
######### Analyses presented in the manuscript #########
########################################################

library(jScore)
library(aricode)
library(clevr)
library(ggplot2)
library(gridExtra)
library(rgl)
library(dendextend)

####### custom functions
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

get.h.f.j = function(truth, pred) c(h=h_score(truth, pred), f=f_score(truth, pred), j=my.jscore(truth, pred))

plot.by.method = function(score.wk, xlim=c(1, 50), abline.v=10, legend.x=6, legend.y=0.4) {
	ccs = c(ARI='orange', F='lightgoldenrod1', H='turquoise1', J='black', NMI='gray60', NVI='purple', RI='limegreen', V='blue')
	ltys = c(ARI=1, F=1, H=2, J=1, NMI=1, NVI=2, RI=1, V=1)
	score.wk$color = ccs[score.wk$method]
	score.wk$lty = ltys[score.wk$method]
	score.wk.by.method = split(score.wk, score.wk$method)
	
	plot(0, 0, xlim=xlim, ylim=c(0, 1), xlab='Number of Clusters', ylab='Score', type='n')
	xxx = lapply(score.wk.by.method, function(x) lines(x$k, x$score.1, col=x$color, lty=x$lty, lwd=2))
	abline(v=abline.v)
	if(!is.na(legend.x) & ! is.na(legend.y))
		legend(legend.x, legend.y, legend=toupper(names(score.wk.by.method)), col=ccs[toupper(names(score.wk.by.method))], lty=ltys[toupper(names(score.wk.by.method))], lwd=2, bty='n')
}

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

####### simulation data
{
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
	
	random.split = function(data.df, new.k) {
		data.df$label.new = data.df$label
		while(length(unique(data.df$label.new)) < new.k) {
			label.wk = sample(unique(data.df$label.new), 1)
			data.1 = data.df[which(data.df$label.new == label.wk), ]
			data.2 = data.df[-which(data.df$label.new == label.wk), ]
			# j = sample(1:nrow(data.1), 1)
			j = round(nrow(data.1)/2)
			idx = sample(1:nrow(data.1), j)
			data.1[idx, 'label.new'] = max(data.df$label.new) + 1
			data.df = rbind(data.1, data.2)	
		}
		table(data.df$label.new, data.df$label)
		return(data.df)
	}

	random.merge = function(data.df, new.k) {
		data.df$label.new = data.df$label
		while(length(unique(data.df$label.new)) > new.k) {
			label.wk = sample(unique(data.df$label.new), 2)
			data.1 = data.df[which(data.df$label.new %in% label.wk), ]
			data.1$label.new = max(data.df$label.new) + 1
			data.2 = data.df[-which(data.df$label.new %in% label.wk), ]
			data.df = rbind(data.1, data.2)	
		}
		table(data.df$label.new, data.df$label)
		return(data.df)
	}
	
	random.split.merge = function(total, k, new.k) {
		sizes = random.size(total, k)
		data.df = random.simulate.truth(sizes)
		df.list = split(data.df, data.df$label)
		for(i in 1:length(df.list)) {
			this.df = df.list[[i]]
			this.n = nrow(this.df)
			if(this.n > new.k) {
				ss = random.size(this.n, new.k)
				clusters = unlist(lapply(1:new.k, function(x) rep(x, ss[x])))
				clusters = sample(clusters, length(clusters))
				this.df$new.label = clusters
			} else {
				this.df$new.label = 1:nrow(this.df)
			}
			df.list[[i]] = this.df
		}
		data.df.new = do.call(rbind, df.list)
		return(data.df.new)
	}
	
	## Figure 1
	{
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
	}
		
	## Figure 2
	{	
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
	}

	## Figure 4A
	{
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

	save(ss.all, file='ss.all.baseline.sim.RData')	
	}

}

####### read-world data
{ 
	## iris data
	library(datasets)
	dim(iris)  ## 150   5

	## wine data
	wine.fl <- "http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data"
	wine <- read.csv(wine.fl,header = F)
	# Names of the variables
	wine.names=c("Alcohol", "Malic acid", "Ash", "Alcalinity of ash", "Magnesium", "Total phenols", "Flavanoids", "Nonflavanoid phenols", "Proanthocyanins", "Color intensity", "Hue", "OD280/OD315 of diluted wines", "Proline")
	colnames(wine)[2:14]=wine.names
	colnames(wine)[1]="Class"
	dim(wine)  ## 178  14

	data.list = list()
	data.list[[1]] = list(features.wk=iris[, 1:4], class.wk = iris[, 5])
	data.list[[2]] = list(features.wk = wine[, 2:13], class.wk = wine[, 1])
	names(data.list) = c('iris', 'wine')

	## cluster complete data sets -- Figure 3A & 3D
	algorithms = c('kmeans', 'hclust')  
	par(mfrow=c(1, 2))
	scores.list = list()
	for(i in 1:length(data.list)) {
		set.seed(1)
		features.wk = data.list[[i]]$features.wk; 
		class.wk = data.list[[i]]$class.wk; 
		features.wk = scale(features.wk)

		scores = c()
		if(algorithms[i] == 'kmeans') {
			for(k in 2:10) {
				kclust = kmeans(features.wk, centers=k, nstart=nstart)
				cluster = kclust$cluster
				ss = all_measures(class.wk, cluster)				
				ss = cbind(k, ss)
				scores = rbind(scores, ss)
			}
		} else {
			dist.wk = dist(features.wk)
			hc = hclust(dist.wk, method='ward.D')
			clusters.cut = cutree(hc, k=2:10)
			colnames(clusters.cut) = paste0('k', 2:10)
			clusters.cut = data.frame(class=class.wk, clusters.cut)			
			for(j in 2:10) {
				k = j
				cluster = clusters.cut[, j]				
				ss = all_measures(class.wk, cluster)				
				ss = cbind(k, ss)
				scores = rbind(scores, ss)
			}
		}
		scores = data.frame(method=rep(c('J', 'H', 'ARI', 'RI', 'NMI', 'V', 'NVI', 'F'), nrow(scores)/8), k=scores[, 1], score.1=scores[, 2], score.2=scores[, 3], stringsAsFactors=F)
		plot.by.method(scores, xlim=c(2, 10), abline.v=3)
		scores.list[[i]] = scores
	}
	reshape(scores.list[[1]][, c('method', 'k', 'score.1')], idvar='k', timevar='method', direction='wide')
		# k score.1.J score.1.H score.1.ARI score.1.RI score.1.NMI score.1.V score.1.NVI score.1.F
	# 1   2 0.6666667 1.0000000   0.5681159  0.7762864   0.5793802 0.7336804   0.5793802 0.7777778
	# 9   3 0.7333717 0.8333333   0.6201352  0.8322148   0.6591265 0.6594869   0.4919660 0.8331832
	# 17  4 0.5667050 0.6666667   0.4705678  0.7762864   0.5450854 0.5967060   0.4252181 0.7220721
	# 25  5 0.4728118 0.5466667   0.4204810  0.7668904   0.5020087 0.5898618   0.4183007 0.6589452
	# 33  6 0.4228007 0.5000000   0.3671626  0.7492617   0.4637074 0.5625626   0.3913649 0.6168120
	# 41  7 0.3548185 0.4200000   0.2965167  0.7299329   0.4083717 0.5141758   0.3460542 0.5443632
	# 49  8 0.3610387 0.4000000   0.3311016  0.7508725   0.4265020 0.5549993   0.3840824 0.5610329
	# 57  9 0.3996762 0.4600000   0.3608573  0.7604474   0.4201334 0.5471510   0.3766055 0.6092628
	# 65 10 0.3335616 0.3733333   0.3284325  0.7561521   0.4223077 0.5692282   0.3978469 0.5423741
	reshape(scores.list[[2]][, c('method', 'k', 'score.1')], idvar='k', timevar='method', direction='wide')
		# k score.1.J score.1.H score.1.ARI score.1.RI score.1.NMI score.1.V score.1.NVI score.1.F
	# 1   2 0.6293717 0.9831461   0.4765925  0.7232273   0.4872605 0.6280956   0.4578275 0.7457485
	# 9   3 0.8745795 0.9325843   0.8000727  0.9104932   0.7774162 0.7775955   0.6361197 0.9321892
	# 17  4 0.8010674 0.8651685   0.7309463  0.8845934   0.6588583 0.7182416   0.5603564 0.9031044
	# 25  5 0.6976052 0.7584270   0.6280876  0.8459976   0.5743232 0.6671946   0.5005942 0.8227116
	# 33  6 0.5765206 0.6404494   0.5363103  0.8140037   0.5141904 0.6273796   0.4570671 0.7479100
	# 41  7 0.5737475 0.6404494   0.5330024  0.8130515   0.5045281 0.6216243   0.4509832 0.7479100
	# 49  8 0.4488966 0.5112360   0.4136259  0.7739478   0.4506699 0.5790146   0.4074740 0.6588570
	# 57  9 0.4134556 0.4775281   0.3841331  0.7653780   0.4265182 0.5591445   0.3880642 0.6363099
	# 65 10 0.4115977 0.4775281   0.3824662  0.7651241   0.4185427 0.5522539   0.3814577 0.6363099
	
	## cluster subsets of samples -- Figure 3B, 3C, 3E & 3F.
	N = 100; pct = 0.9; seeds = 1:N
	hc.list = list()
	scores.list = list()
	clusters.list = list()
	for(i in 1:length(data.list)) {
		cat('data list element:', i, '    algorithm:', algorithms[i], '\n'); flush.console()
		scores.all = c()
		clusters.all = c()
		for(d in 1:length(seeds)) {
			if(d %% 10 == 0) cat('      processing seeds', d, '\n'); flush.console()
			features.wk = data.list[[i]]$features.wk
			class.wk = data.list[[i]]$class.wk
			this.seed = seeds[d]
			set.seed(this.seed)
			M = round(length(class.wk) * pct)
			idx.sub = sample(1:length(class.wk), M)
			class.wk = class.wk[idx.sub]; names(class.wk) = NULL
			features.wk = features.wk[idx.sub, ]; rownames(features.wk) = NULL
			features.wk = scale(features.wk)

			scores = c()
			clusters = c()
			if(algorithms[i] == 'kmeans') {
				for(k in 2:10) {
					kclust = kmeans(features.wk, centers=k, nstart=nstart)
					cluster = kclust$cluster
					ss = all_measures(class.wk, cluster)				
					ss = cbind(k, ss)
					scores = rbind(scores, ss)
					clusters = rbind(clusters, c(k, cluster))
				}
			} else {
				dist.wk = dist(features.wk)
				hc = hclust(dist.wk, method='ward.D')
				clusters.cut = cutree(hc, k=2:10)
				colnames(clusters.cut) = paste0('k', 2:10)
				clusters.cut = data.frame(class=class.wk, clusters.cut)			
				for(j in 2:10) {
					k = j
					cluster = clusters.cut[, j]				
					ss = all_measures(class.wk, cluster)				
					ss = cbind(k, ss)
					scores = rbind(scores, ss)
					clusters = rbind(clusters, c(k, cluster))
				}
			}
			scores = data.frame(method=rep(c('J', 'H', 'ARI', 'RI', 'NMI', 'V', 'NVI', 'F'), nrow(scores)/8), seed=this.seed, k=scores[, 1], score.1=scores[, 2], score.2=scores[, 3], stringsAsFactors=F)
			scores.all = rbind(scores.all, scores)
			clusters = cbind(this.seed, clusters)
			clusters.all = rbind(clusters.all, clusters)				
		}
		scores.list[[length(scores.list) + 1]] = scores.all
		clusters.list[[length(clusters.list) + 1]] = clusters.all
	}
	names(scores.list) = names(clusters.list) = c('iris.kmeans', 'wine.hclust')
	do.call(rbind, lapply(scores.list, dim))
	do.call(rbind, lapply(clusters.list, dim))

	best.k.list = list()
	cnt.eq.3 = c()
	cnt.dev = c()
	for(a in 1:length(scores.list)) {
		scores.all = scores.list[[a]]
		list.by.method = split(scores.all, scores.all$method)
		best.k = c()
		for(i in 1:length(list.by.method)) {
			by.this.method = list.by.method[[i]]
			list.by.seed = split(by.this.method, by.this.method$seed)
			for(j in 1:length(list.by.seed)) {
				by.this.seed = list.by.seed[[j]]
				best = by.this.seed[which.max(by.this.seed$score.1), c('seed', 'k')]
				best.k = rbind(best.k, c(i, j, as.numeric(best)))
			}
		}
		best.k = data.frame(method=names(list.by.method)[best.k[, 1]], j=best.k[, 2], seed=best.k[, 3], k=best.k[, 4], stringsAsFactors=F)
		best.k.list[[a]] = best.k
		cnt.eq.3 = rbind(cnt.eq.3, unlist(lapply(split(best.k, best.k$method), function(x) length(which(x$k == 3)))))
		cnt.dev = rbind(cnt.dev, unlist(lapply(split(best.k, best.k$method), function(x) sum(abs(x$k - 3)))))
	}
	colnames(cnt.eq.3) = colnames(cnt.dev) = names(list.by.method)
	cnt.eq.3 = data.frame(method=names(scores.list), cnt.eq.3, stringsAsFactors=F  )
	cnt.dev = data.frame(method=names(scores.list), cnt.dev, stringsAsFactors=F)
	cnt.eq.3
		   # method ARI   F H   J NMI NVI RI   V
	# 1 iris.kmeans  98 100 0 100 100   0 94   0
	# 2 wine.hclust  98  98 0 100 100 100 95 100
	cnt.dev
		   # method ARI F   H J NMI NVI RI   V
	# 1 iris.kmeans   3 0 100 0   0 100  7 100
	# 2 wine.hclust   2 2 100 0   0   0  5   0

	## Figure 3C & 3F
	best.k = best.k.list[[1]]   ## iris kmeans
	(tb = table(best.k$method, best.k$k))
	    # 2   3   4   5
	  # ARI   0  98   1   1
	  # F     0 100   0   0
	  # H   100   0   0   0
	  # J     0 100   0   0
	  # NMI   0 100   0   0
	  # NVI 100   0   0   0
	  # RI    0  94   5   1
	  # V   100   0   0   0
	best.k = best.k.list[[2]]   ## wine hclust
	(tb = table(best.k$method, best.k$k))
			# 2   3   4
	  # ARI   0  98   2
	  # F     0  98   2
	  # H   100   0   0
	  # J     0 100   0
	  # NMI   0 100   0
	  # NVI   0 100   0
	  # RI    0  95   5
	  # V     0 100   0

	## random clustering -- Figure 4B & 4C
	i = 1;
	class.wk = data.list[[i]]$class.wk
	features.wk = data.list[[i]]$features.wk
	features.wk = scale(features.wk)

	total = length(class.wk)
	k.max = 50
	ss.all = c()
	for(k in 2:k.max) {
		for(iter in 1:50) {
			sizes = random.size(total=total, k=k)
			clusters = unlist(lapply(1:k, function(x) rep(x, sizes[x])))
			clusters = sample(clusters, length(clusters))
			ss = all_measures(class.wk, clusters)
			ss.all = rbind(ss.all, cbind(iter, k, ss))
		}
	}
	ss.all = data.frame(method=rep(c('J', 'H', 'ARI', 'RI', 'NMI', 'V', 'NVI', 'F'), nrow(ss.all)/8), iter=ss.all[, 1], k=ss.all[, 2], score.1=ss.all[, 3], score.2=ss.all[, 4], stringsAsFactors=F)
	save(ss.all, file='ss.all.baseline.iris.RData')
	save(ss.all, file='ss.all.baseline.wine.RData')
	
	ss.files = c('ss.all.baseline.sim.RData', 'ss.all.baseline.iris.RData', 'ss.all.baseline.wine.RData')
	par(mfrow=c(1, 3), mar=c(4, 5, 4, 3))
	ccs = c(ARI='orange', F='lightgoldenrod1', H='turquoise1', J='black', NMI='gray60', NVI='purple', RI='limegreen', V='blue')
	for(f in ss.files) {
		load(f)
		plot.by.set(ss.all, legend.x=NA, legend.y=NA)
	}
	
}
	
