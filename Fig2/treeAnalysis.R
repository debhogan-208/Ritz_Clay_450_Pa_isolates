# Kenneth B. Hoehn
# 2/10/26
# Analyze BEAST output from P aeruginosa isolates
# Perform root-to-tip regression on genetic distance tree

library(dplyr)
library(tidyr)
library(ggtree)
library(ape)
library(ggtree)
library(treeio)
library(qgraph)

print(sessionInfo())

# recurse down tree and tally up parent/child relationships
getSwitchesLoc = function(node, data, edges, min, output=list()){
	children = edges[edges[,1] == node,2]
	if(length(children) > 0){
		for(ch in children){
			datap = filter(data, !!node==node)
			datach = filter(data, node==ch)
			output = bind_rows(output, tibble(p=node, ch=ch,
				p_loc=datap$location, 
				ch_loc=datach$location))
			output = bind_rows(getSwitchesLoc(ch, data, edges, min, output))
		}
	}
	return(output)
}

# read in TreeAnnotator output
t = read.beast("relaxed_2.6/lobe_tree_with_trait_combined_100000.tree")

# map cluster IDs to tips on tree
t@data$cluster = ""
timekey = c("T1"="128","T2"="279","T3"="613")

clusters = c()
for(cs in 1:5){
	ct = readLines(paste0("clusters/cluster",cs,".txt"))
	ct = gsub("_"," ",ct)
	tt = sapply(strsplit(ct, split=" "), function(x)x[1])
	times = timekey[tt]
	lobe = substr(ct, 5, 5)
	ct = gsub(" ","_",ct)
	ct = paste0(ct, "-", times)
	ct = paste0(ct, "-", lobe)
	temp = rep(cs, length=length(ct))
	names(temp) = ct

	m = match(ct, t@phylo$tip.label)
	for(node in m){
		t@data$cluster[t@data$node == node] = cs
	}
	clusters = c(clusters, temp)
}
clusterkey = 21:26
names(clusterkey) = as.character(1:5)

# Time tree figure
pdf("relaxed_2.6/results/tree_figure_new.pdf", height=9,width=3)
print(ggtree(t) + 
geom_nodepoint(aes(fill=location,alpha=as.numeric(location.prob)),pch=21,size=1.5)+
geom_tippoint(aes(fill=location,alpha=as.numeric(location.prob)),pch=21,size=1.5) + 
labs(alpha="Prob",fill="Lobe")+
scale_fill_brewer(palette="Dark2")+
geom_treescale(width=365))

print(ggtree(t) +
geom_tippoint(aes(fill=cluster),pch=21,size=1.5) + 
labs(fill="Lobe")+
geom_treescale(width=365))
dev.off()

# Tally switches across the posterior
t = readLines("relaxed_2.6/lobe_tree_with_trait_combined_100000.trees")
t = t[grepl("tree STATE_",t)]
print(paste("States in T1T2T3:",length(t)))

#Tally parent/child locations across all trees sampled from posterior
#10% burnin already done by treeannotator
sample = 1:length(t)
switchl = parallel::mclapply(sample, function(x){
	if(x %% 10 == 0){print(x)}
	ts = strsplit(t[x], split="\\s+")[[1]]
	ts = ts[length(ts)]
	ts = read.beast.newick(
             textConnection(ts))

	edges = ts@phylo$edge
	data = ts@data
	node = ape::getMRCA(ts@phylo, tip=ts@phylo$tip.label)
	spl = getSwitchesLoc(node, data, edges, 0)
	spl$rep = x
	spl = filter(filter(spl, p_loc != ch_loc))
	spl
}, mc.cores=8)
switches = bind_rows(switchl)

write.csv(switches, "relaxed_2.6/results/switches_posterior.csv")
switches = read.csv("relaxed_2.6/results/switches_posterior.csv")

# Make graph showing mean switches among lobes
reps = n_distinct(switches$rep)
network = switches %>%
	filter(p_loc != ch_loc) %>%
	mutate(type=paste0(p_loc,"->",ch_loc)) %>%
	group_by(type) %>%
	summarize(thickness = n()/reps)

nrow(filter(switches, p_loc != ch_loc))/reps

# mean number of switches per tree
mean_switches =	switches %>% 
	group_by(rep) %>% 
	summarize(switches = sum(p_loc != ch_loc)) %>%
	ungroup() %>%
	summarize(mean = mean(switches))
# 32.9

network$from = sapply(strsplit(network$type, split="-"),function(x)x[1])
network$to = sapply(strsplit(network$type, split="-"),function(x)gsub(">","",x[2]))

lm = as.matrix(rbind(
	c(0,3,0),
	c(1,0,2)))

fn = network[,c("from","to","thickness")]
pdf("relaxed_2.6/results/relaxed_network.pdf",width=2,height=2)
qgraph(fn,esize=30,asize=15,fade=FALSE, label.cex=3,
	node.width=2, node.height=2, arrowAngle=pi/6, layout=lm)
dev.off()

# Quantify 95% HPD intervals and output table of switches
network_ci = switches %>%
	filter(p_loc != ch_loc) %>%
	mutate(type=paste0(p_loc,"->",ch_loc)) %>%
	group_by(type, rep, p_loc, ch_loc) %>%
	summarize(switches = n()) %>%
	group_by(type, p_loc, ch_loc) %>%
	summarize(mean=mean(switches),
		lci = quantile(switches, probs=c(0.025)),
		uci = quantile(switches, probs=c(0.975)),
		)
write.csv(network_ci, "relaxed_2.6/results/switches_posterior_ci.csv")

# years before T1 MRCA (numbers from tracer)
abs(613-1986.236-128)/365
abs(613-1131.447-128)/365
abs(613-3105.54-128)/365


# Genetic distance tree analysis
# read in Newick tree
n = read.tree("genetic_distance/AZ4X3umPimnSn-KnIqDNGQ_newick.txt")

# get divergence and time in a big dataframe
tips = n$tip.label
times = as.numeric(gsub("T","",sapply(strsplit(tips, split="_"), function(x)x[1])))
days = c(128, 279, 613)[times]
lobes = gsub("^R|L$","",sapply(strsplit(tips, split="_"), function(x)x[2]))
mrca = getMRCA(n, tip=tips[tips != "ERR364087"])

dists = dist.nodes(n)

info = bind_cols(node=1:length(tips),taxa=tips,
	time=days,lobe=lobes)

info$div = dists[mrca,info$node]
dinfo = filter(info, taxa != "ERR364087")

l = summary(lm(div ~ time, data=dinfo))$coefficients

obs = cor(dinfo$div, dinfo$time)
perm = c()
for(i in 1:10000){
	temp = dinfo
	temp$time = sample(temp$time, replace=FALSE)
	perm = c(perm, cor(temp$div, temp$time))
}

# P value
sum(perm >= obs)
# 0

# label each tip by cluster
clusters = c()
for(cs in 1:5){
	ct = readLines(paste0("clusters/cluster",cs,".txt"))
	ct = gsub(" ","_",ct)
	temp = rep(cs, length=length(ct))
	names(temp) = ct
	clusters = c(clusters, temp)
}

dinfo$cluster = as.character(clusters[dinfo$taxa])
clusterkey = 21:25
names(clusterkey) = as.character(1:5)

shapes = c("U"=24,"M"=21,"L"=25)

# root to tip regression plot
pdf("genetic_distance/scatterplot.pdf", width=4,height=3)
ggplot(dinfo, aes(x=time, y=div, pch=cluster, fill=lobe)) +
geom_jitter(width=10,height=0) + theme_bw() + 
scale_fill_brewer(palette="Dark2") +
scale_shape_manual(values=clusterkey) +
geom_abline(intercept=l[1,1], slope=l[2,1]) +
ylab("Divergence from MRCA") + xlab("Sample day") +
ggtitle(paste("Slope =", signif(l[2,1],digits=2), "mutations/day, P =",
	mean(perm >= obs))) +
guides(fill=guide_legend(override.aes = list(pch = 21)))
dev.off()




