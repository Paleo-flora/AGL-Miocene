
# COEXISTENCE FUNCTIONS

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
# LOAD PACKAGES
library(BIEN)
library(raster)
library(KernSmooth)
library(readxl)
library(CommEcol)
library(cluster)
library(rioja)
library(geodata)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
# just some colors for plots 
jet.colors <- colorRampPalette(c("#0C0C38", "#0000FF", "#00CED1", "#CAFF70", "#FFD700", "#FFA500", "#FF0000", "#CD2626"))
clrs <- colorRampPalette(c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","greenyellow"),alpha=T)
# transparent colors
t_col <- function(color, percent = 50, name = NULL) {
	rgb.val <- col2rgb(color)
	t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255, alpha = (100 - percent) * 255 / 100, names = name)
	invisible(t.col)
}   

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
# Coexistencce analysis
      
# Fossil = a matrix of fossil data [samples, species]
# x = two columns affinities: fosil_taxa, modern affinity 
# quant = quantile of the surface to be taken into account; 
# ci = confidence interval for the estimates

coex <- function(
  Fossil, 
  x, 
  quant = 0.95,
  ci = c(0.025,0.975),
  col = colorRampPalette(c("#0C0C38", "#0000FF", "#00CED1", "#CAFF70", "#FFD700", "#FFA500", "#FF0000", "#CD2626"))(20),

# ... BIEN Arguments
  new.world = TRUE,
  cultivated = FALSE,
  all.taxonomy = FALSE,
  native.status = FALSE,
  natives.only = TRUE,
  observation.type = FALSE,
  political.boundaries = FALSE,
  collection.info = F,

# ... worldclim Arguments
  var = "bio", 	# c( tmin', 'tmax', 'prec' and 'bio')
  res = 10,		# c(0.5, 2.5, 5, 10) minutes of a degree
  subset_wclim = c(1,12),
  wclim_names = c("MAT","MAP"),
  range_x = c(-10,30),
  range_y = c(0,3000),
# 
  plot_fossil_tx = TRUE,
  fosil_tx_file = "PDF_fossil_tx.pdf",
  plot_sample2D = TRUE,
  sample2D_file = "PDF_sample_2D.pdf"
  )
 {
# Check data
   if(any(is.na(match(colnames(Fossil),x[ ,1])))){
	stop("Taxa in the affinities do not match those in the fossil database")
   }

# modern occurrences
   mtx <- as.character(unique(x[ ,2]))	# Modern taxa affinities names

# check modern library
   if(is.na(match("BIEN_occurrences.RData",dir()))== FALSE){
	load("BIEN_occurrences.RData")
	newtx <- setdiff(mtx,names(BIEN_occs))
  	  bien_ls <- vector("list", length(newtx))
  	  names(bien_ls) <- newtx
   } else {
	BIEN_occs <- list()
  	bien_ls <- vector("list", length(mtx))
  	names(bien_ls) <- mtx
	newtx <- mtx	
   }

   if(length(bien_ls) != 0){

# by species
     for (i in which(lapply(bien_ls,length)== 0)){
	a <- BIEN_occurrence_species(newtx[i],
  		cultivated = cultivated ,
  		new.world = new.world,
  		all.taxonomy = all.taxonomy,
  		native.status = native.status,
  		natives.only = natives.only,
  		observation.type = observation.type,
  		political.boundaries = political.boundaries ,
  		collection.info = collection.info
		)
	 if(nrow(a) > 0 ){	
		bien_ls[[i]] <- a[,2:3]
	 } 
     }#i
   }  

   if (any(lapply(bien_ls,length)== 0)){
# by by genus
    for (j in which(lapply(bien_ls,length)== 0)){

	a <- BIEN_occurrence_genus(newtx[j],
  		cultivated = cultivated ,
  		new.world = new.world,
  		all.taxonomy = all.taxonomy,
  		native.status = native.status,
  		natives.only = natives.only,
  		observation.type = observation.type,
  		political.boundaries = political.boundaries ,
  		collection.info = collection.info
		)
	if(nrow(a) > 0 ){	
		bien_ls[[j]] <- a[,3:4]
	} 
    } #j
    save(bien_ls, file = "BIEN_NEW_occurrences.RData")
   }

  if (any(lapply(bien_ls,length)== 0)){
# By family
    for (h in which(lapply(bien_ls,length)== 0)){
	a <- BIEN_occurrence_family(newtx[h],
  		cultivated = cultivated ,
  		new.world = new.world,
  		all.taxonomy = all.taxonomy,
  		native.status = native.status,
  		natives.only = natives.only,
  		observation.type = observation.type,
  		political.boundaries = political.boundaries ,
  		collection.info = collection.info
		)
	if(nrow(a) > 0 ){	
		bien_ls[[h]] <- a[,3:4]
	} 
    }#h
    save(bien_ls, file="BIEN_NEW_occurrences.RData")
  }
# 
# Check taxa without enough records

  if (any(lapply(bien_ls,length) < 3)){
	warning("Some taxa has been excluded because they cant be found on the BIEN database or have less than 3 occurences")
	print(newtx[which(lapply(bien_ls,length)< 3)])	
  }

# Eliminate empty taxa
  bien_ls[which(unlist(lapply(bien_ls,is.null)))] <- NULL

# Keep only taxa with 3 o more records
  bien_ls <- bien_ls[lapply(bien_ls,nrow) >= 3]

# save occurences
  BIEN_occs <- c(BIEN_occs,bien_ls)
  save(BIEN_occs, file="BIEN_occurrences.RData")

# fossil affinities

  aff <- match(x[,2],names(BIEN_occs))
  fos_affs <- aggregate(aff ~ x[,1], FUN = c)
  foss <- vector("list", nrow(fos_affs))
  names(foss) <- fos_affs[,1]
  for(f in 1:length(foss)) {
	occs <- BIEN_occs[unlist(fos_affs[f,2])]
	foss[[f]] <- na.omit(Reduce(rbind,occs))
  }

# Extract Environmental variables

  wc1 <- worldclim_global(var = "bio", res = res, path = getwd())
  wc <- wc1[[subset_wclim]]
  names(wc) <- wclim_names
  k2d <- vector("list", length(foss))
  names(k2d)<- names(foss)

  for(fi in 1:length(foss)){
	wcs <- extract(wc,foss[[fi]][,2:1])
	wcs <- na.omit(wcs)
	bw1 <- dpik(wcs[,2],scalest="stdev")
	bw2 <- dpik(wcs[,3],scalest="stdev")
	k2d[[fi]] <- bkde2D(wcs[ ,c(2,3)],bandwidth=c(bw1,bw2),gridsize = c(101L, 101L),
		range.x = list(range_x, range_y))
	names(k2d[[fi]]) <-  c(wclim_names,"Prob")
  }
  if(plot_fossil_tx == TRUE){
	pdf(fosil_tx_file)
	for(ki in 1:length(k2d)){
	image(k2d[[ki]][[1]],k2d[[ki]][[2]],k2d[[ki]][[3]], main= names(k2d)[ki],
	xlab = wclim_names[1], ylab = wclim_names[2], col=jet.colors(50))
	box()
	}
	dev.off()
  }

# Fossil estimation per sample

  fossil.pdf <- vector("list", nrow(Fossil))
  names(fossil.pdf) <- rownames(Fossil)
  for(fsi in 1:nrow(Fossil)){
    mixture <- colnames(Fossil)[which(Fossil[fsi,] != 0)]
    k2d_mix <- lapply(k2d[mixture],function(x) x[[3]])
    fossil.pdf[[fsi]] <- list(k2d[[1]][[1]],k2d[[1]][[2]],Reduce('+',k2d_mix))
    names(fossil.pdf[[fsi]])<- names(k2d[[1]])	
  }
  if(plot_sample2D==TRUE){
    pdf(sample2D_file)
    for(i in 1:length(fossil.pdf)){
      image(fossil.pdf[[i]][[1]],fossil.pdf[[i]][[2]],
            fossil.pdf[[i]][[3]],col = col, main= names(fossil.pdf)[i],
	      xlab=names(fossil.pdf[[i]])[1],ylab=names(fossil.pdf[[i]])[2])
    }
    dev.off()
  }

  estimates <-  vector("list", length(fossil.pdf))
        for(i in 1:length(estimates)){
                estimates[[i]] <- which(fossil.pdf[[i]][[3]] > quantile(
                        fossil.pdf[[i]][[3]],quant),arr.ind=TRUE)
                estimates[[i]][,1] <- fossil.pdf[[i]][[1]][estimates[[i]][,1]]
                estimates[[i]][,2] <- fossil.pdf[[i]][[2]][estimates[[i]][,2]]
                }
        results <- matrix(nrow = nrow(Fossil),ncol = 10)
        colnames(results) <- c("Mean var1","Median var1","SD var1","LL var1","UL var1",
                "Mean var2","Median var2","SD var2","LL var2","UL var2")
        rownames(results) <- rownames(Fossil)
        
        results[,1] <- unlist(lapply(estimates,function(x) mean(x[,1])))
        results[,2] <- unlist(lapply(estimates,function(x) quantile(x[,1],0.5)))
        results[,3] <- unlist(lapply(estimates,function(x) sd(x[,1]))) 
 	  results[,4] <- unlist(lapply(estimates,function(x) quantile(x[,1],ci[1])))
        results[,5] <- unlist(lapply(estimates,function(x) quantile(x[,1],ci[2])))
        results[,6] <- unlist(lapply(estimates,function(x) mean(x[,2])))
        results[,7] <- unlist(lapply(estimates,function(x) quantile(x[,2],0.5)))
        results[,8] <- unlist(lapply(estimates,function(x) sd(x[,2]))) 
        results[,9] <- unlist(lapply(estimates,function(x) quantile(x[,2],ci[1])))
        results[,10] <- unlist(lapply(estimates,function(x) quantile(x[,2],ci[2])))
 
  output <- list(results,fossil.pdf)
  names(output) <- c("Stats","SamplePDFs")
  return(output)
  }

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   

group <- function(x, group){
 if(nrow(x)!=length(group)){
	error("nrow(x) does not match length(group)")} else {
 g <- unique(group)
 gx <- matrix(ncol=ncol(x),nrow=length(g))
 colnames(gx) <- colnames(x)
 rownames(gx) <- g
 for(i in g){ gx[i,] <- colSums(x[group == i,])}
 return(gx)
}}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   

estimate2Dseq <- function(x, lower_tail = 0.75, col= jet.colors(20)){

 pdfs <- x$SamplePDFs 
 Estimations_X <- matrix(nrow=length(pdfs),ncol=512)
 Estimations_Y <- matrix(nrow=length(pdfs),ncol=512)

 for(i in 1:length(pdfs)){	  

        b <- diff(pdfs[[i]][[1]])[1] * diff(pdfs[[i]][[2]])[1]
        z.pdf1 <- pdfs[[i]][[3]]/sum(pdfs[[i]][[3]]*b)
        z.pdf <- list(pdfs[[i]][[1]],pdfs[[i]][[2]],z.pdf1)
        probs <- cbind(as.vector(z.pdf1),as.vector(z.pdf1)*b)
        probs <- probs[order(probs[,1]),]
        probs[,2] <-cumsum(probs[,2])
        colnames(probs) <- c("density","probability")
        densities <- probs[which(probs[,2] > lower_tail),1]
        Estimation <- data.frame(density = densities, x = NA, y = NA)
        for(ei in 1:nrow(Estimation)){
                Estimation[ei,2:3] <- which(z.pdf1 == Estimation[ei,1],arr.ind = T)
        }
        Estimation$x <- pdfs[[i]][[1]][Estimation[,2]]
        Estimation$y <- pdfs[[i]][[2]][Estimation[,3]]

	  Estimations_X[i,] <- (density(Estimation[,2],bw=sd(Estimation[,2])/2,from=pdfs[[i]][[1]][1],to=pdfs[[i]][[1]][101]))$y
	  Estimations_Y[i,] <- (density(Estimation[,3],bw=sd(Estimation[,3])/2,from=pdfs[[i]][[2]][1],to=pdfs[[i]][[2]][101]))$y

	  Estimations_X[i,] <- Estimations_X[i,]/max(Estimations_X[i,])
	  Estimations_Y[i,] <- Estimations_Y[i,]/max(Estimations_Y[i,])

	}

  par(mfrow=c(1,2))
  par(mar=c(4,3,2,0))
	image(t(Estimations_X),col= col, axes=F)
	axis(side=1,at = pretty(pdfs[[i]][[1]])/max(pretty(pdfs[[i]][[1]])), labels=pretty(pdfs[[i]][[1]]))
	axis(side = 2, at= 1:length(pdfs)/length(pdfs),labels= names(pdfs),cex.axis= 0.5, las= 2)
	mtext(names(pdfs[[1]][1]), side = 1, line = 2.5)
  par(mar=c(4,0.5,2,2.5))
	image(t(Estimations_Y),col= col, axes=F)
	axis(side=1,at = pretty(pdfs[[i]][[2]])/max(pretty(pdfs[[1]][[2]])), labels=pretty(pdfs[[1]][[2]]))
	mtext(names(pdfs[[1]][2]), side = 1, line = 2.5)

  }
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   
#funtion "estimate1" extracts specific apparts of the bidimensional pdf. 
# Arguments: 
# 'x' and 'y' are the enviroenmntal axes; 
# 'z.pdf' matrix of densities; 
# quant, lower tail of probability to be ignored

estimate1 <- function(x,y,z.pdf,quant){
        b<-(x[2]-x[1])*(y[2]-y[1])
        z.pdf1 <- z.pdf/sum(z.pdf*b)
        z.pdf<-list(x,y,z.pdf1)
        probs <- cbind(as.vector(z.pdf1),as.vector(z.pdf1)*b)
        probs <- probs[order(probs[,1]),]
        probs[,2] <-cumsum(probs[,2])
        colnames(probs) <- c("density","probability")
        densities <- probs[which(probs[,2]>quant),1]
        estimates <- matrix(ncol=3,nrow=length(densities))
        colnames(estimates) <- c("density","x","y")
        estimates[,1] <- densities
        for(i in 1:nrow(estimates)){
                estimates[i,2:3] <- which(z.pdf1 == estimates[i,1],arr.ind = T)
        }
        estimates[,2] <- x[estimates[,2]]
        estimates[,3] <- x[estimates[,3]]
        return(estimates)
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   


mixture <- function(pdf1,pdf2,lambda1=seq(0,1,0.01),quant=c(0.25,0.75),main,
        ylim){
        results<-matrix(nrow=length(lambda1),ncol=5)
        colnames(results)<-c("proportion","median","mode","ul","ll")
        results[,1]<-lambda1
        for(i in 1:nrow(results)){
                model<-data.frame(cbind(pdf1$x,
                        pdf1$y*results[i,1]+pdf2$y*(1-results[i,1])))
                model$prob<-cumsum((model[2,1]-model[1,1])*model[,2])
                results[i,2]<-approx(model[,3],model[,1],xout=0.5)$y
                results[i,3]<-model[which(model[,2]==max(model[,2])),1]
                results[i,4:5]<-approx(model[,3],model[,1],xout=quant)$y
        }
        dy.dx<-matrix(nrow=(length(lambda1)-1),ncol=2)
        colnames(dy.dx)<-c("x","dy.dx")
        dy.dx[,1]<-(results[-1,1]+results[-nrow(results),1])/2
        dy.dx[,2]<-(results[-1,2]-results[-nrow(results),2])/lambda1[3]
        dy.dx2<-matrix(nrow=(length(lambda1)-2),ncol=2)
        colnames(dy.dx2)<-c("x","dy.dx2")
        dy.dx2[,1]<-(dy.dx[-1,1]+dy.dx[-nrow(dy.dx),1])/2
        dy.dx2[,2]<-(dy.dx[-1,2]-dy.dx[-nrow(dy.dx),2])/lambda1[3]
        plot(results[,c(1,2)],ylim=ylim,t="n",main=main,
                xlab="Lamdba 1",ylab=main)
        for(i in 1:nrow(results)){
                arrows(results[i,1],results[i,4],results[i,1],results[i,5],
                        length=0,col="grey")
        }
        points(results[,c(1,2)],pch=18,cex=0.8)
return(list(results,dy.dx,dy.dx2))
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   


plot_coex_2D <- function(x, lower_tail=0.75, central= mean, col=c("coral", "royalblue"))
{
 pdfs <- x[[2]]
 rx <- list()
 ry <- list()

 Estimations_X <-  matrix(nrow=length(pdfs),ncol=512)
 Estimations_Y <-  matrix(nrow=length(pdfs),ncol=512)

 for(i in 1:length(pdfs)){	  

        b <- diff(pdfs[[i]][[1]])[1] * diff(pdfs[[i]][[2]])[1]
        z.pdf1 <- pdfs[[i]][[3]]/sum(pdfs[[i]][[3]]*b)
        z.pdf <- list(pdfs[[i]][[1]],pdfs[[i]][[2]],z.pdf1)
        probs <- cbind(as.vector(z.pdf1),as.vector(z.pdf1)*b)
        probs <- probs[order(probs[,1]),]
        probs[,2] <-cumsum(probs[,2])
        colnames(probs) <- c("density","probability")
	  
        densities <- probs[which(probs[,2] > lower_tail),1]
	  densities_array <- arrayInd(which(z.pdf1%in% densities, arr.ind =TRUE),.dim=dim(z.pdf1))
	  densities_range <- apply(densities_array,2,range)
	  rx[[i]] <- pdfs[[i]][[1]][densities_range[,1]]
	  ry[[i]] <- pdfs[[i]][[2]][densities_range[,2]]

	  Estimation <- data.frame(density = densities, x = NA, y = NA)
        for(ei in 1:nrow(Estimation)){
                Estimation[ei,2:3] <- which(z.pdf1 == Estimation[ei,1],arr.ind = T)
        }
        Estimation$x <- pdfs[[i]][[1]][Estimation[,2]]
        Estimation$y <- pdfs[[i]][[2]][Estimation[,3]]

	  Estimations_X[i,] <- (density(Estimation[,2],bw=sd(Estimation[,2])/2,from=pdfs[[i]][[1]][1],to=pdfs[[i]][[1]][101]))$y
	  Estimations_Y[i,] <- (density(Estimation[,3],bw=sd(Estimation[,3])/2,from=pdfs[[i]][[2]][1],to=pdfs[[i]][[2]][101]))$y
  }
	
# PLOT ==============
 par(mfrow=c(2,2))
 par(mar=c(0, 4, 10, 0))
# pdfs var 1
 plot(0,t="n", xlim= range(pdfs[[1]][[1]]),ylim= c(0,max(Estimations_X)), xaxt="n",
	xlab="",ylab="Density of probability", las=1, cex.axis=0.7)
 abline(v= pretty(pdfs[[1]][[1]]), lty=3, col="grey80")

 for(i in 1:length(pdfs)){	
  xs <- seq(pdfs[[i]][[1]][[1]],pdfs[[i]][[1]][[101]],length.out=512)
  polygon(c(xs[1],xs,xs[512]),c(0,Estimations_X[i,],0), col=t_col(col[i]),border= col[i])
 }
box()
 mtext(names(pdfs[[1]])[1], side = 3, line = 3, cex = 0.8)
 axis(side = 3, cex.axis= 0.7)

 par(mar=c(0, 0, 10, 10))

 plot(0,t="n", axes=F,xlab="",ylab="")
 legend("center", legend=names(pdfs), pch= 15, col=col, bg= "white",bty="n")

 par(mar=c(5, 4, 0, 0))

 plot(0,t="n", xlim= range(pdfs[[1]][[1]]),ylim= range(pdfs[[1]][[2]]),
	xlab= names(pdfs[[1]])[1],ylab= names(pdfs[[1]])[2], las=1, cex.axis=0.7)
 abline(v= pretty(pdfs[[1]][[1]]), lty=3, col="grey80")
 abline(h= pretty(pdfs[[1]][[2]]), lty=3, col="grey80")
 for(i in 1:length(pdfs)){	  

	  rect(rx[[i]][1],ry[[i]][1],rx[[i]][2],ry[[i]][2], col= t_col(col[i]),border=NA)
	  points(x[[1]][i,2],x[[1]][i,6], col= col[i], pch= 3, lwd=3) 
}


 par(mar=c(5, 0, 0, 10))
# pdfs var 2
 plot(0,t="n", ylim= range(pdfs[[1]][[2]]),xlim= c(0,max(Estimations_Y)), yaxt= "n",
	ylab="",xlab="Density of probability", las=1, cex.axis=0.7)
 abline(h= pretty(pdfs[[1]][[2]]), lty=3, col="grey80")

 for(i in 1:length(pdfs)){	
  ys <- seq(pdfs[[i]][[2]][[1]],pdfs[[i]][[2]][[101]],length.out=512)
  polygon(c(0,Estimations_Y[i,],0),c(ys[1],ys,ys[512]), col=t_col(col[i]),border= col[i])
 }
box()
  mtext(names(pdfs[[1]])[2], side = 4, line = 3, cex = 0.8)
  axis(side = 4, cex.axis= 0.7, las=1)

}


#===================================================================================
  
