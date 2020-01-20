############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:
####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.
####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.
####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.
#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

library(gridExtra)
library(ggplot2)
library(gplots)
library(reshape2)
library(grid)

#################################################################################################
###### SignatureAnalyzer employs a Bayesian variant of non-negative matrix factorization algorithm, and
###### enables optimal inferences for the number of signatures through the automatic relevance determination technique, and 
###### delivers highly interpretable and sparse representations for both signature profiles and attributions 
###### at a balance between data fitting and model complexity
#################################################################################################

###############################
###### This function is used in both PRIMARY and SECONDARY steps of the signature extraction in 2780 PCAWG samples.
###### Bayesian non-negative matrix factoriztion algorithm with an exponential prior for W and a half-normal prior for H.
###### DEFALUT PARAMETERS: a0=10, b0=5, phi=1, n.iter=2000000, tol=1.e-05
###### INPUT: V (mutation counts matrix), K0 (maximum num of signatures usually 96)
###### OUTPUT: W (signature-loading matrix), H (activity-loading matrix)
###############################

BayesNMF.L1W.L2H <- function(V0,n.iter,a0,b0,tol,K,K0,phi) {
        eps <- 1.e-50
        del <- 1.0
        V <- V0-min(V0)
        Vmax <- max(V)
        N <- dim(V)[1]
        M <- dim(V)[2]
        W <- matrix(runif(N * K)*sqrt(mean(V)),ncol=K)
        H <- matrix(runif(M * K)*sqrt(mean(V)),ncol=M)
        V.ap <- W %*% H + eps
        I_NM <- array(1,dim=c(N,M))
        I_KM <- array(1,dim=c(K,M))
        I_NK <- array(1,dim=c(N,K))

        C <- N+M/2+a0-1
        beta <- C/(colSums(W)+0.5*rowSums(H*H)+b0)
        beta.bound <- C/b0
        beta.cut <- beta.bound/1.25

        n.beta <- list()
        n.beta[[1]] <- beta
        iter <- 2
        count <- 1
        while (del >= tol & iter < n.iter) {
                B <- diag(beta)
                H <- H*(t(W)%*%(V/V.ap)/(t(W)%*%I_NM+B%*%H+eps))
                V.ap <- W%*%H+eps
                W <- W*((V/V.ap)%*%t(H)/(I_NM%*%t(H)+I_NK%*%B+eps))
                beta <- C/(colSums(W)+0.5*rowSums(H*H)+b0)
                V.ap <- W%*%H+eps
                n.beta[[iter]] <- beta
                error <- sum((V-V.ap)^2)
                if (iter %% 100 == 0) {
                        del <- max(abs(beta-n.beta[[iter-1]])/n.beta[[iter-1]])
                        like <- sum(V * log((V+eps)/(V.ap+eps)) + V.ap - V)
                        evid <- like + sum((colSums(W)+0.5*rowSums(H*H)+b0)*beta-C*log(beta))
                        cat(iter,evid,like,error,del,sum(colSums(W)>1.e-05),sum(beta<=beta.cut),'\n')
                        res <- list(W,H,like,evid,beta,error)
                        save(res,file=paste(TEMPORARY,paste("Bayes.L1W.L2H.temp",iter,"RData",sep="."),sep=""))
                }
                iter <- iter+1
        }
        lambda <- 1/beta
        return(list(W,H,like,evid,lambda,error))
}

###############################
###### This function is used for the first step of signature attribution to determine an optimal set of signatures out of input signatures (W), given V.
###### DEFALUT PARAMETERS: a0=10, tol=1.e-07, phi = 1.50 
###### INPUTS: V - mutation counts matrix, W - columnwise-normalized signature matrix, Z - binary matrix to indicate a signature availability in each sample.
###### OUTPUTS: H - activity-loading matrix
###############################
BayesNMF.L1.KL.fixed_W.Z <- function(V,W,H,Z,lambda,n.iter,a0,tol,phi) {
        eps <- 1.e-50
        del <- 1.0
        N <- dim(V)[1]
        M <- dim(V)[2]
        V <- V + eps
        K <- ncol(W)
        K0 <- K
        V.ap <- W %*% H + eps
        I <- array(1,dim=c(N,M))
        C <- N + M + a0 + 1
        b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
        lambda.bound <- b0/C
        lambda.cut <- 1.1 * lambda.bound

        n.lambda <- list()
        n.lambda[[1]] <- lambda
        iter <- 2
        count <- 1
        while (del >= tol & iter < n.iter) {
                H <- H * (t(W) %*% (V/V.ap))/(matrix(rep(colSums(W)+phi/lambda,M),ncol=M) + eps)
                H <- H * Z
                V.ap <- W %*% H + eps
                lambda <- (colSums(W) + rowSums(H) + b0)/C
                del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
                n.lambda[[iter]] <- lambda
                if (iter %% 100 == 0) {
                        error <- sum((V-V.ap)^2)
                        like <- sum(V*log((V+1)/(V.ap+1))+V.ap-V)
                        evid <- like/phi+sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
                        cat(iter,evid,like,error,del,sum(rowSums(H)!=0),sum(lambda>=lambda.cut),'\n')
                }
                iter <- iter+1
        }
        error <- sum((V-V.ap)^2)
        like <- sum(V*log((V+1)/(V.ap+1))+V.ap-V)
        evid <- like/phi+sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
        cat("***********************************************************",'\n')
        return(list(W,H,like,evid,lambda,error))
}

###############################
###### This function is used for the second step of signature attribution to determine a sample-level attribution using selected signatures in the first step.
###### DEFALUT PARAMETERS: a0=10, tol=1.e-07, phi=1.0
###### INPUTS: V - mutation counts matrix, W - columnwise-normalized signature matrix, Z - binary matrix to constraint a signature availability.
###### OUTPUTS: H - activity-loading matrix
###############################
BayesNMF.L1.KL.fixed_W.Z.sample <- function(V,W,H,Z,lambda,n.iter,a0,tol,phi) {
        eps <- 1.e-50
        del <- 1.0
        N <- dim(V)[1]
        M <- dim(V)[2]
        V <- V + eps
        K <- ncol(W)
        K0 <- K
        V.ap <- W %*% H + eps
        I <- array(1,dim=c(N,M))
        C <- N + M + a0 + 1
        b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
        lambda.bound <- b0/C
        lambda.cut <- 1.1 * lambda.bound

        n.lambda <- list()
        n.lambda[[1]] <- lambda
        iter <- 2
        count <- 1
        Z0 <- Z
        while (del >= tol & iter < n.iter) {
                H <- H * (t(W) %*% (V/V.ap))/(matrix(rep(colSums(W)+0*phi/lambda,M),ncol=M) + eps)
                H <- H * (Z0)
                V.ap <- W %*% H + eps
                lambda <- (colSums(W) + rowSums(H) + b0)/C
                del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
                n.lambda[[iter]] <- lambda
                if (iter %% 100 == 0) {
                        error <- sum((V-V.ap)^2)
                        like <- sum(V*log((V+1)/(V.ap+1))+V.ap-V)
                        evid <- like/phi+sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
                        cat(iter,evid,like,error,del,sum(rowSums(H)!=0),sum(lambda>=lambda.cut),'\n')
                        res <- list(W,H,like,evid,lambda,error)
                }
                iter <- iter+1
        }
        error <- sum((V-V.ap)^2)
        like <- sum(V*log((V+1)/(V.ap+1))+V.ap-V)
        evid <- like/phi+sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
        cat("***********************************************************",'\n')
        return(list(W,H,like,evid,lambda,error))
}

######################################
######################################
scale <- 0.8
.theme_ss <- theme_bw(base_size=14) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12*scale, family="mono"),
        axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
        axis.title.x = element_text(face="bold",colour="black",size=14*scale),
        axis.title.y = element_text(face="bold",colour="black",size=14*scale),
        axis.text = element_text(size = 16*scale, family = "mono"),
        strip.text = element_text(lineheight=0.5),
        strip.text.x = element_text(size=10*scale,face='bold',angle=00),
        strip.text.y = element_text(size=10*scale,face="bold"),
        strip.background = element_rect(colour="black",fill="gray85"),
        panel.margin = unit(0.20,"lines"),
        plot.title=element_text(lineheight=1.0,face="bold",size=12*scale))
lego.colors <- c("cyan","red","yellow","purple","green","blue")

context96 <- read.delim(paste(INPUT,"context96.txt",sep=""),header=T,sep='\t',as.is=T,skip=0)
context96.label <- context96[,2]

###### This function is to get a lego96-representation from a lego1536-representation.
get.lego96.from.lego1536 <- function(W0) {
        W <- W0[1:1536,]
        index1536 <- rownames(W)
        contig.1536 <- paste(substring(index1536,1,2),substring(index1536,4,4),substring(index1536,5,5),sep="")
        contig.96 <- contig.1536
        for (i in 1:96) {
                contig.96[contig.1536%in%context96.label[i]] <- i
        }
        mapping <- data.frame(index1536,contig.1536,contig.96)
        W96 <- array(0,dim=c(96,ncol(W)))
        for (i in 1:96) {
                W96[i,] <- colSums(W[mapping$contig.96==i,])
        }
        rownames(W96) <- context96.label
        colnames(W96) <- colnames(W0)
        if (nrow(W0) > 1536) {
                W96 <- rbind(W96,W0[1537:nrow(W0),])
        }
        return(W96)
}

###### This function is to compute a cosine similarity between two signatures sets W1 and W2
plot.W.correlation <- function(W1,W2) {
        K1 <- ncol(W1)
        K2 <- ncol(W2)
        x <- array(0,dim=c(K1,K2))
        for (i in 1:K1) {
        for (j in 1:K2) {
                if (sum(W1[,i])!=0 & sum(W2[,j])!=0) {
                        x[i,j] <- W1[,i]%*%W2[,j]/sqrt(sum(W1[,i]^2))/sqrt(sum(W2[,j]^2))
                }
        }
        }
        rownames(x) <- colnames(W1)
        colnames(x) <- colnames(W2)
        return(x)
}

###### This function is for plotting ID (INDELs) signatures.
plot.signature.indel <- function(W,title) {
        df1 <- data.frame(W)
        df1[df1 < 1.e-10] <- 0
        df1[,"feature"] <- rownames(W)
        df1 <- melt(df1,id.var="feature")
        colnames(df1) <- c("feature","signature","activity")
        x1 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][1])
        x2 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][2])
        x3 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][3])
        x4 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][4])
        x5 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][5])
        df1[,"type1"] <- paste(x1,x2,sep="_")
        df1[,"type2"] <- paste(x4,x5,sep="_")
        df1$feature <- factor(df1$feature,levels=rownames(W))
        K <- ncol(W)
        p = ggplot(df1)
        p = p + geom_bar(aes_string(x="feature",y="activity",fill="type1"),stat="identity",position="identity")
        p = p + facet_grid(signature ~ ., scale = "free_y")
        p = p + .theme_ss
        p = p + scale_fill_manual(values=c("cyan","red","orange","purple","green","blue","black")) #p = p + scale_fill_brewer(palette = "Set1")
        p = p + guides(fill=FALSE) #p = p + theme(legend.position = "none")
        p = p + ggtitle(title)
        p = p + xlab("Features") + ylab("Contributions")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=10*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=10*scale))
        p = p + theme(strip.text.x = element_text(size = 10*scale, colour = "black", angle = 0))
        p = p + theme(strip.text.y = element_text(size = 10*scale, colour = "black", angle = 270))
        return(p)
}

###### This function is for plotting DBS (DNPs) signatures.
plot.signature.DNP <- function(W,title) {
        df1 <- data.frame(W)
        df1[df1 < 1.e-10] <- 0
        df1[,"feature"] <- rownames(W)
        df1 <- melt(df1,id.var="feature")
        colnames(df1) <- c("feature","signature","activity")
        x1 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][1])
        x2 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][2])
        x3 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][3])
        x4 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][4])
        x5 <- sapply(df1$feature,function(x) strsplit(x,"_")[[1]][5])
        df1[,"type1"] <- paste(x1,x2,sep="_")
        df1[,"type2"] <- paste(x4,x5,sep="_")
        #df1$feature <- factor(df1$feature,levels=rownames(W))
        K <- ncol(W)
        p = ggplot(df1)
        p = p + geom_bar(aes_string(x="feature",y="activity",fill="type1"),stat="identity",position="identity")
        p = p + facet_grid(signature ~ ., scale = "free_y")
        p = p + .theme_ss
        p = p + guides(fill=FALSE) #p = p + theme(legend.position = "none")
        p = p + ggtitle(title)
        p = p + xlab("Features") + ylab("Contributions")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=10*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=10*scale))
        p = p + theme(strip.text.x = element_text(size = 10*scale, colour = "black", angle = 0))
        p = p + theme(strip.text.y = element_text(size = 10*scale, colour = "black", angle = 270))
        return(p)
}

###### This function is for plotting SBS (SNVs) signatures.
plot.signature.SNV <- function(lego,title) {
        n.sample <- ncol(lego)
        sample.order <- colnames(lego)
        W.SNP <- as.vector(unlist(lego))
        context4 <- rep(rownames(lego),n.sample)
        context3 <- paste(substring(context4,3,3),"-",substring(context4,4,4),sep="")
        base.sub <- paste(substring(context4,1,1),"->",substring(context4,2,2),sep="")
        signature.class <- as.vector(t(matrix(t(rep(colnames(lego),nrow(lego))),nrow=n.sample)))
        df.SNP <- data.frame(W.SNP,context4,context3,base.sub,signature.class)
        colnames(df.SNP) <- c("signature","context4","context3","base.sub","class")
        #sample.order <- colnames(lego)[order(colSums(lego),decreasing=T)]
        df.SNP$class <- factor(df.SNP$class,sample.order)
        df.SNP$context3 <- paste("              ",df.SNP$context3,sep="")
        p = ggplot(df.SNP)
        p = p + geom_bar(aes_string(x="context3",y="signature",fill="base.sub"),stat="identity",position="identity",colour="gray50")
        p = p + facet_grid(class ~ base.sub, scale = "free_y")
        p = p + .theme_ss
        p = p + scale_fill_manual(values=lego.colors) #p = p + scale_fill_brewer(palette = "Set1")
        p = p + guides(fill=FALSE) #p = p + theme(legend.position = "none")
        p = p + ggtitle(title)
        p = p + xlab("Motifs") + ylab("Contributions")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=10*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=10*scale))
        return(p)
}

###### This function is to compute separate attributions of SNVs, DNPs, and INDELs.
get.SNV.INDEL.DNP <- function(tmpW,tmpH) {
        SNV <- array(0,dim=c(nrow(tmpH),ncol(tmpH)))
        INDEL <- array(0,dim=c(nrow(tmpH),ncol(tmpH)))
        DNP <- array(0,dim=c(nrow(tmpH),ncol(tmpH)))
        for (i in 1:nrow(tmpH)) {
                x <- as.matrix(tmpW[,i])%*%t(as.matrix(tmpH[i,]))
                SNV[i,] <- SNV[i,] + colSums(x[1:1536,])
                DNP[i,] <- DNP[i,] + colSums(x[1537:1614,])
                INDEL[i,] <- INDEL[i,] + colSums(x[1615:nrow(x),])
        }
        rownames(SNV) <- rownames(tmpH)
        rownames(DNP) <- rownames(tmpH)
        rownames(INDEL) <- rownames(tmpH)
        colnames(SNV) <- colnames(tmpH)
        colnames(DNP) <- colnames(tmpH)
        colnames(INDEL) <- colnames(tmpH)
        return(list(SNV,INDEL,DNP))
}
