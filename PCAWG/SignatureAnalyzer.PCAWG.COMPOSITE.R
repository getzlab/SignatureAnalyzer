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

#################################################################################################
###### COMPOSITE SIGNATURE EXTRACTIO AND ATTRIBUTION
#################################################################################################

#################################################################################################
###### SIGNATURE EXTRACTION
#################################################################################################
###### SignatureAnalyzer applied two-step signature extraction strategy in 2,780 PCAWG samples. 
###### First, the global signature extraction was performed for the low mutation burden samples (n=2624) without putative POLE and MSI samples, skin tumours, and one TMZ sample. 
###### Second, additional signatures unique to hyper-mutated samples were extracted with maintaining all signatures discovered in the low mutation burden-samples 
###### This approach is expected to minimize a "signature bleeding" or biase of hyper- or ultra-mutated samples on the signature extraction process and 
###### also enalbes key information on the signature availability i.e., the signatures found only in hyper-mutated samples are not allowed in low mutation burden-samples, 
###### while all signatures are available in hyper-mutated samples with no constraints.
#################################################################################################

CURRENT <- paste(getwd(),"/",sep="")
INPUT <- paste(CURRENT,"INPUT_SignatureAnalyzer/",sep="") ### Directory for INPUT data
OUTPUT <- paste(CURRENT,"OUTPUT_SignatureAnalyzer/",sep="") #### Directory for OUTPUT data
TEMPORARY <- paste(CURRENT,"TEMPORARY_SignatureAnalzyer/",sep="")
system(paste("mkdir",OUTPUT,sep=" "))
system(paste("mkdir",TEMPORARY,sep=" "))

source("SignatureAnalyzer.PCAWG.function.R")

###### loading input mutation counts matrix along 1536 SNV + 78 DNP + 83 INDEL features across 2780 PCAWG samples 
###### SNV features in 1536 contexts were defined as six-letters such as "CTAATT", in which first two letters refers to the reference and altered bases, and
###### next four letters represent a base at -2, -1, +1, and +2 position at mutated pyrimidines (C or T), e.g, C>T at AACTT penta-nucleotide sequence motif.
###### SNV features in 96 contexts were defined as four-letters such as "CTAT", in which first two letters refers to the reference and altered bases, and
###### next two letters represent a base at -1 and +1 position at mutated pyrimidines (C or T), e.g, C>T at ACT tri-nucleotide sequence motif.

load(file=paste(INPUT,"lego1536.PAN.SNV.091217.RData",sep="")) ### SNV lego-matrix in 1536 contexts
load(file=paste(INPUT,"lego96.PAN.SNV.091217.RData",sep="")) ### SNV lego-matrix in 96 contexts
load(file=paste(INPUT,"lego1536.DNP.091217.RData",sep="")) ### DNP lego-matrix
load(file=paste(INPUT,"lego1536.INDEL.091217.RData",sep="")) ### INDEL lego-matrix

lego1536.PAN <- rbind(lego1536.SNV,lego1536.DNP,lego1536.INDEL) ### COMPOSITE lego-matrix in 1697 features (1536 SNV + 78 DNP + 83 INDEL)

###### loading information on hyper- or ultra-mutated samples; POLE (n=9), MSI (n=39), all SKIN (n=107), and single TMZ (CNS_GBM__SP25494)
load(file=paste(INPUT,"sample.POLE.RData",sep=""))
load(file=paste(INPUT,"sample.MSI1.RData",sep=""))
load(file=paste(INPUT,"sample.remove.RData",sep=""))
sample.TMZ <- c("CNS_GBM__SP25494")

###### COMPOSITE lego matrix for PRIMARY and SECONDARY datasets
lego.PRIMARY <- lego1536.PAN[,!colnames(lego1536.PAN)%in%sample.remove]
lego.SECONDARY <- lego1536.PAN[,colnames(lego1536.PAN)%in%sample.remove]

######################
####### First step: Signature Extraction for the PRIMARY 2624 sample set 
####### Execute several Bayesian NMF runs and select the solution with the lowest value of res[[4]] (-log(posterior)) at the lowest number of signatures.
####### Bayesian NMF begins with the maximum number of signatures (Kcol) and irrelevant columns and rows in both W and H, respectively,  are iteratively removed through 
####### the automatic relevance determination techinique and the final number of signatures is determined 
####### by counting non-zero columns in W (res[[1]]). 
######################
	if (FALSE) {
		Kcol <- 50
		method <- "L1W.L2H.LEGO1536.COMPOSITE_FIRST"
		for (i in 1:10) {
	        	res <- BayesNMF.L1W.L2H(as.matrix(lego.PRIMARY),200000,10,5,1.e-05,Kcol,Kcol,1)
		        save(res,file=paste(OUTPUT,paste(method,i,"RData",sep="."),sep=""))
		}
	}

######################
####### Second step: Signature Extraction for the SECONDARY 156 sample set
####### Execute several Bayesian NMF runs and select the solution with the lowest value of res[[4]] (-log(posterior)) at the lowest number of signatures
######################
	load(file=paste(OUTPUT,"L1W.L2H.LEGO1536.COMPOSITE_FIRST.RData",sep="")) ### loading BayesNMF solution for the primary data set (example RData file).
	W <- res[[1]]
        H <- res[[2]]
        index <- colSums(W) > 1 ### sum(index) is the number of signatures in the primary data set

	###### All mutation burdens of the primary dataset are now included in the W matrix (signature-loading) and the H matrix (activity-loading).
        W <- W[,index]
        H <- H[index,]
	W1536.PRIMARY <- W
	H1536.PRIMARY <- H
        for (i in 1:ncol(W)) {
                H1536.PRIMARY[i,] <- H1536.PRIMARY[i,]*colSums(W)[i]
                W1536.PRIMARY[,i] <- W1536.PRIMARY[,i]*rowSums(H)[i]
        }
	colnames(W1536.PRIMARY) <- paste("Primary.W",seq(1:ncol(W1536.PRIMARY)),sep="")
	rownames(H1536.PRIMARY) <- colnames(W1536.PRIMARY)
	W1536.PRIMARY.norm <- apply(W1536.PRIMARY,2,function(x) x/sum(x)) ## Normalized primary signatures

	###### PRIMARY signatures (W1536.PRIMARY) with all mutation burdens were added to the muation count matrix of the secondary data set.
	###### Thus all PRIMARY signatures are now treated as extra fake samples with keeping mutation burdens in the PRIMARY 2624 samples.
	###### This process enforces the primary signatures to be reatined in the secondary signature extraction, 
	###### while enablling a discovery of novel signatures unique to the second primary data set.
        lego.PRIMARY.SECONDARY <- cbind(W1536.PRIMARY,lego.SECONDARY)
        method <- "L1W.L2H.LEGO1536.COMPOSITE_SECOND"
	if (FALSE) {
        Kcol <- 96
        for (i in 1:10) {
                res <- BayesNMF.L1W.L2H(as.matrix(lego.PRIMARY.SECONDARY),200000,10,5,5.e-06,Kcol,Kcol,1)
                save(res,file=paste(OUTPUT,paste(method,i,"RData",sep="."),sep=""))
	}
	}

	load(file=paste(OUTPUT,"L1W.L2H.LEGO1536.COMPOSITE_SECOND.RData",sep="")) ### loading BayesNMF solution for the primary data set (example RData file)
	W <- res[[1]]
        H <- res[[2]]
        index <- colSums(W) > 1 ### sum(index) is the number of extracted signatures in the primary data set
        W <- W[,index]
        H <- H[index,]
	W1536.SECONDARY <- W
	H1536.SECONDARY <- H
        for (i in 1:ncol(W)) {
                H1536.SECONDARY[i,] <- H1536.SECONDARY[i,]*colSums(W)[i]
                W1536.SECONDARY[,i] <- W1536.SECONDARY[,i]*rowSums(H)[i]
        }
	colnames(W1536.SECONDARY) <- paste("Secondary.W",seq(1:ncol(W1536.SECONDARY)),sep="")
	rownames(H1536.SECONDARY) <- colnames(W1536.SECONDARY)
	W1536.SECONDARY.norm <- apply(W1536.SECONDARY,2,function(x) x/sum(x))
	
	###### cosine similairy between W1536.SECONDARY (signatures in the second step)  and W1536.PRIMARY (signatures in the first step)
	corr <- plot.W.correlation(W1536.SECONDARY.norm[1:1536,],W1536.PRIMARY.norm[1:1536,])
	corr.max <- apply(corr,1,function(x) max(x)) ### maximum cosine similarity of secondary signatures to primary signatures
	corr.id <- apply(corr,1,function(x) which.max(x))  ### primary signature ID with corr.max
	
	###### summary data-frame mapping the signatures in the second extration step to the signatures in the frist step 
	df.summary <- data.frame(colnames(W1536.SECONDARY),colnames(W1536.PRIMARY)[corr.id],corr.id,corr.max) ### summary for the comparison of secondary signatures to primary ones
	colnames(df.summary) <- c("secondary","primary","primary.id","CS.max") 
	df.summary[,"CS_0.99"] <- df.summary$CS.max > 0.99 

	###### identify primary signatures retained in the secondary signature extraction (W1536.1) and those attributions in 2624 samples (H1536.1)
	W1536.1 <- W1536.PRIMARY[,match(df.summary$primary[df.summary$CS_0.99],colnames(W1536.PRIMARY),nomatch=0)]
	H1536.1 <- H1536.PRIMARY[match(df.summary$primary[df.summary$CS_0.99],rownames(H1536.PRIMARY),nomatch=0),]
	colnames(W1536.1) <- df.summary$secondary[df.summary$CS_0.99]
	rownames(H1536.1) <- colnames(W1536.1)

	###### identify signatures unique to the hyper-mutated samples (W1536.2) in the secondary signature extraction and those attributions in 156 samples (H1536.2.2)
	###### H1536.2.1 is a signature attribution of the primary signatures in hyper-mutated samples.
	W1536.2 <- W1536.SECONDARY[,!df.summary$CS_0.99]
	H1536.2.1 <- H1536.SECONDARY[df.summary$CS_0.99,(ncol(W1536.PRIMARY)+1):ncol(H1536.SECONDARY)]
	H1536.2.2 <- H1536.SECONDARY[!df.summary$CS_0.99,(ncol(W1536.PRIMARY)+1):ncol(H1536.SECONDARY)]
	n.hyper <- ncol(W1536.2) ### # of signatures unique to hyper-mutated samples

	W.tmp1 <- W1536.1
	lego.tmp1 <- get.lego96.from.lego1536(W.tmp1[1:1536,]) ### get lego in 96 contexts from 1536 contexts
	colnames(lego.tmp1) <- gsub("Secondary.","",colnames(lego.tmp1))
	p1 <- plot.signature.SNV(lego.tmp1,"validation")
	pdf(file=paste(OUTPUT,paste("signature.SNV.primary.re_ordered.pdf",sep="."),sep=""),width=(12),height=0.75*ncol(lego.tmp1))
		plot(p1) ### plotting primary signatures
	dev.off()

	W.tmp2 <- cbind(W1536.1,W1536.2)
	lego.tmp2 <- get.lego96.from.lego1536(W.tmp2[1:1536,])
	colnames(lego.tmp2) <- gsub("Secondary.","",colnames(lego.tmp2))
	p2 <- plot.signature.SNV(lego.tmp2,"validation")
	pdf(file=paste(OUTPUT,paste("signature.SNV.secondary.re_ordered.pdf",sep="."),sep=""),width=(12),height=0.75*ncol(lego.tmp2))
		plot(p2) ### plotting all signatures (primary + secondary)
	dev.off()

#################################################################################################
###### SIGNATURE ATTRIBUTION
#################################################################################################
###### SignatureAnalyzer applied a separate process for low-mutation burden and hyper-mutated samples in all COMPOSITE, SBS, DBS, and ID signature attributions. 
###### The signature attribution in low-mutation burden samples was performed in each tumour type, 
###### while the attribution was separately performed in MSI (n=39), POLE (n=9), skin-melanoma cohort (n=107), and a single TMZ sample. 
###### In each cohort, the signature availability (e.g., which signatures are active or not in each cohort) was automatically inferred 
###### through the automatic relevance determination process applied to the activity matrix H only, while fixing the signature W. 
###### To prevent signatures found in hyper-mutated samples in low-mutation burden samples 
###### the H matrix was element-wise-multiplied by a binary, signature indicator matrix Z (signature by samples) in every multiplication update of H.
###### Every elements of Z corresponding to the hyper-mutated signatures and low-mutation burden samples set to zero. 
###### Three additional rules were applied in the SBS signature attribution to enforce biological plausibility and minimize a signature bleeding
###### (i) allow signature SBS4 (smoking signature) only in lung and head and neck cases; 
###### (ii) allow signature SBS11 (TMZ signature) in a single GBM sample;   
###### (iii) allow the 25 signatures found in hyper-mutated samples to be used only in hyper-mutated samples. 
###### On the other hand, no other rules except that signatures found in hyper-mutated samples were not allowed in low-mutation burden-samples
###### were applied in both ID and DBS signature attributions 
#################################################################################################

######################
####### Assigning attributions in the PRIMARY 2624 sample set
######################
	W1536.1.norm <- apply(W1536.1,2,function(x) x/sum(x))
	W1536.2.norm <- apply(W1536.2,2,function(x) x/sum(x))

	W0 <- W1536.1.norm #### fixed PRIMARY W
	H0 <- H1536.1      #### initialization for PRIMARY H
	V0 <- lego1536.PAN[,match(colnames(H0),colnames(lego1536.PAN),nomatch=0)] #### COMPOSITE lego matrix for 2624 samples
	K0 <- ncol(W0)     #### # of signatures

	x1 <- sapply(colnames(V0),function(x) strsplit(x,"__")[[1]][1])
	x2 <- sapply(colnames(V0),function(x) strsplit(x,"__")[[1]][2])
	ttype <- x1
	ttype[ttype%in%c("Breast_AdenoCA","Breast_DCIS","Breast_LobularCA")] <- "Breast"
	ttype[ttype%in%c("Cervix_AdenoCA","Cervix_SCC")] <- "Cervix"
	ttype[ttype%in%c("Myeloid_MDS","Myeloid_MPN")] <- "Myeloid_MDS/MPN"
	ttype[ttype%in%c("Bone_Benign","Bone_Epith")] <- "Bone-Other"
	ttype.unique <- unique(ttype)
	n.ttype.unique <- length(ttype.unique)

	Z0 <- array(1,dim=c(K0,ncol(V0))) ### signature indicator matrix Z (0 = not allowed, 1 = allowed); all elements initialized by one.
	colnames(Z0) <- colnames(H0)
	rownames(Z0) <- rownames(H0)
	Z0[11,!ttype%in%c("Head_SCC","Lung_AdenoCA","Lung_SCC")] <- 0 ### Smoking signature is only allowed in Lung and Head-neck

	a0 <- 10 ### default parameter 
	phi <- 1.5 ### default parameter
	for (i in 1:n.ttype.unique) {
		#### First step: determine a set of optimal signatures best explaining the observed mutations in each cohort.
		#### The automatic relevance determination technique was applied to the H matrix only, while keeping signatures (W0) frozen.
        	cohort <- ttype.unique[i]
	        W1 <- W0
	        V1 <- as.matrix(V0[,ttype==cohort])
	        Z1 <- Z0[,match(colnames(V1),colnames(Z0),nomatch=0)]
	        H1 <- Z1*H0[,match(colnames(V1),colnames(Z0),nomatch=0)]
	        lambda <- rep(1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2,K0)
	        res0 <- BayesNMF.L1.KL.fixed_W.Z(as.matrix(V1),as.matrix(W1),as.matrix(H1),as.matrix(Z1),lambda,2000000,a0,1.e-07,1/phi)
	        H2 <- res0[[2]]
	        colnames(H2) <- colnames(V1)
	        rownames(H2) <- colnames(W0)

		#### Second step: determine a sample-level signature attribution using selected sigantures in the first step.
	        index.H2 <- rowSums(H2)>1 ### identify only active signatures in the cohort
	        Z2 <- Z1
	        Z2[!index.H2,] <- 0 ### only selected signatures in the first step are allowed + the original contraints on the signature availability from Z1.
	        for (j in 1:ncol(H2)) {
	                tmpH <- rep(0,ncol(W0))
	                if (sum(V1[,j])>=5) {
	                        lambda <- 1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2
	                        res <- BayesNMF.L1.KL.fixed_W.Z.sample(as.matrix(V1[,j]),W0,as.matrix(H2[,j]),as.matrix(Z2[,j]),lambda,1000000,a0,1.e-07,1)
	                        tmpH <- res[[2]]
	                }
	                if (j==1) {
	                        H3 <- tmpH
	                } else {
	                        H3 <- cbind(H3,tmpH)
	                        cat(j,'\n')
	                }
	        }
	        colnames(H3) <- colnames(V1)
	        rownames(H3) <- colnames(W0)
	        if (i==1) {
	                H2.all <- H2
	                H3.all <- H3
	        } else {
	                H2.all <- cbind(H2.all,H2)
	        	H3.all <- cbind(H3.all,H3)
	        }
	}
	H.COMPOSITE.nohyper <- H3.all ### attributions of COMPOSITE signatures in the PRIMARY sample set

######################
####### Assigning attributions in the SECONDARY 156 sample set
######################
	W0 <- cbind(W1536.1.norm,W1536.2.norm) ### fixed (PRIMARY + SECONDARY) W
	H0 <- rbind(H1536.2.1,H1536.2.2)       ### initialization for (PRIMARY + SECONDARY) H for 156 samples
	V0 <- lego1536.PAN[,match(colnames(H0),colnames(lego1536.PAN),nomatch=0)] #### COMPOSITE lego matrix for 156 samples
	K0 <- ncol(W0)

	x0 <- colnames(H0)
	ttype <- rep("Skin_Melanoma",length(x0))
	ttype[x0%in%sample.MSI1] <- "MSI"
	ttype[x0%in%sample.POLE] <- "POLE"
	ttype[x0=="CNS_GBM__SP25494"] <- "TMZ"
	ttype.unique <- unique(ttype)
	n.ttype.unique <- length(ttype.unique)

	Z0 <- array(1,dim=c(K0,ncol(V0))) ### signature indicator matrix Z (0 = not allowed, 1 = allowed); all elements initialized by one.
	colnames(Z0) <- colnames(H0)
	rownames(Z0) <- rownames(H0)
	Z0[11,!ttype%in%c("Head_SCC","Lung_AdenoCA","Lung_SCC")] <- 0 ### Smoking signature is only allowed in Lung and Head-neck
	Z0[59,!ttype=="TMZ"] <- 0 ### TMZ signature is only allowed in a single TMZ sample.
	
        for (i in 1:n.ttype.unique) {
                #### First step: determine a set of optimal signatures best explaining the observed mutations in each cohort.
                #### The automatic relevance determination technique was applied to the H matrix only, while keeping signatures (W0) frozen.
                cohort <- ttype.unique[i]
                W1 <- W0
                V1 <- as.matrix(V0[,which(ttype==cohort)])
		if (cohort=="TMZ") {
			colnames(V1) <- sample.TMZ
                	Z1 <- as.matrix(Z0[,match(colnames(V1),colnames(Z0),nomatch=0)])
                	H1 <- Z1*as.matrix(H0[,match(colnames(V1),colnames(Z0),nomatch=0)])
			colnames(Z1) <- sample.TMZ
			colnames(H1) <- sample.TMZ
		}
               	Z1 <- as.matrix(Z0[,match(colnames(V1),colnames(Z0),nomatch=0)])
               	H1 <- Z1*as.matrix(H0[,match(colnames(V1),colnames(Z0),nomatch=0)])
                lambda <- rep(1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2,K0)
                res0 <- BayesNMF.L1.KL.fixed_W.Z(as.matrix(V1),as.matrix(W1),as.matrix(H1),as.matrix(Z1),lambda,2000000,a0,1.e-07,1/phi)
                H2 <- res0[[2]]
                colnames(H2) <- colnames(V1)
                rownames(H2) <- colnames(W0)

                #### Second step: determine a sample-level signature attribution using selected sigantures in the first step.
                index.H2 <- rowSums(H2)>1 ### identify only active signatures in the cohort
                Z2 <- Z1
                Z2[!index.H2,] <- 0 ### only selected signatures in the first step are allowed + the original contraints on the signature availability from Z1.
                for (j in 1:ncol(H2)) {
                        tmpH <- rep(0,ncol(W0))
                        if (sum(V1[,j])>=5) {
                                lambda <- 1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2
                                res <- BayesNMF.L1.KL.fixed_W.Z.sample(as.matrix(V1[,j]),W0,as.matrix(H2[,j]),as.matrix(Z2[,j]),lambda,1000000,a0,1.e-07,1)
                                tmpH <- res[[2]]
                        }
                        if (j==1) {
                                H3 <- tmpH
                        } else {
                                H3 <- cbind(H3,tmpH)
                                cat(j,'\n')
                        }
                }
                colnames(H3) <- colnames(V1)
                rownames(H3) <- colnames(W0)
                if (i==1) {
                        H2.all <- H2
                        H3.all <- H3
                } else {
                        H2.all <- cbind(H2.all,H2)
                        H3.all <- cbind(H3.all,H3)
                }
        }
	H.COMPOSITE.hyper <- H3.all ### attributions of COMPOSITE signatures in the SECONDARY sample set

######################
####### Combining attributions of PIRMARY and SECONARY datasets 
######################
	tmp1 <- array(0,dim=c(n.hyper,ncol(H.COMPOSITE.nohyper)))
	rownames(tmp1) <- colnames(W1536.2)
	colnames(tmp1) <- colnames(H.COMPOSITE.nohyper)
	H.COMPOSITE <- cbind(rbind(H.COMPOSITE.nohyper,tmp1),H.COMPOSITE.hyper) ### attributions of COMPOSITE signatures across 2780 samples

	###### save all results: lego.COMPOSITE - original lego-matrix, W.COMPOSITE - signatures, H.COMPOSITE - attributions
        W.COMPOSITE <- W0
        lego.COMPOSITE <- lego1536.PAN[,match(colnames(H.COMPOSITE),colnames(lego1536.PAN),nomatch=0)]
        summary.COMPOSITE <- list(lego.COMPOSITE,W.COMPOSITE,H.COMPOSITE)
        save(summary.COMPOSITE,file=paste(OUTPUT,"summary.COMPOSITE.RData",sep=""))

	###### calculate a separate attribution of SNVs, INDELs, and DNPs
	x <- get.SNV.INDEL.DNP(apply(W.COMPOSITE,2,function(x) x/sum(x)),H.COMPOSITE)
        H.COMPOSITE.SNV <- x[[1]]
        H.COMPOSITE.INDEL <- x[[2]]
        H.COMPOSITE.DNP <- x[[3]]
