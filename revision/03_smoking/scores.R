library(data.table)
root <- Sys.getenv("OSCC_DATA",Sys.getenv("OSCC_DATA"))
src <- c("Integration_Manuscript/Figure5_SharedNetworkArchitecture/Figure5_MasterData.tsv",
         "Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv")
paths <- c(file.path(root,src),"../01_baseline/results/sample_manifest.tsv")
m <- fread(paths[1],na.strings="NA")
a <- fread(paths[2],na.strings="NA")
b <- fread(paths[3],na.strings=c("NA",""))[sample_type=="Primary Tumor"]
stopifnot(!anyDuplicated(m$PatientID),!anyDuplicated(a$patient_id),!anyDuplicated(b$patient_barcode))
x <- merge(m,a,by.x="PatientID",by.y="patient_id")
x <- merge(x,b[,.(patient_barcode,primary_site,remove_base_and_mixed,p16,ISH)],by.x="PatientID",by.y="patient_barcode")
x[smoking_bin.x=="",smoking_bin.x:=NA_character_]
stopifnot(nrow(x)==341L, identical(x$smoking_bin.x,x$smoking_bin.y))
x[,smoke:=factor(smoking_bin.x,levels=c("Never","Ever"))]
x[,site:=factor(primary_site)]
x[,sex:=factor(Sex)]
x[,stage:=factor(Stage)]
x[,alcohol:=factor(alcohol_bin)]
x[,hpv_status:=factor(fifelse(p16=="Positive"|ISH=="Positive","Positive",
                              fifelse(p16=="Negative"|ISH=="Negative","Negative",NA_character_)),
                      levels=c("Negative","Positive"))]
stopifnot(sum(x$smoke=="Ever",na.rm=TRUE)==242L,sum(x$smoke=="Never",na.rm=TRUE)==91L)
sig <- c("Exposure2","Proliferation2","Stromal5")
base <- c("Age","sex","stage","site")
run <- function(rows,vars,label){
 z <- droplevels(x[rows])
 out <- rbindlist(lapply(sig,function(y){
   fit <- lm(reformulate(c("smoke",vars),y),z,na.action=na.fail)
   X <- model.matrix(fit); h <- hatvalues(fit)
   stopifnot(qr(X)$rank==ncol(X),max(h)<1-1e-8)
   w <- residuals(fit)/(1-h)
   bread <- solve(crossprod(X))
   V <- bread %*% crossprod(X*w) %*% bread
   se <- sqrt(V["smokeEver","smokeEver"])
   beta <- unname(coef(fit)["smokeEver"])
   df <- df.residual(fit)
   data.table(model=label,signature=y,N=nrow(z),ever=sum(z$smoke=="Ever"),
     never=sum(z$smoke=="Never"),difference=beta,lower95=beta-qt(.975,df)*se,
     upper95=beta+qt(.975,df)*se,p=2*pt(-abs(beta/se),df))
 }))
 out[,FDR_3_scores:=p.adjust(p,"BH")]
 out
}
cc <- function(vars) complete.cases(x[,c("smoke",sig,vars),with=FALSE])
all <- cc(NULL); clinical <- cc(base); alc <- cc(c(base,"alcohol"))
p16idx <- cc(c(base,"p16"))
hpvidx <- cc(c(base,"hpv_status"))

x[,site:=factor(fcase(grepl("Base of tongue",primary_site),"base_tongue",
                      grepl("pharynx",primary_site),"mixed",
                      grepl("tongue",primary_site),"other_tongue",
                      grepl("Floor of mouth",primary_site),"floor",
                      default="other_oral"))]
out <- rbindlist(list(run(all,NULL,"all unadjusted"),
 run(clinical,NULL,"clinical cases unadjusted"),run(clinical,base,"clinical adjusted"),
 run(alc,base,"alcohol cases clinical adjusted"),run(alc,c(base,"alcohol"),"clinical and alcohol adjusted"),
 run(clinical & x$remove_base_and_mixed,base,"exclude base and mixed adjusted"),
 run(p16idx,c(base,"p16"),"p16 complete cases additionally adjusted"),
 run(hpvidx,c(base,"hpv_status"),"any available HPV assay complete cases additionally adjusted"),
 run(hpvidx & x$remove_base_and_mixed,c(base,"hpv_status"),
     "exclude base and mixed, HPV assay complete cases adjusted")))
fwrite(out,"score_contrasts.tsv",sep="\t")
fwrite(data.table(input=c(paths,"../01_baseline/results/TCGA_HNSC_HPV_cBioPortal.tsv"),
  md5=unname(tools::md5sum(c(paths,"../01_baseline/results/TCGA_HNSC_HPV_cBioPortal.tsv")))),
  "input_hashes.tsv",sep="\t")
capture.output(sessionInfo(),file="session.txt")
