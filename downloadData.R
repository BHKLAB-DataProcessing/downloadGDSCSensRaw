library(downloader)
library(PharmacoGxPrivate)

options(stringsAsFactors=FALSE)


getGDSCrawData <-
function(path.data="/pfs/out", result.type=c("array", "list")){
  require(stringr) || stop("Library stringr is not available!")
  if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
  ## download gdsc raw sensitivity data
  archivn <- "gdsc_drug_sensitivity_raw_data_w5"
  if(!file.exists(file.path(path.data, "dwl"))){dir.create(file.path(path.data, "dwl"), showWarnings=FALSE)}
  myfn <- file.path(path.data, "dwl", sprintf("%s.csv", archivn))
  if (!file.exists(myfn)) {
    dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/%s", sprintf("%s.zip", archivn)), destfile=file.path(path.data, "dwl", sprintf("%s.zip", archivn)), quiet=TRUE)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  }
  res <- unzip(zipfile=file.path(path.data, "dwl", sprintf("%s.zip", archivn)),exdir=file.path(path.data, "dwl"))
  gdsc.raw.drug.sensitivity <- read.csv(myfn, stringsAsFactors=FALSE)
  
  concentrations.no <- length(grep("raw", names(gdsc.raw.drug.sensitivity)))
  
  if(result.type == "array"){
    ## create the gdsc.drug.response object including information viablilities and concentrations for each cell/drug pair
    obj <- array(NA,  dim=c(length(unique(gdsc.raw.drug.sensitivity[ , "CELL_LINE_NAME"])),  length(unique(gdsc.raw.drug.sensitivity[ , "DRUG_ID"])),  2,  concentrations.no), dimnames=list(unique(gdsc.raw.drug.sensitivity[ , "CELL_LINE_NAME"]), unique(gdsc.raw.drug.sensitivity[ , "DRUG_ID"]), c("concentration", "viability"), 1:concentrations.no))
  }
  
  fnexperiment <- 
  function(values)  {
    cellline <- values["CELL_LINE_NAME"]
    drug <- stringr::str_trim(values["DRUG_ID"])
    fold_dillution <- as.numeric(values["FOLD_DILUTION"])
    max_dose <- as.numeric(values["MAX_CONC"])
      doses <- rev(as.numeric(c(max_dose, unlist(lapply(1:8, function(x){max_dose/(fold_dillution ^ x)}))))) # micro molar
      if(concentrations.no > length(doses)) {doses <- c(doses, rep(NA, concentrations.no - length(doses)))}
      
      responses <- as.numeric(values[grep("raw",names(values))])
      if(concentrations.no > length(responses)) {responses <- c(responses, rep(NA, concentrations.no - length(responses)))}
      controls <- values[grep("control", names(values))]#qc_fail
      controls <- as.numeric(controls[which(!is.na(controls) & controls != "qc_fail")])
      
      blanks <- values[grep("blank", names(values))]#qc_fail
      blanks <- as.numeric(blanks[which(!is.na(blanks) & blanks != "qc_fail")])
      
      responses <- rev((responses - mean(blanks))/(mean(controls) - mean(blanks))) * 100 #mean can be replaced by median         
      
      if(result.type == "array"){
        obj[cellline,drug, "concentration", 1:length(doses)] <<- doses
        obj[cellline,drug, "viability", 1:length(responses)] <<- responses
        }else{
        return(list(cell=cellline, drug=drug, doses=doses, responses=responses)) #paste(doses, collapse=", "), responses=paste(responses, collapse=", ")))
      }
    }

    gdsc.raw.drug.sensitivity.list <- do.call(c, apply(gdsc.raw.drug.sensitivity, 1, list))
    gdsc.raw.drug.sensitivity.res <- mapply(fnexperiment, values=gdsc.raw.drug.sensitivity.list)
    if(result.type == "array"){
      return(list(data=obj, concentrations.no=concentrations.no))
    } else {
        return(list(data =gdsc.raw.drug.sensitivity.res, concentrations.no=concentrations.no))
    }
}

raw.sensitivity <- getGDSCrawData(result.type="list")
save(raw.sensitivity, file="GDSC_sens_raw.RData")

con_tested <- raw.sensitivity$concentrations.no
raw.sensitivity <- t(raw.sensitivity$data)
raw.sensitivity <- t(apply(raw.sensitivity, 1, function(x){unlist(x)}))


raw.sensitivity[ , 2] <- paste("drugid", raw.sensitivity[ , 2], sep="_")
rownames(raw.sensitivity)  <- paste(as.character(raw.sensitivity[ , 2]),raw.sensitivity[ , 1], sep="_")
raw.sensitivity <- raw.sensitivity[, -c(1, 2)]
raw.sensitivity <- array(c(as.matrix(raw.sensitivity[ , 1:con_tested]), as.matrix(raw.sensitivity[,(con_tested+1):(2*con_tested)])), c(nrow(raw.sensitivity), con_tested, 2),
                         dimnames=list(rownames(raw.sensitivity), colnames(raw.sensitivity[ , 1:con_tested]), c("Dose", "Viability")))

## TODO: make sure to add these lines to the makePSet:
# drugconc <- drugconc[rownames(raw.sensitivity),]
# duration <- rep(x=72, length=nrow(drugconc))
# sensitivityInfo <- cbind(drugconc, "duration_h"=duration)
# profiles <- profiles[rownames(sensitivityInfo),]
# # recomputed <- .calculateFromRaw(raw.sensitivity,dose.range=c(log10(2), log10(1000)), cap=100)
myfn <- file.path("/pfs/out/", "GDSC_sens_recomputed.RData")
if(!file.exists(myfn)){
  recomputed <- .calculateFromRaw(raw.sensitivity, cap=100)
  save(recomputed, file=myfn)
} else {
  load(myfn, verbose=TRUE)
}
