# This file creates datasets from NMEG Ameriflux files.
#
# Greg Maurer - Nov 23, 2014

setwd('~/current/NMEG_fluxdata/data_analysis/R/')

source('printfigs.r')
library(ggplot2)
library(plyr)

proc_path <- '../processed_data/'
raw_AF_path <- '~/sftp/eddyflux/Ameriflux_files_GM/2009-2013_Reichstein'
files <- list.files(raw_AF_path, full.names=TRUE)

PJC_files_gf <- files[grepl("gapfilled", files) & grepl("US-Mpj", files)]
PJG_files_gf <- files[grepl("gapfilled", files) & grepl("US-Mpg", files)]

readAF <- function(newfile) {
  newheader <- read.csv(newfile, skip=3, nrows=1)
  new <- read.csv(newfile, skip=5, header=FALSE)
  new[new==-9999] <- NA
  colnames(new) <- colnames(newheader)
  return(new)
}

par(mfrow=c(5,2))

for (i in 1:length(PJC_files_gf)) {
    
  # Read in files
  PJCfile <- PJC_files_gf[i]
  PJGfile <- PJG_files_gf[i]
  
  pjc <- readAF(PJCfile)
  pjg <- readAF(PJGfile)
  
  toks1 <- strsplit(PJCfile, '/')
  fname <- toks1[[1]][[8]]
  toks2 <- strsplit(fname, '_')
  year <- toks2[[1]][2]
  
  in_model <- lm(pjc$Rlong_in ~ pjg$Rlong_in)
  in_coef <- round(coef(in_model), 3) # extract coefficients 
  plot(pjc$Rlong_in ~ pjg$Rlong_in)
  abline(in_model, col='red')
  abline(0,1, col='blue')
  mtext(bquote(slope == .(lm_coef[2]), intercept == .(lm_coef[1])),
        adj=1, padj=0)
  
  out_model <- lm(pjc$Rlong_out ~ pjg$Rlong_out)
  out_coef <- round(coef(out_model), 3) # extract coefficients 
  plot(pjc$Rlong_out ~ pjg$Rlong_out)
  abline(out_model, col='red')
  abline(0,1, col='blue')
}
  
  # Calculate annual sums
  GPPgCm_30min <- new$GPP * (12.011/1e+06) * 1800
  df_reich$GPPsum[i] <- sum(GPPgCm_30min, na.rm = TRUE)
  
  REgCm_30min <- new$RE * (12.011/1e+06) * 1800
  df_reich$REsum[i] <- sum(REgCm_30min, na.rm = TRUE)
  
  NEEgCm_30min <- new$FC * (12.011/1e+06) * 1800
  df_reich$NEEsum[i] <- sum(NEEgCm_30min, na.rm = TRUE)
  
  ETmm_30min <- new$FH2O * (1800/1e+06)
  df_reich$ETsum[i] <- sum(ETmm_30min, na.rm=TRUE)
  
  df_reich$Psum[i] <- sum(new$PRECIP, na.rm=TRUE)
}

for (i in 1:length(lass_files_gf)) {
  
  # Read in files
  file <- lass_files_gf[i]
  header <- read.csv(file, skip=3, nrows=1)
  new <- read.csv(file, skip=5, header=FALSE)
  new[new==-9999] <- NA
  colnames(new) <- colnames(header)
  
  toks1 <- strsplit(file, '/')
  fname <- toks1[[1]][[8]]
  toks2 <- strsplit(fname, '_')
  df_lass$year[i] <- toks2[[1]][2]
  df_lass$site[i] <- toks2[[1]][1]
  # Calculate annual sums
  GPPgCm_30min <- new$GPP * (12.011/1e+06) * 1800
  df_lass$GPPsum[i] <- sum(GPPgCm_30min, na.rm = TRUE)
  
  REgCm_30min <- new$RE * (12.011/1e+06) * 1800
  df_lass$REsum[i] <- sum(REgCm_30min, na.rm = TRUE)
  
  NEEgCm_30min <- new$FC * (12.011/1e+06) * 1800
  df_lass$NEEsum[i] <- sum(NEEgCm_30min, na.rm = TRUE)
  
  ETmm_30min <- new$FH2O * (1800/1e+06)
  df_lass$ETsum[i] <- sum(ETmm_30min, na.rm = TRUE)
  
  df_lass$Psum[i] <- sum(new$PRECIP, na.rm=TRUE)
}

df_reich$Part <- 'Reichstein'
df_lass$Part <- 'Lasslop'

allsums <- rbind(df_reich, df_lass)

allsums$site2 <- 
  revalue(allsums$site, c("US-Mpg"="PJ_girdle",
                                 "US-Mpj"="PJ",
                                 "US-Seg"="GLand",
                                 "US-Sen"="New_GLand",
                                 "US-Ses"="SLand",
                                 "US-Wjs"="JSav",
                                 "US-Vcp"="PPine",
                                 "US-Vcm"="MCon"))

reorder <- c('GLand', 'New_GLand', 'SLand', 'PJ', 'PJ_girdle', 'JSav',
             'PPine', 'MCon')

# Make a new column that includes the citation (made for x axis labeling)
allsums$site2 <- factor(allsums$site2, levels=reorder)

write.table(subset(allsums, Part=='Reichstein'), paste(proc_path, 'export.txt'),
            row.names=F)


REplot <- ggplot(allsums, aes(x=year, y=REsum)) +
  geom_point(aes(colour = factor(Part)), alpha=0.6, size=4) + 
  facet_wrap(~site2) + guides(colour=guide_legend(title="Partitioning"))

GPPplot <- ggplot(allsums, aes(x=year, y=GPPsum)) +
  geom_point(aes(colour = factor(Part)), alpha=0.6, size=4) + 
  facet_wrap(~site2) + guides(colour=guide_legend(title="Partitioning"))

ETplot <- ggplot(allsums, aes(x=year, y=ETsum)) +
  geom_point(aes(colour = factor(Part)), alpha=0.6, size=4) + 
  facet_wrap(~site2) + guides(colour=guide_legend(title="Partitioning"))

Pplot <- ggplot(allsums, aes(x=year, y=Psum)) +
  geom_point(aes(colour = factor(Part)), alpha=0.6, size=4) + 
  facet_wrap(~site2) + guides(colour=guide_legend(title="Partitioning"))

NEEplot <- ggplot(allsums, aes(x=year, y=NEEsum)) +
  geom_point(aes(colour = factor(Part)), alpha=0.6, size=4) + 
  facet_wrap(~site2) + guides(colour=guide_legend(title="Partitioning"))

printfigs(REplot, "RE_part.eps", 8, 6)
printfigs(REplot, "RE_part.svg", 8, 6)
printfigs(GPPplot, "GPP_part.eps", 8, 6)
printfigs(GPPplot, "GPP_part.svg", 8, 6)
printfigs(ETplot, "ET_part.eps", 8, 6)
printfigs(ETplot, "ET_part.svg", 8, 6)
printfigs(Pplot, "P_part.eps", 8, 6)
printfigs(Pplot, "P_part.svg", 8, 6)
printfigs(NEEplot, "NEE_part.eps", 8, 6)
printfigs(NEEplot, "NEE_part.svg", 8, 6)
# Function for mean centering an array by a factor
#source('~/data/code_resources/r_common/ctr.byfactor.r')

# Load plot, tag, and disturbance year data
new <- read.csv(reich_files[1], skip=3)
new <- new[2:nrow(new),]


# Load comparison 13C data for niwot
nwt13c  <- read.csv(paste(datapath, 'niwot_13c.dat', sep=''),header=TRUE)

# Load surface chamber data
surfchamb  <- read.csv(paste(datapath, 'surf_chamb.dat', sep=''),header=TRUE)

# Load ecosystem measurements
ecosflux  <- read.csv(paste(datapath, 'ecosys_flux.dat', sep=''),header=TRUE)

# Load soil environment data and convert dates
soilenv  <- read.csv(paste(datapath, 'soilenv.dat', sep=''),
                        header=TRUE)
soilenv$samp_date <- as.POSIXlt(strptime(soilenv$samp_date, "%y%m%d"))
soilenv$SAMP_DATE <- as.factor(as.character(soilenv$samp_date))

# Load extract data and convert dates
extractdat  <- read.csv(paste(datapath, 'soilextracts.dat', sep=''),skip=1,
                        header=TRUE)
extractdat$samp_date <- as.POSIXlt(strptime(extractdat$samp_date, "%y%m%d"))
extractdat$SAMP_DATE <- as.factor(as.character(extractdat$samp_date))
# Remove the blanks
extractdat <- subset(extractdat, type!='B')
# Should we subtract the background values from blanks?

# Merge in the tag and disturbance data
extractdat <- merge(extractdat, tagdata, by='plotnum', all.x=TRUE,
                    incomparables=NA)

# change distyear column name, then calculate years since disturbance
extractdat$tdist <- extractdat$samp_date$year - extractdat$distyear + 1900
extractdat$tdist[extractdat$tdist < -1] <- 0 # Fix controls
extractdat$tdist_f <- as.factor(extractdat$tdist) # Time since disturbance factor
levels(extractdat$tdist_f)[levels(extractdat$tdist_f)=='0.3'] <- '<1'

extractdat$TOC = ((extractdat$TOC/1000) * 25)/5 * 1000
extractdat$TN = ((extractdat$TN/1000) * 25)/5 * 1000
#write.csv(extractdat, '~/Desktop/extractTOCTN.dat', col.names=T, row.names=F)
# Now create a dataset with microbial C and N
microbe_dat <- merge(subset(extractdat, trmnt=='U'), 
                      subset(extractdat, trmnt=='F'),
                      by=c('samp_date', 'plotnum', 'type'), all.x=T)
                      #incomparables=NA)
microbe_dat <- within(microbe_dat, {TOC_m <- TOC.y - TOC.x; TN_m <- TN.y-TN.x})
microbe_dat <- subset(microbe_dat, select=c('samp_date', 'SAMP_DATE.y', 'plotnum', 'TOC_m',
                                            'TN_m', 'distyear.x'))
# change distyear column name, then calculate years since disturbance
colnames(microbe_dat)[2] <- 'SAMP_DATE'
colnames(microbe_dat)[6] <- 'distyear'
microbe_dat$tdist <- microbe_dat$samp_date$year - microbe_dat$distyear + 1900
microbe_dat$tdist[microbe_dat$tdist < -1] <- 0 # Fix controls
microbe_dat$tdist_f <- as.factor(microbe_dat$tdist) # Time since disturbance factor
levels(microbe_dat$tdist_f)[levels(microbe_dat$tdist_f)=='0.3'] <- '<1'

# Load flask data and bind into one dataframe
nwt11_flasks  <- read.csv(paste(datapath, 'nwt_2011_flaskdata_all.csv', sep=''))
nwt12_flasks  <- read.csv(paste(datapath, 'nwt_2012_flaskdata_all.csv', sep=''))
nwt_flasks <- rbind(nwt11_flasks, nwt12_flasks)
# Change a few columns
colnames(nwt_flasks)[3] <- 'horizon'
nwt_flasks$horizon <- as.factor(nwt_flasks$horizon)
levels(nwt_flasks$horizon) <- c('Organic', 'Mineral')
rm(nwt11_flasks, nwt12_flasks)
# Create date vectors
nwt_flasks$samp_date <- as.POSIXlt(strptime(nwt_flasks$samp_date, "%y%m%d"))
nwt_flasks$SAMP_DATE <- as.factor(as.character(nwt_flasks$samp_date))
nwt_flasks$samp_year <- nwt_flasks$samp_date$year + 1900
# Now add tag data and calculate time since disturbance
nwt_flasks <- merge(nwt_flasks, tagdata, by='plotnum', all.x=TRUE)
                    #incomparables=NA)
nwt_flasks$tdist <- nwt_flasks$samp_date$year - nwt_flasks$distyear + 1900
nwt_flasks$tdist[nwt_flasks$tdist == 0] <- 0.3 # Fix sites girdled/sampled in 2011 
nwt_flasks$tdist[nwt_flasks$tdist < -1] <- 0 # Fix controls
nwt_flasks$tdist_f <- as.factor(nwt_flasks$tdist) # Time since disturbance factor
levels(nwt_flasks$tdist_f)[levels(nwt_flasks$tdist_f)=='0.3'] <- '<1'


# Load bulk soil C and 13C data
bulksoil  <- read.csv(paste(datapath, 'nwt_soilCiso_all.csv', sep=''))
bulksoil$samp_date <- as.POSIXlt(strptime(bulksoil$samp_date, "%y%m%d"))
bulksoil$SAMP_DATE <- as.factor(as.character(bulksoil$samp_date))
bulksoil$samp_year <- bulksoil$samp_date$year + 1900
# Change a few columns
bulksoil$horizon <- as.factor(bulksoil$horizon)
levels(bulksoil$horizon) <- c('Organic', 'Mineral')
# Now add tag data and calculate time since disturbance
bulksoil <- merge(bulksoil, tagdata, by='plotnum', all.x=TRUE)
                  #  incomparables=NA)
bulksoil$tdist <- bulksoil$samp_date$year - bulksoil$distyear + 1900
bulksoil$tdist[bulksoil$tdist == 0] <- 0.3 # Fix sites girdled/sampled in 2011 
bulksoil$tdist[bulksoil$tdist < -1] <- 0 # Fix controls
bulksoil$tdist_f <- as.factor(bulksoil$tdist) # Time since disturbance factor
levels(bulksoil$tdist_f)[levels(bulksoil$tdist_f)=='0.3'] <- '<1'

# Load soil CO2 data and make proper date columns
nwt_dogbone  <- read.csv(paste(datapath, 'NWT_dogbones_clean.dat', sep=''))
#nwt_dogbone$SAMP_DATE <- as.factor(nwt_dogbone$samp_date)
nwt_dogbone$samp_date <- as.POSIXlt(strptime(nwt_dogbone$samp_date, "%y%m%d"))
nwt_dogbone$SAMP_DATE <- as.factor(as.character(nwt_dogbone$samp_date))
nwt_dogbone$samp_year <- nwt_dogbone$samp_date$year + 1900

fef_dogbone  <- read.csv(paste(datapath, 'FEF_dogbones_clean.dat', sep=''))
#fef_dogbone$SAMP_DATE <- as.factor(fef_dogbone$samp_date)
fef_dogbone$samp_date <- as.POSIXlt(strptime(fef_dogbone$samp_date, "%y%m%d"))
fef_dogbone$SAMP_DATE <- as.factor(as.character(fef_dogbone$samp_date))
fef_dogbone$samp_year <- fef_dogbone$samp_date$year + 1900

# Remove 2013 data
nwt_dogbone <- subset(nwt_dogbone, samp_year!=2013)
fef_dogbone <- subset(fef_dogbone, samp_year!=2013)

# Reshape the tagdata to make tag numbers accessible
library(reshape2)
mtags <- melt(tagdata, id=c("site","plotnum", "distyear", "distyear2"))
colnames(mtags)[5:6] <- c('depth', 'tag')
mtags$depth <- revalue(mtags$depth,
                       c("tag_oa"="0", "tag_10"="10", "tag_30"="30"))
# Merge in these tag numbers
nwt_dogbone <- merge(nwt_dogbone, mtags, by=c('plotnum','tag'),
                     all.x=TRUE)#,incomparables=NA)
fef_dogbone <- merge(fef_dogbone, mtags, by=c('plotnum','tag'),
                     all.x=TRUE)#,incomparables=NA)

# change distyear column name, then calculate years since disturbance
nwt_dogbone$tdist <- nwt_dogbone$samp_date$year - nwt_dogbone$distyear + 1900
nwt_dogbone$tdist[nwt_dogbone$tdist == 0] <- 0.3 # Fix sites girdled/sampled in 2011 
nwt_dogbone$tdist[nwt_dogbone$tdist < -1] <- 0 # Fix controls
nwt_dogbone$tdist_f <- as.factor(nwt_dogbone$tdist) # Time since disturbance factor
# Rename factors
levels(nwt_dogbone$tdist_f)[levels(nwt_dogbone$tdist_f)=='0.3'] <- '<1'
nwt_dogbone$depth <- revalue(nwt_dogbone$depth, c("0"="0 cm", "10"="10 cm", "30"="30 cm"))

fef_dogbone$tdist <- fef_dogbone$samp_date$year - fef_dogbone$distyear + 1900
fef_dogbone$tdist[fef_dogbone$tdist < -1] <- 0 # Fix controls
fef_dogbone$tdist_f <- as.factor(fef_dogbone$tdist) # Time since disturbance factor
fef_dogbone$Method <- as.factor("Degradation")
# Rename factors
levels(fef_dogbone$tdist_f)[levels(fef_dogbone$tdist_f)=='0.3'] <- '<1'
fef_dogbone$depth <- revalue(fef_dogbone$depth, c("0"="0 cm", "10"="10 cm", "30"="30 cm"))

# The time since death at FEF was reassessed
fef_dogbone2 <- fef_dogbone
fef_dogbone2$tdist <- fef_dogbone2$samp_date$year - fef_dogbone2$distyear2 + 1900
fef_dogbone2$tdist[fef_dogbone2$tdist < -1] <- 0 # Fix controls
fef_dogbone2$tdist_f <- as.factor(fef_dogbone2$tdist) # Time since disturbance factor
fef_dogbone2$Method <- as.factor("Dendrochronology")

# Calculate deltaJ
air_CO2 = 420; #ppm
air_d13C = -9.5; #permil
nwt_dogbone <- within(nwt_dogbone, delJ <- (co2_ppm*(d13c - 4.4) - air_CO2*(air_d13C - 4.4))
                       /(1.0044*(co2_ppm - air_CO2)))
fef_dogbone <- within(fef_dogbone, delJ <- (co2_ppm*(d13c - 4.4) - air_CO2*(air_d13C - 4.4))
                       /(1.0044*(co2_ppm - air_CO2)))

fef_dogbone2 <- within(fef_dogbone2, delJ <- (co2_ppm*(d13c - 4.4) - air_CO2*(air_d13C - 4.4))
                      /(1.0044*(co2_ppm - air_CO2)))

fef_dogbone_meth <- rbind(fef_dogbone, fef_dogbone2)
#del_j_fef = (fef_CO2.*(fef_d13C - 4.4) - air_CO2.*(air_d13C - 4.4))...
#./(1.0044.*(fef_CO2 - air_CO2));