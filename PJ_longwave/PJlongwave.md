Examination of Radiation at PJ and PJG


```r
setwd("~/current/NMEG_fluxdata/data_analysis/R/")

source("printfigs.r")

proc_path <- "../processed_data/"
AF_path <- "~/sftp/eddyflux/Ameriflux_files_GM/2007-2013_Reichstein"
AFfiles <- list.files(AF_path, full.names = TRUE)

PJC_AF_gf <- AFfiles[grepl("gapfilled", AFfiles) & grepl("US-Mpj", AFfiles)]
PJG_AF_gf <- AFfiles[grepl("gapfilled", AFfiles) & grepl("US-Mpg", AFfiles)]

PJC_fall_path <- "~/sftp/eddyflux/Flux_all_files/PJ"
PJG_fall_path <- "~/sftp/eddyflux/Flux_all_files/PJ_girdle"
PJC_fallfiles <- list.files(PJC_fall_path, full.names = TRUE)
PJG_fallfiles <- list.files(PJG_fall_path, full.names = TRUE)
# Trim to 2009-2013
PJC_fall <- PJC_fallfiles[grepl("20[0-1][9,0-3]", PJC_fallfiles)]
PJG_fall <- PJG_fallfiles[grepl("20[0-1][9,0-3]", PJG_fallfiles)]

readAF <- function(newfile) {
    newheader <- read.csv(newfile, skip = 3, nrows = 1)
    new <- read.csv(newfile, skip = 5, header = FALSE)
    new[new == -9999] <- NA
    colnames(new) <- colnames(newheader)
    return(new)
}

readfall <- function(newfile) {
    # newheader <- read.csv(newfile, skip=0, nrows=1)
    new <- read.table(newfile, header = TRUE, sep = "\t")
    new[new == -9999] <- NA
    # colnames(new) <- colnames(newheader)
    return(new)
}
```


Now plot regression lines for each year - (Rlong_in and Rlong_out)


```r

pdf("../figures/PJ_longwavecompare.pdf", width = 7, height = 12.5)
par(mfrow = c(5, 2))

for (i in 1:length(PJG_AF_gf)) {
    
    # Read in files for now there are 2 more PJC files
    PJCfile <- PJC_AF_gf[i + 2]
    PJGfile <- PJG_AF_gf[i]
    
    pjc <- readAF(PJCfile)
    pjg <- readAF(PJGfile)
    
    toks1 <- strsplit(PJCfile, "/")
    fname <- toks1[[1]][[8]]
    toks2 <- strsplit(fname, "_")
    year <- toks2[[1]][2]
    
    in_model <- lm(pjc$Rlong_in ~ pjg$Rlong_in)
    in_coef <- round(coef(in_model), 3)  # extract coefficients 
    plot(pjc$Rlong_in ~ pjg$Rlong_in, main = paste(year, "Rlong incoming"), 
        ylab = "Control (W/m^2)", xlab = "Girdle (W/m^2)", pch = 18, col = "gray")
    # legend('bottomright', legend='obs', pch=18, col='gray')
    abline(in_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    legend("bottomright", legend = c("Fit", "1:1 line"), lwd = 2, col = c("red", 
        "blue"))
    text(180, 370, paste("Slope =", in_coef[[2]], "\n Int = ", in_coef[[1]]))
    
    out_model <- lm(pjc$Rlong_out ~ pjg$Rlong_out)
    out_coef <- round(coef(out_model), 3)  # extract coefficients 
    plot(pjc$Rlong_out ~ pjg$Rlong_out, main = paste(year, "Rlong outgoing"), 
        ylab = "Control (W/m^2)", xlab = "Girdle (W/m^2)", pch = 18, col = "gray")
    abline(out_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    text(285, 550, paste("Slope =", out_coef[[2]], "\n Int = ", out_coef[[1]]))
}

dev.off()
```

```
## pdf 
##   2
```


Same thing for shortwave


```r
pdf("../figures/PJ_shortwavecompare.pdf", width = 7, height = 12.5)
par(mfrow = c(5, 2))

for (i in 1:length(PJG_AF_gf)) {
    
    # Read in files for now there are 2 more PJC files
    PJCfile <- PJC_AF_gf[i + 2]
    PJGfile <- PJG_AF_gf[i]
    
    pjc <- readAF(PJCfile)
    pjg <- readAF(PJGfile)
    
    toks1 <- strsplit(PJCfile, "/")
    fname <- toks1[[1]][[8]]
    toks2 <- strsplit(fname, "_")
    year <- toks2[[1]][2]
    
    in_model <- lm(pjc$Rg ~ pjg$Rg)
    in_coef <- round(coef(in_model), 3)  # extract coefficients 
    plot(pjc$Rg ~ pjg$Rg, main = paste(year, "Rg incoming"), ylab = "Control (W/m^2)", 
        xlab = "Girdle (W/m^2)", pch = 18, col = "gray")
    # legend('bottomright', legend='obs', pch=18, col='gray')
    abline(in_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    legend("bottomright", legend = c("Fit", "1:1 line"), lwd = 2, col = c("red", 
        "blue"))
    text(180, 1000, paste("Slope =", in_coef[[2]], "\n Int = ", in_coef[[1]]))
    
    out_model <- lm(pjc$Rg_out ~ pjg$Rg_out)
    out_coef <- round(coef(out_model), 3)  # extract coefficients 
    plot(pjc$Rg_out ~ pjg$Rg_out, main = paste(year, "Rg outgoing"), ylab = "Control (W/m^2)", 
        xlab = "Girdle (W/m^2)", pch = 18, col = "gray")
    abline(out_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    text(55, 200, paste("Slope =", out_coef[[2]], "\n Int = ", out_coef[[1]]))
}

dev.off()
```

```
## pdf 
##   2
```



Now plot regression lines for Raw (fluxall) data in each year - (Rlong_in and Rlong_out)


```r
pdf("../figures/PJ_longwavecompare_raw.pdf", width = 9, height = 12.5)
par(mfrow = c(5, 3))

for (i in 1:length(PJG_fall)) {
    
    # Read in files for now there are 2 more PJC files
    PJCfile <- PJC_fall[i]
    PJGfile <- PJG_fall[i]
    
    pjc <- readfall(PJCfile)
    pjg <- readfall(PJGfile)
    
    # Some fluxall files are bigger than others if (nrow(pjc) > nrow(pjg)){ ds
    # <- merge(pjc, pjg, all.y=T, by=jday) } else if (nrow(pjg) > nrow(pjc)) {
    # ds <- merge(pjc, pjg, all=T, by=jday) }
    
    toks1 <- strsplit(PJCfile, "/")
    fname <- toks1[[1]][[8]]
    toks2 <- strsplit(fname, "_")
    toks3 <- strsplit(toks2[[1]][[4]], "\\.")
    year <- toks3[[1]][1]
    
    sh_jdays <- intersect(pjc$jday, pjg$jday)
    
    df <- data.frame(jday = sh_jdays)
    cind <- match(df$jday, pjc$jday)
    gind <- match(df$jday, pjg$jday)
    
    if (year > 2010 & year < 2013) {
        df$in_c <- pjc$Rad_long_Up__Avg[cind]
        df$out_c <- pjc$Rad_long_Dn__Avg[cind]
        df$t_c <- pjc$CNR1TC_Avg[cind]
    } else {
        df$in_c <- pjc$Rad_long_Up_Avg[cind]
        df$out_c <- pjc$Rad_long_Dn_Avg[cind]
        df$t_c <- pjc$CNR1TC_Avg[cind]
    }
    df$in_g <- pjg$Rad_long_Up_Avg[gind]
    df$out_g <- pjg$Rad_long_Dn_Avg[gind]
    if (year > 2012) {
        df$t_g <- pjg$CNR1TC_Avg[gind]
    } else {
        df$t_g <- pjg$Temp_C_Avg[gind]
    }
    
    
    in_model <- lm(df$in_c ~ df$in_g)
    in_coef <- round(coef(in_model), 3)  # extract coefficients 
    plot(df$in_c ~ df$in_g, main = paste(year, "Rlong incoming (uncorrected)"), 
        ylab = "Control (W/m^2)", xlab = "Girdle (W/m^2)", pch = 18, col = "gray")
    # legend('bottomright', legend='obs', pch=18, col='gray')
    abline(in_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    legend("bottomright", legend = c("Fit", "1:1 line"), lwd = 2, col = c("red", 
        "blue"))
    text(-180, -20, paste("Slope =", in_coef[[2]], "\n Int = ", in_coef[[1]]))
    
    out_model <- lm(df$out_c ~ df$out_g)
    out_coef <- round(coef(out_model), 3)  # extract coefficients 
    plot(df$out_c ~ df$out_g, main = paste(year, "Rlong outgoing (uncorrected)"), 
        ylab = "Control (W/m^2)", xlab = "Girdle (W/m^2)", pch = 18, col = "gray")
    abline(out_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    text(-20, 80, paste("Slope =", out_coef[[2]], "\n Int = ", out_coef[[1]]))
    
    T_model <- lm(df$t_c ~ df$t_g)
    T_coef <- round(coef(T_model), 3)  # extract coefficients 
    plot(df$t_c ~ df$t_g, main = paste(year, "CNR1 body temp"), ylab = "Control (^oC)", 
        xlab = "Girdle (^oC)", pch = 18, col = "gray")
    abline(T_model, col = "red", lwd = 2)
    abline(0, 1, col = "blue", lwd = 2)
    text(0, 20, paste("Slope =", T_coef[[2]], "\n Int = ", T_coef[[1]]))
}

dev.off()
```

```
## pdf 
##   2
```


