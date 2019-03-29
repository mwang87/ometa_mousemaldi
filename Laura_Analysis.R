library(MALDIquant)
library(MALDIquantForeign)
library(ggfortify)
library(rbokeh)
library(parallel)
library(data.table)
library(lsa)
prepare_meta_data <- function(){
  files <- dir('~/DataNAS/LauraData/March_mzML/Mouse_lavages_OVCAR8RFP_02222018_valya/', '*.mzML')
}

#analysis using MALDI-quant
run_analysis_with_maldi_quant <- function(){
  #using json_file
  json_file <- '/home/jian/DataNAS/LauraData/Json/Galaxy44.json'
  raw <- fromJSON(json_file)
  name <- names(raw)
  ms_spect_list <- lapply(1:length(raw),
                          function(i){
                            r <- raw[[i]]
                            s <- createMassSpectrum(mass=r[,1], intensity= r[,2], metaData=list(name=name[[i]]))
                          }
  )


  #reading mzml
  data_path <- '~/DataNAS/LauraData/CellLine_Mixtures/MTWT_OVCAR8_mixtures/MTWT_OVCAR8_mix_05092018/'
  meta_data <- read_meta_data('~/DataNAS/LauraData/CellLine_Mixtures/MTWT_OVCAR8_mixtures/MTWT_OVCAR8_mix_05092018/Copy of mtwt_ov8_05092018.txt', F)

  data_path2 <- '~/DataNAS/LauraData/CellLine_Mixtures/MTWT_OVCAR8_mixtures/MTWT_OVCAR8_mix_05022018/'
  meta_data2 <- read_meta_data('~/DataNAS/LauraData/CellLine_Mixtures/MTWT_OVCAR8_mixtures/MTWT_OVCAR8_mix_05022018//Copy of mtwt_ov8_05022018.txt', F)

  ms_spect_list <- import(data_path)
  ms_spect_list2 <- import(data_path2)

  ms_spect_list <- c(ms_spect_list, ms_spect_list2)


  #pre-processing MALDI spectra
  plot(ms_spect_list[[1]])
  spectra <- transformIntensity(ms_spect_list, method="sqrt")
  plot(spectra[[1]])
  spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
  plot(spectra[[1]])

  spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
  plot(spectra[[1]])
  spectra <- calibrateIntensity(spectra, method="TIC")
  plot(spectra[[1]])
  spectra <- alignSpectra(spectra, halfWindowSize=20, SNR=2, tolerance=0.002, warpingMethod="lowess")
  plot(spectra[[1]])

  noise <- estimateNoise(spectra[[1]])
  plot(spectra[[1]], xlim=c(5000, 20000))
  lines(noise, col="red")
  lines(noise[,1], noise[, 2]*2, col="blue")

  peaks <- detectPeaks(spectra, method="MAD", halfWindowSize=20, SNR=2)
  plot(spectra[[1]])
  points(peaks[[1]], col="red", pch=4)

  peaks <- binPeaks(peaks, tolerance=0.002)

  filenames <- basename(unlist(lapply(peaks, function(p){attributes(p)$metaData$file})))
  featureMatrix <- intensityMatrix(peaks, spectra)
  #perform total intensity normalization
  s <- rowSums(featureMatrix)
  Ncol <- ncol(featureMatrix)
  a <- matrix(rep(s, Ncol), ncol = Ncol)
  featureMatrix <- featureMatrix/a

  attributes(featureMatrix)$dimnames[[1]] <- basename(filenames)

  df <- as.data.frame(featureMatrix)
  feature_names <- paste0("FEATURE_", seq(1,nrow(df)), "_", colnames(df), "_1.0")
  colnames(df) <- feature_names

  dt <- data.table(filename=filenames)
  dt <- cbind(dt, df)
  feature_mat <- merge(dt, meta_data, by='filename')
  dat_ind <- max(which(grepl('FEATURE_', colnames(feature_mat))))


  #performing PCA
  df <- as.data.frame(featureMatrix)

  pdf('~/workspace/Laura_data_analysis/PCA_plots.pdf')
  ms.pca <- prcomp(featureMatrix,
                   center = TRUE,
                   scale. = TRUE)
  plot(ms.pca)
  feature_mat$ATTRIBUTE_SampleGroupID <- as.character(feature_mat$ATTRIBUTE_SampleGroupID)
  feature_mat$ATTRIBUTE_LAV <- as.character(feature_mat$ATTRIBUTE_LAV)
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_LAV')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_DAY')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_PLATE_ID')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_SampleGroupID')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_OV8')

  #compute reproducibility
  cv <- compute_feature_abudnance_cv_by_group(feature_mat, groups = feature_mat$ATTRIBUTE_SampleGroupID, min.int = 0.000001)
  cv$Group <- as.character(cv$Group)
  hist(cv$presence)
  dat <- as.matrix(feature_mat[,2:dat_ind])

  #reproducibility per group
  g <- ggplot(data=cv)
  g <- g + stat_boxplot(aes(x=Group, y=CV))
  g

  #distribution of feature intensity
  hist(log10(featureMatrix), 100)


  #separating two type of cell lines
  #identifying features that are abundant in OV8
  int_MTWT <- feature_mat[ATTRIBUTE_OV8==0, colSums(.SD), .SDcols=2:dat_ind]
  int_OV8 <- feature_mat[ATTRIBUTE_OV8==100, colSums(.SD), .SDcols=2:dat_ind]
  feature_mass <- colnames(df)
  ratio <- int_OV8 / int_MTWT
  d <- data.frame(x=as.factor(colnames(df)), y=ratio)
  p <- figure(width = 1500) %>%
    ly_points(x = x, y=y, data = d) %>%
    tool_wheel_zoom(dimensions = c("width"))

  inds_OV8 <- which(ratio > 5)

  ind <- inds_OV8[[1]]
  d <- data.frame(x=feature_mat$ATTRIBUTE_OV8, y=log10(feature_mat[[ind+1]]))

  p <- figure(width = 1500) %>%
    ly_points(x = x, y=y, data = d) %>%
    tool_wheel_zoom(dimensions = c("width"))

  #testing for normality of data
  OV8 <- unique(feature_mat$ATTRIBUTE_OV8)
  OV8 <- sort(OV8[OV8 > 0])
  nt <- mclapply(2:dat_ind,
                 function(i){
                   p <- lapply(OV8,
                               function(o){
                                 I2 <- feature_mat[ATTRIBUTE_OV8==o][[i]]
                                 I2 <- log10(I2+1e-11)
                                 m1 <- median(I2)
                                 m2 <- mad(I2)*3
                                 I2 <- I2[I2 > m1-m2 | I2 < m1+m2]
                                 t <- shapiro.test(I2)
                                 return(list(p_value=t$p.value, m1=median(I2)))
                               })
                   dt <- data.table(rbindlist(p))
                   dt$OV8=OV8
                   dt$mass=feature_mass[i-1]
                   return(dt)
                 }, mc.cores = 20
  )
  nt <- rbindlist(nt)


  #perform wilcoxon rank sum test on per feature basis
  OV8 <- unique(feature_mat$ATTRIBUTE_OV8)
  OV8 <- sort(OV8[OV8 > 0])
  nF <- 24
  wt <- mclapply(2:dat_ind,
                  function(i){
                    I1 <- feature_mat[ATTRIBUTE_OV8==0][[i]]
                    p <- lapply(OV8,
                           function(o){
                             I2 <- feature_mat[ATTRIBUTE_OV8==o][[i]]
                             #t <- wilcox.test(I1[1:nF], I2[1:nF])
                             t <- t.test(log10(I1[1:nF]+1e-11), log10(I2[1:nF]+1e-11))
                             return(list(p_value=t$p.value, m1=median(I1), m2=median(I2)))
                           })
                    dt <- data.table(rbindlist(p))
                    dt$OV8=OV8
                    dt$mass=feature_mass[i-1]
                    return(dt)
                  }, mc.cores = 20
               )

  wt <- rbindlist(wt)
  wt$direction <- ifelse(wt$m1 < wt$m2, 1, -1)

  #performing multiple hypothesis correction
  methods <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY")

  ret_list <- lapply(1:length(methods),
                     function(i){
                        wt[, p_value := p.adjust(p_value_orig, method = methods[[i]]), by=OV8]
                        #getting feature that are differetially expressed features
                        OV8_ft <- wt[OV8==100]
                        OV8_limit <- wt[, max(which(p_value > 0.01), 0), by=mass]
                        OV8_limit <- cbind(OV8_ft, OV8_limit)
                        histc <- hist(OV8_limit$V1, breaks=seq(-1,100,1))
                        xlabs <- as.integer(histc$mids + 1.5)
                        counts <- cumsum(histc$counts)
                        data.table(percent=1:length(counts), counts=counts, method=methods[[i]])
                     })
  ret_list <- rbindlist(ret_list)

  pdf('/home/jian/NAS/workspace/Laura_data_analysis/Differentially_detected_feature_vs_percent_OV8_multiple_correction.pdf')

  png('/home/jian/NAS/workspace/Laura_data_analysis/Differentially_detected_feature_vs_percent_OV8_multiple_correction.png', width = 1000, height = 1000)
  g <- ggplot(ret_list) + geom_bar(aes(x=percent, y=counts, color=method, group=method), stat='identity')
  g <- ggplot(ret_list) + geom_line(aes(x=percent, y=counts, color=method, group=method), stat='identity')
  g <- g + geom_point(aes(x=percent, y=counts, color=method, group=method), stat='identity') + scale_x_continuous(limit=c(0.5,10), breaks=seq(0,10,1)) + scale_y_continuous(limit=c(0,150))
  g <- g + theme(text = element_text(size=15))
  g <- g + xlab("% of OV8 Cancerl cell") + ylab("Number of differetially expressed features") + ggtitle("")
  plot(g)
  dev.off()

  #counting features detection limits
  c_n24 <- hist(OV8_limit_n24$V1, seq(-1,19,1))
  c_n18 <- hist(OV8_limit_n18$V1, seq(-1,19,1))
  c_n16 <- hist(OV8_limit_n16$V1, seq(-1,19,1))
  c_n12 <- hist(OV8_limit_n12$V1, seq(-1,19,1))

  bins <- c(1:10, seq(20,110,10))
  d <- data.frame(counts=cumsum(c_n24$counts), bin=bins, rep="24 replicates")
  #d <- rbind(d, data.frame(counts=cumsum(c_n18$counts), bin=bins, rep="18 replicates"))
  d <- rbind(d, data.frame(counts=cumsum(c_n16$counts), bin=bins, rep="16 replicates"))
  d <- rbind(d, data.frame(counts=cumsum(c_n12$counts), bin=bins, rep="12 replicates"))

  g <- ggplot(data = d, aes(x=bin, y=counts, group=rep, color=rep)) + geom_line() + geom_point()
  g <- g + scale_x_continuous(breaks =seq(1,10,1),limit=c(1,10)) + scale_y_continuous(limit=c(0,50))
  g <- g + xlab("Minimum % of OVCAR8 detected") + ylab("Number of feature differentially expressed")
  plot(g)

  #counting features detection limits for t-test
  c_n24 <- hist(OV8_limit_ttest_n24$V1, seq(-1,19,1))
  c_n16 <- hist(OV8_limit_ttest_n16$V1, seq(-1,19,1))
  c_n12 <- hist(OV8_limit_ttest_n12$V1, seq(-1,19,1))
  c_n10 <- hist(OV8_limit_ttest_n10$V1, seq(-1,19,1))
  c_n5 <- hist(OV8_limit_ttest_n5$V1, seq(-1,19,1))

  bins <- c(1:10, seq(20,110,10))
  d <- data.frame(counts=cumsum(c_n24$counts), bin=bins, rep="24 replicates")
  #d <- rbind(d, data.frame(counts=cumsum(c_n18$counts), bin=bins, rep="18 replicates"))
  d <- rbind(d, data.frame(counts=cumsum(c_n16$counts), bin=bins, rep="16 replicates"))
  d <- rbind(d, data.frame(counts=cumsum(c_n12$counts), bin=bins, rep="12 replicates"))
  d <- rbind(d, data.frame(counts=cumsum(c_n10$counts), bin=bins, rep="10 replicates"))
  d <- rbind(d, data.frame(counts=cumsum(c_n5$counts), bin=bins, rep="5 replicates"))

  g <- ggplot(data = d, aes(x=bin, y=counts, group=rep, color=rep)) + geom_line() + geom_point()
  g <- g + scale_x_continuous(breaks =seq(1,10,1),limit=c(1,10)) + scale_y_continuous(limit=c(0,50))
  g <- g + xlab("Minimum % of OVCAR8 detected") + ylab("Number of feature differentially expressed")
  plot(g)

  #compute profile similarity
  OV8_inds <- which(wt[OV8==100, p_value < 1e-5])

  N <- nrow(feature_mat)
  sims <- rep(0,1e6)
  cnt <- 1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(feature_mat$ATTRIBUTE_OV8[[i]] == 0 & feature_mat$ATTRIBUTE_OV8[[j]] == 5){
        #browser()
        p1 <- feature_mat[i, OV8_inds+1, with=F]
        p2 <- feature_mat[j, OV8_inds+1, with=F]
        p1 <- p1/sum(p1^2)^0.5
        p2 <- p2/sum(p2^2)^0.5
        sims[cnt] <- sum(p1*p2)
        cnt <- cnt + 1
      }
    }
  }

  sims2 <- rep(0,1e6)
  cnt <- 1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(feature_mat$ATTRIBUTE_OV8[[i]] == 0 & feature_mat$ATTRIBUTE_OV8[[j]] == 0){
        #browser()
        p1 <- feature_mat[i, OV8_inds+1, with=F]
        p2 <- feature_mat[j, OV8_inds+1, with=F]
        p1 <- p1/sum(p1^2)^0.5
        p2 <- p2/sum(p2^2)^0.5
        sims2[cnt] <- sum(p1*p2)
        cnt <- cnt + 1
      }
    }
  }

  df <- data.frame(sim=sims2[sims2 > 0], OV8=0)
  df <- rbind(df, data.frame(sim=sims0_1, OV8="1% OV8"))
  df <- rbind(df, data.frame(sim=sims0_5, OV8="5% OV8"))
  df <- rbind(df, data.frame(sim=sims0_10, OV8="10% OV8"))
  df <- rbind(df, data.frame(sim=sims0_20, OV8="20% OV8"))

  g <- ggplot(df) + geom_density(aes(x=sim, group=OV8, color=OV8))
}


import_all_spec <- function(file_list){
  spect_list <- unlist(lapply(file_list, function(d){import(d)}))
}


analyze_TIC <- function(spect_list, meta_data, title=""){
  TIC <- unlist(lapply(spect_list, function(s){sum(intensity(s))}))
  filenames <- basename(unlist(lapply(spect_list, function(p){attributes(p)$metaData$file})))
  dt <- data.table(TIC=TIC, filename=filenames)
  dt_merge <- merge(dt, meta_data, by='filename')
  g <- ggplot(dt_merge)
  g <- g + geom_boxplot(aes(x=as.factor(dt_merge$ATTRIBUTE_DAY), y=dt_merge$TIC)) + ggtitle(title)
  plot(g)

  g <- ggplot(dt_merge)
  g <- g + geom_boxplot(aes(x=as.factor(dt_merge$ATTRIBUTE_Mouse_Number), y=dt_merge$TIC)) + ggtitle(title)
  plot(g)

  g <- ggplot(dt_merge)
  g <- g + geom_boxplot(aes(x=as.factor(paste0(dt_merge$ATTRIBUTE_Mouse_Number, "_", dt_merge$ATTRIBUTE_DAY)), y=dt_merge$TIC)) + ggtitle(title)
  plot(g)
}

analyze_mouse_study <- function(){
  #reading mzml
  data_path <- '~/DataNAS/LauraData/MouseStudy/Mice_901_902_valya/'
  data_path2 <- '~/DataNAS/LauraData/MouseStudy/Mouse Study_1_Mice 903_904_905_Day0_56_Melissa/'
  data_path3 <- '~/DataNAS/LauraData/MouseStudy/Pooled_Lavages_for_Tap_the_Pap_Mouse_Study/'
  data_path_study2 <- '~/DataNAS/LauraData/MouseStudy/Mouse Study 2 - Mice 77-79/'
  data_list <- list(data_path, data_path2, data_path3)

  meta_data <- read_meta_data('~/DataNAS/LauraData/MouseStudy/Meta_data_merged.txt', F)
  ms_spect_list1 <- import_all_spec(data_list)
  ms_spect_list2 <- import_all_spec(list(data_path_study2))


  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/TIC_analysis_with_calibrated_intensity.pdf')
  #checking TIC stats for raw spectra
  analyze_TIC(ms_spect_list, meta_data, "TIC Raw spectra")

  #pre-processing MALDI spectra
  plot(ms_spect_list[[1]])
  spectra <- transformIntensity(ms_spect_list, method="sqrt")
  plot(spectra[[1]])

  TIC <- unlist(lapply(ms_spect_list, function(s){sum(intensity(s))}))

  spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
  plot(spectra[[1]])
  analyze_TIC(spectra, meta_data, "TIC signal-processed spectra")

  spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
  plot(spectra[[1]])
  analyze_TIC(spectra, meta_data, "TIC baseline-removed spectra")


  spectra <- calibrateIntensity(spectra, method="TIC")
  plot(spectra[[1]])
  spectra <- alignSpectra(spectra[1:1000], halfWindowSize=20, SNR=2, tolerance=0.002, warpingMethod="lowess")
  plot(spectra[[1]])

  #preprocessing spectra for study1 and study2 separately; code above this mainly when processing study 1 only
  spectra1<- transformIntensity(ms_spect_list, method='sqrt')
  spectra1 <- removeBaseline(spectra1, method="SNIP", iteration=100)
  spectra2 <- transformIntensity(ms_spect_list2, method='sqrt')
  spectra2 <- removeBaseline(spectra2, method="SNIP", iteration=100)


  #check TIC difference to see if there are any anormality
  TIC <- unlist(lapply(spectra1, function(s){sum(intensity(s))}))
  TIC2 <- unlist(lapply(spectra2, function(s){sum(intensity(s))}))
  hist(log10(TIC),50)
  hist(log10(TIC2), 50)

  spectra1 <- calibrateIntensity(spectra, method="TIC")
  spectra2 <- calibrateIntensity(spectra2, method="TIC")


  #we try investigate the reference peaks manually since it fail to align automatically (complaining only 3 peaks in reference peaks)
  peaks <- detectPeaks(spectra1, halfWindowSize=20, SNR=3)

  freqs <- seq(0.7, 0.9, 0.05)

  ref_peak_cnt <- unlist(lapply(freqs,
                          function(freq){
                            reference_peaks <- referencePeaks(peaks, minFrequency = freq, tolerance = 0.002)
                            length(reference_peaks)
                          }
                         ))

  peaks2 <- detectPeaks(spectra2, halfWindowSize=20, SNR=3)
  ref_peak_cnt2 <- unlist(lapply(freqs,
                         function(freq){
                           reference_peaks2 <- referencePeaks(peaks2, minFrequency = freq, tolerance = 0.002)
                           length(reference_peaks2)
                         })
  )

  peaks_combined <- detectPeaks(c(spectra1, spectra2),halfWindowSize=20, SNR=4)
  reference_combine <- referencePeaks(peaks_combined, minFrequency = 0.75, tolerance = 0.002)
  length(reference_combine)
  #using a lower minFrequency of 0.8 seems to return reasonable number of reference peaks for both set together: 17 peaks found

  #try to visualize some spectra to see why they fail to align or if replicate spectra align well
  filenames1 <- unlist(lapply(spectra1, function(s){basename(attributes(s)$metaData$file)}))
  filenames2 <- unlist(lapply(spectra2, function(s){basename(attributes(s)$metaData$file)}))


  inds <- which(grepl('Study 2', lapply(spectra, function(s){attributes(s)$metaData$file})))
  i <- 2
  j <- 3
  pair_spect <- rbind(data.table(x=mass(peaks[[i]]), y=intensity(peaks[[i]]), label=paste0('spectrum ', i)),
                      data.table(x=mass(peaks[[j]]), y=intensity(peaks[[j]]), label = paste0('spectrum', j)))
  compare_spectrum_interactive(pair_spect$x, pair_spect$y, pair_spect$label)

  pair_spect2 <- rbind(data.table(x=mass(peaks2[[i]]), y=intensity(peaks2[[i]]), label=paste0('spectrum ', i)),
                      data.table(x=mass(peaks2[[j]]), y=intensity(peaks2[[j]]), label = paste0('spectrum', j)))

  compare_spectrum_interactive(pair_spect2$x, pair_spect2$y, pair_spect2$label)


  #aligning spectra from study 1 and 2 together
  #we are supplying our custom reference peaks, from above
  spectra <- alignSpectra(c(spectra1), halfWindowSize=20, SNR=4, tolerance=0.002, warpingMethod="lowess", reference = reference_combine)
  spectra <- alignSpectra(c(spectra1, spectra2), halfWindowSize=20, SNR=4, tolerance=0.002, warpingMethod="lowess", reference = reference_combine)

  #####
  noise <- estimateNoise(spectra[[1]])
  plot(spectra[[1]], xlim=c(5000, 20000))
  lines(noise, col="red")
  lines(noise[,1], noise[, 2]*2, col="blue")

  peaks <- detectPeaks(spectra, method="MAD", halfWindowSize=20, SNR=4)
  plot(spectra[[1]])
  points(peaks[[1]], col="red", pch=4)

  peaks <- binPeaks(peaks, tolerance=0.002)
  analyze_TIC(peaks, meta_data, "TIC peak-picked spectra")
  dev.off()

  filenames <- basename(unlist(lapply(peaks, function(p){attributes(p)$metaData$file})))
  featureMatrix <- intensityMatrix(peaks, spectra)

  #do peak intensity normalization
  #featurMatrix <- featureMatrix^0.5  #we do not need to do this here since this is already done during spectra preprocessing
  s <- rowSums(featureMatrix)
  Ncol <- ncol(featureMatrix)
  a <- matrix(rep(s, Ncol), ncol = Ncol)
  featureMatrix <- featureMatrix/a


  #adding filename to feature matrix
  attributes(featureMatrix)$dimnames[[1]] <- basename(filenames)

  df <- as.data.frame(featureMatrix)
  feature_names <- paste0("FEATURE_", seq(1,nrow(df)), "_", colnames(df), "_1.0")
  colnames(df) <- feature_names

  dt <- data.table(filename=filenames)
  dt <- cbind(dt, df)
  feature_mat <- merge(dt, meta_data, by='filename')
  attr(feature_mat, 'dataColInd') <- max(which(grepl('FEATURE', colnames(feature_mat))))
  dat_ind <- max(which(grepl('FEATURE_', colnames(feature_mat))))

  #performing PCA
  df <- as.data.frame(featureMatrix)

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/PCA_stats_withoaut_sqrt_transform.pdf')
  ms.pca <- prcomp(featureMatrix,
                   center = TRUE,
                   scale. = TRUE)
  plot(ms.pca)
  feature_mat$ATTRIBUTE_SampleGroupID <- as.character(feature_mat$ATTRIBUTE_SampleGroupID)
  feature_mat$ATTRIBUTE_DAY <- as.character(feature_mat$ATTRIBUTE_DAY)
  feature_mat$AvgRadiantEff <- as.character(feature_mat$AvgRadiantEff)
  feature_mat$ATTRIBUTE_Mouse_Number <- as.character(feature_mat$ATTRIBUTE_Mouse_Number)
  feature_mat$ATTRIBUTE_AvgRadiantEff <- as.character(feature_mat$ATTRIBUTE_AvgRadiantEff)
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_DAY')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_SampleGroupID')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_Mouse_Number')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_AvgRadiantEff')

  #looking at each mouse
  mouse_id <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  colfunc <- colorRampPalette(c("blue", "red", "green"))
  col <- colfunc(max(as.numeric(feature_mat$ATTRIBUTE_DAY))+1)
  cols <- col[(as.numeric(feature_mat$ATTRIBUTE_DAY)+1)]
  plot(rep(1,length(col)), col=col)
  xlim <- c(min(ms.pca$x[,1]), max(ms.pca$x[,1]))
  ylim <- c(min(ms.pca$x[,2]), max(ms.pca$x[,2]))
  plots <- lapply(mouse_id,
          function(id){
            inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == id)
            g <- ggplot() + geom_point(aes(x=ms.pca$x[inds,1], y=ms.pca$x[inds, 2], color=as.numeric(feature_mat$ATTRIBUTE_DAY[inds])))
            g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
            g <- g + ggtitle(paste0("Mouse Id ", id)) + scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim)
            plot(g)
          })
  dev.off()

  #compute reproducibility
  cv <- compute_feature_abudnance_cv_by_group(feature_mat, groups = as.numeric(feature_mat$ATTRIBUTE_SampleGroupID), min.int = 0.00000000001)
  cv$logInt <- as.integer(log10(cv$Mean.V1+0.00000000001))
  hist(cv$presence)
  dat <- as.matrix(feature_mat[,2:dat_ind])

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/CV_stats_without_sqrt_transform_Study_1_and_2_SNR4.pdf')
  #reproducibility per group
  g <- ggplot(data=cv)
  g <- g + stat_boxplot(aes(x=as.factor(Group), y=CV)) + scale_y_continuous(limit=c(0,2))
  g

  g <- ggplot(data=cv[logInt > -6])
  g <- g + stat_boxplot(aes(x=as.factor(logInt), y=CV)) + scale_y_continuous(limit=c(0,2))
  g

  dev.off()

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/Similarity_stats_without_sqrt_transform_Study1_and_2_SNR4.pdf')

  tech_rep_sim <- computeGroupSim(feature_mat, groups = feature_mat$ATTRIBUTE_SampleGroupID)
  tech_rep_sim <- merge(tech_rep_sim, meta_data, by='filename')
  g <- ggplot(tech_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_DAY), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("DAY") + ylab("Avg sim across technical replicates")
  plot(g)

  g <- ggplot(tech_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_Mouse_Number), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Mouse Number") + ylab("Avg sim across technical replicates")
  plot(g)

  bio_rep_sim <- computeGroupSim(feature_mat, groups = feature_mat$ATTRIBUTE_DAY)
  bio_rep_sim <- merge(bio_rep_sim, meta_data, by='filename')
  g <- ggplot(bio_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_DAY), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("DAY") + ylab("Avg sim across samples from same time point")
  plot(g)

  bio_rep_sim2 <- computeGroupSim(feature_mat, groups = feature_mat$ATTRIBUTE_Mouse_Number)
  bio_rep_sim2 <- merge(bio_rep_sim2, meta_data, by='filename')

  g <- ggplot(bio_rep_sim2) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_Mouse_Number), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Mouse Number") + ylab("Avg sim across samples from same mouse")
  plot(g)

  g <- ggplot(bio_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_Mouse_Number), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Mouse Number") + ylab("Avg sim across samples from same time point")
  plot(g)

  g <- ggplot(bio_rep_sim2) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_DAY), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Day") + ylab("Avg sim across samples from same mouse")
  plot(g)

  dev.off()

  dat <- t(as.matrix(feature_mat[,2:dat_ind, with=F]))
  sims <- cosine(dat)

  days <- unique(feature_mat$ATTRIBUTE_DAY)
  out_sims <- rbindlist(lapply(1:(length(days)-1),
                function(i){
                  d <- days[[i]]
                  inds <- which(feature_mat$ATTRIBUTE_DAY == d)
                  rbindlist(lapply(1:length(days),
                      function(j){
                          d2 <- days[j]
                           inds2 <- which(feature_mat$ATTRIBUTE_DAY == d2)
                           pairInds <- expand.grid(inds, inds2)
                           pairInds <- (pairInds[,1])*nrow(sims)+pairInds[,2]
                           data.table(ind1 = d, ind2 = d2, sim=sims[pairInds])
                      }))
                }
      ))
  avg_out_sims <- out_sims[,mean(sim), by=.(ind1, ind2)]
  g <- ggplot(avg_out_sims) + geom_point(aes(x=ind1, y=ind2, color=V1, size = 2))
  g <-  g + geom_point(aes(x=ind2, y=ind1, color=V1, size = 2))
  g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
  pdf('/home/jian/DataNAS/LauraData/MouseStudy/Analysis/Avg_matrix_in_time_point_without_sqrt_transform_study_1_and_2_SNR4.pdf')
  plot(g)
  dev.off()

  #doing MDS on cosine similarity
  mds <- cmdscale(1-sims, k=2)

  #trying other distance measure
  d <- dist(t(dat), method = 'man')
  mds <- cmdscale(d, k=3)

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/MDS_plot_without_sqrt_transform_study_1_and_2_SNR4.pdf')
  mouse_id <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  colfunc <- colorRampPalette(c("blue", "red", "green"))
  col <- colfunc(max(as.numeric(feature_mat$ATTRIBUTE_DAY)+1))
  cols <- col[as.numeric(feature_mat$ATTRIBUTE_DAY)+1]
  plot(rep(1,length(col)), col=col)
  xlim <- c(min(mds[,1]), max(mds[,1]))
  ylim <- c(min(mds[,2]), max(mds[,2]))
  plots <- lapply(mouse_id,
                  function(id){
                    inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == id)
                    g <- ggplot() + geom_point(aes(x=mds[inds,1], y=mds[inds, 2], color=as.numeric(feature_mat$ATTRIBUTE_DAY[inds])))
                    g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
                    g <- g + ggtitle(paste0("Mouse Id ", id)) + scale_x_continuous(limits = xlim) + scale_y_continuous(limits=ylim)
                    plot(g)
                  })
  dev.off()

  #trying out t-SNE embedding
  library(Rtsne)

  rtsne_embedd <- Rtsne(t(dat))
  rtsne_embedd <- Rtsne(1-sims, is_distance = T)

  xlim <- c(min(rtsne_embedd$Y[,1]), max(rtsne_embedd$Y[,1]))
  ylim <- c(min(rtsne_embedd$Y[,2]), max(rtsne_embedd$Y[,2]))
  plots <- lapply(mouse_id,
                  function(id){
                    inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == id)
                    g <- ggplot() + geom_point(aes(x=rtsne_embedd$Y[inds,1], y=rtsne_embedd$Y[inds, 2], color=as.numeric(feature_mat$ATTRIBUTE_DAY[inds])))
                    g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
                    g <- g + ggtitle(paste0("Mouse Id ", id)) + scale_x_continuous(limits = xlim) + scale_y_continuous(limits=ylim)
                    plot(g)
                  })
  #distribution of intensity
  hist(log10(featureMatrix))


  #checking features that correlate with tumor burden
  cor_list <- unlist(lapply(2:dat_ind,
                     function(i){
                        corr <- cor(feature_mat[,i,with=F], as.numeric(feature_mat$ATTRIBUTE_AvgRadiantEff, method='kendall'))
                     }))

  hist(cor_list)
  which(abs(cor_list) > 0.6)
  plot(feature_mat[[1459]], as.numeric(feature_mat$ATTRIBUTE_AvgRadiantEff))

  #it seems last tumor burden separate from no turmor burden searching for features
  p_list <- unlist(lapply(2:dat_ind,
                            function(i){
                              a <- as.numeric(feature_mat$ATTRIBUTE_AvgRadiantEff)
                              p <- wilcox.test(feature_mat[a ==0, i, with=F][[1]],
                                               feature_mat[a == 2.7e9, i, with=F][[1]])
                              p$p.value
                            }))

  #compute statistically significant features across time-points
  base_line_inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == -2)
  days <- unique(feature_mat$ATTRIBUTE_DAY)
  days <- sort(days[days >1])
  mouses <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  mouses <- mouses[mouses > 10]
  ft1 <- feature_mat[base_line_inds,2:dat_ind]
  diff_ft_allmouses <- rbindlist(lapply(days,
                        function(d){
                          rbindlist(lapply(days,
                            function(d2){
                             if(d <d2){
                               rbindlist(
                                  lapply(mouses,
                                        function(m){
                                           base_inds <- which(feature_mat$ATTRIBUTE_DAY == d & feature_mat$ATTRIBUTE_Mouse_Number == m)
                                           inds <- which(feature_mat$ATTRIBUTE_DAY == d2 & feature_mat$ATTRIBUTE_Mouse_Number == m)
                                              ft1 <- feature_mat[base_inds,2:dat_ind]
                                              ft2 <- feature_mat[inds,2:dat_ind]
                                              if(nrow(ft2) > 5 & nrow(ft1) > 5){
                                                    p <- sapply(1:ncol(ft1),
                                                    function(i){
                                                      wilcox.test(ft1[,i, with=F][[1]], ft2[,i, with=F][[1]])$p.value
                                                    })
                                                    list(day=rep(d, length(p)), day2= rep(d2, length(p)), mouse=rep(m, length(p)), pval=p, ft_id=1:length(p))
                                              }
                                         }
                               ))
                              }
                            }
                          ))
                        }))

  mouses <- mouses[mouses > 100]
  diff_ft <- rbindlist(lapply(days,
                              function(d){
                                rbindlist(lapply(days,
                                                 function(d2){
                                                   if(d <d2){
                                                     rbindlist(
                                                       lapply(mouses,
                                                              function(m){
                                                                base_inds <- which(feature_mat$ATTRIBUTE_DAY == d & feature_mat$ATTRIBUTE_Mouse_Number == m)
                                                                inds <- which(feature_mat$ATTRIBUTE_DAY == d2 & feature_mat$ATTRIBUTE_Mouse_Number == m)
                                                                ft1 <- feature_mat[base_inds,2:dat_ind]
                                                                ft2 <- feature_mat[inds,2:dat_ind]
                                                                if(nrow(ft2) > 5 & nrow(ft1) > 5){
                                                                  p <- sapply(1:ncol(ft1),
                                                                              function(i){
                                                                                wilcox.test(ft1[,i, with=F][[1]], ft2[,i, with=F][[1]])$p.value
                                                                              })
                                                                  dir <- sapply(1:ncol(ft1),
                                                                              function(i){
                                                                                ifelse(mean(ft1[,i, with=F][[1]]) > mean(ft2[,i,with=F][[1]]), 1, -1)
                                                                              })

                                                                  list(day=rep(d, length(p)), day2= rep(d2, length(p)), mouse=rep(m, length(p)), pval=p, dir=dir, ft_id=1:length(p))
                                                                }
                                                              }
                                                       ))
                                                   }
                                                 }
                                ))
                              }))

  ft1_mean <- colSums(ft1)/12
  diff_ft$base_mean <- ft1_mean

  cnt_mouses <- diff_ft[pval*3175 < 0.01, list(cnt=length(unique(mouse))),by=c('ft_id', 'day','day2', 'dir')]
  cnt_days <- diff_ft[pval*3175 < 0.01, list(cnt=length(unique(day2))),by=c('ft_id', 'mouse', 'day', 'dir')]

  cnt_mouses[cnt==5 & dir==1, table(day,day2)]
  cnt_mouses[cnt==5 & dir==-1, table(day,day2)]

  head(cnt_mouses[cnt==5 & day==7 & day2==49])

  #plotting features
  days <- unique(feature_mat$ATTRIBUTE_DAY)
  mouses <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  days_label <- sort(as.integer(days))
  mouses_label <- sort(as.numeric(mouses))
  label <- NULL
  for(m in mouses_label){
    for(d in days_label){
      label <- c(label, paste0(m, "_", d))
    }
  }

  dat_labels <- paste0(feature_mat$ATTRIBUTE_Mouse_Number, "_", feature_mat$ATTRIBUTE_DAY)
  dat_labels <- factor(dat_labels, levels=label)

  feature_mat$ATTRIBUTE_LABEL <- dat_labels

  g <- ggplot() + geom_boxplot(aes(x=feature_mat$ATTRIBUTE_LABEL, y=feature_mat[,46][[1]]))
  plot(g)

  ind <- 96
  ind <- ind + 1
  day <- 7
  day2 <- 49
  mouse <- 910

  d <- feature_mat[(ATTRIBUTE_DAY==day | ATTRIBUTE_DAY==day | ATTRIBUTE_DAY ==day2) & ATTRIBUTE_Mouse_Number > 0 & ATTRIBUTE_Mouse_Number > 100]
  g <-ggplot() + geom_boxplot(aes(x=as.factor(d$ATTRIBUTE_LABEL), y=d[,ind, with=F][[1]]))
  plot(g)

  d <- feature_mat[ATTRIBUTE_Mouse_Number == 901]
  g <-ggplot() + geom_boxplot(aes(x=as.factor(d$ATTRIBUTE_LABEL), y=d[,ind, with=F][[1]]))
  plot(g)



  plot_list<- lapply(mouses,
                     function(m){
                       d <- feature_mat[ATTRIBUTE_Mouse_Number==m,]
                       g <-ggplot() + geom_boxplot(aes(x=as.factor(d$ATTRIBUTE_DAY), y=d[,ind, with=F][[1]]))
                       return(g)
                     }
               )
  marrangeGrob(plot_list, ncol=1, nrow=length(plot_list))

}

#the following analysis is only for study 1
analyze_mouse_study1 <- function(){
  #reading mzml
  data_path <- '~/DataNAS/LauraData/MouseStudy/Mice_901_902_valya/'
  data_path2 <- '~/DataNAS/LauraData/MouseStudy/Mouse Study_1_Mice 903_904_905_Day0_56_Melissa/'
  data_path3 <- '~/DataNAS/LauraData/MouseStudy/Pooled_Lavages_for_Tap_the_Pap_Mouse_Study/'
  data_list <- list(data_path, data_path2, data_path3)

  meta_data <- read_meta_data('~/DataNAS/LauraData/MouseStudy/Meta_data_merged.txt', F)
  ms_spect_list1 <- import_all_spec(data_list)

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/TIC_analysis_with_calibrated_intensity.pdf')
  #checking TIC stats for raw spectra
  analyze_TIC(ms_spect_list, meta_data, "TIC Raw spectra")

  #pre-processing MALDI spectra

  #####
  noise <- estimateNoise(spectra[[1]])
  plot(spectra[[1]], xlim=c(5000, 20000))
  lines(noise, col="red")
  lines(noise[,1], noise[, 2]*2, col="blue")

  #lines(noise[,1], noise[, 2]*4, col="blue")

  peaks <- detectPeaks(spectra, method="MAD", halfWindowSize=20, SNR=4)
  plot(spectra[[1]])
  points(peaks[[1]], col="red", pch=4)

  peaks <- binPeaks(peaks, tolerance=0.002)
  analyze_TIC(peaks, meta_data, "TIC peak-picked spectra")
  dev.off()

  filenames <- basename(unlist(lapply(peaks, function(p){attributes(p)$metaData$file})))
  featureMatrix <- intensityMatrix(peaks, spectra)

  #do peak intensity normalization
  #featurMatrix <- featureMatrix^0.5  #we do not need to do this here since this is already done during spectra preprocessing
  #Do we need to redo normalization here? Since we already does this in TIC normalization
  s <- rowSums(featureMatrix)
  Ncol <- ncol(featureMatrix)
  a <- matrix(rep(s, Ncol), ncol = Ncol)
  featureMatrix <- featureMatrix/a


  #adding filename to feature matrix
  attributes(featureMatrix)$dimnames[[1]] <- basename(filenames)

  df <- as.data.frame(featureMatrix)
  feature_names <- paste0("FEATURE_", seq(1,ncol(df)), "_", colnames(df), "_1.0")
  colnames(df) <- feature_names

  dt <- data.table(filename=filenames)
  dt <- cbind(dt, df)
  feature_mat <- merge(dt, meta_data, by='filename')
  attr(feature_mat, 'dataColInd') <- max(which(grepl('FEATURE', colnames(feature_mat))))
  dat_ind <- max(which(grepl('FEATURE_', colnames(feature_mat))))

  #performing PCA
  df <- as.data.frame(featureMatrix)

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/PCA_stats_withoaut_sqrt_transform.pdf')
  ms.pca <- prcomp(featureMatrix,
                   center = TRUE,
                   scale. = TRUE)
  plot(ms.pca)
  feature_mat$ATTRIBUTE_SampleGroupID <- as.character(feature_mat$ATTRIBUTE_SampleGroupID)
  feature_mat$ATTRIBUTE_DAY <- as.character(feature_mat$ATTRIBUTE_DAY)
  feature_mat$AvgRadiantEff <- as.character(feature_mat$AvgRadiantEff)
  feature_mat$ATTRIBUTE_Mouse_Number <- as.character(feature_mat$ATTRIBUTE_Mouse_Number)
  feature_mat$ATTRIBUTE_AvgRadiantEff <- as.character(feature_mat$ATTRIBUTE_AvgRadiantEff)
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_DAY')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_SampleGroupID')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_Mouse_Number')
  autoplot(ms.pca, data = feature_mat, colour='ATTRIBUTE_AvgRadiantEff')

  #looking at each mouse
  mouse_id <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  colfunc <- colorRampPalette(c("blue", "red", "green"))
  col <- colfunc(max(as.numeric(feature_mat$ATTRIBUTE_DAY))+1)
  cols <- col[(as.numeric(feature_mat$ATTRIBUTE_DAY)+1)]
  plot(rep(1,length(col)), col=col)
  xlim <- c(min(ms.pca$x[,1]), max(ms.pca$x[,1]))
  ylim <- c(min(ms.pca$x[,2]), max(ms.pca$x[,2]))
  plots <- lapply(mouse_id,
                  function(id){
                    inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == id)
                    g <- ggplot() + geom_point(aes(x=ms.pca$x[inds,1], y=ms.pca$x[inds, 2], color=as.numeric(feature_mat$ATTRIBUTE_DAY[inds])))
                    g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
                    g <- g + ggtitle(paste0("Mouse Id ", id)) + scale_x_continuous(limits = xlim) + scale_y_continuous(limits = ylim)
                    plot(g)
                  })
  dev.off()

  #compute reproducibility
  cv <- compute_feature_abudnance_cv_by_group(feature_mat, groups = as.numeric(feature_mat$ATTRIBUTE_SampleGroupID), min.int = 0.00000000001)
  cv$logInt <- as.integer(log10(cv$Mean.V1+0.00000000001))
  hist(cv$presence)
  dat <- as.matrix(feature_mat[,2:dat_ind])

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/CV_stats_without_sqrt_transform_Study_1_and_2_SNR4.pdf')
  #reproducibility per group
  g <- ggplot(data=cv)
  g <- g + stat_boxplot(aes(x=as.factor(Group), y=CV)) + scale_y_continuous(limit=c(0,2))
  g

  g <- ggplot(data=cv[logInt > -6])
  g <- g + stat_boxplot(aes(x=as.factor(logInt), y=CV)) + scale_y_continuous(limit=c(0,2))
  g

  dev.off()

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/Similarity_stats_without_sqrt_transform_Study1_and_2_SNR4.pdf')

  tech_rep_sim <- computeGroupSim(feature_mat, groups = feature_mat$ATTRIBUTE_SampleGroupID)
  tech_rep_sim <- merge(tech_rep_sim, meta_data, by='filename')
  g <- ggplot(tech_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_DAY), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("DAY") + ylab("Avg sim across technical replicates")
  plot(g)

  g <- ggplot(tech_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_Mouse_Number), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Mouse Number") + ylab("Avg sim across technical replicates")
  plot(g)

  bio_rep_sim <- computeGroupSim(feature_mat, groups = feature_mat$ATTRIBUTE_DAY)
  bio_rep_sim <- merge(bio_rep_sim, meta_data, by='filename')
  g <- ggplot(bio_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_DAY), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("DAY") + ylab("Avg sim across samples from same time point")
  plot(g)

  bio_rep_sim2 <- computeGroupSim(feature_mat, groups = feature_mat$ATTRIBUTE_Mouse_Number)
  bio_rep_sim2 <- merge(bio_rep_sim2, meta_data, by='filename')

  g <- ggplot(bio_rep_sim2) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_Mouse_Number), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Mouse Number") + ylab("Avg sim across samples from same mouse")
  plot(g)

  g <- ggplot(bio_rep_sim) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_Mouse_Number), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Mouse Number") + ylab("Avg sim across samples from same time point")
  plot(g)

  g <- ggplot(bio_rep_sim2) + geom_boxplot(aes(x=as.factor(ATTRIBUTE_DAY), y=avgSim)) + scale_y_continuous(limit=c(0.6,1))
  g <- g + xlab("Day") + ylab("Avg sim across samples from same mouse")
  plot(g)

  dev.off()

  feature_mat <- feature_mat[ATTRIBUTE_DAY > 5]
  dat <- t(as.matrix(feature_mat[,2:dat_ind, with=F]))
  sims <- cosine(dat)

  days <- unique(feature_mat$ATTRIBUTE_DAY)
  out_sims <- rbindlist(lapply(1:(length(days)-1),
                               function(i){
                                 d <- days[[i]]
                                 inds <- which(feature_mat$ATTRIBUTE_DAY == d)
                                 rbindlist(lapply(1:length(days),
                                                  function(j){
                                                    d2 <- days[j]
                                                    inds2 <- which(feature_mat$ATTRIBUTE_DAY == d2)
                                                    pairInds <- expand.grid(inds, inds2)
                                                    pairInds <- (pairInds[,1])*nrow(sims)+pairInds[,2]
                                                    data.table(ind1 = d, ind2 = d2, sim=sims[pairInds])
                                                  }))
                               }
  ))
  avg_out_sims <- out_sims[,mean(sim, na.rm=T), by=.(ind1, ind2)]
  avg_out_sims$week1 <- avg_out_sims$ind1/7
  avg_out_sims$week2 <- avg_out_sims$ind2/7
  avg_out_sims$cosine_sim <- avg_out_sims$V1
  g <- ggplot(avg_out_sims)
  g <-  g + geom_point(aes(x=week1, y=week2, color=cosine_sim, size=2))
  g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
  ggsave('/home/jian/DataNAS/LauraData/MouseStudy/Analysis/Avg_matrix_in_time_point_without_sqrt_transform_study_1_and_2_SNR4.svg', g)
  #pdf('/home/jian/DataNAS/LauraData/MouseStudy/Analysis/Avg_matrix_in_time_point_without_sqrt_transform_study_1_and_2_SNR4.pdf')
  #plot(g)
  #dev.off()

  #doing MDS on cosine similarity
  mds <- cmdscale(1-sims, k=2)

  #trying other distance measure
  d <- dist(t(dat), method = 'man')
  mds <- cmdscale(d, k=3)

  pdf('~/DataNAS/LauraData/MouseStudy/Analysis/MDS_plot_without_sqrt_transform_study_1_and_2_SNR2_day_7_after.pdf')
  mouse_id <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  colfunc <- colorRampPalette(c("blue", "red", "green"))
  col <- colfunc(max(as.numeric(feature_mat$ATTRIBUTE_DAY)+1))
  cols <- col[as.numeric(feature_mat$ATTRIBUTE_DAY)+1]
  plot(rep(1,length(col)), col=col)
  xlim <- c(min(mds[,1]), max(mds[,1]))
  ylim <- c(min(mds[,2]), max(mds[,2]))
  plots <- lapply(mouse_id,
                  function(id){
                    inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == id)
                    g <- ggplot() + geom_point(aes(x=mds[inds,1], y=mds[inds, 2], color=as.numeric(feature_mat$ATTRIBUTE_DAY[inds])))
                    g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
                    g <- g + ggtitle(paste0("Mouse Id ", id)) + scale_x_continuous(limits = xlim) + scale_y_continuous(limits=ylim)
                    plot(g)
                  })
  dev.off()

  #trying out t-SNE embedding
  library(Rtsne)

  rtsne_embedd <- Rtsne(t(dat))
  rtsne_embedd <- Rtsne(1-sims, is_distance = T)

  xlim <- c(min(rtsne_embedd$Y[,1]), max(rtsne_embedd$Y[,1]))
  ylim <- c(min(rtsne_embedd$Y[,2]), max(rtsne_embedd$Y[,2]))
  plots <- lapply(mouse_id[mouse_id > 100],
                  function(id){
                    inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == id)
                    g <- ggplot() + geom_point(aes(x=rtsne_embedd$Y[inds,1], y=rtsne_embedd$Y[inds, 2], color=as.numeric(feature_mat$ATTRIBUTE_DAY[inds])))
                    g <- g  + scale_color_gradientn(colours = c("red", "blue", "green"))
                    g <- g + ggtitle(paste0("Mouse Id ", id)) + scale_x_continuous(limits = xlim) + scale_y_continuous(limits=ylim)
                    plot(g)
                  })
  #distribution of intensity
  hist(log10(featureMatrix))


  #checking features that correlate with tumor burden
  cor_list <- unlist(lapply(2:dat_ind,
                            function(i){
                              corr <- cor(feature_mat[,i,with=F], as.numeric(feature_mat$ATTRIBUTE_AvgRadiantEff, method='kendall'))
                            }))

  hist(cor_list)
  which(abs(cor_list) > 0.6)
  plot(feature_mat[[1459]], as.numeric(feature_mat$ATTRIBUTE_AvgRadiantEff))

  #it seems last tumor burden separate from no turmor burden searching for features
  p_list <- unlist(lapply(2:dat_ind,
                          function(i){
                            a <- as.numeric(feature_mat$ATTRIBUTE_AvgRadiantEff)
                            p <- wilcox.test(feature_mat[a ==0, i, with=F][[1]],
                                             feature_mat[a == 2.7e9, i, with=F][[1]])
                            p$p.value
                          }))

  #compute statistically significant features across time-points
  base_line_inds <- which(feature_mat$ATTRIBUTE_Mouse_Number == -2)
  days <- unique(feature_mat$ATTRIBUTE_DAY)
  days <- sort(days[days >1])
  mouses <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  mouses <- mouses[mouses > 10]
  mouses <- mouses[mouses > 100]
  diff_ft <- rbindlist(lapply(days,
                              function(d){
                                rbindlist(lapply(days,
                                                 function(d2){
                                                   if(d <d2){
                                                     rbindlist(
                                                       lapply(mouses,
                                                              function(m){
                                                                base_inds <- which(feature_mat$ATTRIBUTE_DAY == d & feature_mat$ATTRIBUTE_Mouse_Number == m)
                                                                inds <- which(feature_mat$ATTRIBUTE_DAY == d2 & feature_mat$ATTRIBUTE_Mouse_Number == m)
                                                                ft1 <- feature_mat[base_inds,2:dat_ind]
                                                                ft2 <- feature_mat[inds,2:dat_ind]
                                                                if(nrow(ft2) > 5 & nrow(ft1) > 5){
                                                                  p <- sapply(1:ncol(ft1),
                                                                              function(i){
                                                                                wilcox.test(ft1[,i, with=F][[1]], ft2[,i, with=F][[1]])$p.value
                                                                              })
                                                                  dir <- sapply(1:ncol(ft1),
                                                                                function(i){
                                                                                  ifelse(mean(ft1[,i, with=F][[1]]) > mean(ft2[,i,with=F][[1]]), 1, -1)
                                                                                })

                                                                  list(day=rep(d, length(p)), day2= rep(d2, length(p)), mouse=rep(m, length(p)), pval=p, dir=dir, ft_id=1:length(p))
                                                                }
                                                              }
                                                       ))
                                                   }
                                                 }
                                ))
                              }))


  feature_names <- colnames(feature_mat)[2:dat_ind]
  feature_mass <- strsplit(feature_names, split = "_")
  feature_mass <- sapply(feature_mass, function(m){m[[3]]})
  diff_ft$featue_mass <- feature_mass[diff_ft$ft_id]


  cnt_mouses <- diff_ft[pval*3175 < 0.01, list(cnt=length(unique(mouse))),by=c('ft_id', 'day','day2', 'dir')]
  cnt_days <- diff_ft[pval*3175 < 0.01, list(cnt=length(unique(day2))),by=c('ft_id', 'mouse', 'day', 'dir')]


  cnt_mouses$ft_mass <- feature_mass[cnt_mouses$ft_id]

  cnt_mouses[cnt==5 & dir==1, table(day,day2)]
  cnt_mouses[cnt==5 & dir==-1, table(day,day2)]




  saveRDS(diff_ft, '~/DataNAS/LauraData/MouseStudy/Analysis/diff_features_withtout_total_peak_intensity_normalization.RDS')
  write.table(diff_ft, '~/DataNAS/LauraData/MouseStudy/Analysis/diff_features_SNR3.txt', sep='\t', row.names = F, quote = F)
  write.table(cnt_mouses, '~/DataNAS/LauraData/MouseStudy/Analysis/diff_features_mousesCnt.txt', sep='\t', row.names = F, quote = F)

  head(cnt_mouses[cnt==5 & day==7 & day2==49])

  #plotting features
  days <- unique(feature_mat$ATTRIBUTE_DAY)
  mouses <- unique(feature_mat$ATTRIBUTE_Mouse_Number)
  days_label <- sort(as.integer(days))
  mouses_label <- sort(as.numeric(mouses))
  label <- NULL
  for(m in mouses_label){
    for(d in days_label){
      label <- c(label, paste0(m, "_", d))
    }
  }

  dat_labels <- paste0(feature_mat$ATTRIBUTE_Mouse_Number, "_", feature_mat$ATTRIBUTE_DAY)
  dat_labels <- factor(dat_labels, levels=label)

  feature_mat$ATTRIBUTE_LABEL <- dat_labels

  g <- ggplot() + geom_boxplot(aes(x=feature_mat$ATTRIBUTE_LABEL, y=feature_mat[,46][[1]]))
  plot(g)

  ind <- 223
  ind <- ind + 1
  day <- 7
  day2 <- 49
  mouse <- 910

  d <- feature_mat[(ATTRIBUTE_DAY==day | ATTRIBUTE_DAY ==day2) & ATTRIBUTE_Mouse_Number > 0 & ATTRIBUTE_Mouse_Number > 100]
  g <-ggplot() + geom_boxplot(aes(x=as.factor(d$ATTRIBUTE_LABEL), y=d[,ind, with=F][[1]]))
  plot(g)

  d <- feature_mat[ATTRIBUTE_Mouse_Number == 901]
  g <-ggplot() + geom_boxplot(aes(x=as.factor(d$ATTRIBUTE_LABEL), y=d[,ind, with=F][[1]]))
  plot(g)



  plot_list<- lapply(mouses,
                     function(m){
                       d <- feature_mat[ATTRIBUTE_Mouse_Number==m,]
                       g <-ggplot() + geom_boxplot(aes(x=as.factor(d$ATTRIBUTE_DAY), y=d[,ind, with=F][[1]]))
                       return(g)
                     }
  )
  marrangeGrob(plot_list, ncol=1, nrow=length(plot_list))

}

get_sample_name <- function(){
  df$sample_name <- as.double(unlist(lapply(peaks,
                                            function(s){
                                              n <- metaData(s)$name
                                              s <- strsplit(n, 'MTWT')[[1]]
                                              print(s)
                                              if(length(s) == 1){
                                                return("0")
                                              }
                                              s <- strsplit(s[[2]], '_')[[1]]
                                              print(s)
                                              if(s[[1]]==""){
                                                return("100")
                                              }
                                              return(s[[1]])
                                            })))
}

feature_selection <- function(){
  #feature selection
  pc1_loading <- ms.pca$rotation[,1]
  inds <- order(abs(pc1_loading), decreasing = T)

  for(i in head(inds,20)){
    plot(df$sample_name, df[[i]], main=paste0("Peak mass: ", names(df)[[i]]))
  }

  dev.off()

  #computing linear correlation for each mass
  corrs <- unlist(lapply(inds, function(i){cor(df$sample_name, df[[i]])}))

  corrs <- unlist(lapply(1:ncol(df), function(i){cor(df$sample_name, df[[i]])}))

  df_out <- data.frame(mass=names(df), corr=corrs)
  df_out <- df_out[order(abs(df_out$corr), decreasing=T),]

  write.table(df_out, '~/workspace/Laura_data_analysis/peak_correlation_table.txt', quote = F, row.names = F)

  hist(corrs, seq(-1,1,0.05))
  plot(corrs)

  inds_rev <- order(abs(pc1_loading))

  for(i in head(inds_rev,20)){
    plot(df$sample_name, df[[i]])
  }

}

#import all spectrum file from the path
importAll <- function(root_data_path, ext='mzML'){
  p <- '*.mzML'
  #root_data_path <- '~/DataNAS/LauraData/March_mzML/'
  files <- dir(root_data_path, pattern=p, recursive = T, full.names = T)
  ms_spec_list <- lapply(files[1:10], function(f){import(f)})
  return(ms_spec_list)
}






