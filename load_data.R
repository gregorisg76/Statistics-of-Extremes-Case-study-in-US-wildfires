#############################
#### DATA COMPETITION #######
#############################
library(INLA)
library(tidyverse)
library(Hmisc)
load("data/data_train.RData")
source("utils.R")
## data_train_DF
## -------------
data_train_DF <- data_train_DF %>% mutate(key=1:nrow(data_train_DF))
df <- data_train_DF %>% group_by(lon, lat) %>% nest()
## sort by latitude
## df <- data_train_DF %>% arrange(lat) %>% group_by(lon, lat) %>% nest()
df.id <- data.frame(cell_ID=1:nrow(df))
df    <- bind_cols(df, df.id) %>% unnest(data) %>% arrange(cell_ID) %>%
    mutate(
        date       = as.Date(paste0(year, paste0(0,month), "01" ),  "%Y%m%d"),
        no_of_days = monthDays(date),
        RH         = exp(17.502*(clim3 - 273.16)/(clim3 - 32.19)) / exp(17.502*(clim4 - 273.16)/(clim4 - 32.19)),
        WS         = sqrt(clim1^2 + clim2^2),
        Wdir       = atan2(clim1/WS, clim2/WS),
        month      = month.abb[month],
        gntr       = clim6-clim7,
        month_year = paste0(month, year),
        month_cell_ID = paste0(month, cell_ID),
        year_cell_ID = paste0(year, cell_ID),
        year_month_cell_ID = paste0(year, month, cell_ID))  %>% ungroup %>%
    arrange(date) %>% group_by(cell_ID) %>% nest() %>%
    mutate(data=map(data, function(x) {
        x.nest <- x %>% group_by(year) %>% nest %>%
            mutate(data=map(data, function(y){
                ## -------------------                
                ## **Lags and Leads**
                ## -------------------
                ## A problem here is the missing of lag 1 and lag nrow(y)
                ## values for all March and September months recorded respectively.
                ## Approach below
                ## imputes by respective months' values                
                y <- y %>% mutate(clim10.lag = dplyr::lag(clim10),
                                  clim5.lead = dplyr::lead(clim5),
                                  clim7.lead = dplyr::lead(clim7),
                                  clim9.lead = dplyr::lead(clim9))
                y$clim10.lag[1]       <- y$clim10[1]
                y$clim5.lead[nrow(y)] <- y$clim5[nrow(y)]
                y$clim7.lead[nrow(y)] <- y$clim7[nrow(y)]
                y$clim9.lead[nrow(y)] <- y$clim9[nrow(y)]
                y
            }))
        x <- x.nest %>% unnest(data)
        x
    })) %>% unnest(data)


if(FALSE){
    load("mclust.RData")
    df <- df %>% mutate(land.label=mod.land$classification)
    
}



coords <- df %>% ungroup %>% dplyr::select(lon, lat) %>% unique %>% as.matrix
## ---------------------------------------------------
## split domain in regions
## ---------------------------------------------------
J <- 4
qlon <- unname(quantile(coords[,1], seq(0,1,len=J)))
qlat <- unname(quantile(coords[,2], seq(0,1,len=J)))
mat.lon <- do.call(rbind, replicate(length(qlon)-1, cbind(lag(qlon), qlon)[-1,], simplify=FALSE))
mat.lat <- do.call(rbind, replicate(length(qlat)-1, cbind(lag(qlat), qlat)[-1,], simplify=FALSE))
## store intervals in matrix. col1 and col2 are lon.low and lon.high,
## col3 and col4 are lat.low and lat high for a box (row of matrix)
mat.lon.lat.intervals <- as.matrix(data.frame(lon.low = mat.lon[,1], lon.high = mat.lon[,2]) %>%
                                   arrange(lon.low) %>%
                                   mutate(lat.low=mat.lat[,1], lat.high=mat.lat[,2],
                                          region.index = 1:((length(qlon)-1)*(length(qlat)-1))))
df$region.index <- unlist(lapply(1:nrow(df), function(k) which.region.index(df$lon[k],df$lat[k],mat.lon.lat.intervals=mat.lon.lat.intervals)))


if(FALSE){ #plot regions
    plot(df %>% ungroup %>%dplyr::select(lon,lat) %>% unique, asp=1)
    for(i in 1:(J^2)) points(df %>% ungroup %>%filter(region.index ==i) %>% dplyr::select(lon,lat) %>% unique, col=i+1, pch=16)
}


## ------------------------------------------------
## data frames for zeroinflated0 INLA likelihood
## ------------------------------------------------

## small data frame
df.small      <- df %>% filter(cell_ID <=100)

df.small.prob <- df.small %>% mutate(CNT_binary = ifelse(CNT>0, 1, 0), CNT = NA)
df.small.rate <- df.small %>% filter(CNT > 0) %>% mutate(CNT_binary=NA)
df.small.zeroinflated <- rbind(df.small.prob, df.small.rate)
n.prob.small <- nrow(df.small.prob)
n.rate.small <- nrow(df.small.rate)

ldf.small.zeroinflated <- list(
    CNT = df.small.zeroinflated %>% ungroup %>% dplyr::select(CNT_binary, CNT) %>% as.matrix,
    Intercept.prob = c(rep(1, n.prob.small), rep(NA,n.rate.small)), Intercept.rate = c(rep(NA,n.prob.small), rep(1,n.rate.small)),
    lc1.prob = c(df.small.prob$lc1, rep(NA,n.rate.small)), lc1.rate = c(rep(NA,n.prob.small), df.small.rate$lc1),
    lc2.prob = c(df.small.prob$lc2, rep(NA,n.rate.small)), lc2.rate = c(rep(NA,n.prob.small), df.small.rate$lc2),
    lc3.prob = c(df.small.prob$lc3, rep(NA,n.rate.small)), lc3.rate = c(rep(NA,n.prob.small), df.small.rate$lc3),
    lc4.prob = c(df.small.prob$lc4, rep(NA,n.rate.small)), lc4.rate = c(rep(NA,n.prob.small), df.small.rate$lc4),
    lc5.prob = c(df.small.prob$lc5, rep(NA,n.rate.small)), lc5.rate = c(rep(NA,n.prob.small), df.small.rate$lc5),
    lc6.prob = c(df.small.prob$lc6, rep(NA,n.rate.small)), lc6.rate = c(rep(NA,n.prob.small), df.small.rate$lc6),
    lc7.prob = c(df.small.prob$lc7, rep(NA,n.rate.small)), lc7.rate = c(rep(NA,n.prob.small), df.small.rate$lc7),
    lc8.prob = c(df.small.prob$lc8, rep(NA,n.rate.small)), lc8.rate = c(rep(NA,n.prob.small), df.small.rate$lc8),
    lc9.prob = c(df.small.prob$lc9, rep(NA,n.rate.small)), lc9.rate = c(rep(NA,n.prob.small), df.small.rate$lc9),
    lc10.prob = c(df.small.prob$lc10, rep(NA,n.rate.small)), lc10.rate = c(rep(NA,n.prob.small), df.small.rate$lc10),
    lc11.prob = c(df.small.prob$lc11, rep(NA,n.rate.small)), lc11.rate = c(rep(NA,n.prob.small), df.small.rate$lc11),
    lc12.prob = c(df.small.prob$lc12, rep(NA,n.rate.small)), lc12.rate = c(rep(NA,n.prob.small), df.small.rate$lc12),
    lc13.prob = c(df.small.prob$lc13, rep(NA,n.rate.small)), lc13.rate = c(rep(NA,n.prob.small), df.small.rate$lc13),
    lc14.prob = c(df.small.prob$lc14, rep(NA,n.rate.small)), lc14.rate = c(rep(NA,n.prob.small), df.small.rate$lc14),
    lc15.prob = c(df.small.prob$lc15, rep(NA,n.rate.small)), lc15.rate = c(rep(NA,n.prob.small), df.small.rate$lc15),
    lc16.prob = c(df.small.prob$lc16, rep(NA,n.rate.small)), lc16.rate = c(rep(NA,n.prob.small), df.small.rate$lc16),
    lc17.prob = c(df.small.prob$lc17, rep(NA,n.rate.small)), lc17.rate = c(rep(NA,n.prob.small), df.small.rate$lc17),
    lc18.prob = c(df.small.prob$lc18, rep(NA,n.rate.small)), lc18.rate = c(rep(NA,n.prob.small), df.small.rate$lc18),
    clim1.prob = c(df.small.prob$clim1, rep(NA,n.rate.small)), clim1.rate = c(rep(NA,n.prob.small), df.small.rate$clim1),
    clim2.prob = c(df.small.prob$clim2, rep(NA,n.rate.small)), clim2.rate = c(rep(NA,n.prob.small), df.small.rate$clim2),
    clim3.prob = c(df.small.prob$clim3, rep(NA,n.rate.small)), clim3.rate = c(rep(NA,n.prob.small), df.small.rate$clim3),
    clim4.prob = c(df.small.prob$clim4, rep(NA,n.rate.small)), clim4.rate = c(rep(NA,n.prob.small), df.small.rate$clim4),
    clim5.prob = c(df.small.prob$clim5, rep(NA,n.rate.small)), clim5.rate = c(rep(NA,n.prob.small), df.small.rate$clim5),
    clim5.lead.prob = c(df.small.prob$clim5.lead, rep(NA,n.rate.small)), clim5.lead.rate = c(rep(NA,n.prob.small), df.small.rate$clim5.lead),
    clim6.prob = c(df.small.prob$clim6, rep(NA,n.rate.small)), clim6.rate = c(rep(NA,n.prob.small), df.small.rate$clim6),
    clim7.prob = c(df.small.prob$clim7, rep(NA,n.rate.small)), clim7.rate = c(rep(NA,n.prob.small), df.small.rate$clim7),
    clim8.prob = c(df.small.prob$clim8, rep(NA,n.rate.small)), clim8.rate = c(rep(NA,n.prob.small), df.small.rate$clim8),
    clim9.prob = c(df.small.prob$clim9, rep(NA,n.rate.small)), clim9.rate = c(rep(NA,n.prob.small), df.small.rate$clim9),
    clim9.lead.prob = c(df.small.prob$clim9.lead, rep(NA,n.rate.small)), clim9.lead.rate = c(rep(NA,n.prob.small), df.small.rate$clim9.lead),
    clim10.prob = c(df.small.prob$clim10, rep(NA,n.rate.small)), clim10.rate = c(rep(NA,n.prob.small), df.small.rate$clim10),
    clim10.lag.prob = c(df.small.prob$clim10.lag, rep(NA,n.rate.small)), clim10.lag.rate = c(rep(NA,n.prob.small), df.small.rate$clim10.lag),
    altiMean.prob.small = c(df.small.prob$altiMean, rep(NA,n.rate.small)), altiMean.rate.small = c(rep(NA,n.prob.small), df.small.rate$altiMean),
    altiSD.prob = c(df.small.prob$altiSD, rep(NA,n.rate.small)), altiSD.rate = c(rep(NA,n.prob.small), df.small.rate$altiSD),
    RH.prob = c(df.small.prob$RH, rep(NA,n.rate.small)), RH.rate = c(rep(NA,n.prob.small), df.small.rate$RH),
    WS.prob = c(df.small.prob$WS, rep(NA,n.rate.small)), WS.rate = c(rep(NA,n.prob.small), df.small.rate$WS),
    cell_ID.prob = c(df.small.prob$cell_ID, rep(NA,n.rate.small)), cell_ID.rate = c(rep(NA,n.prob.small), df.small.rate$cell_ID),
    year.prob = c(df.small.prob$year, rep(NA,n.rate.small)), year.rate = c(rep(NA,n.prob.small), df.small.rate$year),
    year_cell_ID.prob = c(df.small.prob$year_cell_ID, rep(NA,n.rate.small)), year_cell_ID.rate = c(rep(NA,n.prob.small), df.small.rate$year_cell_ID),
    year_month_cell_ID.prob = c(df.small.prob$year_month_cell_ID, rep(NA,n.rate.small)), year_month_cell_ID.rate = c(rep(NA,n.prob.small), df.small.rate$year_month_cell_ID),
    month.prob = c(df.small.prob$month, rep(NA,n.rate.small)), month.rate = c(rep(NA,n.prob.small), df.small.rate$month))




## for full data set
df.prob <- df %>% mutate(CNT_binary = ifelse(CNT>0, 1, 0), CNT = NA)
df.rate <- df %>% filter(CNT > 0) %>% mutate(CNT_binary=NA)
df.zeroinflated <- rbind(df.prob, df.rate)
n.prob <- nrow(df.prob)
n.rate <- nrow(df.rate)
ldf.zeroinflated <- list(
    CNT = df.zeroinflated %>% ungroup %>% dplyr::select(CNT_binary, CNT) %>% as.matrix,
    Intercept.prob = c(rep(1, n.prob), rep(NA,n.rate)), Intercept.rate = c(rep(NA,n.prob), rep(1,n.rate)),
    lc1.prob = c(df.prob$lc1, rep(NA,n.rate)), lc1.rate = c(rep(NA,n.prob), df.rate$lc1),
    lc2.prob = c(df.prob$lc2, rep(NA,n.rate)), lc2.rate = c(rep(NA,n.prob), df.rate$lc2),
    lc3.prob = c(df.prob$lc3, rep(NA,n.rate)), lc3.rate = c(rep(NA,n.prob), df.rate$lc3),
    lc4.prob = c(df.prob$lc4, rep(NA,n.rate)), lc4.rate = c(rep(NA,n.prob), df.rate$lc4),
    lc5.prob = c(df.prob$lc5, rep(NA,n.rate)), lc5.rate = c(rep(NA,n.prob), df.rate$lc5),
    lc6.prob = c(df.prob$lc6, rep(NA,n.rate)), lc6.rate = c(rep(NA,n.prob), df.rate$lc6),
    lc7.prob = c(df.prob$lc7, rep(NA,n.rate)), lc7.rate = c(rep(NA,n.prob), df.rate$lc7),
    lc8.prob = c(df.prob$lc8, rep(NA,n.rate)), lc8.rate = c(rep(NA,n.prob), df.rate$lc8),
    lc9.prob = c(df.prob$lc9, rep(NA,n.rate)), lc9.rate = c(rep(NA,n.prob), df.rate$lc9),
    lc10.prob = c(df.prob$lc10, rep(NA,n.rate)), lc10.rate = c(rep(NA,n.prob), df.rate$lc10),
    lc11.prob = c(df.prob$lc11, rep(NA,n.rate)), lc11.rate = c(rep(NA,n.prob), df.rate$lc11),
    lc12.prob = c(df.prob$lc12, rep(NA,n.rate)), lc12.rate = c(rep(NA,n.prob), df.rate$lc12),
    lc13.prob = c(df.prob$lc13, rep(NA,n.rate)), lc13.rate = c(rep(NA,n.prob), df.rate$lc13),
    lc14.prob = c(df.prob$lc14, rep(NA,n.rate)), lc14.rate = c(rep(NA,n.prob), df.rate$lc14),
    lc15.prob = c(df.prob$lc15, rep(NA,n.rate)), lc15.rate = c(rep(NA,n.prob), df.rate$lc15),
    lc16.prob = c(df.prob$lc16, rep(NA,n.rate)), lc16.rate = c(rep(NA,n.prob), df.rate$lc16),
    lc17.prob = c(df.prob$lc17, rep(NA,n.rate)), lc17.rate = c(rep(NA,n.prob), df.rate$lc17),
    lc18.prob = c(df.prob$lc18, rep(NA,n.rate)), lc18.rate = c(rep(NA,n.prob), df.rate$lc18),
    clim1.prob = c(df.prob$clim1, rep(NA,n.rate)), clim1.rate = c(rep(NA,n.prob), df.rate$clim1),
    clim2.prob = c(df.prob$clim2, rep(NA,n.rate)), clim2.rate = c(rep(NA,n.prob), df.rate$clim2),
    clim3.prob = c(df.prob$clim3, rep(NA,n.rate)), clim3.rate = c(rep(NA,n.prob), df.rate$clim3),
    clim4.prob = c(df.prob$clim4, rep(NA,n.rate)), clim4.rate = c(rep(NA,n.prob), df.rate$clim4),
    clim5.prob = c(df.prob$clim5, rep(NA,n.rate)), clim5.rate = c(rep(NA,n.prob), df.rate$clim5),
    clim5.lead.prob = c(df.prob$clim5.lead, rep(NA,n.rate)), clim5.lead.rate = c(rep(NA,n.prob), df.rate$clim5.lead),
    clim6.prob = c(df.prob$clim6, rep(NA,n.rate)), clim6.rate = c(rep(NA,n.prob), df.rate$clim6),
    clim7.prob = c(df.prob$clim7, rep(NA,n.rate)), clim7.rate = c(rep(NA,n.prob), df.rate$clim7),
    clim8.prob = c(df.prob$clim8, rep(NA,n.rate)), clim8.rate = c(rep(NA,n.prob), df.rate$clim8),
    clim9.prob = c(df.prob$clim9, rep(NA,n.rate)), clim9.rate = c(rep(NA,n.prob), df.rate$clim9),
    clim9.lead.prob = c(df.prob$clim9.lead, rep(NA,n.rate)), clim9.lead.rate = c(rep(NA,n.prob), df.rate$clim9.lead),
    clim10.prob = c(df.prob$clim10, rep(NA,n.rate)), clim10.rate = c(rep(NA,n.prob), df.rate$clim10),
    clim10.lag.prob = c(df.prob$clim10.lag, rep(NA,n.rate)), clim10.lag.rate = c(rep(NA,n.prob), df.rate$clim10.lag),
    altiMean.prob = c(df.prob$altiMean, rep(NA,n.rate)), altiMean.rate = c(rep(NA,n.prob), df.rate$altiMean),
    altiSD.prob = c(df.prob$altiSD, rep(NA,n.rate)), altiSD.rate = c(rep(NA,n.prob), df.rate$altiSD),
    RH.prob = c(df.prob$RH, rep(NA,n.rate)), RH.rate = c(rep(NA,n.prob), df.rate$RH),
    WS.prob = c(df.prob$WS, rep(NA,n.rate)), WS.rate = c(rep(NA,n.prob), df.rate$WS),
    cell_ID.prob = c(df.prob$cell_ID, rep(NA,n.rate)), cell_ID.rate = c(rep(NA,n.prob), df.rate$cell_ID),
    year.prob = c(df.prob$year, rep(NA,n.rate)), year.rate = c(rep(NA,n.prob), df.rate$year),
    year_cell_ID.prob = c(df.prob$year_cell_ID, rep(NA,n.rate)), year_cell_ID.rate = c(rep(NA,n.prob), df.rate$year_cell_ID),
    year_month_cell_ID.prob = c(df.prob$year_month_cell_ID, rep(NA,n.rate)), year_month_cell_ID.rate = c(rep(NA,n.prob), df.rate$year_month_cell_ID),
    month.prob = c(df.prob$month, rep(NA,n.rate)), month.rate = c(rep(NA,n.prob), df.rate$month))





## set link function for predictions in INLA

## link = rep(NA, n.prob + n.rate)
## link[which(is.na(ldf.zeroinflated$CNT[,1]))] = 1
## link[which(is.na(ldf.zeroinflated$CNT[,2]))] = 2
