# This file creates datasets from NMEG Ameriflux files and calculates
# Annual sums of fluxes.
#
# Greg Maurer - Nov 23, 2014

setwd('~/current/NMEG_seasonality/')

source('~/current/NMEG_misc/r_functions/printfigs.r')
source('~/current/NMEG_misc/r_functions/get_daily_data.r')

proc_path <- '../processed_data/'

library(ggplot2)
library(plyr)
library(reshape2)

theme_set(theme_bw())

# Look at cumulative NEE over a 5 year period, divided into the 3 seasons
FC.m <- melt(FC_daily, id.vars=c('season', 'year_w'))
FC.m <- ddply(FC.m, .(variable, season, year_w), summarise, sum = sum(value, na.rm = T))

NEE_seas_by_yr <- ggplot(FC.m, aes(x=year_w, y=sum, fill=factor(season))) +
  geom_bar(stat='identity', position='dodge') + facet_wrap(~variable, nrow=1) + 
  xlab('Year') + ylab(bquote('NEE (cumulative '~ g/m^2 ~')')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

# Same for GPP
GPP.m <- melt(GPP_daily, id.vars=c('season', 'year_w'))
GPP.m <- ddply(GPP.m, .(variable, season, year_w), summarise, sum = sum(value, na.rm = T))

GPP_seas_by_yr <- ggplot(GPP.m, aes(x=year_w, y=sum, fill=factor(season))) +
  geom_bar(stat='identity') + facet_wrap(~variable, nrow=1) + 
  xlab('Year') + ylab(bquote('GPP (cumulative '~ g/m^2 ~')')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

# Same for RE
RE.m <- melt(RE_daily, id.vars=c('season', 'year_w'))
RE.m <- ddply(RE.m, .(variable, season, year_w), summarise, sum = sum(value, na.rm = T))

RE_seas_by_yr <- ggplot(RE.m, aes(x=year_w, y=sum, fill=factor(season))) +
  geom_bar(stat='identity') + facet_wrap(~variable, nrow=1) + 
  xlab('Year') + ylab(bquote('RE (cumulative '~ g/m^2 ~')')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")


# Look at 5 year cumulative GPP, divided into the 3 seasons
GPP.m <- melt(GPP_daily, id.vars='season')
GPP.m <- subset(GPP.m, subset=(variable != 'year_w'))
GPP.m <- ddply(GPP.m, .(variable, season), summarise, sum = sum(value, na.rm = T))

GPP_seas_multiyr <- ggplot(GPP.m, aes(x=variable, y=sum, fill=factor(season))) +
  geom_bar(stat='identity') + xlab('Site') + 
  ylab(bquote('GPP ('~ g/m^2 ~', 5 year sum)')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

# Look at cumulative RE over a 5 year period, divided into the 3 seasons
RE.m <- melt(RE_daily, id.vars='season')
RE.m <- subset(RE.m, subset=(variable != 'year_w'))
RE.m <- ddply(RE.m, .(variable, season), summarise, sum = sum(value, na.rm = T))

RE_seas_multiyr <- ggplot(RE.m, aes(x=variable, y=sum, fill=factor(season))) +
  geom_bar(stat='identity') + xlab('Site') + 
  ylab(bquote('RE ('~ g/m^2 ~', 5 year sum)')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

# Look at cumulative Precip over a 6 year period, divided into the 3 seasons
P.m <- melt(P_daily, id.vars=c('season', 'year_w'))
P.m <- ddply(P.m, .(variable, season, year_w), summarise, sum = sum(value, na.rm = T))

P_seas_by_yr <- ggplot(P.m, aes(x=year_w, y=sum, fill=factor(season))) +
  geom_bar(stat='identity') + facet_wrap(~ variable, nrow=1) +
  xlab('Year') + ylab('Precipitation (mm)') + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

P.m <- melt(P_daily, id.vars='season')
P.m <- subset(P.m, subset=(variable != 'year_w'))
P.m <- ddply(P.m, .(variable, season), summarise, sum = sum(value, na.rm = T))

P_seas_multiyr <- ggplot(P.m, aes(x=variable, y=sum, fill=factor(season))) +
  geom_bar(stat='identity') + xlab('Site') + ylab('Precipitation (mm, 6 year sum)') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

printfigs(P_seas_multiyr, 'P_seas_multiyr.svg', 7, 6)
printfigs(GPP_seas_by_yr, 'GPP_seas_by_yr.svg', 12, 3)
printfigs(RE_seas_by_yr, 'RE_seas_by_yr.svg', 12, 3)
printfigs(NEE_seas_by_yr, 'NEE_seas_by_yr.svg', 12, 3)
printfigs(P_seas_by_yr, 'P_seas_by_yr.svg', 12, 3)

# Get seasonal FC totals for each season at all sites and yearly total
get_precip_flux_dep <- function(df_daily){
  flux.m <- melt(df_daily, id.vars=c('season', 'year_w'))
  flux_tot <- ddply(flux.m, .(variable, year_w),
                    summarise, year_flux_tot = sum(value, na.rm = T))
  flux_wyr_seas.m <- ddply(flux.m, .(variable, season, year_w),
                           summarise, year_seas_tot = sum(value, na.rm = T))
  
  # subset spring flux
  flux_spring <- subset(flux_wyr_seas.m, subset=(season=='Spring'))
  flux_spring <- rename(flux_spring, c("year_seas_tot"="Spring_flux_tot"))
  flux_spring$season <- NULL
  
  # subset monsoon flux
  flux_monsoon <- subset(flux_wyr_seas.m, subset=(season=='Monsoon'))
  flux_monsoon <- rename(flux_monsoon, c("year_seas_tot"="Monsoon_flux_tot"))
  flux_monsoon$season <- NULL
  
  # Get seasonal Precip totals for each season at all sites
  P.m <- melt(P_daily, id.vars=c('season', 'year_w'))
  P_wyr_seas.m <- ddply(P.m, .(variable, season, year_w),
                        summarise, P_tot = sum(value, na.rm = T))
  
  # Merge in flux_spring, flux_monsoon, and flux_tot by site/year
  P_dep <- P_wyr_seas.m
  P_dep <- merge(P_dep, flux_tot, by=c('variable', 'year_w'))
  P_dep <- merge(P_dep, flux_spring, by=c('variable', 'year_w'))
  P_dep <- merge(P_dep, flux_monsoon, by=c('variable', 'year_w'))
  
  return(P_dep)
}

NEE_P_dep <- get_precip_flux_dep(FC_daily)
GPP_P_dep <- get_precip_flux_dep(GPP_daily)
#  Yearly, Spring, and Monsoon fluxes vs cold season precip
subs <- subset(GPP_P_dep, subset=(season=='Cold'))
names(subs)[1]<-"Site"
subs = melt(subs, measure.vars=c("year_flux_tot", "Spring_flux_tot", "Monsoon_flux_tot"))
GPP_on_coldP <- ggplot(subs, aes(x=P_tot, y=value, col=variable)) +
  geom_point(size=3) + stat_smooth(method=lm, alpha = 0.1) +
  facet_wrap(~ Site, scales = 'free') + xlab('Cold season precip (cumulative mm)') + 
  ylab(bquote('GPP (cumulative '~ g/m^2 ~')')) +
  scale_colour_manual(values=c("black", "#999999", "#E69F00"),
                      labels=c("Full year", "Spring", "Monsoon"),
                      name="Season")

#  Yearly, Spring, and Monsoon fluxes vs spring season precip
subs <- subset(GPP_P_dep, subset=(season=='Spring'))
names(subs)[1]<-"Site"
subs = melt(subs, measure.vars=c("year_flux_tot", "Spring_flux_tot", "Monsoon_flux_tot"))
GPP_on_springP <- ggplot(subs, aes(x=P_tot, y=value, col=variable)) +
  geom_point(size=3) + stat_smooth(method=lm, alpha = 0.1) +
  facet_wrap(~ Site, scales = 'free') + xlab('Spring precip (cumulative mm)') + 
  ylab(bquote('GPP (cumulative '~ g/m^2 ~')')) +
  scale_colour_manual(values=c("black", "#999999", "#E69F00"),
                      labels=c("Full year", "Spring", "Monsoon"),
                      name="Season")

#  Yearly, Spring, and Monsoon fluxes vs monsoon season precip
subs <- subset(GPP_P_dep, subset=(season=='Monsoon'))
names(subs)[1]<-"Site"
subs = melt(subs, measure.vars=c("year_flux_tot", "Monsoon_flux_tot"))
GPP_on_monsoonP <- ggplot(subs, aes(x=P_tot, y=value, col=variable)) +
  geom_point(size=3) + stat_smooth(method=lm, alpha = 0.1) +
  facet_wrap(~ Site, scales = 'free') + xlab('Monsoon precip (cumulative mm)') + 
  ylab(bquote('GPP (cumulative '~ g/m^2 ~')')) +
  scale_colour_manual(values=c("black", "#E69F00"),
                      labels=c("Full year", "Monsoon"),
                      name="Season")


printfigs(GPP_on_monsoonP, 'GPP_on_monsoonP.svg', 8.5, 5)
printfigs(GPP_on_coldP, 'GPP_on_coldP.svg', 8.5, 5)


# On to soil water content
deep <- grep('30_AVG', colnames(VWC_daily), value=T)
shallow <- grep('5_AVG', colnames(VWC_daily), value=T)

VWC_daily['Deep'] <- rowMeans(VWC_daily[deep], na.rm=T)
VWC_daily['Shallow'] <- rowMeans(VWC_daily[shallow], na.rm=T)


# Make a smaller dataset to melt
VWC <- VWC_daily[,c('year_month_mday','year_w','doy_w','season','Deep','Shallow')]
VWC.m <- melt(VWC, id.vars=c('year_month_mday', 'year_w', 'doy_w', 'season'))

VWC_ts <- ggplot(VWC.m, aes(x=as.Date(year_month_mday), y=value, color=variable)) +
  geom_line() + scale_color_manual(name="Soil depth",values=c("blue", "green")) +
  ylab('Volumetric Wate Content') + xlab('Date')

printfigs(VWC_ts, 'VWC_ts.svg', 10, 2)

VWC.m <- ddply(VWC.m, .(variable, season, year_w),
               summarise, VWC_med = median(value, na.rm = T))

VWC_seas_by_yr <- ggplot(VWC.m, aes(x=year_w, y=VWC_med, fill=factor(season))) +
  geom_bar(stat='identity', position='dodge') + facet_wrap(~ variable, nrow=1) +
  xlab('Year') + ylab('Mean VWC') + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
  scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                    name="Season")

printfigs(VWC_seas_by_yr, 'VWC_seas_by_yr.svg', 6.5, 4)

GPP_sub <- subset(GPP_P_dep, subset=(year_w > 2008 & variable=='PJ woodland'))
VWC.m <- merge(VWC.m, GPP_sub, by=c('season', 'year_w'))

GPP_vs_VWC <- ggplot(VWC.m, aes(x=VWC_med, y=year_flux_tot, col=factor(season))) +
  geom_point(size=3) + stat_smooth(method=lm, alpha=0.1, fullrange=T) +
  facet_wrap(~variable.x, nrow=1) + xlab('Median VWC') + 
  ylab(bquote('Yearly GPP (cumulative '~ g/m^2 ~')')) +
  scale_colour_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
                      name="Season")

printfigs(GPP_vs_VWC, 'GPP_vs_VWC.svg', 6.5, 4)

# #  Yearly GPP vs total precip in each of 3 seasons
# ggplot(GPP_P_dep, aes(x=P_tot, y=year_flux_tot, col=factor(season))) +
#   geom_point() + facet_wrap(~ variable)
# 
# #  Spring NEE vs total precip in each of 3 seasons
# ggplot(NEE_P_dep, aes(x=P_tot, y=Spring_flux_tot, col=factor(season))) +
#   geom_point() + facet_wrap(~ variable)
# #  Spring GPP vs total precip in each of 3 seasons
# ggplot(GPP_P_dep, aes(x=P_tot, y=Spring_flux_tot, col=factor(season))) +
#   geom_point() + facet_wrap(~ variable)
# 
# #  Monsoon NEE vs total precip in each of 3 seasons
# ggplot(NEE_P_dep, aes(x=P_tot, y=Monsoon_flux_tot, col=factor(season))) +
#   geom_point() + facet_wrap(~ variable)
# #  Monsoon GPP vs total precip in each of 3 seasons
# ggplot(GPP_P_dep, aes(x=P_tot, y=Monsoon_flux_tot, col=factor(season))) +
#   geom_point() + facet_wrap(~ variable)
# 
# + xlab('Site') + ylab('Precipitation (mm, 6 year sum)') +
#   theme(axis.text.x=element_text(angle=45, hjust=1, size=11)) + 
#   scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00"), 
#                     name="Season")



