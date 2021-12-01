### INFO: tutorial from https://www.highcharts.com/blog/data-science/highcharts-for-r-users/
### DATE: Tue May 07 07:53:03 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test/other")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(timetk)
library(kableExtra)
library(highcharter)
library(quantmod)
library(PerformanceAnalytics)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# get and clean data
symbols <- c("SPY","EFA", "IJS", "EEM","AGG")

prices <- 
  getSymbols(symbols, 
             src = 'yahoo', 
             from = "2013-01-01",
             to = "2017-12-31",
             auto.assign = TRUE, 
             warnings = FALSE) %>% 
  map(~Ad(get(.))) %>%
  reduce(merge) %>% 
  `colnames<-`(symbols)

prices_monthly <- to.monthly(prices, 
                             indexAt = "last", 
                             OHLC = FALSE)

asset_returns_xts <- na.omit(Return.calculate(prices_monthly, method = "log"))

asset_returns_xts <- asset_returns_xts * 100

asset_returns_long <-  
  prices %>% 
  to.monthly(indexAt = "last", OHLC = FALSE) %>% 
  tk_tbl(preserve_index = TRUE, rename_index = "date") %>%
  gather(asset, returns, -date) %>% 
  group_by(asset) %>%  
  mutate(returns = (log(returns) - log(lag(returns))) *100) %>% 
  na.omit()


### plot 1
p <- 
  highchart(type = "stock") %>% 
  hc_add_series(asset_returns_xts$SPY, type = "line")

htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "plot1.html"),
                        selfcontained = T)

### plot 2
p <- 
  highchart(type = "stock") %>% 
  hc_add_series(asset_returns_xts$SPY, type = "column")

htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "plot2.html"),
                        selfcontained = T)

### plot 3
p <- 
  highchart(type = "stock") %>% 
  hc_add_series(asset_returns_xts$SPY, type = "scatter", name = "SPY") %>% 
  hc_add_series(asset_returns_xts$EFA, type = "scatter", name = "EFA") %>% 
  hc_tooltip(pointFormat = '{point.x: %Y-%m-%d}
             
             
             {point.y:.4f}%')

htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "plot3.html"),
                        selfcontained = T)

### plot 4
p <- 
  asset_returns_long %>% 
  filter(asset == "SPY") %>% 
  hchart(., 
         type = "line", 
         hcaes(x = date, 
               y = returns)) %>% 
  hc_yAxis(opposite = TRUE,
           labels = list(format = "{value}%"))

htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "plot4.html"),
                        selfcontained = T)

### plot 5
p <- 
  asset_returns_long %>% 
  hchart(., 
         type = "scatter", 
         hcaes(x = date, 
               y = returns, 
               group = asset)) %>% 
  hc_yAxis(opposite = TRUE,
           labels = list(format = "{value}%")) %>% 
  hc_tooltip(pointFormat = '{point.x:%Y-%m-%d}
                             

                             {point.y: .4f}%')

htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "plot5.html"),
                        selfcontained = T)

### plot 6
p <- 
  asset_returns_long %>% 
  filter(date >= "2017-01-01" & date < "2018-01-01") %>% 
  hchart(., type = "column", 
         hcaes(x = date, 
               y = returns, 
               group = asset)) %>% 
  hc_yAxis(opposite = FALSE,
           labels = list(format = "{value}%")) %>% 
  hc_tooltip(pointFormat = '{point.x: %Y-%m-%d}
             
             
             {point.y:.4f}% ')

htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "plot6.html"),
                        selfcontained = T)

