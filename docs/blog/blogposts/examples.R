library(tradestatistics)
`%>%` <- magrittr::`%>%`
us_trade <- ots_create_tidy_data(years = 2000:2018, reporters = "usa", table = "yrpc") %>%
  ots_inflation_adjustment

fdim(us_trade)

fndistinct(us_trade)
qsu(with(us_trade, fnobs(year, list(partner_iso, product_code))))
qsu(us_trade, by = export_value_usd + import_value_usd ~ year, pid = ~ partner_iso + year)
us_trade_agg <- collap(us_trade, ~ partner_iso + group_code + year, fsum)
qsu(with(us_trade_agg, fnobs(year, list(partner_iso, group_code))))

L(us_trade_agg, 1, ~ partner_iso + group_code)

settfm(us_trade_agg, id = finteraction(partner_iso, group_code))
settfm(us_trade_agg, id = seqid(year, radixorder(partner_iso, group_code, year)))

View(slt(us_trade_agg, id, year, t))

View(L(us_trade_agg, 1, ~id, ~year))

settransform(us_trade, t = seqid(year, radixorder(partner_iso, product_code, year)))
View(slt(us_trade, year, t))

us_trade <- ftransform(us_trade,
             export_value_usd = replace_NA(export_value_usd, 0),
             import_value_usd = replace_NA(import_value_usd, 0)) %>%
             tfm(value_usd = export_value_usd + import_value_usd) %>%
             tfm(partner_value_usd = fsum(value_usd, partner_iso, TRA = "replace_fill"))



plot(fdiff(EuStockMarkets, 0:1))

library(dplyr)
library(data.table) # Running with 2 threads
library(microbenchmark)
suppressMessages(
                microbenchmark(collapse = collap(us_trade, export_value_usd + import_value_usd ~ partner_iso + group_code + year, fsum),
                               data.table = us_trade[, list(export_value_usd = sum(export_value_usd, na.rm = TRUE),
                                                             import_value_usd = sum(import_value_usd, na.rm = TRUE)),
                                                      by = c("partner_iso", "group_code", "year")],
                               dplyr = group_by(us_trade, partner_iso, group_code, year) %>%
                                       dplyr::select(export_value_usd, import_value_usd) %>% summarise_all(sum, na.rm = TRUE), times = 10))

suppressMessages(
  microbenchmark(collapse = collap(us_trade, export_value_usd + import_value_usd ~ partner_iso + group_code + year, fmean),
                 data.table = us_trade[, list(export_value_usd = mean(export_value_usd, na.rm = TRUE),
                                              import_value_usd = mean(import_value_usd, na.rm = TRUE)),
                                       by = c("partner_iso", "group_code", "year")],
                 dplyr = group_by(us_trade, partner_iso, group_code, year) %>%
                   dplyr::select(export_value_usd, import_value_usd) %>% summarise_all(mean, na.rm = TRUE), times = 10))

suppressMessages(
  microbenchmark(collapse = fgroup_by(us_trade, partner_iso, group_code, year) %>%
                   fselect(export_value_usd, import_value_usd) %>% fsum(TRA = "replace_fill"),
                 data.table = us_trade[, `:=`(export_value_usd2 = sum(export_value_usd, na.rm = TRUE),
                                               import_value_usd2 = sum(import_value_usd, na.rm = TRUE)),
                                        by = c("partner_iso", "group_code", "year")],
                 dplyr = group_by(us_trade, partner_iso, group_code, year) %>%
                   dplyr::select(export_value_usd, import_value_usd) %>% mutate_all(sum, na.rm = TRUE), times = 10))

# centering, partner-product
suppressMessages(
  microbenchmark(collapse = fgroup_by(us_trade, partner_iso, product_code) %>%
                   fselect(export_value_usd, import_value_usd) %>% fwithin,
                 data.table = us_trade[, `:=`(export_value_usd2 = export_value_usd - mean(export_value_usd, na.rm = TRUE),
                                              import_value_usd2 = import_value_usd - mean(import_value_usd, na.rm = TRUE)),
                                       by = c("partner_iso", "group_code", "year")],
                 dplyr = group_by(us_trade, partner_iso, group_code, year) %>%
                   dplyr::select(export_value_usd, import_value_usd) %>% mutate_all(function(x) x - mean(x, na.rm = TRUE)), times = 1))

# Lag
setorder(us_trade, partner_iso, product_code, year)
suppressMessages(
  microbenchmark(collapse = L(us_trade, 1, export_value_usd + import_value_usd ~ partner_iso + product_code),
                 data.table = us_trade[, shift(.SD), keyby = c("partner_iso", "product_code"),
                                       .SDcols = c("export_value_usd","import_value_usd")],
                 dplyr = group_by(us_trade, partner_iso, product_code) %>%
                   dplyr::select(export_value_usd, import_value_usd) %>% mutate_all(lag), times = 1))





library(collapse)
             impo
             value_usd = export_value_usd + import_value_usd,
             ctry_exp_value = fsum(export_value_usd, partner_iso, TRA = "replace"))
View(collap(us_trade, ~ product_code, fsum, fmode, w = ~))


library(rjson)

PRIO <- fromJSON(file = "https://grid.prio.org/api/data/pop_hyd_sum")

PRIO <- jsonlite::fromJSON("https://grid.prio.org/api/data/prec_gpcc,temp,pop_hyd_sum,nlights_calib_mean/")
