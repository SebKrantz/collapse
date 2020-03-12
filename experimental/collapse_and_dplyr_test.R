res <- GGDC10S %>%
  group_by(Variable,Country) %>%
  select_at(6:16) %>% list(fmedian(.), fmean(.))
View(res)

GGDC10S %>%
  group_by(Variable,Country) %>% {
    cbind(group_keys(.),
          add_stub(fmean(get_vars(., 6:8)), "mean_"),
          add_stub(fmedian(get_vars(., 9:12)), "median_"),
          add_stub(fmin(get_vars(., 9:10)), "min_"))
  } %>% head

data <- qDT(list(g1 = sample.int(1e5, 1e6, TRUE), g2 = sample.int(1e4, 1e6, TRUE)))
add_vars(data) <- replicate(5, rnorm(1e6), FALSE)

str(data)

system.time(collapv(data, 1:2, fsum))

system.time(data[, lapply(.SD, sum, na.rm = TRUE), by = "g1,g2"])

system.time(summarise_all(group_by(data,g1,g2), sum, na.rm = TRUE))

system.time(fsum(group_by(data,g1,g2)))

system.time(group_by(data,g1,g2))

data = group_by(data,g1,g2)

system.time(summarise_all(data, sum, na.rm = TRUE))
system.time(fsum(data))

system.time(GRP.default(data, 1:2))

system.time(qG(data$g2))

data <- data.frame()
for (i in 1:200) {
  temp <- GGDC10S
  get_vars(temp, c(1,4)) <- lapply(get_vars(temp, c(1,4)), paste0, i)
  data <- rbind(data, temp)
}
data <- data %>% group_by(Variable,Country) %>% select_at(6:16)


