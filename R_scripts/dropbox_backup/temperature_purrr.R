library(tidyverse)

setwd(dir = "C:/Users/fhorvat/Dropbox/Praksa bioinfo/Projekti/test/")

temperature_wide <- 
  read_csv("temperature.csv") %>%
  print()

temperature_tall <-
  temperature_wide %>%
  gather(key = "id_sensor", value = "temperature", starts_with("temp")) %>% 
  mutate(id_sensor = str_replace(id_sensor, "temperature_", "")) %>%
  print()

temperature_tall %>%
  ggplot(aes(x = instant, y = temperature, color = id_sensor)) +
  geom_line()

delta <- 
  temperature_tall %>%
  arrange(id_sensor, instant) %>%
  group_by(id_sensor) %>%
  mutate(delta_time = as.numeric(instant) - as.numeric(instant[[1]]),
         delta_temperature = temperature - temperature[[1]]) %>%
  select(id_sensor, delta_time, delta_temperature)

delta %>%
  ggplot(aes(x = delta_time, y = delta_temperature, color = id_sensor)) +
  geom_line()  

# definitions
# reference: http://stackoverflow.com/questions/29067916/r-error-function-erfz
# (see Abramowitz and Stegun 29.2.29)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

newton_cooling <- function(x) {
  nls(
    delta_temperature ~ delta_temperature_0*(1 - exp(-delta_time/tau_0)),
    start = list(delta_temperature_0 = -10, tau_0 = 50),
    data = x
  )
}

semi_infinite_simple <- function(x) {
  nls(
    delta_temperature ~ delta_temperature_0*erfc(sqrt(tau_0/delta_time)),
    start = list(delta_temperature_0 = -10, tau_0 = 50),
    data = x
  )    
}

semi_infinite_convection <- function(x){
  nls(
    delta_temperature ~
      delta_temperature_0*(
        erfc(sqrt(tau_0/delta_time)) -
          exp(Bi_0 + (Bi_0/2)^2*delta_time/tau_0)*
          erfc(sqrt(tau_0/delta_time) + 
                 (Bi_0/2)*sqrt(delta_time/tau_0))
      ),
    start = list(delta_temperature_0 = -5, tau_0 = 50, Bi_0 = 1.e6),
    data = x
  )
}

tmp_data <- 
  delta %>% 
  filter(id_sensor == "a")
tmp_model <- newton_cooling(tmp_data)
summary(tmp_model)


delta_nested <- 
  delta %>%
  nest(-id_sensor) %>%
  print()

model_nested <-
  delta_nested %>%
  mutate(model = map(data, newton_cooling)) %>%
  print()

predict_nested <-
  model_nested %>%
  mutate(pred = map2(model, data, predict)) %>%
  print()

predict_unnested <- 
  predict_nested %>%
  unnest(data, pred) %>% 
  print()

predict_tall <- 
  predict_unnested %>%
  rename(modeled = pred, measured = delta_temperature) %>%
  gather("type", "delta_temperature", modeled, measured) %>%
  print()

predict_tall %>%
  ggplot(aes(x = delta_time, y = delta_temperature)) +
  geom_line(aes(color = id_sensor, linetype = type))

list_model <-
  list(
    newton_cooling = newton_cooling,
    semi_infinite_simple = semi_infinite_simple,
    semi_infinite_convection = semi_infinite_convection
  )

fn_model <- function(.model, df){
  # safer to avoid non-standard evaluation
  # df %>% mutate(model = map(data, .model)) 
  
  df$model <- map(df$data, possibly(.model, NULL))
  df
}

model_nested_new <-
  list_model %>%
  map_df(fn_model, delta_nested, .id = "id_model") %>%
  mutate(is_null = map_lgl(model, is.null)) %>%
  filter(!is_null) %>%
  select(-is_null) %>%
  print()

predict_nested <- 
  model_nested_new %>%
  mutate(pred = map2(model, data, predict)) %>%
  print()

predict_tall <-
  predict_nested %>%
  unnest(data, pred) %>% 
  rename(modeled = pred, measured = delta_temperature) %>%
  gather("type", "delta_temperature", modeled, measured) %>%
  print()

predict_tall %>%
  ggplot(aes(x = delta_time, y = delta_temperature)) +
  geom_line(aes(color = id_sensor, linetype = type)) +
  facet_grid(id_model ~ .)

resid <-
  model_nested_new %>%
  mutate(resid = map(model, resid)) %>%
  unnest(data, resid) %>%
  print()

resid %>%
  ggplot(aes(x = delta_time, y = resid)) +
  geom_line(aes(color = id_sensor)) +
  facet_grid(id_model ~ .)

model_parameters <- 
  model_nested_new %>%
  select(id_model, id_sensor, model) %>%
  mutate(tidy = map(model, tidy)) %>%
  select(-model) %>%
  unnest() %>%
  print()

model_summary <-
  model_parameters %>%
  select(id_model, id_sensor, term, estimate) %>%
  spread(key = "term", value = "estimate") %>%
  print()
