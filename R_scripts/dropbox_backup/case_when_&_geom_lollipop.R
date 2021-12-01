library(ggalt)
library(hrbrthemes)
library(tidyverse)
library(magrittr)

os <- read_lines("C:/Users/Filip/Dropbox/Praksa bioinfo/functions_R_UNIX/os.txt", na = "NA")
str(table(os, useNA = "always"))
sort(unique(os))

case_when(
  grepl("Windows", os) ~ "Windows-ish",
  grepl("Red Hat|CentOS|Fedora", os) ~ "Fedora-ish",
  grepl("Ubuntu|Debian", os) ~ "Debian-ish",
  grepl("CoreOS|Amazon", os) ~ "Amazon-ish",
  is.na(os) ~ "Unknown",
  TRUE ~ "Other") %>%
  table() %>%
  as_data_frame() %>%
  set_names(c("os", "Node Count")) %>%
  arrange(`Node Count`) %>%
  mutate(os = factor(os, os)) %>%
  ggplot(aes(`Node Count`, os)) +
  geom_lollipop(horizontal = TRUE, size=1.5, color="#54278f") +
  scale_x_comma(limits=c(0,300)) +
  labs(y=NULL, title="OS Types") +
  theme_ipsum_rc(grid="X")
