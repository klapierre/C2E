### Plot synchrony and asynchrony simulations

#file.choose()
library(tidyverse)
library(ggthemes)
asynchrony_simulation_data <- read.csv("C:\\Users\\kwilcox4\\Dropbox\\Shared working groups\\C2E\\GCD asynchrony\\midi files\\asynchrony_sim_data.csv") %>%
  gather(key=Plot, value=Cover, -time)
synchrony_simulation_data <- read.csv("C:\\Users\\kwilcox4\\Dropbox\\Shared working groups\\C2E\\GCD asynchrony\\midi files\\synchrony_sim_data.csv") %>%
  gather(key=Plot, value=Cover, -time)

ggplot(asynchrony_simulation_data, aes(x=time, y=Cover, col=Plot)) +
  geom_line() +
  theme_few()

ggplot(synchrony_simulation_data, aes(x=time, y=Cover, col=Plot)) +
  geom_line() +
  theme_few()



