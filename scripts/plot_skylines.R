#Bethany Allen
#Code to plot exponential coalescent and fossilised-birth-death skylines

#Install the `tidyverse` if you haven't already
install.packages("tidyverse")

library(tidyverse)

###Plotting the exponential coalescent skyline

# Navigate to Session > Set Working Directory > Choose Directory (on RStudio)
# or change file name to the full path to the log file
#(Use "dinosaur_coal_final.log" if you used our pre-cooked XML)
coal_file <- "dinosaur_coal.log"

#Read in coalescent log and trim 10% burn-in
coalescent <- read.table(coal_file, header = T) %>% slice_tail(prop = 0.9)

#Pivot the table to stack the rate estimates into a single column
coalescent <- pivot_longer(coalescent, c(growthRate1, growthRate2, growthRate3,
                                   growthRate4),
                         names_to = "time_bin", names_prefix = "growthRate")

#Summarise the rates in the log
coalescent_summary <- coalescent %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add the interval names
coalescent_summary$interval <- c("Late Cretaceous", "Early Cretaceous",
                                 "Jurassic", "Triassic")

#Ensure that the time intervals plot in the correct order
coalescent_summary$interval <- factor(coalescent_summary$interval,
                                      levels = c("Triassic", "Jurassic",
                                                 "Early Cretaceous",
                                                 "Late Cretaceous"))
#Plot diversification skyline as error bars
ggplot(data = coalescent_summary, aes(x = interval, y = median, ymin = lowCI,
                             ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Diversification rate") +
  theme_classic(base_size = 17)

#Plot diversification skyline as a ribbon plot
ages <- seq.int(252, 66)
interval <- c(rep("Triassic", ((252 - 202) + 1)),
              rep("Jurassic", ((201 - 146) + 1)),
              rep("Early Cretaceous", ((145 - 101) + 1)),
              rep("Late Cretaceous", ((100 - 66) + 1)))
age_table <- as.data.frame(cbind(ages, interval))
to_plot <- left_join(coalescent_summary, age_table, by = "interval")
to_plot$ages <- as.numeric(to_plot$ages)

ggplot(to_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  xlab("Age (Ma)") + ylab("Diversification rate") +
  theme_classic(base_size = 17)

#Extract estimated diversity
pop_data <- pull(coalescent, "ePopSize")

#Summarise the diversity estimates in the log
pop_data <- as.data.frame(rbind(c(median(pop_data), quantile(pop_data, 0.025),
                                  quantile(pop_data, 0.975))))
colnames(pop_data) <- c("median", "lowCI", "highCI")
print(pop_data)


###Plotting the fossilised-birth-death skylines

# Navigate to Session > Set Working Directory > Choose Directory (on RStudio)
# or change file name to the full path to the log file
#(Use "dinosaur_BDSKY_final.log" if you used our pre-cooked XML)
fbd_file <- "dinosaur_BDSKY.log"

#Read in coalescent log and trim 10% burn-in
fbd <- read.table(fbd_file, header = T) %>% slice_tail(prop = 0.9)

#Calculate diversification and turnover
birth_rates <- select(fbd, starts_with("birthRate"))
death_rates <- select(fbd, starts_with("deathRate"))

div_rates <- birth_rates - death_rates
colnames(div_rates) <- paste0("divRate.",
                              seq(1:ncol(div_rates)))

TO_rates <- birth_rates / death_rates
colnames(TO_rates) <- paste0("TORate.",
                             seq(1:ncol(TO_rates)))

#Pivot the table to stack the rate estimates into a single column
div_data <- pivot_longer(div_rates,
                         c(divRate.1, divRate.2, divRate.3, divRate.4),
                         names_to = "time_bin", names_prefix = "divRate.")

#Summarise the diversification estimates
div_data <- div_data %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add interval names
div_data$interval <- c("Triassic", "Jurassic", "Early Cretaceous",
                       "Late Cretaceous")

#Pivot the table to stack the rate estimates into a single column
turn_data <- pivot_longer(TO_rates, c(TORate.1, TORate.2, TORate.3, TORate.4),
                          names_to = "time_bin", names_prefix = "TORate.")

#Summarise the turnover estimates
turn_data <- turn_data %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add interval names
turn_data$interval <- c("Triassic", "Jurassic", "Early Cretaceous",
                        "Late Cretaceous")

#Pivot the table to stack the rate estimates into a single column
samp_data <- select(fbd, starts_with("sampling"))
samp_data <- pivot_longer(samp_data, c(samplingBDS.1, samplingBDS.2,
                                       samplingBDS.3, samplingBDS.4),
                          names_to = "time_bin", names_prefix = "samplingBDS.")

#Summarise log
samp_data <- samp_data %>% group_by(time_bin) %>%
  summarise(median = median(value),
            lowCI = quantile(value, 0.025),
            highCI = quantile(value, 0.975))

#Add interval names
samp_data$interval <- c("Triassic", "Jurassic", "Early Cretaceous",
                        "Late Cretaceous")

#Ensure that the time intervals plot in the correct order
div_data$interval <- factor(div_data$interval,
                            levels = c("Triassic", "Jurassic",
                                       "Early Cretaceous", "Late Cretaceous"))

turn_data$interval <- factor(turn_data$interval,
                             levels = c("Triassic", "Jurassic",
                                        "Early Cretaceous", "Late Cretaceous"))

samp_data$interval <- factor(samp_data$interval,
                             levels = c("Triassic", "Jurassic",
                                        "Early Cretaceous", "Late Cretaceous"))

#Plot skylines as error bars
ggplot(data = div_data, aes(x = interval, y = median, ymin = lowCI,
                                      ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Diversification rate") +
  theme_classic(base_size = 17)

ggplot(data = turn_data, aes(x = interval, y = median, ymin = lowCI,
                            ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 1), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Turnover rate") +
  theme_classic(base_size = 17)

ggplot(data = samp_data, aes(x = interval, y = median, ymin = lowCI,
                            ymax = highCI)) +
  geom_point(size = 1.5) +
  geom_errorbar(size = 1, width = 0.5) +
  scale_colour_manual(values = c("black")) +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Interval", y = "Sampling rate") +
  theme_classic(base_size = 17)

#Plot skylines as a ribbon plot
ages <- seq.int(252, 66)
interval <- c(rep("Triassic", ((252 - 202) + 1)),
              rep("Jurassic", ((201 - 146) + 1)),
              rep("Early Cretaceous", ((145 - 101) + 1)),
              rep("Late Cretaceous", ((100 - 66) + 1)))
age_table <- as.data.frame(cbind(ages, interval))

div_plot <- left_join(div_data, age_table, by = "interval")
div_plot$ages <- as.numeric(div_plot$ages)
ggplot(div_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  xlab("Age (Ma)") + ylab("Diversification rate") +
  theme_classic(base_size = 17)

turn_plot <- left_join(turn_data, age_table, by = "interval")
turn_plot$ages <- as.numeric(turn_plot$ages)
ggplot(turn_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 1), colour = "black") +
  xlab("Age (Ma)") + ylab("Turnover rate") +
  theme_classic(base_size = 17)

samp_plot <- left_join(samp_data, age_table, by = "interval")
samp_plot$ages <- as.numeric(samp_plot$ages)
ggplot(samp_plot) +
  geom_ribbon(aes(x = ages, ymin = lowCI, ymax = highCI), alpha = 0.5) +
  geom_line(aes(x = ages, y = median)) +
  scale_x_reverse() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  xlab("Age (Ma)") + ylab("Sampling rate") +
  theme_classic(base_size = 17)

#Extract origin data
origin_data <- pull(fbd, "origin")

#Summarise the origin estimates in the log
origin_data <- as.data.frame(rbind(c((median(origin_data) + 66),
                                     (quantile(origin_data, 0.025) + 66),
                                     (quantile(origin_data, 0.975) + 66))))
colnames(origin_data) <- c("median", "lowCI", "highCI")
print(origin_data)
