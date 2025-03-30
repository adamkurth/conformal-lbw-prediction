# boot.R
rm(list=ls())
library(tidymodels,quietly = TRUE)
library(tidyverse,quietly = TRUE)
library(dbplyr,quietly = TRUE)
library(ggplot2,quietly = TRUE)

historical_tuition <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-03-10/historical_tuition.csv')
tuition_df <- historical_tuition %>% 
  pivot_wider(names_from = type,
              values_from = tuition_cost
  ) %>%
  na.omit() %>% 
  janitor::clean_names()

# plot 
tuition_df %>% 
  ggplot(aes(public, private)) + 
  geom_point() +
  scale_y_continuous(labels=scales::dollar_format()) + 
  scale_x_continuous(labels=scales::dollar_format()) +
  ggtitle("Private vs Public Tuition Costs") +
  labs(x = "Public Tuition Cost",
       y = "Private Tuition Cost") +
  theme_linedraw() + 
  theme(axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 20, face = "bold"))

tuition_fit <- lm(formula = private ~ 0 + public, data=tuition_df)
tidy(tuition_fit)
set.seed(123)

tuition_boot <- bootstraps(tuition_df, times=1e3, apparent = TRUE)
tuition_models <- tuition_boot %>% 
  mutate(model = map(splits, ~lm(private ~ 0 + public, data = .)), 
         coef_inf = map(model, tidy))

tuition_coefs <- tuition_models %>% 
  unnest(coef_inf)

tuition_coefs %>% 
  ggplot(aes(estimate))+
  geom_histogram(alpha = .7)+
  ggtitle("Distribution of Estimated Slope")+
  xlab("Estimate")+
  ylab("Frequency")+
  theme_linedraw()+
  theme(axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size = 20, face = "bold"))
