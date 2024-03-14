rm(list=ls())
library(tidyverse)
library(cowplot)
library(readxl)
library(gt)
library(ggplotify)
library(gridExtra)
library(grid)
library(flextable)
library(plotrix)
library(zoo)

'%not_in%' = Negate("%in%")


dat = readRDS("3. Interventions/scenario_point_estimates.rds") %>% 
  mutate(scenario=factor(scenario, levels=c("Scenario A","Scenario B", "Scenario C", "Scenario D",
                                            "Scenario E","Scenario F", "Scenario G", "Scenario H",
                                            "Counterfactual")),
         date = as.Date(date),
         `<1yr`=`<6m`+`>6m`) %>% 
  pivot_longer(cols=c(`<6m`:`All ages`,`<1yr`),names_to = "age",values_to="hospitalizations")



plot1 = ggplot(data=dat %>% filter(age %in% c("<1yr","60+yrs","All ages")))+
  theme_bw()+
  geom_line(aes(x=date, y=hospitalizations,color=scenario))+
  facet_wrap(~age, ncol=3)+
  labs(x=NULL,y="Weekly Hospitalizations",title="Weekly Projected Hospitalizations")



dat2 = readRDS("3. Interventions/scenarios_with_PI.rds") %>% 
  mutate(scenario=factor(scenario, levels=c("Scenario A","Scenario B", "Scenario C", "Scenario D",
                                            "Scenario E","Scenario F", "Scenario G", "Scenario H",
                                            "Counterfactual"))) %>% 
  mutate(`<1yr`=`<6m`+`>6m`) %>% 
  pivot_longer(cols=c(`<1yr`,`60+yrs`,`All ages`),names_to = "age",values_to="hospitalizations") %>% 
  group_by(scenario, age,sample) %>% 
  summarize(total = sum(hospitalizations)) %>% 
  ungroup() %>% 
  group_by(scenario, age) %>% 
  summarize(mean = mean(total),
            sderr=std.error(total),
            median = median(total),
            lower = quantile(total, probs=0.025),
            upper = quantile(total, probs=0.975)) %>% 
  mutate(var = ((upper-lower)/3.92)^2,
         mean.upper=mean+1.96*sderr,
         mean.lower=mean-1.96*sderr)


counter =dat2 %>%ungroup() %>%  filter(scenario=="Counterfactual") %>% select("Counterfactual"=median, age)

diffs = dat2 %>% 
  full_join(counter, by="age") %>% 
  mutate(percent_diff = round((1-(median/Counterfactual))*100,1))


plot2 = ggplot(data=diffs)+
  theme_bw()+
  geom_bar(aes(y=scenario, x=median,fill=scenario),stat="identity",position="dodge")+
  geom_errorbar(aes(y=scenario, xmin=lower, xmax=upper),width=.1)+
  facet_grid(cols=vars(age),scales="free")+
  geom_text(aes(x=50, y=scenario,label=paste0("-",round(percent_diff),"%")))+
  labs(y=NULL, x="Total Hospitalizations",title="% Difference with Counterfactual")

all_plots = plot_grid(plot1, plot2, nrow=2)
all_plots
ggsave(plot=all_plots, "intervention_plots.png",height=8,width=13,units="in")
