library(shiny)
library(bslib)
library(tidyverse)
library(stringr)
library(plotrix)
library(pander)
library(readxl)


data1 = readRDS("scenario_point_estimates.rds") %>% #mutate(sample="point estimate") %>% 
  pivot_longer(cols=`<6m`:`All ages`,names_to="age_grp",values_to="hosp") %>%
  mutate(scenario = factor(scenario),
         Age=factor(age_grp,levels=c("All ages","<6m",">6m","1-4yr","5-59yrs","60+yrs"),
                    labels=c("All ages","<6m",">6m","1-4yrs","5-59yrs","60+yrs"))) %>% 
  mutate(lower_point = NA, upper_point=NA)

data2 = readRDS("scenarios_with_PI.rds")%>% 
  pivot_longer(cols=`<6m`:`All ages`,names_to="age_grp",values_to="hosp") %>%
  mutate(scenario = factor(scenario),
         Age=factor(age_grp,levels=c("All ages","<6m",">6m","1-4yr","5-59yrs","60+yrs"),
                    labels=c("All ages","<6m",">6m","1-4yrs","5-59yrs","60+yrs"))) %>% 
  group_by(Age, scenario, date) %>% 
  summarize(mean = mean(hosp),
            sderr=std.error(hosp),
            lower_50 = quantile(hosp, probs=0.25),
            upper_50 = quantile(hosp, probs=0.75),
            lower_95 = quantile(hosp, probs=0.025),
            upper_95 = quantile(hosp, probs=0.975)) %>% 
  ungroup() %>% 
  mutate(date=as.Date(date)) %>% 
  left_join(data1, by=c("Age","scenario","date")) %>% 
  pivot_longer(cols=c(lower_50,upper_50,lower_95,upper_95,lower_point,upper_point),
               names_to="PI",values_to="value") %>% 
  mutate(PI_level = case_when(grepl("50",PI)~"50%",
                              grepl("95",PI)~"95%",
                              grepl("point",PI)~"Point Estimate"),
         PI = substr(PI,1,5)) 


data3 = readRDS("scenarios_with_PI.rds")%>% 
  pivot_longer(cols=`<6m`:`All ages`,names_to="age_grp",values_to="hosp") %>%
  mutate(scenario = factor(scenario),
         Age=factor(age_grp,levels=c("All ages","<6m",">6m","1-4yr","5-59yrs","60+yrs"),
                    labels=c("All ages","<6m",">6m","1-4yrs","5-59yrs","60+yrs"))) %>% 
  group_by(Age, scenario,sample) %>% 
  summarize(total = sum(hosp)) %>% 
  ungroup() %>% 
  group_by(Age, scenario) %>% 
  summarize(mean = mean(total),
            sderr=std.error(total)) 

counter =data3 %>% filter(scenario=="Counterfactual") %>% select(Age,"Counterfactual_mean"=mean,"Counterfactual_se"=sderr)

diffs = data3 %>% 
  left_join(counter, by=c("Age")) %>% 
  mutate(diff = (1-(mean/Counterfactual_mean))*100,
         err = sqrt(Counterfactual_se^2+Counterfactual_se^2),
         lower_50 = diff-0.67*err,
         upper_50 = diff+0.67*err,
         lower_95 = diff-1.96*err,
         upper_95 = diff+1.96*err,
         lower_point=NA,upper_point=NA) %>% 
  pivot_longer(cols=c(lower_50,upper_50,lower_95,upper_95,lower_point,upper_point),
               names_to="PI",values_to="value") %>% 
  mutate(value = ifelse(scenario=="Counterfactual",0,value),
         PI_level = case_when(grepl("50",PI)~"50%",
                              grepl("95",PI)~"95%",
                              grepl("point",PI)~"Point Estimate"),
         PI = substr(PI,1,5)) 
scenarios = read_excel("scenarios.xlsx")
effectiveness_parameters = read_excel("effectiveness_parameters.xlsx")

curves = readRDS("coverage_curves_2023_24.rds") %>% 
  mutate(mon_opt = monoclonal_cum*40,
         mon_pes = monoclonal_cum*15,
         mat_opt = maternal_cum*25,
         mat_pes = maternal_cum*10,
         sen_opt = senior_cum*30,
         sen_pes = senior_cum*15) %>% 
  filter(date>='2023-08-01') %>% 
  pivot_longer(cols=mon_opt:sen_pes) %>% 
  mutate(intervention = case_when(grepl("mon",name)~"Monoclonal Antibodies",
                                  grepl("mat",name)~"Maternal Vaccination",
                                  grepl("sen",name)~"Senior Vaccination"),
         scenario = case_when(grepl("opt",name)~"Optimistic",
                              grepl("pes",name)~"Pessimistic"))



# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "RSV Scenario Projections",
  sidebar = sidebar(
    checkboxGroupInput(
      inputId ="checkGroup",label=tags$h3("Select Scenarios"),
      choices = levels(data1$scenario),
      selected = "Counterfactual"
    ),
    radioButtons(
      inputId="radio",label=tags$h3("Projection Interval"),
      choices = list("Point Estimate","50%","95%")
    )
    
  ),
  
  navset_card_underline(
    nav_panel("Scenario Descriptions", tableOutput("scenarios")),
    nav_panel("Immunization Effectiveness", tableOutput("ve")),
    
    nav_panel("Coverage Curves", plotOutput("cov_curve",hover="plot_hover"),
              verbatimTextOutput("info")),
    nav_panel("Time Series", plotOutput("timeseries")),
    
    nav_panel("Summary Table", card(strong("Cumulative Seasonal RSV Hospitalizations by Scenario",size=20),tableOutput("table"))),
    
    nav_panel("% Difference with Counterfactual", plotOutput("barplot"))
    
    
    
    
  ) )
# card(
# plotOutput(outputId = "timeseries"),
#  ),
# card(tableOutput(outputId = "table")
#   ),
#card(plotOutput(outputId = "barplot")
#)
#)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  dat_new = reactive({filter(data2, scenario %in% input$checkGroup & PI_level %in% input$radio) %>% 
      pivot_wider(names_from=PI, values_from=value)})
  diffs_new = reactive({filter(diffs, scenario %in% input$checkGroup & PI_level %in% input$radio) %>% 
      pivot_wider(names_from=PI, values_from=value)})
  
  
  output$timeseries <- renderPlot({
    req(input$checkGroup)
    ggplot(data = dat_new())+
      theme_bw()+
      geom_line(aes(x = date, y = hosp,group=scenario,color=scenario)) +
      geom_ribbon(aes(x = date, ymin=lower,ymax=upper,group=scenario,fill=scenario),alpha=0.5) +
      facet_wrap(~Age,ncol=3,scales="free")+
      theme(legend.position="top",
            axis.text=element_text(size=15),
            strip.text = element_text(size=15),
            axis.title = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=20),
            title=element_text(size=20))+
      guides(fill=FALSE)+
      labs(x=NULL, y="Hospitalizations",color="Scenario",title="Weekly RSV hospitalizations by scenario")
  })
  
  
  
  output$table <- renderTable(
    {
      req(input$checkGroup)
      dat_new() %>% 
        group_by(Age, scenario) %>%
        summarise(Total = round(sum(hosp)),
                  lower = round(sum(lower)),
                  upper = round(sum(upper))) %>%
        ungroup() %>% 
        mutate(combine =paste0(Total," (",lower,"-",upper,")"),
               combine = ifelse(grepl("NA",combine),str_sub(combine,end=-8),combine)) %>% 
        select(-Total,-lower,-upper) %>% 
        pivot_wider(names_from=scenario,values_from=combine)
    },
    striped = TRUE,
    spacing = "s",
    align = "c",
    digits = 0,
    width = "auto"
  ) 
  
  output$barplot <- renderPlot({
    req(input$checkGroup)
    ggplot(data = diffs_new())+
      theme_bw()+
      geom_errorbar(aes(y=reorder(scenario,desc(scenario)), xmin=lower*-1, xmax=upper*-1,color=scenario),width=0,size=5,alpha=0.5)+
      geom_point(aes(y=reorder(scenario,desc(scenario)), x=diff*-1,color=scenario),position="dodge",size=8)+
      geom_vline(xintercept=0)+
      facet_wrap(~Age,ncol=3,scales="free")+
      geom_text(aes(x=diff*-1, y=reorder(scenario,desc(scenario)),label=paste0("-",round(diff))),size=4)+
      labs(y=NULL, x=NULL,color="Scenario",title="% Difference between Scenario and Counterfactual")+
      theme(#legend.position="top",
        axis.text.y=element_text(size=15),
        axis.text.x =element_text(size=15),
        strip.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        title=element_text(size=20))
  })
  
  
  
  output$scenarios <- renderTable(
    {
      scenarios
    },
    striped = TRUE,
    spacing = "s",
    align = "c",
    digits = 0,
    width = "auto",
    caption = 'See "Immunization Effectiveness" and "Coverage Curves" panels for scenario parameter values'
  ) 
  
  
  output$ve <- renderTable(
    {
      effectiveness_parameters
    },
    striped = TRUE,
    spacing = "s",
    align = "c",
    digits = 0,
    width = "auto",
    caption = 'Scenario values based on data from clinical trials'
  ) 
  
  output$cov_curve <- renderPlot({
    #req(input$checkGroup)
    ggplot(curves)+
      theme_bw()+
      geom_line(aes(x=date, y=value,color=scenario),size=2)+
      facet_grid(cols=vars(intervention))+
      labs(x=NULL, y="Cumulative % Immunized",title=paste0("Cumulative RSV Immunization Coverage Scenarios, ",min(year(curves$date)),"-",max(year(curves$date))))+
      scale_color_manual(name="Scenario",values=c("navy","maroon3"))+
      theme(axis.text=element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=20),
            strip.text.x=element_text(size=15),
            axis.title.y=element_text(size=20),
            title = element_text(size=20))
  })
  
  
  output$info <- renderText({
    req(input$checkGroup)
    paste0("Date: ",as.Date(input$plot_hover$x)," Coverage: ",round(input$plot_hover$y,1),"%")
  })
  
  
}




# Run the application 
shinyApp(ui = ui, server = server)


