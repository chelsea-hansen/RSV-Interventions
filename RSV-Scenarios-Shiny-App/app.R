library(shiny)
library(bslib)
library(tidyverse)
library(stringr)
library(plotrix)
library(pander)
library(readxl)



data2 = readRDS("results_rates_2024-25.rds")%>% 
  mutate(scenario = factor(scenario, levels=c("scenarioE","scenarioA","scenarioB","scenarioC","scenarioD"),
                           labels=c("Counterfactual","Scenario A","Scenario B","Scenario C","Scenario D")),
         Age=factor(Age,levels=c("All","<6m","6-11m","1-4yrs","5-59yrs","60+yrs"),
                    labels=c("All ages","<6m","6-11m","1-4yrs","5-59yrs","60+yrs")),
         lower_median=NA,upper_median=NA,media_95=median,media_50=median) %>% 
  pivot_longer(cols=c(median,lower_50,upper_50,lower_95,upper_95,lower_median,upper_median,media_95,media_50),
               names_to="PI",values_to="value") %>% 
  mutate(PI_level = case_when(grepl("50",PI)~"50%",
                              grepl("95",PI)~"95%",
                              grepl("med",PI)~"Median"),
         PI = substr(PI,1,5)) 


diffs= readRDS("percentage_reduction_2024-25.rds")%>% 
  mutate(scenario = factor(scenario,levels=c("scenarioE","scenarioA","scenarioB","scenarioC","scenarioD"),
                           labels=c("Counterfactual","Scenario A","Scenario B","Scenario C","Scenario D")),
         Age=factor(Age,levels=c("All","<6m","6-11m","1-4yrs","5-59yrs","60+yrs"),
                    labels=c("All ages","<6m","6-11m","1-4yrs","5-59yrs","60+yrs")),
         lower_median=NA,upper_median=NA,media_95=median,media_50=median) %>% 
  pivot_longer(cols=c(median,lower_50,upper_50,lower_95,upper_95,lower_median,upper_median,media_95,media_50),
               names_to="PI",values_to="value") %>% 
  mutate(PI_level = case_when(grepl("50",PI)~"50%",
                              grepl("95",PI)~"95%",
                              grepl("med",PI)~"Median"),
         PI = substr(PI,1,5)) 
  


scenarios = read_excel("scenarios.xlsx")
effectiveness_parameters = read_excel("effectiveness_parameters.xlsx")

curves = readRDS("coverage for figures.rds") 


# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "RSV Scenario Projections for the 2024-25 Season in King County, WA",
  sidebar = sidebar(
    checkboxGroupInput(
      inputId ="checkGroup",label=tags$h3("Select Scenarios"),
      choices = levels(data2$scenario),
      selected = "Counterfactual"
    ),
    radioButtons(
      inputId="radio",label=tags$h3("Projection Interval"),
      choices = list("Median","50%","95%")
    )
    
  ),
  
  navset_card_underline(
    nav_panel("Scenario Descriptions", tableOutput("scenarios")),
    nav_panel("Immunization Effectiveness", tableOutput("ve")),
    
    nav_panel("Coverage Curves", plotOutput("cov_curve",hover="plot_hover"),
              verbatimTextOutput("info")),
    nav_panel("Time Series", plotOutput("timeseries")),
    
    nav_panel("Summary Table", card(strong("Cumulative Seasonal RSV Hosp Rate per 100,000 by Scenario",size=20),tableOutput("table"))),
    
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
  
  check = diffs %>% filter(scenario %in% c("scenarioA"),PI_level %in% c("95%")) %>% 
    pivot_wider(names_from=PI, values_from=value)
  
  output$timeseries <- renderPlot({
    req(input$checkGroup)
    ggplot(data = dat_new())+
      theme_bw()+
      geom_line(aes(x = date, y = media,group=scenario,color=scenario)) +
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
      labs(x=NULL, y="RSV Hosp Rate",color="Scenario",title="Weekly RSV hosp rate per 100,000 by scenario")
  })
  
  
  
  output$table <- renderTable(
    {
      req(input$checkGroup)
      dat_new() %>% 
        group_by(Age, scenario) %>%
        summarise(Total = round(sum(media),1),
                  lower = round(sum(lower),1),
                  upper = round(sum(upper),1)) %>%
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
      geom_point(aes(y=reorder(scenario,desc(scenario)), x=media*-1,color=scenario),position="dodge",size=8)+
      geom_vline(xintercept=0)+
      facet_wrap(~Age,ncol=3,scales="free")+
      geom_text(aes(x=media*-1, y=reorder(scenario,desc(scenario)),label=paste0("-",round(media))),size=4)+
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
      geom_line(aes(x=date, y=doses,color=scenario),size=2)+
      facet_wrap(~population,ncol=3,scales="free")+
      labs(x=NULL, y="Cumulative Doses Administered",title=paste0("Cumulative RSV Immunization Coverage Scenarios, ",min(year(curves$date)),"-",max(year(curves$date))))+
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
    paste0("Date: ",as.Date(input$plot_hover$x)," Coverage: ",round(input$plot_hover$y),"Doses")
  })
  
  
}




# Run the application 
shinyApp(ui = ui, server = server)


