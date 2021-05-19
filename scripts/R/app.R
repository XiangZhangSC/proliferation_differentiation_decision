rm(list=ls(all=t))

library(deSolve)
library(shiny)
library(ggplot2)
source("./scripts/R/rate_equations.R")
source("./scripts/R/forcing_functions.R")
source("./scripts/R/plot_functions.R")
source("./scripts/R/proliferation-differentiation-model.R")


ui <- fluidPage(
  fluidRow(
    column(6, 
           wellPanel(
             fluidRow(
               column(4, checkboxInput("display_dat", "Display gene expression data?", value = TRUE)), 
               column(4, radioButtons("scenario", "Scenario",  
                                      choices = list("Molly", "Mitogen ON and mls-2 OFF", "Mitogen ON and mls-2 ON"))),
               column(4, downloadButton("download", "Download the plot"))
             )
           ),
           wellPanel(
             h5("HLH-1"), 
             fluidRow(
               column(3, sliderInput("tau_hlh1", "tau_hlh1", min = 0, max = 1, value = 0.5, step = 0.05)),
               column(3, sliderInput("k_mls2_hlh1", "k_mls2_hlh1", min = 0, max = 1, value = 0.3, step = 0.05)),
               column(3, sliderInput("k_hlh1_hlh1", "k_hlh1_hlh1", min = 0, max = 1, value = 0.1, step = 0.05)), 
               column(3, sliderInput("k_fos1_hlh1", "k_fos1_hlh1", min = 0, max = 1, value = 0.4, step = 0.05))
             ), 
             fluidRow(
               column(6, sliderInput("k_cye1_hlh1", "k_cye1_hlh1", min = 0, max = 1, value = 0.8, step = 0.05)), 
               column(6, sliderInput("k_cki1_hlh1", "k_cki1_hlh1", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("FOS-1"), 
             fluidRow(
               column(4, sliderInput("tau_fos1", "tau_fos1", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(4, sliderInput("k_lin1_fos1", "k_lin1_fos1", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(4, sliderInput("k_hlh1_fos1", "k_hlh1_fos1", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ), 
           wellPanel(
             h5("CYD-1"), 
             fluidRow(
               column(6, sliderInput("tau_cyd1", "tau_cyd1", min = 0, max = 1, value = 1, step = 0.05)), 
               column(6, sliderInput("k_fos1_cyd1", "k_fos1_cyd1", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("CYE-1"), 
             fluidRow(
               column(6, sliderInput("tau_cye1", "tau_cye1", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(6, sliderInput("k_e2f_cye1", "k_e2f_cye1", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("CKI-1"), 
             fluidRow(
               column(6, sliderInput("tau_cki1", "tau_cki1", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(6, sliderInput("k_hlh1_cki1", "k_hlh1_cki1", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("LIN-35"), 
             fluidRow(
               column(6, sliderInput("tau_lin35", "tau_lin35", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(6, sliderInput("k_hlh1_lin35", "k_hlh1_lin35", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("UNC-120"), 
             fluidRow(
               column(4, sliderInput("tau_unc120", "tau_unc120", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(4, sliderInput("k_hlh1_unc120", "k_hlh1_unc120", min = 0, max = 1, value = 0.8, step = 0.05)), 
               column(4, sliderInput("k_fos1_unc120", "k_fos1_unc120", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("Protein complex"), 
             fluidRow(
               column(4, sliderInput("tau_HLH1LIN35", "tau_HLH1LIN35", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(4, sliderInput("kd_HLH1LIN35", "kd_HLH1LIN35", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(4, sliderInput("tau_E2F", "tau_E2F", min = 0, max = 1, value = 0.5, step = 0.05))
             )
           ),
           wellPanel(
             h5("RNR-1"), 
             fluidRow(
               column(6, sliderInput("tau_rnr1", "tau_rnr1", min = 0, max = 1, value = 1, step = 0.05)), 
               column(6, sliderInput("k_e2f_rnr1", "k_e2f_rnr1", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ), 
           wellPanel(
             h5("UNC-15"), 
             fluidRow(
               column(3, sliderInput("tau_unc15", "tau_unc15", min = 0, max = 1, value = 1, step = 0.05)), 
               column(3, sliderInput("k_hlh1_unc15", "k_hlh1_unc15", min = 0, max = 1, value = 0.8, step = 0.05)), 
               column(3, sliderInput("k_unc120_unc15", "k_unc120_unc15", min = 0, max = 1, value = 0.8, step = 0.05)), 
               column(3, sliderInput("k_fos1_unc15", "k_fos1_unc15", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           )
    ), 
    column(6, 
           plotOutput("plot", height = "950px")
    )
  )
)
server <- function(input, output, session) {
  
  params_chosen <- reactive(c(tau_hlh1 = input$tau_hlh1, tau_fos1 = input$tau_fos1, tau_cyd1 = input$tau_cyd1, tau_cye1 = input$tau_cye1, 
                     tau_cki1 = input$tau_cki1, tau_lin35 = input$tau_lin35, tau_unc120 = input$tau_unc120, 
                     tau_HLH1LIN35 = input$tau_HLH1LIN35, 
                     tau_E2F = input$tau_E2F, 
                     tau_rnr1 = input$tau_rnr1, tau_unc15 = input$tau_unc15, 
                     k_mls2_hlh1 = input$k_mls2_hlh1, k_hlh1_hlh1 = input$k_hlh1_hlh1, k_fos1_hlh1 = input$k_fos1_hlh1, 
                     k_lin1_fos1 = input$k_lin1_fos1, k_hlh1_fos1 = input$k_hlh1_fos1, 
                     k_fos1_cyd1 = input$k_fos1_cyd1, k_e2f_cye1 = input$k_e2f_cye1, 
                     k_hlh1_cki1 = input$k_hlh1_cki1, k_hlh1_lin35 = input$k_hlh1_lin35, 
                     k_hlh1_unc120 = input$k_hlh1_unc120, k_fos1_unc120 = input$k_fos1_unc120, 
                     k_e2f_rnr1 = input$k_e2f_rnr1, 
                     k_hlh1_unc15 = input$k_hlh1_unc15, k_unc120_unc15 = input$k_unc120_unc15, k_fos1_unc15 = input$k_fos1_unc15, 
                     k_cye1_hlh1 = input$k_cye1_hlh1, k_cki1_hlh1 = input$k_cki1_hlh1, 
                     kd_HLH1LIN35 = input$kd_HLH1LIN35, k_mls2_cye1 = input$k_mls2_cye1))
  
  display_dat <- reactive(input$display_dat)
  scenario <- reactive(input$scenario)
  
  plot_sim <- reactive(chart_dynamics(solve_prodiff(params = params_chosen(), scenario = scenario()), display_dat()))
  
  output$plot <- renderPlot({
    print(plot_sim())
  })
  
  output$download <- downloadHandler(
    filename = function() {paste(input$scenario, '.png', sep = '')}, 
    content = function(file) {
      ggsave(file, plot_sim())
    }
  )
}

shinyApp(ui, server)
