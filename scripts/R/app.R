library(deSolve)
library(shiny)
library(ggplot2)
source("./scripts/R/rate_equations.R")
source("./scripts/R/forcing_functions.R")
source("./scripts/R/plot_functions.R")
source("./scripts/R/proliferation-differentiation-model.R")


ui <- fluidPage(
  titlePanel("Proliferation-differentiation decision in M lineage"), 
  fluidRow(
    column(6, 
           wellPanel(
             h5("HLH-1"), 
             fluidRow(
               column(3, sliderInput("tau_hlh1", "tau_hlh1", min = 0, max = 1, value = 0.5, step = 0.05)),
               column(3, sliderInput("k_mls2_hlh1", "k_mls2_hlh1", min = 0, max = 1, value = 0.3, step = 0.05)),
               column(3, sliderInput("k_hlh1_hlh1", "k_hlh1_hlh1", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(3, sliderInput("k_fos1_hlh1", "k_fos1_hlh1", min = 0, max = 1, value = 0.3, step = 0.05))
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
               column(6, sliderInput("tau_cyd1", "tau_cyd1", min = 0, max = 1, value = 0.5, step = 0.05)), 
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
             h5("E2F"), 
             fluidRow(
               column(12, sliderInput("tau_E2F", "tau_E2F", min = 0, max = 1, value = 0.5, step = 0.05))
             )
           ),
           wellPanel(
             h5("MEF-2"), 
             fluidRow(
               column(3, sliderInput("tau_mef2", "tau_mef2", min = 0, max = 1, value = 0.5, step = 0.05)), 
               column(3, sliderInput("k_hlh1_mef2", "k_hlh1_mef2", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(3, sliderInput("k_mef2_mef2", "k_mef2_mef2", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(3, sliderInput("k_fos1_mef2", "k_fos1_mef2", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ),
           wellPanel(
             h5("Proliferation"), 
             fluidRow(
               column(12, sliderInput("k_e2f_proliferation", "k_e2f_proliferation", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ), 
           wellPanel(
             h5("Differentiation"), 
             fluidRow(
               column(4, sliderInput("k_hlh1_differentiation", "k_hlh1_differentiation", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(4, sliderInput("k_mef2_differentiation", "k_mef2_differentiation", min = 0, max = 1, value = 0.3, step = 0.05)), 
               column(4, sliderInput("k_fos1_differentiation", "k_fos1_differentiation", min = 0, max = 1, value = 0.3, step = 0.05))
             )
           ), 
           wellPanel(
             fluidRow(
               column(6, checkboxInput("display_dat", "Display gene expression data?", value = TRUE))
             )
           )
    ), 
    column(6, 
           plotOutput("plot", height = "950px")
    )
  )
)
server <- function(input, output, session) {
  
  output$plot <- renderPlot({
    params_chosen <- c(tau_hlh1 = input$tau_hlh1, tau_fos1 = input$tau_fos1, tau_cyd1 = input$tau_cyd1, tau_cye1 = input$tau_cye1, 
                       tau_cki1 = input$tau_cki1, tau_lin35 = input$tau_lin35, tau_mef2 = input$tau_mef2, tau_E2F = input$tau_E2F, 
                       k_mls2_hlh1 = input$k_mls2_hlh1, k_hlh1_hlh1 = input$k_hlh1_hlh1, k_fos1_hlh1 = input$k_fos1_hlh1, 
                       k_lin1_fos1 = input$k_lin1_fos1, k_hlh1_fos1 = input$k_hlh1_fos1, 
                       k_fos1_cyd1 = input$k_fos1_cyd1, k_e2f_cye1 = input$k_e2f_cye1, 
                       k_hlh1_cki1 = input$k_hlh1_cki1, k_hlh1_lin35 = input$k_hlh1_lin35, 
                       k_hlh1_mef2 = input$k_hlh1_mef2, k_mef2_mef2 = input$k_mef2_mef2, k_fos1_mef2 = input$k_fos1_mef2, 
                       k_e2f_proliferation = input$k_e2f_proliferation, 
                       k_hlh1_differentiation = input$k_hlh1_differentiation, k_mef2_differentiation = input$k_mef2_differentiation, k_fos1_differentiation = input$k_fos1_differentiation)
    display_dat = input$display_dat
    chart_dynamics(solve_prodiff(params_chosen), display_dat)
  })
}

shinyApp(ui, server)