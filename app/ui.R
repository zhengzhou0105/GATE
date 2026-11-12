library(shiny)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  
  titlePanel("Genotoxicity Adaptive Testing & Evaluation (GATE)"),
  
  # ==========================================
  # Page 1: Welcome Page
  # ==========================================
  div(id = "page1",
      h3("Welcome to the GATE app!"),
      h4("Version 0.2.1"),
      br(),
      p("This version includes five in vitro genotoxicity assays: Ames, MLA, HPRT, in vitro CA (CAvt), and in vitro MN (MNvt)."),
      
      wellPanel(
        h4("Included In Vitro Assays:"),
        tags$ul(
          tags$li(strong("Ames:"), " Bacterial reverse mutation test (OECD TG 471)"),
          tags$li(strong("MLA:"), " In vitro mammalian cell gene mutation test using the mouse lymphoma tk gene (OECD TG 490)"),
          tags$li(strong("HPRT:"), " In vitro mammalian cell gene mutation test using the Hprt and xprt genes (OECD TG 476)"),
          tags$li(strong("CAvt:"), " In vitro mammalian chromosome aberration test (OECD TG 473)"),
          tags$li(strong("MNvt:"), " In vitro mammalian cell micronucleus test (OECD TG 487)")
        )
      ),
      br(),
      fluidRow(
        column(12, align = "center",
               actionButton("start_btn", "Start", class = "btn-primary btn-lg"),
               br(), br(),
               p("Results will not be saved.")
        )
      )
  ),
  
  # ==========================================
  # Page 2: First User Action & Thresholds
  # ==========================================
  hidden(
    div(id = "page2",
        h3("For the compound of interest:"),
        radioButtons("app_intent", label = NULL,
                     choices = list(
                       "a) No data has been collected from these five assays" = "a",
                       "b) Input existing data" = "b"
                     ),
                     width = "100%"),
        hr(),
        h4("Conclusive Decision Threshold"),
        radioButtons("keep_thresh_p2", 
                     "Default conclusive decision threshold: Positive if probability > 0.975 and Negative if probability < 0.025. Do you want to keep the default threshold?",
                     choices = c("Yes (Keep default)" = "yes", "No (Customize)" = "no"),
                     selected = "yes", inline = TRUE),
        hidden(
          div(id = "custom_thresh_p2_div",
              fluidRow(
                column(4, numericInput("thresh_upper_p2", "Positive Threshold (Upper):", value = 0.975, min = 0.50, max = 0.9999, step = 0.01)),
                column(4, numericInput("thresh_lower_p2", "Negative Threshold (Lower):", value = 0.025, min = 0.0001, max = 0.50, step = 0.01))
              )
          )
        ),
        br(),
        actionButton("next_btn", "Next", class = "btn-primary")
    )
  ),
  
  # ==========================================
  # Page 3a: Path A - Recommended Starting Assays
  # ==========================================
  hidden(
    div(id = "page3a",
        h3("Recommended Starting Assays"),
        p("Based on the results from our paper, the recommended starting in vitro assays are:"),
        tags$ul(
          tags$li(strong("Gene Mutation (GM):"), " HPRT (In vitro mammalian cell gene mutation test, OECD TG 476)"),
          tags$li(strong("Chromosome Aberration (CA):"), " CAvt (In vitro mammalian chromosome aberration test, OECD TG 473)")
        ),
        hr(),
        h4("Evaluate Starting Assay Outcomes"),
        p("Select the outcomes for the recommended starting assays (or keep 'Positive' / 'Negative' to preview specific scenarios) to calculate the corresponding P(GM) and P(CA):"),
        fluidRow(
          column(4, selectInput("hprt_3a", "HPRT (Gene Mutation):", choices = c("Positive", "Negative", "NA"), selected = "Positive")),
          column(4, selectInput("cavt_3a", "CAvt (Chromosome Aberration):", choices = c("Positive", "Negative", "NA"), selected = "Positive"))
        ),
        actionButton("calc_3a_btn", "Calculate P(GM) and P(CA)", class = "btn-success"),
        br(), br(),
        
        uiOutput("results_3a_ui"),
        
        hidden(
          div(id = "more_testing_3a_div",
              hr(),
              h4("Would you like to suggest a follow-up assay?"),
              actionButton("yes_to_3b_btn", "Yes", class = "btn-info"),
              actionButton("no_to_end_3a_btn", "No, finish", class = "btn-secondary")
          )
        ),
        
        br(), hr(),
        p("Reference: See our paper for full details on optimal testing sequences: ",
          tags$a(href = "https://doi.org/", "Link to Paper", target = "_blank")),
        actionButton("back_btn_3a", label = "Back", icon = icon("arrow-left"))
    )
  ),
  
  # ==========================================
  # Page 3b: Path B - Input Existing Data & Follow-up
  # ==========================================
  hidden(
    div(id = "page3b",
        h3("Input Existing Data"),
        fluidRow(
          column(4, selectInput("ames", "Ames (GM):", choices = c("NA", "Positive", "Negative"))),
          column(4, selectInput("mla", "MLA (GM):", choices = c("NA", "Positive", "Negative"))),
          column(4, selectInput("hprt", "HPRT (GM):", choices = c("NA", "Positive", "Negative")))
        ),
        fluidRow(
          column(4, selectInput("mnvt", "MNvt (CA):", choices = c("NA", "Positive", "Negative"))),
          column(4, selectInput("cavt", "CAvt (CA):", choices = c("NA", "Positive", "Negative")))
        ),
        
        hr(),
        h4("Conclusive Decision Threshold"),
        radioButtons("keep_thresh_3b", 
                     "Default conclusive decision threshold: Positive if probability > 0.975 and Negative if probability < 0.025. Do you want to keep the default threshold?",
                     choices = c("Yes (Keep default)" = "yes", "No (Customize)" = "no"),
                     selected = "yes", inline = TRUE),
        hidden(
          div(id = "custom_thresh_3b_div",
              fluidRow(
                column(4, numericInput("thresh_upper_3b", "Positive Threshold (Upper):", value = 0.975, min = 0.500, max = 0.999, step = 0.001)),
                column(4, numericInput("thresh_lower_3b", "Negative Threshold (Lower):", value = 0.025, min = 0.001, max = 0.500, step = 0.001))
              )
          )
        ),
        
        hr(),
        h4("Prior & Rodent Carcinogenicity (RC) Settings"),
        radioButtons("provide_prior", 
                     tags$span("A default prior belief of 0.5 is used for the probability of the compound being a rodent carcinogen. Do you want to update it?", 
                               actionLink("help_prior", icon("question-circle"))),
                     choices = c("No" = "no", "Yes" = "yes"), inline = TRUE),
        hidden(
          div(id = "prior_input_div",
              numericInput("pi_RC_input", "Prior Probability (pi_RC):", value = 0.500, min = 0.000, max = 1.000, step = 0.001)
          )
        ),
        
        radioButtons("show_rc",
                     "Do you want to see P(RC) (Probability of Rodent Carcinogenicity)?",
                     choices = c("No" = "no", "Yes" = "yes"), selected = "no", inline = TRUE),
        p(em("Note: P(RC) is not stable in the current version as discussed in the paper."), style = "color: #666;"),
        
        br(),
        actionButton("calc_btn", "Calculate", class = "btn-success"),
        br(), br(),
        
        # Results Output Section
        uiOutput("results_ui"),
        
        # Ask if user wants follow-up when more testing is needed
        hidden(
          div(id = "ask_followup_prompt_div",
              hr(),
              h4("Would you like to suggest a follow-up assay?"),
              actionButton("ask_followup", "Yes", class = "btn-info"),
              actionButton("no_to_end_3b_btn", "No, finish", class = "btn-secondary")
          )
        ),
        
        # Follow-up Recommendation Output Section
        hidden(
          div(id = "followup_div",
              hr(),
              h4("Follow-up Assay Recommendations"),
              uiOutput("followup_ui"),
              br(),
              actionButton("to_end_from_followup_btn", "Finish", class = "btn-primary")
          )
        ),
        
        br(), hr(),
        actionButton("back_btn_3b", label = "Back", icon = icon("arrow-left"))
    )
  ),
  
  # ==========================================
  # Page 4: Ending Page
  # ==========================================
  hidden(
    div(id = "page4",
        h3("Thank you for using the GATE app!"),
        p("Version 0.2.1"),
        hr(),
        h4("Contact Information"),
        p("Send comments, questions, and feedback to:"),
        p(strong("Zheng Zhou:"), " ", tags$a(href = "mailto:zheng.zhou@ruokavirasto.fi", "zheng.zhou@ruokavirasto.fi")),
        p("Reference Paper: ", tags$a(href = "https://doi.org/", "Link to Paper", target = "_blank")),
        p("Github repository: ", tags$a(href = "https://github.com/zhengzhou0105/GATE", "G.A.T.E.")),
        br(),
        actionButton("restart_btn", "Return to Start", class = "btn-primary", icon = icon("redo"))
    )
  )
)
