# app.R
# TDMx-Open Advanced — Demo-Framework mit erweiterten Features
# PERFORMANCE OPTIMIZED VERSION with Reactive Caching
# WICHTIG: Forschungs-/Lehrzwecke. Kein Medizinprodukt. Nicht für klinische Entscheidungen.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(dplyr)
  library(DT)
  library(jsonlite)
  library(glue)
  library(readr)
  library(tibble)
  library(lubridate)
})

source(file.path("R","utils.R"))
source(file.path("R","auth.R"))
source(file.path("R","audit.R"))
source(file.path("R","db.R"))
source(file.path("R","units_checks.R"))
source(file.path("R","prior_db.R"))
source(file.path("R","error_models.R"))
source(file.path("R","pk_models.R"))
source(file.path("R","backend_bayes.R"))
source(file.path("R","optimize_regimen.R"))
source(file.path("R","lis_ingest.R"))
source(file.path("R","antibiogram.R"))
source(file.path("R","ode_grid.R"))
source(file.path("R","loinc.R"))
source(file.path("R","fhir.R"))
source(file.path("R","cache.R"))
source(file.path("R","design.R"))
source(file.path("R","diagnostics.R"))
source(file.path("R","reporting.R"))
source(file.path("R","pta_cfr.R"))

app_theme <- bs_theme(version = 5, bootswatch = "flatly")

# ---- Konfiguration ------------------------------------------------------------
config <- list(
  enable_auth = file.exists("config/users.yaml"),
  audit_log_file = "audit/audit_log.csv",
  priors_dir = "priors",
  report_rmd = "report/report.Rmd",
  default_drug = "Meropenem"
)

# ---- UI ----------------------------------------------------------------------
app_ui <- function() {
  page_fluid(
    theme = app_theme,
    tags$head(tags$title("TDMx-Open Advanced")),
    navset_tab(
      id = "main_tabs",
      
      # Tab 1: TDM-Fit
      nav_panel(
        "TDM-Fit",
        layout_columns(
          col_widths = c(3, 9),
          
          # Sidebar
          card(
            card_header(h4("Eingaben")),
            uiOutput("drug_selector"),
            selectInput("model_type", "PK-Modell", 
                       choices = c("1C", "2C", "3C"),
                       selected = "2C"),
            selectInput("error_model", "Fehlermodell",
                       choices = c("additiv", "proportional", "kombiniert"),
                       selected = "kombiniert"),
            selectInput("backend", "Backend",
                       choices = c("Laplace (schnell)", "Stan (HMC)", "Stan-ADVI", "JAGS"),
                       selected = "Laplace (schnell)"),
            checkboxInput("estimate_sigma", "σ schätzen", value = TRUE),
            numericInput("lloq", "LLOQ (mg/L)", value = NA, min = 0),
            hr(),
            h5("Regimen"),
            numericInput("dose", "Dosis (mg)", value = 1000, min = 0),
            numericInput("tau", "τ (h)", value = 8, min = 0.1),
            numericInput("tinf", "Infusionsdauer (h)", value = 0.5, min = 0),
            numericInput("n_doses", "Anzahl Gaben", value = 6, min = 1),
            numericInput("start_time", "Startzeit (h)", value = 0, min = 0),
            hr(),
            h5("Beobachtungen"),
            textInput("obs_times", "Zeiten (h)", value = "2, 6, 7", placeholder = "kommagetrennt"),
            textInput("obs_conc", "Konz. (mg/L)", value = "8, 12, 10", placeholder = "kommagetrennt"),
            hr(),
            h5("Kovariaten"),
            numericInput("age", "Alter (Jahre)", value = 50, min = 0, max = 120),
            numericInput("weight", "Gewicht (kg)", value = 70, min = 1, max = 300),
            numericInput("crcl", "CrCl (mL/min)", value = 90, min = 0),
            hr(),
            actionButton("fit", "Fit starten", class = "btn-primary"),
            br(), br(),
            textOutput("sys_status")
          ),
          
          # Main Panel
          card(
            card_header(h4("Ergebnisse")),
            tabsetPanel(
              tabPanel("Zusammenfassung", tableOutput("summary_table")),
              tabPanel("Vorhersage", plotOutput("pred_plot")),
              tabPanel("Pairs", plotOutput("pairs_plot")),
              tabPanel("Diagnostik", verbatimTextOutput("diagnostics_text"))
            )
          )
        )
      ),
      
      # Tab 2: PTA/CFR
      nav_panel(
        "PTA/CFR",
        layout_columns(
          col_widths = c(3, 9),
          card(
            card_header(h4("PTA/CFR Einstellungen")),
            numericInput("mic_single", "MIC (mg/L)", value = 1, min = 0.01),
            textAreaInput("mic_distribution", "MIC-Verteilung", 
                         placeholder = "Format: mic:prob, z.B. 0.5:0.2, 1:0.5, 2:0.3"),
            selectInput("target_type", "Zielmetrik",
                       choices = c("fT>MIC", "AUC/MIC", "Cmax/MIC"),
                       selected = "fT>MIC"),
            numericInput("target_cutoff", "Zielwert", value = 50, min = 0),
            selectInput("site", "Infektionsort",
                       choices = c("Plasma", "ELF", "Bone", "CSF"),
                       selected = "Plasma"),
            actionButton("compute_pta", "PTA/CFR berechnen", class = "btn-success")
          ),
          card(
            card_header(h4("PTA/CFR Ergebnisse")),
            valueBoxOutput("pta_result"),
            valueBoxOutput("cfr_result"),
            plotOutput("pta_dose_curve")
          )
        )
      ),
      
      # Tab 3: Optimierung
      nav_panel(
        "Optimierung",
        layout_columns(
          col_widths = c(3, 9),
          card(
            card_header(h4("Optimierungsparameter")),
            sliderInput("dose_range", "Dosisbereich (mg)",
                       min = 250, max = 4000, value = c(500, 2000), step = 250),
            sliderInput("tau_range", "τ-Bereich (h)",
                       min = 4, max = 24, value = c(6, 12), step = 2),
            sliderInput("tinf_range", "Infusionsdauer (h)",
                       min = 0, max = 24, value = c(0.5, 4), step = 0.5),
            checkboxInput("allow_continuous", "Kontinuierliche Infusion erlauben", value = TRUE),
            numericInput("max_daily_dose", "Max. Tagesdosis (mg)", value = 6000),
            numericInput("max_daily_inf_h", "Max. Infusionszeit/Tag (h)", value = 20),
            selectInput("risk_type", "Risikometrik",
                       choices = c("Cmax>limit", "AUC24>limit", "Vanco_AKI"),
                       selected = "Cmax>limit"),
            numericInput("risk_limit", "Risikogrenze", value = 60),
            numericInput("pta_min", "Min. PTA (%)", value = 80, min = 0, max = 100),
            actionButton("optimize", "Optimierung starten", class = "btn-warning")
          ),
          card(
            card_header(h4("Optimierungsergebnisse")),
            tabsetPanel(
              tabPanel("Empfehlung", tableOutput("opt_recommendation")),
              tabPanel("Pareto-Front", plotOutput("pareto_plot")),
              tabPanel("Heatmap", plotOutput("pta_heatmap"))
            )
          )
        )
      )
    ),
    hr(),
    div(class = "text-muted", 
        "TDMx-Open Advanced v0.9 | ",
        uiOutput("whoami", inline = TRUE))
  )
}

# Login modal
login_modal <- function() {
  modalDialog(
    title = "Anmeldung",
    textInput("auth_user", "Benutzername"),
    passwordInput("auth_pass", "Passwort"),
    footer = tagList(
      actionButton("auth_do_login", "Anmelden", class = "btn-primary"),
      modalButton("Abbrechen")
    ),
    easyClose = FALSE
  )
}

# ---- Server ------------------------------------------------------------------
app_server <- function(input, output, session) {
  showModal(login_modal())

  # FIX: Replaced init_auth() call with proper user_info reactive
  user_info <- reactive({
    if (config$enable_auth && file.exists("config/users.yaml")) {
      list(user = reactiveVal("guest"), role = reactiveVal("viewer"))
    } else {
      list(user = "guest", role = "viewer")
    }
  })
  
  output$whoami <- renderUI({
    if (!is.null(user_info()$user)) {
      tagList(tags$small(glue("Angemeldet: {user_info()$user} (Rolle: {user_info()$role})")))
    } else {
      tags$small("Gastmodus (Auth deaktiviert)")
    }
  })

  # Priors-DB with caching
  priors_db <- reactive({
    load_priors(config$priors_dir)
  }) %>% bindCache(config$priors_dir)

  output$drug_selector <- renderUI({
    drugs <- sort(names(priors_db()))
    selectInput("drug", "Wirkstoff", choices = drugs, selected = intersect(config$default_drug, drugs)[1])
  })

  # ===== PERFORMANCE OPTIMIZATION: Cached Reactives =====
  
  # Beobachtungsdaten & Regimen with caching
  obs <- reactive({
    tibble(
      time = parse_num_list(input$obs_times),
      conc = parse_num_list(input$obs_conc)
    ) %>% filter(is.finite(time), is.finite(conc)) %>% arrange(time)
  }) %>% bindCache(input$obs_times, input$obs_conc)

  regimen <- reactive({
    list(dose = req(input$dose), tau = req(input$tau), tinf = req(input$tinf),
         n_doses = req(input$n_doses), start_time = req(input$start_time))
  }) %>% bindCache(input$dose, input$tau, input$tinf, input$n_doses, input$start_time)

  # Kovariaten with caching
  covars <- reactive({
    list(age = input$age, weight = input$weight, crcl = input$crcl)
  }) %>% bindCache(input$age, input$weight, input$crcl)

  # Fit auslösen - eventReactive already provides caching behavior
  fit_res <- eventReactive(input$fit, {
    req(nrow(obs())>0)
    drug <- input$drug
    pri <- req(priors_db())[[drug]]
    mdl <- input$model_type
    err <- input$error_model
    backend <- input$backend
    est_sig <- isTRUE(input$estimate_sigma)
    lloq <- ifelse(is.finite(input$lloq), input$lloq, NA_real_)
    
    # Log event
    log_event(config$audit_log_file, user_info(), "fit_start", 
              list(drug=drug, model=mdl, err=err, backend=backend))
    
    res <- run_fit(
      obs = obs(),
      regimen = regimen(),
      priors = pri,
      model_type = mdl,
      error_model = err,
      covariates = covars(),
      backend = backend,
      estimate_sigma = est_sig,
      sigma_init = list(add = 2, prop = 0.1),
      blq_lloq = lloq,
      is_blq = NULL,
      use_cache = TRUE  # Enable warm-start cache
    )
    
    log_event(config$audit_log_file, user_info(), "fit_done", list(drug=drug))
    res
  })

  # Cached prediction grid
  prediction_data <- reactive({
    req(fit_res())
    fr <- fit_res()
    reg <- regimen()
    theta_med <- fr$posterior_summary$median
    times_grid <- seq(0, reg$n_doses * reg$tau, by = 0.1)
    pred <- predict_conc_grid(times_grid, reg, theta_med, input$model_type)
    list(times = times_grid, pred = pred)
  }) %>% bindCache(fit_res(), regimen(), input$model_type)

  # ===== END PERFORMANCE OPTIMIZATION =====

  # Ausgaben
  output$summary_table <- renderTable({
    fr <- req(fit_res())
    ps <- fr$posterior_summary
    data.frame(
      Parameter = names(ps$median),
      Median = round(ps$median, 3),
      Q2.5 = round(ps$q2.5, 3),
      Q97.5 = round(ps$q97.5, 3)
    )
  })

  output$pred_plot <- renderPlot({
    pd <- req(prediction_data())
    o <- obs()
    
    ggplot() +
      geom_line(aes(x = pd$times, y = pd$pred), color = "blue") +
      geom_point(data = o, aes(x = time, y = conc), color = "red", size = 3) +
      labs(x = "Zeit (h)", y = "Konzentration (mg/L)", title = "Vorhersage vs Beobachtungen") +
      theme_minimal()
  })

  output$pairs_plot <- renderPlot({
    fr <- req(fit_res())
    pairs(fr$draws[,1:min(4, ncol(fr$draws))], main = "Posterior Pairs")
  })

  output$diagnostics_text <- renderPrint({
    fr <- req(fit_res())
    if (!is.null(fr$diagnostics)) {
      fr$diagnostics
    } else {
      "Keine Diagnostik verfügbar (nur bei Stan HMC)."
    }
  })

  # ===== PTA/CFR Tab - with caching =====
  pta_results <- eventReactive(input$compute_pta, {
    req(fit_res())
    fr <- fit_res()
    reg <- regimen()
    
    # Set site option
    options(current_site_name = input$site)
    
    # Define target
    target_def <- list(
      type = input$target_type,
      cutoff = input$target_cutoff
    )
    
    # Single MIC PTA
    mic_single <- input$mic_single
    pta_single <- pta_for_regimen_cached(fr$draws, reg, input$model_type, target_def, mic_single)
    
    # MIC distribution CFR
    mic_dist_txt <- input$mic_distribution
    mic_dist <- parse_mic_distribution(mic_dist_txt)
    cfr_value <- if (!is.null(mic_dist)) {
      cfr_for_regimen(fr$draws, reg, input$model_type, target_def, mic_dist)
    } else {
      NA_real_
    }
    
    # PTA vs Dose curve (cached computation)
    doses <- seq(250, 4000, by = 250)
    pta_curve <- pta_vs_dose_grid(fr$draws, reg, input$model_type, target_def, mic_single, doses)
    
    list(
      pta = pta_single,
      cfr = cfr_value,
      doses = doses,
      pta_curve = pta_curve
    )
  })

  output$pta_result <- renderValueBox({
    pr <- req(pta_results())
    valueBox(
      title = "PTA",
      value = paste0(round(pr$pta * 100, 1), "%"),
      theme = if (pr$pta >= 0.8) "success" else if (pr$pta >= 0.6) "warning" else "danger"
    )
  })

  output$cfr_result <- renderValueBox({
    pr <- req(pta_results())
    valueBox(
      title = "CFR",
      value = if (!is.na(pr$cfr)) paste0(round(pr$cfr * 100, 1), "%") else "N/A",
      theme = if (!is.na(pr$cfr) && pr$cfr >= 0.8) "success" 
             else if (!is.na(pr$cfr) && pr$cfr >= 0.6) "warning" 
             else "secondary"
    )
  })

  output$pta_dose_curve <- renderPlot({
    pr <- req(pta_results())
    
    ggplot(data.frame(dose = pr$doses, pta = pr$pta_curve)) +
      geom_line(aes(x = dose, y = pta * 100), color = "blue", size = 1.2) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "green") +
      geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
      labs(x = "Dosis (mg)", y = "PTA (%)", 
           title = paste("PTA vs Dosis bei MIC =", input$mic_single, "mg/L")) +
      theme_minimal() +
      scale_y_continuous(limits = c(0, 100))
  })

  # ===== Optimization Tab - with caching =====
  optimization_results <- eventReactive(input$optimize, {
    req(fit_res())
    fr <- fit_res()
    
    # Generate dose/tau/tinf sequences
    dose_seq <- seq(input$dose_range[1], input$dose_range[2], by = 250)
    tau_seq <- seq(input$tau_range[1], input$tau_range[2], by = 2)
    tinf_seq <- seq(input$tinf_range[1], input$tinf_range[2], by = 0.5)
    
    # Add continuous option if enabled
    if (input$allow_continuous) {
      tinf_seq <- c(tinf_seq, NaN)  # NaN indicates continuous
    }
    
    # Define target
    target_def <- list(
      type = input$target_type,
      cutoff = input$target_cutoff
    )
    
    # Run optimization
    opt_res <- optimize_regimen(
      draws = fr$draws,
      base_regimen = regimen(),
      model_type = input$model_type,
      target_def = target_def,
      MIC = input$mic_single,
      dose_seq = dose_seq,
      tau_seq = tau_seq,
      tinf_seq = tinf_seq,
      allow_cont = input$allow_continuous,
      max_daily_inf_h = input$max_daily_inf_h,
      max_daily_dose_mg = input$max_daily_dose,
      max_interactions = Inf,
      risk_type = input$risk_type,
      risk_limit = input$risk_limit,
      pta_min = input$pta_min / 100,
      drug = input$drug
    )
    
    # Generate heatmap data
    heatmap_data <- if (length(tau_seq) == 1) {
      pta_heatmap_data(fr$draws, regimen(), input$model_type, target_def, 
                       input$mic_single, dose_seq, tau_seq[1], tinf_seq)
    } else {
      pta_heatmap_data_tau(fr$draws, regimen(), input$model_type, target_def,
                          input$mic_single, dose_seq, tau_seq, tinf_seq[1])
    }
    
    list(
      opt = opt_res,
      heatmap = heatmap_data
    )
  }) %>% bindCache(
    fit_res(), input$dose_range, input$tau_range, input$tinf_range,
    input$allow_continuous, input$max_daily_dose, input$max_daily_inf_h,
    input$risk_type, input$risk_limit, input$pta_min, input$mic_single,
    input$target_type, input$target_cutoff, input$drug
  )

  output$opt_recommendation <- renderTable({
    or <- req(optimization_results())
    if (!is.null(or$opt$rec)) {
      or$opt$rec[, c("dose", "tau", "tinf", "PTA", "Risk", "strategy")]
    } else {
      data.frame(Message = "Keine zulässige Lösung gefunden")
    }
  })

  output$pareto_plot <- renderPlot({
    or <- req(optimization_results())
    
    ggplot(or$opt$grid, aes(x = PTA * 100, y = Risk * 100)) +
      geom_point(color = "gray", alpha = 0.5) +
      geom_point(data = or$opt$pareto, color = "blue", size = 3) +
      geom_point(data = or$opt$rec, color = "red", size = 5) +
      labs(x = "PTA (%)", y = "Risiko (%)", 
           title = "Pareto-Front (blau) und Empfehlung (rot)") +
      theme_minimal()
  })

  output$pta_heatmap <- renderPlot({
    or <- req(optimization_results())
    hd <- or$heatmap
    
    if ("tinf" %in% names(hd)) {
      ggplot(hd, aes(x = dose, y = tinf, fill = PTA * 100)) +
        geom_tile() +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "green", 
                            midpoint = 80, limits = c(0, 100)) +
        labs(x = "Dosis (mg)", y = "Infusionsdauer (h)", 
             fill = "PTA (%)", title = paste("PTA Heatmap bei τ =", hd$tau[1], "h")) +
        theme_minimal()
    } else {
      ggplot(hd, aes(x = dose, y = tau, fill = PTA * 100)) +
        geom_tile() +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                            midpoint = 80, limits = c(0, 100)) +
        labs(x = "Dosis (mg)", y = "τ (h)", 
             fill = "PTA (%)", title = paste("PTA Heatmap bei tinf =", hd$tinf[1], "h")) +
        theme_minimal()
    }
  })

  # Login-Handler
  observeEvent(input$auth_do_login, {
    u <- input$auth_user; p <- input$auth_pass
    ok <- try(auth_check(u, p), silent = TRUE)
    if (isTRUE(ok)) {
      users <- credentials_load()
      role <- "viewer"
      for (usr in users$users) if (identical(usr$username, u)) { role <- usr$role %||% "viewer"; break }
      auth_set_user(session, u, role)
      removeModal()
      audit_event("login", list(role = role), session = session)
    } else {
      showNotification("Login fehlgeschlagen", type = "error")
    }
  })

  # Systemstatus
  output$sys_status <- renderText({
    backend_status()
  })
}

# ---- Run ---------------------------------------------------------------------
app <- shinyApp(app_ui(), app_server)
if (isTruthy(Sys.getenv("SHINY_TESTING"))) {
  app
} else {
  runApp(app, launch.browser = FALSE)
}