ui <- fluidPage(
  theme = shinytheme("slate"),
  shinyWidgets::chooseSliderSkin(skin = "Flat", color = "#37312c"),
  tags$head(
    tags$style(HTML("
      .myclass pre {
        color: white;
        background-color: #272b30;
        font-family: arial;
        border-color: #272b30;
      }")),
    tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: #272b30;  color:white; border-style: none; font-size: 0.8em}
    .tabbable > .nav > li[class=active]    > a {background-color: #3e444c; color:white; border-style: none;}
  ")),
    tags$style("[type = 'number'] { font-size:0.7em;height:50px;}"),
    tags$style(type = 'text/css', ".irs-grid-text { font-size: 5pt; color: white}"),
    tags$style(type = "text/css", ".irs {max-width: 200px; max-height: 50px; }"),
    tags$style(type = "text/css", ".form-control {max-width: 100px; max-height: 20px; }"),

    tags$style(type = "text/css", ".form-group {max-width: 20px; max-height: 70px}"),
    tags$style(type = "text/css", ".control-label {font-size: 0.8em; color: white}"),
    tags$style(type = "text/css", ".checkbox {font-size: 0.8em; color: white}"),
    tags$style(type = "text/css", ".selectize-control {font-size: 0.8em; max-width: 100px; max-height: 30px}")
  ),

  titlePanel(""),
  sidebarLayout(
    sidebarPanel(position = "left",
                 tabsetPanel(type = "tabs", id = "tabs_settings",
                             tabPanel("Main", value = 1,
                                      conditionalPanel(condition = "input.tabs1==1",
                                                       sliderInput(inputId = 'init_pop_size_1',
                                                                   label = "Initial Population Size",
                                                                   value = 100, min = 2, max = 1000),
                                                       sliderInput(inputId = 'num_gen_simple',
                                                                   label = "Number of Generations",
                                                                   value = 20, min = 2, max = 100),
                                                       sliderInput(inputId = 'put',
                                                                   label = 'Putting individuals',
                                                                   value = 0, min = 0, max = 100),
                                                       sliderInput(inputId = 'pull',
                                                                   label = 'Pulling individuals',
                                                                   value = 0, min = 0, max = 100),
                                                       sliderInput(inputId = 'init_frac_simple',
                                                                   label = 'Starting fraction of focal ancestry',
                                                                   value = 0.5, min = 0, max = 1),
                                                       selectInput("density_model", "Density dependence:",
                                                                   c("Weak", "Strong", "Manual")),
                                                       checkboxInput("model_used_single",
                                                                     "Use explicit recombination",
                                                                     value = FALSE)
                                      ),
                                      conditionalPanel(condition = "input.tabs1==3",
                                                       sliderInput(inputId = 'init_pop_size_3',
                                                                   label = "Initial Population Size",
                                                                   value = 100, min = 2, max = 1000),
                                                       sliderInput(inputId = 'num_gen_optim_s',
                                                                   label = "Number of Generations",
                                                                   value = 20, min = 2, max = 100),
                                                       sliderInput(inputId = 'init_frac_optim',
                                                                   label = 'Starting fraction of Hawaiian Ancestry',
                                                                   value = 0.8, min = 0, max = 1),
                                                       selectInput("density_model_2", "Density dependence:",
                                                                   c("Weak", "Strong", "Manual")),
                                                       checkboxGroupInput("optim_choice",
                                                                          label = "Optimize",
                                                                          choices = list("Put",
                                                                                         "Pull"),
                                                                          selected = "Put"),
                                                       checkboxInput("model_used_s",
                                                                     "Use explicit recombination",
                                                                     value = FALSE)
                                      ),
                                      conditionalPanel(condition = "input.tabs1==4",
                                                       sliderInput(inputId = 'init_pop_size_4',
                                                                   label = "Initial Population Size",
                                                                   value = 100, min = 2, max = 1000),
                                                       sliderInput(inputId = 'num_gen_optim_c',
                                                                   label = "Number of Generations",
                                                                   value = 20, min = 2, max = 100),
                                                       sliderInput(inputId = 'total_put',
                                                                   label = 'Put: Total number of individuals',
                                                                   value = 100, min = 0, max = 1000),
                                                       sliderInput(inputId = 'total_pull',
                                                                   label = 'Pull: Total number of individuals',
                                                                   value = 100, min = 0, max = 1000),
                                                       sliderInput(inputId = 'init_frac_optim_complex',
                                                                   label = 'Starting fraction of focal Ancestry',
                                                                   value = 0.8, min = 0, max = 1),
                                                       selectInput("density_model_3", "Density dependence:",
                                                                   c("Weak", "Strong", "Manual")),
                                                       checkboxInput("model_used_c",
                                                                     "Use explicit recombination",
                                                                     value = FALSE)
                                      )
                             ),
                             tabPanel("Advanced", value = 2,
                                      numericInput(inputId = 'K',
                                                   label = "Carrying Capacity of ecosystem",
                                                   value = 400, min = 2, max = 1000, step = 50),
                                      numericInput(inputId = 'fmd',
                                                   label = "Nesting Risk",
                                                   value = 0.2, min = 0, max = 1, step = 0.01),
                                      numericInput(inputId = 'nest_succes_rate',
                                                   label = "Nest Succes Rate",
                                                   value = 0.387, min = 0, max = 1, step = 0.01),
                                      numericInput(inputId = 'morgan',
                                                   label = "Size of Genome (in Morgan)",
                                                   value = 1.0, min = 0, max = 3, step = 0.1),
                                      numericInput(inputId = 'num_repl',
                                                   label = "Number of replicates",
                                                   value = 1, min = 1, max = 10, step = 1),
                                      numericInput(inputId = 'max_age',
                                                   label = "Maximum Age",
                                                   value = 6, min = 1, max = 20, step = 1),
                                      numericInput(inputId = 'clutch_size',
                                                   label = "Clutch Size",
                                                   value = 6, min = 1, max = 20, step = 1),
                                      numericInput(inputId = 'clutch_sd',
                                                   label = "SD Clutch Size",
                                                   value = 1, min = 0, max = 2, step = 0.1),
                                      numericInput(inputId = 'sex_ratio_put',
                                                   label = "Sex Ratio of Put individuals (males / females)",
                                                   value = 0.5, min = 0, max = 1, step = 0.05),
                                      numericInput(inputId = 'sex_ratio_offspring',
                                                   label = "Sex Ratio of offspring (males / females)",
                                                   value = 0.5, min = 0, max = 1, step = 0.05),
                                      numericInput(inputId = 'target_frequency',
                                                   label = "Target Frequency (used in optimization)",
                                                   value = 0.999, min = 0, max = 1, step = 0.01)
                             ),
                             tabPanel("Density Dependence", value = 3,
                                      numericInput(inputId = 'smin',
                                                  label = "Minimum Survival Rate",
                                                  value = 0.5, min = 0, max = 1, step = 0.05),
                                      numericInput(inputId = 'smax',
                                                  label = "Maximum Survival Rate",
                                                  value = 0.9, min = 0, max = 1, step = 0.05),
                                      numericInput(inputId = 'b',
                                                  label = "Steepness Survival curve",
                                                  value = -2, min = -3, max = 0, step = 0.05),
                                      numericInput(inputId = 'p',
                                                  label = "Density of maximum steepness Survival curve",
                                                  value = 0.5, min = 0, max = 2, step = 0.05)
                             )
                 )
    ),
    mainPanel("",
              tabsetPanel(type = "tabs", id = "tabs1",
                          tabPanel("Simulation", value = 1,
                                   plotOutput("simple_plots")),
                          tabPanel("Simple Optimization", value = 3,
                                   plotOutput('Optim_simple_plots'),
                                   div(class = "myclass",
                                       verbatimTextOutput("selected_var"))
                          ),
                          tabPanel("Complex Optimization", value = 4,
                                   plotOutput('Optim_complex_plots'),
                                   div(class = "myclass",
                                       verbatimTextOutput("complex_text_output")
                                   ))
              )
    )
  )
)

server <- function(input, output) {


  simple_data <- reactive({ simulate_policy(initial_population_size = input$init_pop_size_1,
                                            nest_success_rate = input$nest_succes_rate,
                                            nesting_risk = input$fmd,
                                            num_generations = input$num_gen_simple,
                                            K = input$K,
                                            pull = input$pull,
                                            put = input$put,
                                            num_replicates  = input$num_repl,
                                            starting_freq = input$init_frac_simple,
                                            morgan = input$morgan,
                                            max_age = input$max_age,
                                            mean_clutch_size = input$clutch_size,
                                            sd_clutch_size = input$clutch_sd,
                                            smin = ifelse(input$density_model == "Manual", input$smin, 0.5),
                                            smax = ifelse(input$density_model == "Manual", input$smax, 0.9),
                                            b = ifelse(input$density_model == "Manual", input$b,
                                                       ifelse(input$density_model == "Strong", -5, -2)),
                                            p = ifelse(input$density_model == "Manual", input$p,
                                                       ifelse(input$density_model == "Strong", 0.45, 0.5)),
                                            sex_ratio_put = input$sex_ratio_put,
                                            sex_ratio_offspring = input$sex_ratio_offspring,
                                            use_simplified_model = 1 - input$model_used_single)  })

  output$simple_plots <- renderPlot({
    to_plot <- simple_data()

    p1 <- to_plot$results %>%
      ggplot(aes(x = t, y = freq_hawaii, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Average focal ancestry") +
      ylim(0, 1) +
      theme(legend.position = "top")

    focal_y <- 1.03 * tail(to_plot$results$freq_hawaii, 1)
    if (round(tail(to_plot$results$freq_hawaii, 1), 2) >= 0.99) {
      focal_y <- 0.95 * tail(to_plot$results$freq_hawaii, 1)
    }

    p1 <- p1 +
      annotate("text", x = max(to_plot$results$t), y = focal_y,
               label = round(tail(to_plot$results$freq_hawaii, 1), 2),
               hjust = 1)

    p2 <- ggplot(to_plot$results, aes(x = t, y = Num_individuals, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Total number of individuals")

    p3 <- to_plot$results %>%
      dplyr::mutate("Males" = Num_males) %>%
      dplyr::mutate("Females" = Num_females) %>%
      ggplot(aes(x = t, group = replicate)) +
      geom_line(aes(y = Males, color = "Males")) +
      geom_line(aes(y = Females, group = replicate, color = "Females")) +
      labs(x = "Generation",
           y = "Number of individuals",
           color = "Sex") +
      theme(legend.position = "top")

    ggthemr::ggthemr(palette = "earth",
                     type = "outer",
                     spacing = 2)

    egg::ggarrange(p1, p2, p3, nrow = 1)

  } , bg = "transparent")

  ########  OPTIMIZATION #########################################################

  optim_data_simple <-  reactive({
    get_optim_data_simple(initial_population_size = input$init_pop_size_3,
                          nest_success_rate = input$nest_succes_rate,
                          female_death_rate = input$fmd,
                          num_generations = input$num_gen_optim_s,
                          K = input$K,
                          num_replicates = input$num_repl,
                          target_frequency = input$target_frequency,
                          optim_choice = input$optim_choice,
                          morgan = input$morgan,
                          starting_freq = input$init_frac_optim,
                          use_complex_model = input$model_used_s,
                          max_age = input$max_age,
                          mean_clutch_size = input$clutch_size,
                          sd_clutch_size = input$clutch_sd,
                          smin = ifelse(input$density_model_2 == "Manual", input$smin, 0.5),
                          smax = ifelse(input$density_model_2 == "Manual", input$smax, 0.9),
                          b = ifelse(input$density_model_2 == "Manual", input$b,
                                     ifelse(input$density_model_2 == "Strong", -5, -2)),
                          p = ifelse(input$density_model_2 == "Manual", input$p,
                                     ifelse(input$density_model_2 == "Strong", 0.45, 0.5)),
                          sex_ratio_put = input$sex_ratio_put,
                          sex_ratio_offspring = input$sex_ratio_offspring)
  })

  optim_data_complex <-  reactive({
    get_optim_data_complex(initial_population_size = input$init_pop_size_4,
                           nest_success_rate = input$nest_succes_rate,
                           female_death_rate = input$fmd,
                           num_generations = input$num_gen_optim_c,
                           K = input$K,
                           num_replicates = input$num_repl,
                           target_frequency = input$target_frequency,
                           morgan = input$morgan,
                           total_put = input$total_put,
                           total_pull = input$total_pull,
                           starting_freq = input$init_frac_optim_complex,
                           use_complex_model = input$model_used_c,
                           max_age = input$max_age,
                           mean_clutch_size = input$clutch_size,
                           sd_clutch_size = input$clutch_sd,
                           smin = ifelse(input$density_model_3 == "Manual", input$smin, 0.5),
                           smax = ifelse(input$density_model_3 == "Manual", input$smax, 0.9),
                           b = ifelse(input$density_model_3 == "Manual", input$b,
                                      ifelse(input$density_model_3 == "Strong", -5, -2)),
                           p = ifelse(input$density_model_3 == "Manual", input$p,
                                      ifelse(input$density_model_3 == "Strong", 0.45, 0.5)),
                           sex_ratio_put = input$sex_ratio_put,
                           sex_ratio_offspring = input$sex_ratio_offspring)
  })


  output$Optim_simple_plots <- renderPlot({
    to_plot <- optim_data_simple()


    final_freq <- round(to_plot$final_freq, digits = 3)

    for_render_text <- c()
    for_render_text <- c(for_render_text, " Target frequency was: ", input$target_frequency, "\n")
    for_render_text <- c(for_render_text, " Final frequency was: ", final_freq, "\n")

    if (final_freq >= input$target_frequency) {
      for_render_text <- c(for_render_text, "Target frequency was reached\n")
    } else {
      for_render_text <- c(for_render_text, "Target frequency was NOT reached\n")
    }

    if (length(input$optim_choice) == 1) {
      if (input$optim_choice == "Pull") {
        for_render_text <- c(for_render_text, "\n",
                             "Advice is to pull ", round(to_plot$pull), " individuals per generation")
      }
      if (input$optim_choice == "Put") {
        for_render_text <- c(for_render_text, "\n", "Advice is to put ", round(to_plot$put), " individuals per generation")
      }
    }
    if (length(input$optim_choice) == 2) {

      for_render_text <- c(for_render_text, "\n",
                           "Advice is to put ", round(to_plot$put), " individuals per generation\n", "         and pull ", round(to_plot$pull), " individuals per generation")
    }

    output$selected_var <- renderText({ for_render_text })

    p1 <- to_plot$results %>%
      ggplot(aes(x = t, y = freq_hawaii, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Average focal ancestry") +
      ylim(0, 1) +
      theme(legend.position = "top")

    focal_y <- 1.03 * tail(to_plot$results$freq_hawaii, 1)
    if (round(tail(to_plot$results$freq_hawaii, 1), 2) >= 0.99) {
      focal_y <- 0.95 * tail(to_plot$results$freq_hawaii, 1)
    }

    p1 <- p1 +
      annotate("text", x = max(to_plot$results$t), y = focal_y,
               label = round(tail(to_plot$results$freq_hawaii, 1), 2),
               hjust = 1)

    p2 <- ggplot(to_plot$results, aes(x = t, y = Num_individuals, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Total number of individuals")

    p3 <- to_plot$results %>%
      dplyr::mutate("Males" = Num_males) %>%
      dplyr::mutate("Females" = Num_females) %>%
      ggplot(aes(x = t, group = replicate)) +
      geom_line(aes(y = Males, color = "Males")) +
      geom_line(aes(y = Females, group = replicate, color = "Females")) +
      labs(x = "Generation",
           y = "Number of individuals",
           color = "Sex") +
      theme(legend.position = "top")

    p4 <-  gather(to_plot$curve, key = "type", value = "number", -t) %>%
      ggplot(aes(x = t, y = number, col = type)) +
      geom_line() +
      ylab("Amount") +
      xlab("Number of Generations") +
      theme(legend.position = "top")
    ggthemr::ggthemr(palette = "earth",
                     type = "outer",
                     spacing = 2)
    egg::ggarrange(p1, p2, p3, p4, nrow = 1)

  }, bg = "transparent")

  output$Optim_complex_plots <- renderPlot({
    to_plot <- optim_data_complex()
    for_text <- to_plot$curve
    # tibble with t, pull, put
    final_freq <- round(to_plot$final_freq, 3)
    for_render_text <- c()
    for_render_text <- c(for_render_text, " Target frequency was: ", input$target_frequency, "\n")
    for_render_text <- c(for_render_text, "Final frequency was: ", final_freq, "\n")
    if (final_freq >= input$target_frequency) {
      for_render_text <- c(for_render_text, "Target frequency was reached\n")
    } else {
      for_render_text <- c(for_render_text, "Target frequency was NOT reached\n")
    }

    for_render_text <- c(for_render_text, "Advice:", "\n")
    if (input$total_put > 0 && input$total_pull > 0) {
      for_render_text <- c(for_render_text, c("Generation", "\t", "Put", "\t", "Pull", "\n"))
    }
    if (input$total_put == 0 && input$total_pull > 0) {
      for_render_text <- c(for_render_text, c("Generation", "\t", "Pull", "\n"))
    }

    if (input$total_put > 0 && input$total_pull == 0) {
      for_render_text <- c(for_render_text, c("Generation", "\t", "Put", "\n"))
    }

    for (i in seq_along(for_text$t)) {

      if (input$total_put > 0 && input$total_pull > 0) {
        add_text <- paste(round(for_text$t[i]), "\t\t\t",
                          round(for_text$put[i]), "\t",
                          round(for_text$pull[i]), "\n")
      }
      if (input$total_put == 0 && input$total_pull > 0) {
        add_text <- paste(round(for_text$t[i]), "\t\t\t",
                          round(for_text$pull[i]), "\n")
      }
      if (input$total_put > 0 && input$total_pull == 0) {
        add_text <- paste(round(for_text$t[i]), "\t\t\t",
                          round(for_text$put[i]), "\n")
      }
      for_render_text <- c(for_render_text, add_text)
    }

    output$complex_text_output <- renderText({ for_render_text })



    p1 <- to_plot$results %>%
      ggplot(aes(x = t, y = freq_hawaii, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Average focal ancestry") +
      ylim(0, 1) +
      theme(legend.position = "top")

    focal_y <- 1.03 * tail(to_plot$results$freq_hawaii, 1)
    if (round(tail(to_plot$results$freq_hawaii, 1), 2) >= 0.99) {
      focal_y <- 0.95 * tail(to_plot$results$freq_hawaii, 1)
    }

    p1 <- p1 +
      annotate("text", x = max(to_plot$results$t), y = focal_y,
               label = round(tail(to_plot$results$freq_hawaii, 1), 2),
               hjust = 1)

    p2 <- ggplot(to_plot$results, aes(x = t, y = Num_individuals, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Total number of individuals")

    p3 <- to_plot$results %>%
      dplyr::mutate("Males" = Num_males) %>%
      dplyr::mutate("Females" = Num_females) %>%
      ggplot(aes(x = t, group = replicate)) +
      geom_line(aes(y = Males, color = "Males")) +
      geom_line(aes(y = Females, group = replicate, color = "Females")) +
      labs(x = "Generation",
           y = "Number of individuals",
           color = "Sex") +
      theme(legend.position = "top")

    p4 <-  gather(to_plot$curve, key = "type", value = "number", -t) %>%
      ggplot(aes(x = t, y = number, col = type)) +
      geom_step() +
      ylab("Amount") +
      xlab("Number of Generations") +
      theme(legend.position = "top")
    ggthemr::ggthemr(palette = "earth",
                     type = "outer",
                     spacing = 2)
    egg::ggarrange(p1, p2, p3, p4, nrow = 1)

  }, bg = "transparent")

  output$downloadData_s <- downloadHandler(
    filename = function() {
      paste0("dataset_", Sys.Date(), ".txt")
    },
    content = function(file) {
      stored_data <- read.table(input$data_for_download)
      write.table(stored_data, file, quote = F)
    }
  )

}

get_optim_data_simple <- function(initial_population_size,
                                  nest_success_rate,
                                  female_death_rate,
                                  num_generations,
                                  K,
                                  num_replicates,
                                  target_frequency,
                                  optim_choice,
                                  morgan,
                                  starting_freq,
                                  use_complex_model,
                                  max_age,
                                  mean_clutch_size,
                                  sd_clutch_size,
                                  smin,
                                  smax,
                                  b,
                                  p,
                                  sex_ratio_put,
                                  sex_ratio_offspring) {

  opt_pull = FALSE
  opt_put  = FALSE

  if (length(optim_choice) == 2) {
    opt_pull = TRUE
    opt_put = TRUE
  } else {
    if (length(optim_choice) == 1) {
      if (optim_choice == "Put") {
        opt_put = TRUE
        opt_pull = FALSE
      }
      if (optim_choice == "Pull") {
        opt_pull = TRUE
        opt_put = FALSE
      }
    }
  }

  use_simple_model = TRUE
  if (use_complex_model == TRUE) use_simple_model = FALSE

  return(simRestore::optimize_policy(initial_population_size = initial_population_size,
                                    nest_success_rate = nest_success_rate,
                                    nesting_risk = female_death_rate,
                                    num_generations = num_generations,
                                    K = K,
                                    num_replicates = num_replicates,
                                    use_simplified_model = use_simple_model,
                                    target_frequency = target_frequency,
                                    optimize_pull = opt_pull,
                                    optimize_put = opt_put,
                                    morgan = morgan,
                                    starting_freq = starting_freq,
                                    max_age = max_age,
                                    mean_clutch_size = mean_clutch_size,
                                    sd_clutch_size = sd_clutch_size,
                                    smin = smin,
                                    smax = smax,
                                    b = b,
                                    p = p,
                                    sex_ratio_put = sex_ratio_put,
                                    sex_ratio_offspring = sex_ratio_offspring,
                                    verbose = FALSE))
}

get_optim_data_complex <- function(initial_population_size,
                                   nest_success_rate,
                                   female_death_rate,
                                   num_generations,
                                   K,
                                   num_replicates,
                                   target_frequency,
                                   total_put,
                                   total_pull,
                                   morgan,
                                   starting_freq,
                                   use_complex_model,
                                   max_age,
                                   mean_clutch_size,
                                   sd_clutch_size,
                                   smin,
                                   smax,
                                   b,
                                   p,
                                   sex_ratio_put,
                                   sex_ratio_offspring) {

  use_simple_model = TRUE
  if (use_complex_model == TRUE) use_simple_model = FALSE

  return(simRestore::optimize_policy_beta_curve(
    initial_population_size = initial_population_size,
    nest_success_rate = nest_success_rate,
    nesting_risk = female_death_rate,
    num_generations = num_generations,
    K = K,
    num_replicates = num_replicates,
    use_simplified_model = use_simple_model,
    target_frequency = target_frequency,
    optimize_pull = total_pull,
    optimize_put = total_put,
    morgan = morgan,
    starting_freq = starting_freq,
    max_age = max_age,
    mean_clutch_size = mean_clutch_size,
    sd_clutch_size = sd_clutch_size,
    smin = smin,
    smax = smax,
    b = b,
    p = p,
    sex_ratio_put = sex_ratio_put,
    sex_ratio_offspring = sex_ratio_offspring,
    verbose = FALSE))
}

shinyApp(ui = ui, server = server)
