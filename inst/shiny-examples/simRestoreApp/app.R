require(magrittr)
require(ggplot2)
require(shinyBS)
require(simRestore)

data_storage <- c()

ui <- fluidPage(
  theme = shinythemes::shinytheme("slate"),
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
                                                       shinyBS::bsTooltip('init_pop_size_1',
                                                                 "Population size at the start of the simulation"),

                                                       sliderInput(inputId = 'num_gen_simple',
                                                                   label = "Number of Generations",
                                                                   value = 20, min = 2, max = 100),
                                                       shinyBS::bsTooltip('num_gen_simple',
                                                                          "Total number of generations simulated"),

                                                       sliderInput(inputId = 'put',
                                                                   label = 'Putting individuals',
                                                                   value = 0, min = 0, max = 100),
                                                       shinyBS::bsTooltip('put',
                                                                          "Number of individuals added per generation"),

                                                       sliderInput(inputId = 'pull',
                                                                   label = 'Pulling individuals',
                                                                   value = 0, min = 0, max = 100),
                                                       shinyBS::bsTooltip('pull',
                                                                          "Number of individuals removed per generation"),

                                                       sliderInput(inputId = 'init_frac_simple',
                                                                   label = 'Starting fraction of focal ancestry',
                                                                   value = 0.5, min = 0, max = 1),
                                                       shinyBS::bsTooltip('init_frac_simple',
                                                                          "Initial frequency of focal ancestry"),

                                                       selectInput(inputId = "density_model",
                                                                   label = "Density dependence: ",
                                                                   choices = c("Weak", "Strong", "Manual")),
                                                       shinyBS::bsTooltip("density_model",
                                                                          "The user can pick between two pre-defined parameter sets that implement weak or strong density dependence. Alternatively the user can modify parameters manually in the corresponding tab."),

                                                       checkboxInput(inputId = "model_used_single",
                                                                     label = "Use explicit recombination",
                                                                     value = FALSE),
                                                       shinyBS::bsTooltip("model_used_single",
                                                                          "When unchecked, a simplified genetic model is used. When checked, explicit recombination is modeled"),
                                                       downloadButton("download_gen1", label = "Download Genetics"),
                                                       shinyBS::bsTooltip("download_gen1",
                                                                          "Download local ancestry information of the last generation as a tibble"),
                                                       downloadButton("download_res1", label = "Download results"),
                                                       shinyBS::bsTooltip("download_res1",
                                                                          "Download results as text file")
                                      ),
                                      conditionalPanel(condition = "input.tabs1==3",
                                                       sliderInput(inputId = 'init_pop_size_3',
                                                                   label = "Initial Population Size",
                                                                   value = 100, min = 2, max = 1000),
                                                       shinyBS::bsTooltip("init_pop_size_3",
                                                                          "Population size at the start of the simulation"),

                                                       sliderInput(inputId = 'num_gen_optim_s',
                                                                   label = "Number of Generations",
                                                                   value = 20, min = 2, max = 100),
                                                       shinyBS::bsTooltip("num_gen_optim_s",
                                                                          "Total number of generations simulated"),

                                                       sliderInput(inputId = 'init_frac_optim',
                                                                   label = 'Starting frequency of focal Ancestry',
                                                                   value = 0.8, min = 0, max = 1),
                                                       shinyBS::bsTooltip('init_frac_optim',
                                                                          "Initial frequency of focal ancestry"),

                                                       selectInput("density_model_2", "Density dependence:",
                                                                   c("Weak", "Strong", "Manual")),
                                                       shinyBS::bsTooltip("density_model_2",
                                                                          "The user can pick between two pre-defined parameter sets that implement weak or strong density dependence. Alternatively the user can modify parameters manually in the corresponding tab."),

                                                       checkboxGroupInput("optim_choice",
                                                                          label = "Optimize",
                                                                          choices = list("Put",
                                                                                         "Pull"),
                                                                          selected = "Put"),
                                                       shinyBS::bsTooltip("optim_choice",
                                                                          "Should only putting be optimized, only pulling, or both?"),

                                                       checkboxInput("model_used_s",
                                                                     "Use explicit recombination",
                                                                     value = FALSE),
                                                       shinyBS::bsTooltip("model_used_s",
                                                                          "When unchecked, a simplified genetic model is used. When checked, explicit recombination is modeled"),
                                                       downloadButton("download_gen2", label = "Download Genetics"),
                                                       shinyBS::bsTooltip("download_gen2",
                                                                          "Download local ancestry information of the last generation as a tibble"),
                                                       downloadButton("download_res2", label = "Download results"),
                                                       shinyBS::bsTooltip("download_res2",
                                                                          "Download results as text file")

                                      ),
                                      conditionalPanel(condition = "input.tabs1==4",
                                                       sliderInput(inputId = 'init_pop_size_4',
                                                                   label = "Initial Population Size",
                                                                   value = 100, min = 2, max = 1000),
                                                       shinyBS::bsTooltip('init_pop_size_4',
                                                                          "Population size at the start of the simulation"),

                                                       sliderInput(inputId = 'num_gen_optim_c',
                                                                   label = "Number of Generations",
                                                                   value = 20, min = 2, max = 100),
                                                       shinyBS::bsTooltip("num_gen_optim_c",
                                                                          "Total number of generations simulated"),

                                                       sliderInput(inputId = 'total_put',
                                                                   label = 'Put: Total number of individuals',
                                                                   value = 100, min = 0, max = 1000),
                                                       shinyBS::bsTooltip('total_put',
                                                                          "Total number of individuals added, summed over all generations"),

                                                       sliderInput(inputId = 'total_pull',
                                                                   label = 'Pull: Total number of individuals',
                                                                   value = 100, min = 0, max = 1000),
                                                       shinyBS::bsTooltip('total_pull',
                                                                          "Total number of individuals removed, summed over all generations"),

                                                       sliderInput(inputId = 'init_frac_optim_complex',
                                                                   label = 'Starting fraction of focal Ancestry',
                                                                   value = 0.8, min = 0, max = 1),
                                                       shinyBS::bsTooltip('init_frac_optim_complex',
                                                                          "Initial frequency of focal ancestry"),

                                                       selectInput("density_model_3",
                                                                   "Density dependence:",
                                                                   c("Weak", "Strong", "Manual")),
                                                       shinyBS::bsTooltip("density_model_3",
                                                                          "The user can pick between two pre-defined parameter sets that implement weak or strong density dependence. Alternatively the user can modify parameters manually in the corresponding tab."),

                                                       checkboxInput("model_used_c",
                                                                     "Use explicit recombination",
                                                                     value = FALSE),
                                                       shinyBS::bsTooltip("model_used_c",
                                                                          "When unchecked, a simplified genetic model is used. When checked, explicit recombination is modeled"),
                                                       downloadButton("download_gen3", label = "Download Genetics"),
                                                       shinyBS::bsTooltip("download_gen3",
                                                                          "Download local ancestry information of the last generation as a tibble"),
                                                       downloadButton("download_res3", label = "Download results"),
                                                       shinyBS::bsTooltip("download_res3",
                                                                          "Download results as text file")

                                      )
                             ),
                             tabPanel("Advanced", value = 2,
                                      numericInput(inputId = 'K',
                                                   label = "Carrying Capacity of ecosystem",
                                                   value = 400, min = 2, max = 1000, step = 50),
                                      shinyBS::bsTooltip("K",
                                                         "Carrying Capacity of the ecosystem, e.g. the maximum number of individuals that can be sustained by the ecosystem"),

                                      numericInput(inputId = 'f_n_r',
                                                   label = "Breeding Risk Female",
                                                   value = 0.2, min = 0, max = 1, step = 0.01),
                                      shinyBS::bsTooltip("f_n_r",
                                                         "Breeding risk for females, caused by for instance increased predation in defending offspring"),

                                      numericInput(inputId = 'm_n_r',
                                                   label = "Breeding Risk Male",
                                                   value = 0.0, min = 0, max = 1, step = 0.01),
                                      shinyBS::bsTooltip("m_n_r",
                                                         "Breeding risk for males, caused by for instance increased predation in defending offspring"),

                                      numericInput(inputId = 'nest_succes_rate',
                                                   label = "Reproduction Succes Rate",
                                                   value = 0.387, min = 0, max = 1, step = 0.01),
                                      shinyBS::bsTooltip("nest_succes_rate",
                                                         "Success rate of producing offspring per mating"),

                                      numericInput(inputId = 'morgan',
                                                   label = "Size of Genome (in Morgan)",
                                                   value = 1.0, min = 0, max = 3, step = 0.1),
                                      shinyBS::bsTooltip("morgan",
                                                         "Size of the modeled chromosome in Morgan, this influences the expected number of crossover events per meiosis."),

                                      numericInput(inputId = 'num_repl',
                                                   label = "Number of replicates",
                                                   value = 1, min = 1, max = 10, step = 1),
                                      shinyBS::bsTooltip("num_repl",
                                                         "Multiple replicates using different random seeds are shown"),

                                      numericInput(inputId = 'max_age',
                                                   label = "Maximum Age",
                                                   value = 6, min = 1, max = 20, step = 1),
                                      shinyBS::bsTooltip("max_age",
                                                         "Maximum age an individual can obtain. This is modeled as a hard upper limit that individuals can not exceed. This is mainly usefull to avoid individuals with extreme old-age."),

                                      numericInput(inputId = 'clutch_size',
                                                   label = "Number of offspring",
                                                   value = 6, min = 1, max = 20, step = 1),
                                      shinyBS::bsTooltip("clutch_size",
                                                         "Total number of offspring generated per mated female"),

                                      numericInput(inputId = 'clutch_sd',
                                                   label = "SD Number of Offspring",
                                                   value = 1, min = 0, max = 2, step = 0.1),
                                      shinyBS::bsTooltip("clutch_sd",
                                                         "Standard deviation of number of offspring generated per mated female"),

                                      numericInput(inputId = 'sex_ratio_put',
                                                   label = "Sex Ratio of Put individuals (males / females)",
                                                   value = 0.5, min = 0, max = 1, step = 0.05),
                                      shinyBS::bsTooltip("sex_ratio_put",
                                                         "Sex ratio of individuals added, where values > 0.5 indicate a male biased sex ratio, and values < 0.5 indicate a female biased sex ratio"),

                                      numericInput(inputId = 'sex_ratio_pull',
                                                   label = "Sex Ratio of Pulled individuals (males / females)",
                                                   value = 0.5, min = 0, max = 1, step = 0.05),
                                      shinyBS::bsTooltip("sex_ratio_pull",
                                                         "Sex ratio of individuals removed, where values > 0.5 indicate a male biased sex ratio, and values < 0.5 indicate a female biased sex ratio"),

                                      numericInput(inputId = 'sex_ratio_offspring',
                                                   label = "Sex Ratio of offspring (males / females)",
                                                   value = 0.5, min = 0, max = 1, step = 0.05),
                                      shinyBS::bsTooltip("sex_ratio_offspring",
                                                         "Sex ratio of born offspring, where values > 0.5 indicate a male biased sex ratio, and values < 0.5 indicate a female biased sex ratio"),

                                      numericInput(inputId = 'target_frequency',
                                                   label = "Target Frequency (used in optimization)",
                                                   value = 0.999, min = 0, max = 1, step = 0.01),
                                      shinyBS::bsTooltip("target_frequency",
                                                         "The optimizer tries to optimize pull and or put to reach this frequency after the set number of generations."),
                             ),
                             tabPanel("Density Dependence", value = 3,
                                      numericInput(inputId = 'smin',
                                                   label = "Minimum Survival Rate",
                                                   value = 0.5, min = 0, max = 1, step = 0.05),
                                      shinyBS::bsTooltip("smin",
                                                         "Minimum survival rate, e.g. the survival rate even at extremely high densities does not drop below this value"),

                                      numericInput(inputId = 'smax',
                                                   label = "Maximum Survival Rate",
                                                   value = 0.9, min = 0, max = 1, step = 0.05),
                                      shinyBS::bsTooltip("smax",
                                                         "Maximum survival rate, e.g. the survival rate even at extremely low densities does not exceed this value"),

                                      numericInput(inputId = 'b',
                                                   label = "Steepness Survival curve",
                                                   value = -2, min = -3, max = 0, step = 0.05),
                                      shinyBS::bsTooltip("b",
                                                         "Steepness of the survival curve, where negative values indicate decreasing survival with increasing density, and positive values indicate increasing survival with density (this typically causes the simulation to grind to a halt, because it drives population explosion, and should be avoided)."),


                                      numericInput(inputId = 'p',
                                                   label = "Density of maximum steepness Survival curve",
                                                   value = 0.5, min = 0, max = 2, step = 0.05),
                                      shinyBS::bsTooltip("p",
                                                         "Density at which the survival curve shows maximum steepness, or in other words, density at which survival is exactly (smax + smin) / 2"),

                             )
                 )
    ),
    mainPanel("",
              tabsetPanel(type = "tabs", id = "tabs1",
                          tabPanel("Simulation", value = 1,
                                   plotOutput("simple_plots")),
                          tabPanel("Static Optimization", value = 3,
                                   plotOutput('Optim_simple_plots'),
                                   div(class = "myclass",
                                       verbatimTextOutput("selected_var"))
                          ),
                          tabPanel("Adaptive Optimization", value = 4,
                                   plotOutput('Optim_complex_plots'),
                                   div(class = "myclass",
                                       verbatimTextOutput("complex_text_output")
                                   ))
              )
    )
  )
)

server <- function(input, output, session) {

  simple_data <- reactive({
    simRestore::simulate_policy(initial_population_size = input$init_pop_size_1,
                                reproduction_success_rate = input$nest_succes_rate,
                                reproductive_risk = c(input$f_n_r, input$m_n_r),
                                num_generations = input$num_gen_simple,
                                K = input$K,
                                pull = input$pull,
                                put = input$put,
                                num_replicates  = input$num_repl,
                                starting_freq = input$init_frac_simple,
                                morgan = input$morgan,
                                max_age = input$max_age,
                                mean_number_of_offspring = input$clutch_size,
                                sd_number_of_offspring = input$clutch_sd,
                                genetic_model =
                                  ifelse(input$model_used_single == 1,
                                                       "junctions",
                                                       "simplified"),
                                smin = ifelse(input$density_model == "Manual",
                                              input$smin, 0.5),
                                smax = ifelse(input$density_model == "Manual",
                                              input$smax, 0.9),
                                b = ifelse(input$density_model == "Manual",
                                           input$b,
                                           ifelse(input$density_model == "Strong",
                                                  -5, -2)),
                             p = ifelse(input$density_model == "Manual",
                                        input$p,
                                        ifelse(input$density_model == "Strong",
                                               0.45, 0.5)),
                                sex_ratio_put = input$sex_ratio_put,
                                sex_ratio_pull = input$sex_ratio_pull,
                                sex_ratio_offspring = input$sex_ratio_offspring,
                             return_genetics = TRUE)
    })

  output$simple_plots <- renderPlot({
    to_plot <- simple_data()
    data_storage <<- to_plot

    p1 <- to_plot$results %>%
      ggplot(aes(x = t, y = freq_focal_ancestry, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Average focal ancestry") +
      ylim(0, 1) +
      theme(legend.position = "top")

    focal_y <- 1.03 * tail(to_plot$results$freq_focal_ancestry, 1)
    if (round(tail(to_plot$results$freq_focal_ancestry, 1), 2) >= 0.99) {
      focal_y <- 0.95 * tail(to_plot$results$freq_focal_ancestry, 1)
    }

    p1 <- p1 +
      annotate("text", x = max(to_plot$results$t), y = focal_y,
               label = round(tail(to_plot$results$freq_focal_ancestry, 1), 2),
               hjust = 1)

    p2 <- ggplot(to_plot$results, aes(x = t,
                                      y = num_individuals,
                                      group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Total number of individuals")

    p3 <- to_plot$results %>%
      dplyr::mutate("Males" = num_males) %>%
      dplyr::mutate("Females" = num_females) %>%
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

  ########  OPTIMIZATION #####################################################

  optim_data_static <-  reactive({
    get_optim_data_static(initial_population_size = input$init_pop_size_3,
                          reproduction_success_rate = input$nest_succes_rate,
                          reproductive_risk = c(input$f_n_r, input$m_n_r),
                          num_generations = input$num_gen_optim_s,
                          K = input$K,
                          num_replicates = input$num_repl,
                          target_frequency = input$target_frequency,
                          optim_choice = input$optim_choice,
                          morgan = input$morgan,
                          starting_freq = input$init_frac_optim,
                          use_complex_model = input$model_used_s,
                          max_age = input$max_age,
                          mean_number_of_offspring = input$clutch_size,
                          sd_number_of_offspring = input$clutch_sd,
                          smin = ifelse(input$density_model_2 == "Manual",
                                        input$smin, 0.5),
                          smax = ifelse(input$density_model_2 == "Manual",
                                        input$smax, 0.9),
                          b = ifelse(input$density_model_2 == "Manual",
                                     input$b,
                                     ifelse(input$density_model_2 == "Strong",
                                            -5, -2)),
                          p = ifelse(input$density_model_2 == "Manual",
                                     input$p,
                                     ifelse(input$density_model_2 == "Strong",
                                            0.45, 0.5)),
                          sex_ratio_put = input$sex_ratio_put,
                          sex_ratio_pull = input$sex_ratio_pull,
                          sex_ratio_offspring = input$sex_ratio_offspring)
  })

  optim_data_complex <-  reactive({
    get_optim_data_adaptive(initial_population_size = input$init_pop_size_4,
                           reproduction_success_rate = input$nest_succes_rate,
                           reproductive_risk = c(input$f_n_r, input$m_n_r),
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
                           mean_number_of_offspring = input$clutch_size,
                           sd_number_of_offspring = input$clutch_sd,
                           smin = ifelse(input$density_model_3 == "Manual",
                                         input$smin, 0.5),
                           smax = ifelse(input$density_model_3 == "Manual",
                                         input$smax, 0.9),
                           b = ifelse(input$density_model_3 == "Manual",
                                      input$b,
                                      ifelse(input$density_model_3 == "Strong",
                                             -5, -2)),
                           p = ifelse(input$density_model_3 == "Manual",
                                      input$p,
                                      ifelse(input$density_model_3 == "Strong",
                                             0.45, 0.5)),
                           sex_ratio_put = input$sex_ratio_put,
                           sex_ratio_pull = input$sex_ratio_pull,
                           sex_ratio_offspring = input$sex_ratio_offspring)
  })


  output$Optim_simple_plots <- renderPlot({
    to_plot <- optim_data_static()
    data_storage <<- to_plot

    final_freq <- round(to_plot$final_freq, digits = 3)

    for_render_text <- c()
    for_render_text <- c(for_render_text,
                         " Target frequency was: ",
                         input$target_frequency, "\n")
    for_render_text <- c(for_render_text,
                         " Final frequency was: ",
                         final_freq, "\n")

    if (final_freq >= input$target_frequency) {
      for_render_text <- c(for_render_text,
                           "Target frequency was reached\n")
    } else {
      for_render_text <- c(for_render_text,
                           "Target frequency was NOT reached\n")
    }

    if (length(input$optim_choice) == 1) {
      if (input$optim_choice == "Pull") {
        for_render_text <- c(for_render_text, "\n",
                             "Advice is to pull ",
                             round(to_plot$pull),
                             " individuals per generation")
      }
      if (input$optim_choice == "Put") {
        for_render_text <- c(for_render_text, "\n",
                             "Advice is to put ",
                             round(to_plot$put),
                             " individuals per generation")
      }
    }
    if (length(input$optim_choice) == 2) {

      for_render_text <- c(for_render_text, "\n",
                           "Advice is to put ",
                           round(to_plot$put),
                           " individuals per generation\n",
                           "         and pull ", round(to_plot$pull),
                           " individuals per generation")
    }

    output$selected_var <- renderText({ for_render_text })

    p1 <- to_plot$results %>%
      ggplot(aes(x = t, y = freq_focal_ancestry, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Average focal ancestry") +
      ylim(0, 1) +
      theme(legend.position = "top")

    focal_y <- 1.03 * tail(to_plot$results$freq_focal_ancestry, 1)
    if (round(tail(to_plot$results$freq_focal_ancestry, 1), 2) >= 0.99) {
      focal_y <- 0.95 * tail(to_plot$results$freq_focal_ancestry, 1)
    }

    p1 <- p1 +
      annotate("text", x = max(to_plot$results$t), y = focal_y,
               label = round(tail(to_plot$results$freq_focal_ancestry, 1), 2),
               hjust = 1)

    p2 <- ggplot(to_plot$results,
                 aes(x = t, y = num_individuals, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Total number of individuals")

    p3 <- to_plot$results %>%
      dplyr::mutate("Males" = num_males) %>%
      dplyr::mutate("Females" = num_females) %>%
      ggplot(aes(x = t, group = replicate)) +
      geom_line(aes(y = Males, color = "Males")) +
      geom_line(aes(y = Females, group = replicate, color = "Females")) +
      labs(x = "Generation",
           y = "Number of individuals",
           color = "Sex") +
      theme(legend.position = "top")

    p4 <-  tidyr::gather(to_plot$curve, key = "type", value = "number", -t) %>%
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
    data_storage <<- to_plot
    for_text <- to_plot$curve
    # tibble with t, pull, put
    final_freq <- round(to_plot$final_freq, 3)
    for_render_text <- c()
    for_render_text <- c(for_render_text,
                         " Target frequency was: ",
                         input$target_frequency, "\n")
    for_render_text <- c(for_render_text,
                         "Final frequency was: ", final_freq, "\n")
    if (final_freq >= input$target_frequency) {
      for_render_text <- c(for_render_text,
                           "Target frequency was reached\n")
    } else {
      for_render_text <- c(for_render_text,
                           "Target frequency was NOT reached\n")
    }

    for_render_text <- c(for_render_text, "Advice:", "\n")
    if (input$total_put > 0 && input$total_pull > 0) {
      for_render_text <- c(for_render_text,
                           c("Generation", "\t", "Put", "\t", "Pull", "\n"))
    }
    if (input$total_put == 0 && input$total_pull > 0) {
      for_render_text <- c(for_render_text,
                           c("Generation", "\t", "Pull", "\n"))
    }

    if (input$total_put > 0 && input$total_pull == 0) {
      for_render_text <- c(for_render_text,
                           c("Generation", "\t", "Put", "\n"))
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

    output$complex_text_output <- renderText({for_render_text})



    p1 <- to_plot$results %>%
      ggplot(aes(x = t, y = freq_focal_ancestry, group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Average focal ancestry") +
      ylim(0, 1) +
      theme(legend.position = "top")

    focal_y <- 1.03 * tail(to_plot$results$freq_focal_ancestry, 1)
    if (round(tail(to_plot$results$freq_focal_ancestry, 1), 2) >= 0.99) {
      focal_y <- 0.95 * tail(to_plot$results$freq_focal_ancestry, 1)
    }

    p1 <- p1 +
      annotate("text", x = max(to_plot$results$t), y = focal_y,
               label = round(tail(to_plot$results$freq_focal_ancestry, 1), 2),
               hjust = 1)

    p2 <- ggplot(to_plot$results, aes(x = t,
                                      y = num_individuals,
                                      group = replicate)) +
      geom_line() +
      xlab("Number of generations") +
      ylab("Total number of individuals")

    p3 <- to_plot$results %>%
      dplyr::mutate("Males" = num_males) %>%
      dplyr::mutate("Females" = num_females) %>%
      ggplot(aes(x = t, group = replicate)) +
      geom_line(aes(y = Males, color = "Males")) +
      geom_line(aes(y = Females, group = replicate, color = "Females")) +
      labs(x = "Generation",
           y = "Number of individuals",
           color = "Sex") +
      theme(legend.position = "top")

    p4 <-  tidyr::gather(to_plot$curve, key = "type", value = "number", -t) %>%
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

  output$download_res1 <- downloadHandler(
    filename = function() {
      paste0("dataset_", Sys.Date(), ".txt")
    },
    content = function(file) {
      # stored_data <- read.table(input$data_for_download)
      write.table(data_storage$results, file, quote = FALSE)
    }
  )

  output$download_gen1 <- downloadHandler(
    filename = function() {
      paste0("genetics_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(data_storage$genetics, file, quote = FALSE)
    }
  )

  output$download_res2 <- downloadHandler(
    filename = function() {
      paste0("dataset_", Sys.Date(), ".txt")
    },
    content = function(file) {
      # stored_data <- read.table(input$data_for_download)
      write.table(data_storage$results, file, quote = FALSE)
    }
  )

  output$download_gen2 <- downloadHandler(
    filename = function() {
      paste0("genetics_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(data_storage$genetics, file, quote = FALSE)
    }
  )

  output$download_res3 <- downloadHandler(
    filename = function() {
      paste0("dataset_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(data_storage$results, file, quote = FALSE)
    }
  )

  output$download_gen3 <- downloadHandler(
    filename = function() {
      paste0("genetics_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(data_storage$genetics, file, quote = FALSE)
    }
  )
}

get_optim_data_static <- function(initial_population_size,
                                  reproduction_success_rate,
                                  reproductive_risk,
                                  num_generations,
                                  K,
                                  num_replicates,
                                  target_frequency,
                                  optim_choice,
                                  morgan,
                                  starting_freq,
                                  use_complex_model,
                                  max_age,
                                  mean_number_of_offspring,
                                  sd_number_of_offspring,
                                  smin,
                                  smax,
                                  b,
                                  p,
                                  sex_ratio_put,
                                  sex_ratio_pull,
                                  sex_ratio_offspring) {

  opt_pull <- FALSE
  opt_put  <- FALSE

  if (length(optim_choice) == 2) {
    opt_pull <- TRUE
    opt_put <- TRUE
  } else {
    if (length(optim_choice) == 1) {
      if (optim_choice == "Put") {
        opt_put <- TRUE
        opt_pull <- FALSE
      }
      if (optim_choice == "Pull") {
        opt_pull <- TRUE
        opt_put <- FALSE
      }
    }
  }

  return(simRestore::optimize_static(initial_population_size =
                                       initial_population_size,
                                     reproduction_success_rate = reproduction_success_rate,
                                     reproductive_risk = reproductive_risk,
                                     num_generations = num_generations,
                                     K = K,
                                     num_replicates = num_replicates,
                                     target_frequency = target_frequency,
                                     optimize_pull = opt_pull,
                                     optimize_put = opt_put,
                                     morgan = morgan,
                                     starting_freq = starting_freq,
                                     max_age = max_age,
                                     mean_number_of_offspring = mean_number_of_offspring,
                                     sd_number_of_offspring = sd_number_of_offspring,
                                     genetic_model =
                                       ifelse(use_complex_model == TRUE,
                                              "junctions",
                                              "simplified"),
                                     smin = smin,
                                     smax = smax,
                                     b = b,
                                     p = p,
                                     sex_ratio_put = sex_ratio_put,
                                     sex_ratio_pull = sex_ratio_pull,
                                     sex_ratio_offspring = sex_ratio_offspring,
                                     verbose = FALSE,
                                     return_genetics = TRUE))
}

get_optim_data_adaptive <- function(initial_population_size,
                                   reproduction_success_rate,
                                   reproductive_risk,
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
                                   mean_number_of_offspring,
                                   sd_number_of_offspring,
                                   smin,
                                   smax,
                                   b,
                                   p,
                                   sex_ratio_put,
                                   sex_ratio_pull,
                                   sex_ratio_offspring) {

  return(simRestore::optimize_adaptive(
    initial_population_size = initial_population_size,
    reproduction_success_rate = reproduction_success_rate,
    reproductive_risk = reproductive_risk,
    num_generations = num_generations,
    K = K,
    num_replicates = num_replicates,
    target_frequency = target_frequency,
    optimize_pull = total_pull,
    optimize_put = total_put,
    morgan = morgan,
    starting_freq = starting_freq,
    max_age = max_age,
    mean_number_of_offspring = mean_number_of_offspring,
    sd_number_of_offspring = sd_number_of_offspring,
    genetic_model =
      ifelse(use_complex_model == TRUE,
             "junctions",
             "simplified"),
    smin = smin,
    smax = smax,
    b = b,
    p = p,
    sex_ratio_put = sex_ratio_put,
    sex_ratio_pull = sex_ratio_pull,
    sex_ratio_offspring = sex_ratio_offspring,
    verbose = FALSE,
    return_genetics = TRUE))
}

shinyApp(ui = ui, server = server)
