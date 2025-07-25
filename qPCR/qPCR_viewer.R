# Title: Interactive qPCR Analysis App (Command-Line Version)
# Description: This Shiny app visualizes qPCR data from a single Excel file
# provided as a command-line argument. It displays four interactive plots
# with persistent zoom, advanced filtering, and detailed hover information.
# This version suppresses console output and launches automatically in a browser.

# --- 1. Load Required Libraries (Silently) ---
# Ensure you have these packages installed:
# install.packages(c("shiny", "plotly", "dplyr", "readr", "bslib", "scales", "tidyr", "readxl"))

suppressPackageStartupMessages({
  library(shiny)
  library(plotly)
  library(dplyr)
  library(readr)
  library(bslib)
  library(scales)
  library(tidyr)
  library(readxl)
})

# --- 2. Helper Function for Dynamic Header Finding ---
find_header_row <- function(file, sheet_name) {
  df_preview <- tryCatch({
    suppressMessages(read_excel(file, sheet = sheet_name, n_max = 60, col_names = FALSE))
  }, error = function(e) {
    stop(paste("Could not read sheet:", sheet_name, "in file:", file, ". Error:", e$message))
    return(NULL)
  })
  if (is.null(df_preview)) return(47) # Should not be reached due to stop()
  for (i in 1:nrow(df_preview)) {
    row_values <- as.character(df_preview[i, ])
    if ("Well" %in% row_values || "Well Position" %in% row_values) {
      return(i - 1)
    }
  }
  warning(paste("Header row not found in sheet:", sheet_name, ". Defaulting to 47."))
  return(47)
}


# --- 3. Data Loading and Preprocessing Function ---
load_qpcr_data <- function(excel_file) {
  
  if (!file.exists(excel_file)) {
      stop(paste("CRITICAL ERROR: The file '", excel_file, "' was not found.", sep=""))
  }
  
  # --- Dynamically find the number of rows to skip for each sheet ---
  setup_skip <- find_header_row(excel_file, "Sample Setup")
  results_skip <- find_header_row(excel_file, "Results")
  melt_raw_skip <- find_header_row(excel_file, "Melt Curve Raw Data")
  amp_skip <- find_header_row(excel_file, "Amplification Data")
  
  # --- Load data from Excel worksheets ---
  setup_df <- read_excel(excel_file, sheet = "Sample Setup", skip = setup_skip, col_types = "text")
  results_df <- read_excel(excel_file, sheet = "Results", skip = results_skip, col_types = "text")
  melt_raw_df <- read_excel(excel_file, sheet = "Melt Curve Raw Data", skip = melt_raw_skip, col_types = "text")
  amp_df <- read_excel(excel_file, sheet = "Amplification Data", skip = amp_skip, col_types = "text")
  
  # --- Create Master Information Table ---
  setup_map <- setup_df %>%
    select(`Well`, `Well Position`, `Sample Name`, `Target Name`) %>%
    rename(
      WellNum = `Well`,
      WellPosition = `Well Position`,
      SampleName = `Sample Name`,
      TargetName = `Target Name`
    ) %>%
    mutate(
      WellNum = as.integer(WellNum),
      WellPosition = as.character(WellPosition)
    )

  results_info <- results_df %>%
    select(`Well Position`, `CT`, `Tm1`, `Tm2`) %>%
    rename(WellPosition = `Well Position`) %>%
    mutate(
      WellPosition = as.character(WellPosition),
      CT = as.numeric(CT),
      Tm1 = as.numeric(Tm1),
      Tm2 = as.numeric(Tm2)
    )
    
  well_info <- left_join(setup_map, results_info, by = "WellPosition") %>%
    filter(!is.na(SampleName) & SampleName != "") %>%
    mutate(ID = paste(SampleName, TargetName, WellPosition, sep = "-")) %>%
    mutate(
      hover_text = paste0(
        "<b>Well:</b> ", WellPosition, "<br>",
        "<b>Sample ID:</b> ", SampleName, "<br>",
        "<b>Target:</b> ", TargetName, "<br>",
        "<b>CT:</b> ", ifelse(is.na(CT), "N/A", round(CT, 2)), "<br>",
        "<b>Tm1:</b> ", ifelse(is.na(Tm1), "N/A", round(Tm1, 2)), "<br>",
        "<b>Tm2:</b> ", ifelse(is.na(Tm2), "N/A", round(Tm2, 2))
      )
    )

  # --- Clean and Join Raw Data ---
  melt_data_final <- melt_raw_df %>%
    rename(WellNum = `Well`) %>%
    mutate(
        WellNum = as.integer(WellNum),
        Temperature = as.numeric(Temperature),
        Fluorescence = as.numeric(Fluorescence),
        Derivative = as.numeric(Derivative)
    ) %>%
    left_join(well_info, by = "WellNum") %>%
    filter(!is.na(ID))

  amp_data_final <- amp_df %>%
    rename(WellNum = `Well`) %>%
    mutate(
        WellNum = as.integer(WellNum),
        Cycle = as.numeric(Cycle),
        Rn = as.numeric(Rn),
        DeltaRn = as.numeric(`Delta Rn`)
    ) %>%
    left_join(well_info, by = "WellNum") %>%
    filter(!is.na(ID))

  # --- Prepare data for CT point markers ---
  ct_points <- amp_data_final %>%
    filter(!is.na(CT)) %>%
    mutate(Cycle_match = round(CT)) %>%
    group_by(ID) %>%
    filter(Cycle == Cycle_match) %>%
    slice(1) %>%
    ungroup() %>%
    select(ID, CT, Rn, DeltaRn, hover_text)

  # --- Get unique values for filters and color palette ---
  unique_targets <- sort(unique(well_info$TargetName))
  unique_samples <- sort(unique(well_info$SampleName))
  unique_wells <- sort(unique(well_info$WellPosition))
  
  all_ids <- sort(unique(well_info$ID))
  id_color_palette <- scales::hue_pal()(length(all_ids))
  names(id_color_palette) <- all_ids
  
  # Return all processed data in a list
  list(
    well_info = well_info,
    melt_data_final = melt_data_final,
    amp_data_final = amp_data_final,
    ct_points = ct_points,
    unique_targets = unique_targets,
    unique_samples = unique_samples,
    unique_wells = unique_wells,
    id_color_palette = id_color_palette
  )
}


# --- 4. UI Function ---
# The UI now needs to be a function to accept the loaded data
create_ui <- function(unique_targets, unique_samples) {
  fluidPage(
    theme = bslib::bs_theme(version = 5, bootswatch = "cerulean"),
    titlePanel("Interactive qPCR Analysis"),
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Controls"),
        p("Filter the data to display. Double-click a plot to reset its zoom."),
        
        h5("Filter by Target Name"),
        actionButton("select_all_targets", "Select All"),
        actionButton("deselect_all_targets", "Deselect All"),
        checkboxGroupInput("selected_targets", NULL, choices = unique_targets, selected = unique_targets),
        
        hr(),
        
        h5("Filter by Sample ID"),
        actionButton("select_all_samples", "Select All"),
        actionButton("deselect_all_samples", "Deselect All"),
        checkboxGroupInput("selected_samples", NULL, choices = unique_samples, selected = unique_samples),
        
        hr(),
        
        h5("Filter by Well ID"),
        actionButton("select_all_wells", "Select All"),
        actionButton("deselect_all_wells", "Deselect All"),
        fluidRow(
          column(6, checkboxGroupInput("selected_wells_1", NULL, choices = list())),
          column(6, checkboxGroupInput("selected_wells_2", NULL, choices = list()))
        )
      ),
      
      mainPanel(
        width = 9,
        tabsetPanel(
          type = "tabs",
          tabPanel("Fluorescence vs. Temp", plotlyOutput("fluorescence_plot", height = "800px")),
          tabPanel("Derivative vs. Temp", plotlyOutput("derivative_plot", height = "800px")),
          tabPanel("Rn vs. Cycle", plotlyOutput("rn_plot", height = "800px")),
          tabPanel("Delta Rn vs. Cycle", plotlyOutput("delta_rn_plot", height = "800px"))
        )
      )
    )
  )
}


# --- 5. Server Function ---
# The server logic is also wrapped in a function to receive the data
create_server <- function(data) {
  
  # Unpack the data list for easier access
  well_info <- data$well_info
  melt_data_final <- data$melt_data_final
  amp_data_final <- data$amp_data_final
  ct_points <- data$ct_points
  unique_targets <- data$unique_targets
  unique_samples <- data$unique_samples
  unique_wells <- data$unique_wells
  id_color_palette <- data$id_color_palette
  
  # The actual server function
  function(input, output, session) {
    
    observe({
      req(unique_wells)
      wells <- unique_wells
      n <- length(wells)
      half <- ceiling(n / 2)
      updateCheckboxGroupInput(session, "selected_wells_1", choices = wells[1:half], selected = wells[1:half])
      updateCheckboxGroupInput(session, "selected_wells_2", choices = wells[(half + 1):n], selected = wells[(half + 1):n])
    })
    
    observeEvent(input$select_all_targets, { updateCheckboxGroupInput(session, "selected_targets", selected = unique_targets) })
    observeEvent(input$deselect_all_targets, { updateCheckboxGroupInput(session, "selected_targets", selected = character(0)) })
    observeEvent(input$select_all_samples, { updateCheckboxGroupInput(session, "selected_samples", selected = unique_samples) })
    observeEvent(input$deselect_all_samples, { updateCheckboxGroupInput(session, "selected_samples", selected = character(0)) })
    observeEvent(input$select_all_wells, {
      updateCheckboxGroupInput(session, "selected_wells_1", selected = unique_wells)
      updateCheckboxGroupInput(session, "selected_wells_2", selected = unique_wells)
    })
    observeEvent(input$deselect_all_wells, {
      updateCheckboxGroupInput(session, "selected_wells_1", selected = character(0))
      updateCheckboxGroupInput(session, "selected_wells_2", selected = character(0))
    })

    plot_zoom <- reactiveValues(
      fluorescence_plot = NULL,
      derivative_plot = NULL,
      rn_plot = NULL,
      delta_rn_plot = NULL
    )

    create_zoom_observer <- function(plot_name) {
      observeEvent(event_data("plotly_relayout", source = plot_name), {
        d <- event_data("plotly_relayout", source = plot_name)
        is_reset <- any(grepl("autorange", names(d)))
        if (is_reset) {
          plot_zoom[[plot_name]] <- NULL
        } else {
          plot_zoom[[plot_name]] <- list(
            xaxis = list(range = c(d$`xaxis.range[0]`, d$`xaxis.range[1]`)),
            yaxis = list(range = c(d$`yaxis.range[0]`, d$`yaxis.range[1]`))
          )
        }
      })
    }
    
    create_zoom_observer("fluorescence_plot")
    create_zoom_observer("derivative_plot")
    create_zoom_observer("rn_plot")
    create_zoom_observer("delta_rn_plot")
    
    base_filtered_ids <- reactive({
        validate(need(nrow(well_info) > 0, "Data could not be loaded. Please check R console for errors."))
        selected_wells <- c(input$selected_wells_1, input$selected_wells_2)
        req(input$selected_targets, input$selected_samples, selected_wells)
        well_info %>%
          filter(
            TargetName %in% input$selected_targets, 
            SampleName %in% input$selected_samples,
            WellPosition %in% selected_wells
          ) %>%
          pull(ID)
    })
    
    filtered_melt_data <- reactive({ melt_data_final %>% filter(ID %in% base_filtered_ids()) })
    filtered_amp_data <- reactive({ amp_data_final %>% filter(ID %in% base_filtered_ids()) })
    filtered_ct_points <- reactive({ ct_points %>% filter(ID %in% base_filtered_ids()) })

    render_zoomable_plot <- function(zoom_info, base_plot_expr) {
      p <- base_plot_expr
      if (!is.null(zoom_info)) {
        p <- p %>% layout(xaxis = zoom_info$xaxis, yaxis = zoom_info$yaxis)
      }
      p
    }
    
    output$fluorescence_plot <- renderPlotly({
      validate(need(nrow(filtered_melt_data()) > 0, "No data for current selection."))
      base_plot <- plot_ly(data = filtered_melt_data(), x = ~Temperature, y = ~Fluorescence, color = ~ID, colors = id_color_palette, type = 'scatter', mode = 'lines', 
              source = "fluorescence_plot", hoverinfo = 'text', text = ~hover_text) %>%
        # **FIX**: Register the event to suppress the warning
        event_register('plotly_relayout') %>%
        layout(title = "Melt Curve: Fluorescence vs. Temperature", xaxis = list(title = "Temperature (°C)"), yaxis = list(title = "Fluorescence"), legend = list(title=list(text='<b>ID</b>')))
      render_zoomable_plot(plot_zoom$fluorescence_plot, base_plot)
    })
    
    output$derivative_plot <- renderPlotly({
      validate(need(nrow(filtered_melt_data()) > 0, "No data for current selection."))
      base_plot <- plot_ly(data = filtered_melt_data(), x = ~Temperature, y = ~Derivative, color = ~ID, colors = id_color_palette, type = 'scatter', mode = 'lines', 
              source = "derivative_plot", hoverinfo = 'text', text = ~hover_text) %>%
        # **FIX**: Register the event to suppress the warning
        event_register('plotly_relayout') %>%
        layout(title = "Melt Curve: Derivative vs. Temperature", xaxis = list(title = "Temperature (°C)"), yaxis = list(title = "Derivative (-dF/dT)"), legend = list(title=list(text='<b>ID</b>')))
      render_zoomable_plot(plot_zoom$derivative_plot, base_plot)
    })

    output$rn_plot <- renderPlotly({
      validate(need(nrow(filtered_amp_data()) > 0, "No data for current selection."))
      base_plot <- plot_ly(data = filtered_amp_data(), x = ~Cycle, y = ~Rn, color = ~ID, colors = id_color_palette, type = 'scatter', mode = 'lines', legendgroup = ~ID, 
              source = "rn_plot", hoverinfo = 'text', text = ~hover_text) %>%
        add_markers(data = filtered_ct_points(), x = ~CT, y = ~Rn, color = ~ID, colors = id_color_palette, legendgroup = ~ID, showlegend = FALSE,
                    marker = list(size = 8, symbol = 'diamond', line = list(color = 'black', width = 1)),
                    hoverinfo = 'text', text = ~hover_text) %>%
        # **FIX**: Register the event to suppress the warning
        event_register('plotly_relayout') %>%
        layout(title = "Amplification: Rn vs. Cycle", xaxis = list(title = "Cycle"), yaxis = list(title = "Rn (Normalized Reporter)"), legend = list(title=list(text='<b>ID</b>')))
      render_zoomable_plot(plot_zoom$rn_plot, base_plot)
    })

    output$delta_rn_plot <- renderPlotly({
      validate(need(nrow(filtered_amp_data()) > 0, "No data for current selection."))
      base_plot <- plot_ly(data = filtered_amp_data(), x = ~Cycle, y = ~DeltaRn, color = ~ID, colors = id_color_palette, type = 'scatter', mode = 'lines', legendgroup = ~ID, 
              source = "delta_rn_plot", hoverinfo = 'text', text = ~hover_text) %>%
        add_markers(data = filtered_ct_points(), x = ~CT, y = ~DeltaRn, color = ~ID, colors = id_color_palette, legendgroup = ~ID, showlegend = FALSE,
                    marker = list(size = 8, symbol = 'diamond', line = list(color = 'black', width = 1)),
                    hoverinfo = 'text', text = ~hover_text) %>%
        # **FIX**: Register the event to suppress the warning
        event_register('plotly_relayout') %>%
        layout(title = "Amplification: Delta Rn vs. Cycle", xaxis = list(title = "Cycle"), yaxis = list(title = "Delta Rn"), legend = list(title=list(text='<b>ID</b>')))
      render_zoomable_plot(plot_zoom$delta_rn_plot, base_plot)
    })
  }
}

# --- 6. Main Application Logic ---
# This block runs when the script is executed.
main <- function() {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check if a file path was provided
  if (length(args) == 0) {
    stop("Usage: Rscript app.R <path_to_excel_file>", call. = FALSE)
  }
  excel_file_path <- args[1]
  
  # Load and process the data from the specified file
  loaded_data <- tryCatch({
    load_qpcr_data(excel_file_path)
  }, error = function(e) {
    # If data loading fails, print a message and stop.
    stop(paste("Failed to load or process data:", e$message), call. = FALSE)
  })
  
  # Create the UI and Server using the loaded data
  ui <- create_ui(loaded_data$unique_targets, loaded_data$unique_samples)
  server <- create_server(loaded_data)
  
  # Launch the Shiny app, automatically opening a browser window.
  shinyApp(ui, server, options = list(launch.browser = TRUE))
}

# Run the main function
main()
