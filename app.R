# Initialisation ----------------------------------------------------------

library(shiny)
library(ggplot2)
library(readr)

oscorad_to_easi <- read_delim("oscorad_to_easi.csv", ";")
colnames(oscorad_to_easi) <- c("x", "y")

easi_to_oscorad <- read_delim("easi_to_oscorad.csv", ";")
colnames(easi_to_oscorad) <- c("x", "y")

abacus <- read_delim("abacus.csv", ";", escape_double = FALSE, trim_ws = TRUE)

ref_article <- tags$p("The detailed methods are described in",
                      tags$a(href = "https://doi.org/10.1111/bjd.16916",
                             "Hurault et al. (2018), Br. J. Dermat."))

# Prediction from abacus --------------------------------------------------

pred_abacus <- function(abc, x0) {
  # Prediction (weighted mean) of x0 given an abacus
  if (x0 == 0){
    p <- 0
  } else {
    inf <- tail(which(abc$x <= x0), 1)
    sup <- head(which(abc$x >= x0), 1)
    if (inf >= sup){
      p <- abc$y[inf]
    } else {
      dinf <- x0 - abc$x[inf]
      dsup <- abc$x[sup] - x0
      p <- (abc$y[inf] * dsup + abc$y[sup] * dinf) / (dsup + dinf)
      p <- round(p * 2) / 2
    }
  }
  return(p)
}

pred_abacus2 <- function(abc, X) {
  # Prediction for vector
  res <- rep(NA, length(X))
  for (i in 1:length(X)){
    res[i] <- pred_abacus(abc, X[i])
  }
  return(res)
}

# Plot conversion ---------------------------------------------------------

plot_conv <- function(df, xlab) {
  if (xlab == "EASI") {
    p <- ggplot(data = df, aes(x = EASI, y = oSCORAD)) +
      lims(x = c(0, 72), y = c(0, 83))
  } else {
    p <- ggplot(data = df, aes(x = oSCORAD, y = EASI)) +
      lims(x = c(0 ,83), y = c(0, 72))
  }
  p <- p +
    geom_point(size = 2) +
    theme_bw(base_size = 15)
  if (ncol(df) > 2) {
    p <- p +
      geom_text(aes(label = ID), vjust = -1)
  }
  return(p)
}


# Heatmap -----------------------------------------------------------------

pred <- function(abc, lY, x0) {
  # lY: EASI or oSCORAD
  # x0: iga to predict
  inf <- tail(which(abc$IGA <= x0), 1)
  sup <- head(which(abc$IGA >= x0), 1)
  if (length(inf) == 0) {
    p <- x0 / abc$IGA[sup] * abc[sup, lY][[1]]
  } else {
    if (length(sup) == 0){
      p <- x0 / abc$IGA[inf] * abc[inf, lY][[1]]
    } else {
      if (inf >= sup) {
        p <- abc[inf, lY][[1]]
      } else {
        dinf <- x0 - abc$IGA[inf]
        dsup <- abc$IGA[sup] - x0
        p <- (abc[inf, lY][[1]] * dsup + abc[sup, lY][[1]] * dinf) / (dsup + dinf)
      }
    }
  }
  round( p * 2) / 2
}

heatmap <- function(ct, xi) {
  # ct: cutoffs, xi: IGA values to display
  
  # Position of severity label
  xl <- c(0, ct, 5)
  xl <- (xl[-1] + xl[-length(xl)]) / 2 # Between cutoffs
  
  sl <- c("Clear", "Almost\nclear", "Mild", "Moderate", "Severe", "Very\nsevere") # Severity label
  cl <- c(rep("black", 5), "white") # Colour for severity labels
  cl <- rep(c("black", "white"), each = 3) # Colour for severity labels
  
  xi <- sort(c(xi, ct)) # Position for IGA
  xe <- 0 * xi # Position for EASI
  xo <- xe # Position for oSCORAD
  
  for (i in 1:length(xi)) {
    xe[i] <- pred(abacus, "EASI", xi[i])
    xo[i] <- pred(abacus, "oSCORAD", xi[i])
  }
  
  ds <- data.frame(x = c(rep(-.5, 3), rep(xi, 3)),
                   y = c(1.1, 1.2, 1.3, rep(c(1.1, 1.2, 1.3), each = length(xi))),
                   l = c("IGA", "EASI", "oSCORAD", xi, xe, xo))
  
  ds$x[ds$x == 0] <- -.1
  xl[1] <- -.1 # Repositioning text for IGA = 0
  
  ggplot() +
    geom_tile(data = data.frame(IGA = seq(0, 5, .01)),
              aes(x = IGA, fill = IGA, y = .5)) + # x Colour as a function of IGA
    scale_fill_gradientn(colours = c("#FFFFFF", rev(heat.colors(3)), "#000000")) + # Colour gradient
    geom_text(aes(x = xl, y = .5, label = sl),
              colour = cl, fontface = "bold", size = 8) + # Severity labels
    geom_text(data = ds,
              aes(x = x, y = y, label = l),
              fontface = "bold", size = 7)+ # Scale
    geom_segment(aes(x = ct, xend = ct, y = 0, yend = 1),
                 colour = head(cl, 5), size = 2, linetype = "longdash") + # Vertical line cutoff
    labs(x = "", y = "")+ # Remove axis labels
    ylim(c(0, 2))+ xlim(-.5, 5.5) + # Window size (y,x)
    theme_classic(base_size = 15) +
    theme(legend.position = "none",
          axis.text.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_blank(), axis.line.x = element_blank()) # Remove graphical elements
}

# User Interface ----------------------------------------------------------

ui <- fluidPage(
  titlePanel(""),
  tabsetPanel(
    tabPanel("Conversion",
             sidebarLayout(
               sidebarPanel(
                 # Conversion
                 radioButtons(inputId = "conversion",
                              label = "Conversion",
                              choices = c("from oSCORAD to EASI", "from EASI to oSCORAD"),
                              selected = "from oSCORAD to EASI"),
                 
                 tags$hr(), # Horizontal line
                 
                 # One value to convert
                 tags$h3("Value-wise conversion"),
                 textInput(inputId = "inp1",
                           label = "Severity score to be converted",
                           placeholder = "Enter a value"),
                 
                 # Output conversion
                 textOutput("out_conv"),
                 
                 tags$hr(), # Horizontal line
                 tags$h3("Dataset conversion"),
                 
                 # File Input
                 fileInput(inputId = "up_df",
                           label = "Choose CSV File",
                           accept = c("text/csv","text/comma-separated-values,text/plain", ".csv")
                 ),
                 # Header
                 checkboxInput("header", "Header", TRUE),
                 # Separator
                 radioButtons("sep",
                              "Separator",
                              choices = c(Comma = ",",
                                          Semicolon = ";",
                                          Tab = "\t"),
                              selected = ","),
                 # Missing values
                 textInput("na", "Treat as missing"),
                 "Please separate strings with a comma",
                 tags$hr(), # Horizontal line
                 
                 # Select column
                 uiOutput("id_name"),
                 uiOutput("score"),
                 
                 # Save
                 tags$hr(), # Horizontal line
                 uiOutput("download")
                 
               ),
               mainPanel(
                 ref_article,
                 # Help
                 actionButton("help_button1", label = "Hide this page"),
                 uiOutput("help_conv"),
                 uiOutput("text_data"),
                 column(4, tableOutput("disp_data")), # Preview dataset
                 column(6, plotOutput("plot1")) # Preview plot
               )
             )
    ),
    
    tabPanel("Stratification",
             sidebarLayout(
               sidebarPanel(
                 # Cutoffs
                 radioButtons(inputId = "cutoffs",
                              label = "Cut-off",
                              choiceNames = c("Default", "Same cut-off between all strata", "Personalised cut-offs"),
                              choiceValues = c("Default", "Same", "Personalised"), 
                              selected = "Default"),
                 
                 # Slider inputs
                 uiOutput("sl1"),
                 uiOutput("sl2"),
                 uiOutput("sl3"),
                 uiOutput("sl4"),
                 uiOutput("sl5"),
                 
                 # Options heatmap
                 radioButtons(inputId = "disp_iga",
                              label = "Values on display",
                              choiceNames = c("Cut-offs only", "Cut-offs and IGA=0 and IGA=5", "Cut-offs and IGA"),
                              choiceValues = c("no", "only", "yes"),
                              selected = "only")
                 
               ),
               mainPanel(
                 ref_article,
                 plotOutput("strat"),
                 # Help button
                 actionButton("help_button2", label = "Hide this page"),
                 uiOutput("help_strat")
               )
             )
    )
  )
  
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Conversion one value
  output$out_conv <- renderText({
    req(input$inp1)
    if (input$conversion == "from oSCORAD to EASI"){
      paste("EASI =", pred_abacus(oscorad_to_easi, as.numeric(input$inp1)))}
    else{
      paste("oSCORAD =", pred_abacus(easi_to_oscorad, as.numeric(input$inp1)))
    }
  })
  
  # Read data
  df <- reactive({
    req(input$up_df)
    return(read.csv(input$up_df$datapath,
                    header = input$header,
                    sep = input$sep,
                    na.strings = c(strsplit(input$na, ",")[[1]], "", "NaN")))
  })
  
  # Display input data
  output$disp_data <- renderTable({
    df()
  })
  
  # Process data
  df2 <- reactive({
    req(df())
    req(input$ID_select)
    req(input$score_select)
    
    if (input$score_select == "None") {
      res <- data.frame()
    } else {
      if (input$conversion == "from oSCORAD to EASI") {
        res <- data.frame(oSCORAD = df()[, colnames(df()) == input$score_select])
        res$EASI <- pred_abacus2(oscorad_to_easi, res$oSCORAD)
      } else {
        res <- data.frame(EASI = df()[, colnames(df()) == input$score_select])
        res$oSCORAD <- pred_abacus2(easi_to_oscorad, res$EASI)
      }
      if (input$ID_select != "None") {
        res$ID <- df()[, colnames(df()) == input$ID_select]
        res <- res[, c(3, 1, 2)] # Reorder columns
      }
    }
    
    return(res)
  })
  
  # df2 valid ?
  isvalid <- reactive({
    req(df2())
    ncol(df2()) > 1   
  })
  
  # Display conversion
  output$plot1 <- renderPlot({
    req(isvalid())
    if (isvalid()) {
      if (input$conversion == "from oSCORAD to EASI") {
        plot_conv(df2(), "oSCORAD")
      } else {
        plot_conv(df2(), "EASI")
      }
    }
    
  })
  
  # ID selection
  output$id_name <- renderUI({
    selectInput("ID_select",
                "Select patient ID (optional):",
                c("None", colnames(df())))
  })
  
  # Input selection
  output$score <- renderUI({
    selectInput("score_select",
                "Select score:",
                c("None",colnames(df())))
  })
  
  # Download data button
  output$download <- renderUI({
    req(isvalid())
    if (isvalid()) {
      downloadButton("downloadData", "Save")
    }
  })
  
  # Download data
  output$downloadData <- downloadHandler(
    filename = function() {"converted_data.csv"},
    content = function(file) {
      write.csv(df2(), file, row.names = FALSE)
    }
  )
  
  # Change text help button 1
  observeEvent(input$help_button1, {
    if ((input$help_button1 + 1) %% 2) {
      btn_lbl <- "Hide this page"
    } else {
      btn_lbl <- "Show more"
    }
    updateActionButton(session, "help_button1", label = btn_lbl)
  })
  
  # Help conversion
  output$help_conv <- renderUI({
    if ((input$help_button1 + 1) %% 2) {
      tags$div(
        tags$h3("Conversion between oSCORAD and EASI"),
        tags$p("
From the patients' data that measures both oSCORAD and EASI, we can derive a relationship between oSCORAD and EASI. 
For example, from the data in",
               tags$a(href = "https://doi.org/10.1016/j.jaci.2011.03.024","Schram et al. JACI, 2011"),
               "plotted below (filled black circles, each corresponding to a patient), we can derive a relationship between oSCORAD and EASI (orange solid line).
               "),
        tags$img(src = "fit_EASI_oSCORAD.jpeg", width = "505px", height = "400px"),
        tags$p("
The root-mean-square error for the cross-validation was 4.77 and the accuracy, measured by the ratio of the EASI values predicted from oSCORAD being within a MCID (=6.6 for EASI), was 83% (Hurault et al., 2018).
               "),
        
        tags$h3("User Guide"),
        tags$p("This app executes conversion between oSCORAD and EASI."),
        tags$ol(
          tags$li("Select whether you would like to perform ",
                  tags$strong("Conversion \"from oSCORAD to EASI\""),
                  "(default) or \"",
                  tags$strong("from EASI to oSCORAD\".")),
          tags$li("\"",
                  tags$strong("Value-wise conversion"),
                  "\" provides the value that is converted from the value of either oSCORAD or EASI you put in the \"",
                  tags$strong("Severity score to be converted"),
                  "\" box."),
          tags$li(
            tags$strong("Dataset conversion"),
            " makes a plot that visualises the conversion of the entire dataset you upload, by the following steps:",
            tags$ol(
              tags$li("Upload your dataset as a .csv file, and check the preview of the dataset. 
                      oSCORAD is defined between 0 and 83 and EASI between 0 and 72. 
                      Please double-check that all the values in your dataset is within these ranges."),
              tags$li("Tick the box to indicate whether the dataset has a ",tags$strong("Header")," (cf. variable names)."),
              tags$li("Select the column ",tags$strong("Separator")," used in your dataset."), 
              tags$li("Enter the strings (e.g. NA, NaN, not done) that the program should ",
                      tags$strong("Treat as missing"),
                      ". This step is important for successful conversion."),
              tags$li("Select the column corresponding to the",tags$strong("score"),"oSCORAD or EASI in the dataset",
                      " and the optional column corresponding to ",
                      tags$strong("patient ID."),
                      "These buttons appears when a dataset is uploaded."),
              tags$li("
A plot will appear to visualise the conversion. You can also export it as a .csv file.
If not, please double-check what you did during the previous steps, especially what the strings that should treated as missing.
                      ")
            )
          )
        )
      )
    } else {
      return()
    }
  })
  
  output$text_data <- renderUI({
    req(df())
    tags$h1("Data")
  })
  
  # Slider inputs
  output$sl1 <- renderUI({
    if (input$cutoffs != "Default") {
      if (input$cutoffs == "Same"){
        lbl <- "Choose a cut-off"
      } else {
        lbl <- "\"Clear\" vs \"Almost clear\""
      }
      sliderInput(inputId = "sl1",
                  label = lbl,
                  value = 0.5,
                  min = 0,
                  max = 1)
    }
  })
  
  output$sl2 <- renderUI({
    if (input$cutoffs == "Personalised") {
      sliderInput(inputId = "sl2",
                  label = "\"Almost clear\" vs \"Mild\"",
                  value = 0.5,
                  min = 0,
                  max = 1)
    }
  })
  
  output$sl3 <- renderUI({
    if (input$cutoffs == "Personalised") {
      sliderInput(inputId = "sl3",
                  label = "\"Mild\" vs \"Moderate\"",
                  value = 0.5,
                  min = 0,
                  max = 1)
    }
  })
  
  output$sl4 <- renderUI({
    if (input$cutoffs == "Personalised") {
      sliderInput(inputId = "sl4",
                  label = "\"Moderate\" vs \"Severe\"",
                  value = 0.5,
                  min = 0,
                  max = 1)
    }
  })
  
  output$sl5 <- renderUI({
    if (input$cutoffs == "Personalised") {
      sliderInput(inputId = "sl5",
                  label = "\"Severe\" vs \"Very severe\"",
                  value = 0.5,
                  min = 0,
                  max = 1)
    }
  })
  
  ct <- reactive({
    res <- 0:4
    if (input$cutoffs == "Default"){
      res <- res + 0.5
    } else {
      if (input$cutoffs == "Same") {
        res <- res + input$sl1
      } else {
        res <- res + c(input$sl1, input$sl2, input$sl3, input$sl4, input$sl5)
      }
    }
    return(res)
  })
  
  xi <- reactive({
    if (input$disp_iga == "no"){
      res <- c()
    } else {
      if (input$disp_iga == "only"){
        res <- c(0, 5)
      } else {
        res <- 0:5
      }
    }
    return(res)
  })
  
  # Heatmap
  output$strat <- renderPlot({
    heatmap(ct(), xi())
  })
  
  # Change text help button
  observeEvent(input$help_button2, {
    if ((input$help_button2 + 1) %% 2) {
      btn_lbl <- "Hide this page"
    } else {
      btn_lbl <- "Show User Guide"
    }
    updateActionButton(session, "help_button2", label = btn_lbl)
  })
  
  # Help stratification
  output$help_strat <- renderUI({
    if ((input$help_button2 + 1) %% 2) {
      tags$div(
        withMathJax(),
        tags$h1("User Guide"),
        tags$p("
This app provides a stratification of EASI into 6 severity strata (clear, almost clear, mild, moderate, severe and very severe) defined by IGA. 
The correspondence between oSCORAD and EASI is defined by the conversion shown in the ",tags$strong("Conversion tab"),"."
        ),
        tags$p("
               IGA is normally defined as integers.
               Here we show probabilistic IGA (p-IGA) which is a continuous value ranging from 0 to 5.
               For example, p-IGA=2.2 means that the patient has 20% probability of IGA=3 (moderate) and 80% probability of IGA=2 (mild)."
        ),
        tags$p("The",tags$strong("cut-off values")," between 6 categories are defined as 0.5, 1.5, 2.5 ... by ",tags$strong("Default"),"."),
        tags$p("You can also",tags$strong("customise the cut-offs"), "to set the cut-off values of your choice."),
        
        tags$p("You can also change the ",tags$strong("Values on display")," for the plot to show:"),
        tags$ul(
          tags$li(tags$strong("Cut-offs only")),
          tags$li(tags$strong("Cut-offs and IGA=0 and IGA=5")),
          tags$li(tags$strong("Cut-offs and IGA"))
        )
      )
      
    } else {
      return()
    }
  })
  
}

# App ---------------------------------------------------------------------

shinyApp(ui = ui, server = server)
