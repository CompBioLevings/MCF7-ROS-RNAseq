#
# This is a Shiny web application for looking through the Slattery Lab's published MCF7
# RNA-seq experiment results for ROS (MEN/tBOOH) exposure time course
#

# Load necessary apps
library("shiny", quietly = TRUE)

# Load in the data for the app
load(file="MCF7-ROS-RNAseq.RData")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # App title ----
    titlePanel(title = "MCF-7 ROS RNA-seq experiment (2017)"),
    
    # Add sidebar
    sidebarLayout(position = "left",
                  
        # Sidebar panel for inputs ----
        sidebarPanel(strong("Input options"), width = 3,
            div("-----------------------------------------------------", style = "margin-bottom: 10px; font-weight: bold;"),
          
            # Input: Text gene "symbol" of gene to be plotted  ----
            div("Please input one of the following:", style = "margin-top: 10px; margin-left: 2px, margin-bottom: 0px; font-weight: bold;"),
            div("1) an official gene symbol", style = "margin-top: 1px; margin-left: 3px; margin-bottom: 0px; font-weight: bold;"),
            div("2) a RefSeq transcript identifier (NM_#)", style = "margin-top: 1px; margin-left: 3px; margin-bottom: 1px; font-weight: bold;"),
            textInput(inputId = "symbol", label = div("3) an Ensembl gene ID (ENSG#)", style = "margin-left: 3px; margin-top: 1px; font-weight: bold;"), value = "NQO1"),
            
            # Input: Which treatments to plot
            div("Choose treatments to plot:", style = "margin-top: 15px; margin-bottom: 0px; font-weight: bold;"),
            sliderInput(inputId = "treatment",
                        label = "(1 = both, 2 = MEN, 3 = tBOOH)",
                        min = 1,
                        max = 3,
                        value = 1,
                        round = TRUE),
            
            # Input: Time points to be plotted/x-window  ----
            div("Time point 'window' in hours", style = "margin-top: 15px; margin-bottom: 0px; font-weight: bold;"),
            textInput(inputId = "times", label = "(input 2 comma-separated numbers)", value = "0, 24"),
            
            # Input: Gene abundance y-limits for plotting  ----
            div("Limits for y-axis", style = "margin-top: 15px; margin-bottom: 0px; font-weight: bold;"),
            textInput(inputId = "ylimits", label = "(input 2 comma-separated numbers)", value = ""),
            
            
            # Input: Checkboxes for whether to use batch-corrected values
            checkboxInput(inputId = "correct_batch", label = strong("Use batch-corrected values for rep"), value = TRUE),
            
            # Input: Checkboxes for whether to display replicate label
            checkboxInput(inputId = "include_rep", label = strong("Display replicate on plot"), value = FALSE),
            
            # Input: Number of breaks for y-axis of plot 1 (left)
            sliderInput(inputId = "plot_breaks",
                        label = "Plot y-axis breaks:",
                        min = 5,
                        max = 30,
                        value = 10,
                        round = TRUE),
        ),
       
        # Main panel for displaying outputs ----
        mainPanel(
          tabsetPanel(
              tabPanel(fluidRow("",
                                column(6,plotOutput(outputId="plotgraphs", width="900px", height="650px"))
              ),
              fluidRow("",
                       column(6,tableOutput(outputId = "table"))
              )
              )
          )
        )
))

# Define server logic required to draw time course plot and statistics table
server <- function(input, output) {
    
    # Load necessary packages and data
    library("ggplot2", quietly = TRUE)
    library("ggrepel", quietly = TRUE)
    library("ggtext", quietly = TRUE)
    library("scales", quietly = TRUE)
    library("gridExtra", quietly = TRUE)
    library("xtable", quietly = TRUE)
    library("magrittr", quietly = TRUE, warn.conflicts = FALSE)
    library("tidyverse", quietly = TRUE, warn.conflicts = FALSE)
    library("shiny", quietly = TRUE)
    
    # Function for drawing the plot of gene expression
    gene_expr_plot <- function(genename = NA, treat=1, break.num=10, 
                               batchcorrect = TRUE, includerep = FALSE, timepoints = NA, limits = "") {
      
      if (batchcorrect) {
        expression_df <- batch_corrected_TPM
      } else {
        expression_df <- TPM_df
      }
      
      # Turn off warnings
      oldw <- getOption("warn")
      options(warn = -1)
      
      # Convert gene name/symbol to Ensembl ID to pass to RNA seq data, then import that data
      geneid <- as.character(gene_to_id_map[grep(x=as.character(gene_to_id_map$symbol), pattern = paste0("^",genename,"$"), ignore.case = T, perl = T),1])
      
      # Import Ethanol and either Menadione or tBOOH data- and add 0 hr timepoints for these
      if (treat == 2) {
        genedata <- subset(x = expression_df, subset = (expression_df$gene_id %in% geneid), select = which(colnames(expression_df) %in% 
             c('Ctrl_0hr_B', 'EtOH_1hr_B', 'EtOH_8hr_B', 'EtOH_24hr_B', 'Men_1hr_B', 'Men_8hr_B', 
               'Men_24hr_B', 'Ctrl_0hr_F', 'EtOH_1hr_F', 'EtOH_8hr_F', 'EtOH_24hr_F', 
               'Men_1hr_F', 'Men_8hr_F', 'Men_24hr_F')))
        genedata <- cbind.data.frame(colnames(genedata), t(genedata)) %>%  
          set_colnames(c("sample", "TPM"))
        genecolors <- c("royalblue2", "firebrick2")
      } else if (treat == 3) {
        genedata <- subset(x = expression_df, subset = (expression_df$gene_id %in% geneid), select = which(colnames(expression_df) %in% 
             c('Ctrl_0hr_B', 'EtOH_1hr_B', 'EtOH_8hr_B', 'EtOH_24hr_B', 'TBOOH_1hr_B', 'TBOOH_8hr_B', 
               'TBOOH_24hr_B', 'Ctrl_0hr_F', 'EtOH_1hr_F', 'EtOH_8hr_F', 'EtOH_24hr_F', 
               'TBOOH_1hr_F', 'TBOOH_8hr_F', 'TBOOH_24hr_F')))
        genedata <- cbind.data.frame(colnames(genedata), t(genedata)) %>%  
          set_colnames(c("sample", "TPM"))
        genecolors <- c("royalblue2", "chocolate2")
      } else {
        genedata <- subset(x = expression_df, subset = (expression_df$gene_id %in% geneid), select = -1)
        genedata <- cbind.data.frame(colnames(genedata), t(genedata)) %>%  
          set_colnames(c("sample", "TPM"))
        genecolors <- c("royalblue2", "firebrick2", "chocolate2")
      }
      
      # Now join the data above with the sample characteristics
      genedata %<>% left_join(y=sampleTable[,c("sample", "treatment", "time", "rep")])
      
      # Duplicate control timepoint and rename so each treatment starts at time 0
      if (treat == 2) {
        genedata <- rbind.data.frame(genedata, genedata %>% filter(time == 0) %>% 
             mutate(treatment = "MEN", rep = NA))
      } else if (treat == 3) {
        genedata <- rbind.data.frame(genedata, genedata %>% filter(time == 0) %>% 
             mutate(treatment = "tBOOH", rep = NA))
      } else {
        genedata <- rbind.data.frame(genedata, genedata %>% filter(time == 0) %>% 
             mutate(treatment = "MEN", rep = NA), genedata %>% filter(time == 0) %>% 
            mutate(treatment = "tBOOH", rep = NA))
      }
      
      # Factorize time variable
      genedata$time <- factor(genedata$time, levels = c(0, 1, 8, 24))
      
      # Arrange data order
      genedata %<>% arrange(treatment, time)
      
      # Factorize treatment variables
      genedata$treatment <- factor(genedata$treatment, levels = unique(genedata$treatment))
      
      # Get the breaks I should use based on the max data values
      TPM.breaks <- seq(0, max(genedata$TPM)+(max(genedata$TPM)-min(genedata$TPM)), by = 
                          as.numeric(signif(x = max(genedata$TPM)+(max(genedata$TPM)-min(genedata$TPM)), digits = 1))/break.num)
      TPM.breaks <- round(x = TPM.breaks, digits = ifelse(nchar(trunc(max(TPM.breaks, na.rm=TRUE))) < 3, 2, 0))
      
      #  Now plot it
      q <- ggplot(genedata, mapping = aes(x=as.numeric(as.character(time)), y=TPM, color=treatment, group=treatment, label = rep))
      q <- q + geom_point(size=2) + geom_smooth(se=TRUE, method="loess", alpha=0.15) + scale_y_continuous(breaks=TPM.breaks)
      if (includerep) {
        q <- q + geom_text_repel(color="black", size=5, alpha=0.7)
      } else {
      }
      
      
      # Now set the borders/limits of the plot
      # X axis (time) first
      if (is.na(timepoints)) {
        timepoints <- c(0, 24)
      } else if (is.character(timepoints)) {
        timepoints <- as.numeric(trimws(as.character(unlist(strsplit(timepoints, split = ",")[1:2]))))
        if (all(is.numeric(timepoints))){
        } else {
          stop("Invalid y limits specified.  Please type in only two comma-delimited numbers: eg- 1, 100")
        }
      } else {
        stop("Invalid y limits specified.  Please type in only two comma-delimited numbers: eg- 1, 100")
      }
      
      # Now Y axis (gene abundance/TPM) limits
      if (nchar(as.character(limits)) == 0) {
        q <- q + scale_color_manual(values=genecolors) + coord_cartesian(xlim=timepoints, 
             ylim=c(0+(max(genedata$TPM,na.rm=TRUE)+max(genedata$TPM,na.rm=TRUE)*0.25)*0.1, 
             (max(genedata$TPM,na.rm=TRUE)+max(genedata$TPM,na.rm=TRUE)*0.25)), expand = F)
      } else if (nchar(as.character(limits)) > 0) {
        limits <- as.numeric(trimws(as.character(unlist(strsplit(limits, split = ",")[1:2]))))
        if (all(is.numeric(limits))){
          q <- q + scale_color_manual(values=genecolors) + coord_cartesian(xlim=timepoints, 
                ylim=c(limits[1], limits[2]), expand = F)
        } else {
          stop("Invalid y limits specified.  Please type in only two comma-delimited numbers: eg- 1, 100")
        }
      } else {
        stop("Invalid y limits specified.  Please type in only two comma-delimited numbers: eg- 1, 100")
      }
      
      if (any(geneid %in% all.DEGs.final.stringent)){
        q <- q + labs(y = "Gene expression level (TPM)", x="Time (hrs)", 
            title = paste0("__Expression change of <span style = 'color:#FF0000;'>*", paste0(gene_to_id_map %>% 
           filter(gene_id == geneid, type == "symbol") %>% dplyr::select(symbol) %>% unlist %>% unique(), collapse = ';'),
       "</span>* in ", ifelse(treat == 2, "MEN", ifelse(treat == 3, "tBOOH", "MEN/tBOOH")), 
       "__<br> <span style = 'font-size:12pt;'> \\*Significant DEGs (p < 0.05) are highlighted in </span><span style = 'color:#FF0000; font-size:12pt;'>red</span>")) + 
          theme_bw() 
      } else {
        q <- q + labs(y = "Gene expression level (TPM)", x="Time (hrs)", 
            title = paste0("__Expression change of *", paste0(gene_to_id_map %>% 
          filter(gene_id == geneid, type == "symbol") %>% dplyr::select(symbol) %>% unlist %>% unique(), collapse = ';'),
           "* in  ", ifelse(treat == 2, "MEN", ifelse(treat == 3, "tBOOH", "MEN/tBOOH")), "__")) + 
          theme_bw() 
      }
      q <- q + theme(plot.title = element_textbox_simple(size=18, lineheight = 1, halign = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)), 
                     axis.text = element_text(color="black", size=14),
                     axis.title.x = element_text(face="bold", color="black", size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
                     axis.title.y = element_text(face="bold", color="black", size=16, margin = margin(t = 0, r = 7.5, b = 0, l = 0)), 
                     axis.ticks = element_line(color="black"),
                     panel.border = element_rect(size=0.5, color="black"), 
                     legend.text = element_text(size = 12), 
                     legend.title = element_text(size=14, face="bold")
      )
      
      return(q)
      
      # turn warnings back on
      options(warn = oldw)
    }
    
    # Function for getting and plotting table of statistical values
    get_stats <- function(genename = NULL, treat = 1) {
      
        # Convert gene name/symbol to Ensembl ID to pass to RNA seq data, then import that data
        geneid <- as.character(unlist(gene_to_id_map[grep(x=as.character(unlist(gene_to_id_map$symbol)), 
                                                          pattern = paste0("^",genename,"$"), ignore.case = T, perl = T),1]))
        
        # Now get only the relevant
        results.data <- rbind.data.frame(results.EtOH.1hr %>% mutate(grouping = "EtOH 1hr vs 0hr"), 
             results.EtOH.8hr %>% mutate(grouping = "EtOH 8hr vs 0hr"), results.EtOH.24hr %>% mutate(grouping = "EtOH 24hr vs 0hr"), 
             results.MEN.1hr %>% mutate(grouping = "MEN 1hr vs 0hr"), results.MEN.8hr %>% mutate(grouping = "MEN 8hr vs 0hr"), 
             results.MEN.24hr %>% mutate(grouping = "MEN 24hr vs 0hr"), results.tBOOH.1hr %>% mutate(grouping = "tBOOH 1hr vs 0hr"), 
             results.tBOOH.8hr %>% mutate(grouping = "tBOOH 8hr vs 0hr"), results.tBOOH.24hr %>% mutate(grouping = "tBOOH 24hr vs 0hr")) %>%
          filter(gene_id == geneid) %>% rowwise %>% mutate(log2FC_vs_0 = round(log2FC_vs_0, digits = 8),
              lfcSE = round(lfcSE, digits = 8)) %>% rbind.data.frame(results.LRT %>% filter(gene_id == geneid) %>% 
               mutate(grouping = "LRT - MEN & tBOOH vs EtOH/0hr", log2FC_vs_0 = NA, lfcSE = NA) %>% 
               dplyr::select(c(gene_id, log2FC_vs_0, lfcSE, pval, FDR_p, grouping))) %>% as.data.frame()
        
        if (treat == 2) {
          results.data %<>% filter(grouping %in% as.character(paste0(c("EtOH 1hr", "EtOH 8hr", "EtOH 24hr", "MEN 1hr", 
              "MEN 8hr", "MEN 24hr"), " vs 0hr")))
        } else if (treat == 3) {
          results.data %<>% filter(grouping %in% paste0(c("EtOH 1hr", "EtOH 8hr", "EtOH 24hr", "tBOOH 1hr", 
              "tBOOH 8hr", "tBOOH 24hr"), " vs 0hr"))
        } else {
        }
        
        # rename for display
        colnames(results.data)[6] <- "Comparison"
        
        # Print results data
        return(results.data %>% dplyr::select(Comparison, log2FC_vs_0, lfcSE, pval, FDR_p))
    }
    
    set.seed(123)
    # Make plot of MCF7 ROS RNA-seq experiment results for gene of interest
    pt1 <- reactive({
        plot(gene_expr_plot(genename=input$symbol, break.num=input$plot_breaks, treat=input$treatment,
            batchcorrect=input$correct_batch, includerep=input$include_rep, timepoints=input$times,
            limits=input$ylimits))
    })
    
    # Output the plot
    output$plotgraphs <- renderPlot(res = 100, {
        ptlist <- list(pt1())
        wtlist <- c(6)
        grid.arrange(grobs=ptlist,widths=wtlist,ncol=length(ptlist))
    })
    
    # Make combined datatable to render, with results
    output$table <- renderTable(digits = 3, {
        values <- get_stats(genename=input$symbol, treat = input$treatment) %>% 
            set_colnames(c("Experimental_comparison________", "log2(fold-change)", 
            "SE of log2FC", "pvalue________", "FDR p_________"))
        values$`log2(fold-change)` <- format(values$`log2(fold-change)`, digits = 5)
        values$`SE of log2FC` <- format(values$`SE of log2FC`, digits = 5)
        values$`pvalue________` <- format(values$`pvalue________`, scientific=TRUE)
        values$`FDR p_________` <- format(values$`FDR p_________`, scientific=TRUE)
        return(values) 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
