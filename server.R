
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(GenomicRanges)
source("functions_intersect.R")

shinyFiles2load = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$files), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

shinyFiles2save = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$name), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

shinyServer(function(input, output, session) {
  roots_load <<- dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T)
  names(roots_load) <- basename(roots_load)
  shinyFileChoose(input, 'FilesLoadData', roots= roots_load, filetypes=c("narrowPeak", "broadPeak"))
  roots_output =  c("bed_dir" = "intersectR_beds")
  dir.create(roots_output, showWarnings = F)
  shinyFileSave(input, 'FilesSaveResults', roots= roots_output, filetypes=c("bed"))
  peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")
  
  peak_prev_file = reactiveVal(value = NULL, label = "peak_prev_file")
  peak_prev_name = reactiveVal(value = NULL, label = "peak_prev_name")
  peak_files = reactiveValues() 
  
  observeEvent(input$BtnAddFile, {
    if(is.null(peak_dfr())) return(NULL)
    peak_files[[input$TxtFileName]] = peak_dfr()
    curr_sel = isolate(input$ListIntervalSets)
    updateSelectInput(session, inputId = "ListIntervalSets", choices = isolate(names(reactiveValuesToList(peak_files))), c(curr_sel, input$TxtFileName), label = NULL)
    # print(object.size(isolate(reactiveValuesToList(peak_files))))
    peak_prev_file(NULL)
    peak_prev_name(NULL)
  })
  
  observeEvent(input$BtnCancelFile, {
    peak_prev_file(NULL)
    peak_prev_name(NULL)
  })
  
  observeEvent(input$FilesLoadData, {
    file_path = shinyFiles2load(input$FilesLoadData, roots_load)
    peak_prev_file(file_path)
    peak_prev_name(basename(file_path))
  })
  
  observeEvent(input$uploadPeakfile, {
    peak_prev_file(input$uploadPeakfile$datapath)
    peak_prev_name(input$uploadPeakfile$name)
  })
  
  observe({
    new_val =  peak_prev_name()
    if(is.null(new_val)) new_val = ""
    print(paste('update:', new_val))
    updateTextInput(session, "TxtFileName", value = new_val)
  })
  
  output$chooser = renderUI({
    pf = names(reactiveValuesToList(peak_files))
    print(pf)
    return(chooserInput(inputId = "ChooseIntervalSets", 
                        leftLabel = "Ready", rightLabel = "For Analysis", 
                        leftChoices = character(), rightChoices = pf, size = 8, multiple = T))
  })
  
  peak_dfr = reactive({
    if(is.null(peak_prev_file())) return(NULL)
    df = read.table(peak_prev_file(), stringsAsFactors = F)
    colnames(df) = peak_cn
    return(df)
  })
  
  output$peaksHeader = DT::renderDataTable({
    if(is.null(peak_dfr())){
      m = matrix(0, ncol = length(peak_cn), nrow = 0)
      colnames(m) = peak_cn
      return(as.data.frame(m))
    }
    df = peak_dfr()
    # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
    DT::datatable(df)#, caption = paste0("Preview of", basename(peak_prev_file()), ":"))
  })
  
  observe({
    print(input$StrategyRadio)
  })
  
  output$peaksFilter = renderUI({
    if(is.null(peak_dfr())) 
      return(    
        # fixedRow(
        #   selectInput("FilterColname", label = "", choices = "waiting..."),
        #   selectInput("FilterOperator", label = "", choice = c("==", "!=", ">", "<", ">=", "<=")),
        #   textInput("FilterValue", label = "")
        # )
        NULL
      )
    fixedRow(
      column(selectInput("FilterColname", label = "", choices = colnames(peak_dfr())), width = 4),
      column(selectInput("FilterOperator", label = "", choice = c("==", "!=", ">", "<", ">=", "<=")), width = 2),
      column(textInput("FilterValue", label = ""), width = 6)
    )
  })
  
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
  observeEvent(input$BtnQuickFlat, handlerExpr = {
    load("peak_gr.save")
    for(i in 1:length(peak_gr)){
      peak_files[[names(peak_gr)[i]]] = peak_gr[[i]]
    }
    
    updateSelectInput(session, inputId = "ListIntervalSets", choices = isolate(names(reactiveValuesToList(peak_files))), label = NULL)
    # print(object.size(isolate(reactiveValuesToList(peak_files))))
    # peak_prev_file(NULL)
    # peak_prev_name(NULL)
  })
  
  gr_to_plot = reactiveVal(label = "gr_to_plot")
  
  observeEvent(input$BtnAnalyze, handlerExpr = {
    print(input$StrategyRadio)
    to_anlayze = input$ChooseIntervalSets$right
    names(to_anlayze) = to_anlayze
    peak_df = reactiveValuesToList(peak_files)
    peak_df = peak_df[to_anlayze]
    peak_gr = lapply(peak_df, GRanges)
    
    gr_to_plot(intersectR(peak_gr, use_first = input$StrategyRadio == "serial"))
  })
  
  output$AnalysisPlot = renderPlot({
    if(is.null(gr_to_plot())) return(NULL)
    gr = gr_to_plot()
    nam = names(reactiveValuesToList(peak_files))
    df = as.data.frame(elementMetadata(gr))
    tp = which(colnames(df) != "group")
    
    switch(input$RadioPlotType,
           "pie" = {
             print("plot: pie")
             pie(table(gr$group))
           },
           "venn" = {
             print("plot: venn")
             limma::vennDiagram(df[,tp])
           },
           "heatmap" = {  
             print("plot: heatmap")
             mat = ifelse(as.matrix(df[,tp]), 1, 0)
             mat = mat + runif(length(mat), min = -.02, .02)
             gplots::heatmap.2(mat[sample(nrow(mat), 2000),]+.2, col = c("white", "black"), trace = "n")
           })
  })
  
  observeEvent(input$FilesSaveResults, {
    if(is.null(input$FilesSaveResults)) return(NULL)
    print(input$FilesSaveResults)
    file_path = shinyFiles2save(input$FilesSaveResults, roots_output)
    print(file_path)
    if(!grepl(".bed$", file_path)){
      file_path = paste0(file_path, ".bed")
    }
    write.table(gr_to_plot(), file = file_path, row.names = F, col.names = T, quote = F)
  })
})