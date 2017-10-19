
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(GenomicRanges)
source("functions_intersect.R")

bed_path = "~/ShinyApps/shiny_peak_data/beds"

# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
dataModal <- function(sets, failed = FALSE) {
  sets_html = HTML(paste(sapply(sets$selected[1], function(x){
    as.character(textInput("TxtRename", label = paste("rename", x), value = x))
  }), collapse = "\n"))
  modalDialog(
    sets_html,
    span('(Please rename the selected sample)'),
    if (failed)
      div(tags$b("One or more names no longer unique!", style = "color: red;")),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("BtnConfirmRename", "Confirm")
    )
  )
}

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
  roots_output =  c("bed_dir" = bed_path)
  dir.create(roots_output, showWarnings = F)
  shinyFileSave(input, 'FilesSaveResults', roots= roots_output, filetypes=c("bed"))
  peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")
  peak_prev_file = reactiveVal(value = "", label = "peak_prev_file")
  peak_prev_name = reactiveVal(value = "", label = "peak_prev_name")
  #stores all set data, keyed by set_id
  peak_files = reactiveVal(value = list(), label = "peak_files")
  #stores the order of peak_files as well as any updated names
  peak_names = reactiveVal(value = data.frame(set_id = character(), display_name = character()))
  
  observeEvent(input$BtnAddFile, {
    if(is.null(peak_dfr())) return(NULL)
    tmp = peak_files()
    tmp[[input$TxtFileName]] = peak_dfr()
    peak_files(tmp)
    sets = input$ChooseIntervalSets
    output$SetChooser = renderUI({
      pfl = sets$left
      #add to active
      pfr = c(sets$right, isolate(input$TxtFileName))
      
      # print(pf)
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    peak_prev_file("")
    peak_prev_name("")
  })
  
  observeEvent(peak_files, {
    if(length(peak_files()) == 0){
      output$SetChooser = renderUI({
        pfl = character()
        #add to active
        pfr = character()
        
        # print(pf)
        return(chooserInput(inputId = "ChooseIntervalSets", 
                            leftLabel = "Ready", rightLabel = "For Analysis", 
                            leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
      })
    }
  })
  
  observeEvent(input$BtnCancelFile, {
    peak_prev_file("")
    peak_prev_name("")
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
  
  #whenever peak_prev_name changes (new file selected or file added/cancelled) update text box.
  observeEvent(peak_prev_name(), {
    new_val =  peak_prev_name()
    if(is.null(new_val)) new_val = ""
    print(paste('update:', new_val))
    updateTextInput(session, "TxtFileName", value = new_val)
  })
  
  observeEvent(input$BtnDeleteSet, {
    sets = input$ChooseIntervalSets
    to_del = sets$selected
    
    sets$left = setdiff(sets$left, to_del)
    sets$right = setdiff(sets$right, to_del)
    output$SetChooser = renderUI({
      pfl = sets$left
      pfr = sets$right
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    
  })
  observeEvent(input$BtnRenameSet, {
    sets = input$ChooseIntervalSets
    if(is.null(sets)) return(NULL)
    if(length(sets$selected) == 0) return(NULL)
    # print(sets)
    showModal(dataModal(sets = sets))
  })
  
  # When OK button is pressed, attempt to load the data set. If successful,
  # remove the modal. If not show another modal, but this time with a failure
  # message.
  observeEvent(input$BtnConfirmRename, {
    sets = input$ChooseIntervalSets
    old_name = sets$selected[1]
    new_name = input$TxtRename
    print(paste("rename", old_name, "-->", new_name))
    if(old_name == new_name){
      removeModal()
    }else if(any(c(unique(c(sets$left, sets$right)) == new_name))){
      showModal(dataModal(failed = TRUE))
    }else{
      tmp = peak_files()
      names(tmp)[names(tmp) == old_name] = new_name
      peak_files(tmp)
      
      
      sets$left[sets$left == old_name] = new_name
      sets$right[sets$right == old_name] = new_name
      output$SetChooser = renderUI({
        pfl = sets$left
        pfr = sets$right
        return(chooserInput(inputId = "ChooseIntervalSets", 
                            leftLabel = "Ready", rightLabel = "For Analysis", 
                            leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
      })
      removeModal()
    }
  })
  
  output$NumericMergeExtensionOut = renderUI({
    if(!is.null(input$SliderMergeExtension))
      num = input$SliderMergeExtension
    numericInput(inputId = "NumericMergeExtension", label = "", min = 0, max = Inf, value = num)
  })
  
  peak_dfr = reactive({
    if(is.null(peak_prev_file())) return(NULL)
    if(peak_prev_file() == "") return(NULL)
    df = read.table(peak_prev_file(), stringsAsFactors = F)
    if(ncol(df) == length(peak_cn)){
      showNotification("assuming file is narrowPeak.", type = "warning")
      colnames(df) = peak_cn  
    }else{
      showNotification("file not narrowPeak, loading as minimal bed file.", type = "warning")
      bed_cn = peak_cn[1:4]
      nc = min(ncol(df), length(bed_cn))
      colnames(df)[1:nc] = peak_cn[1:nc]
    }
    
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
    DT::datatable(df, 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    pageLength = 5), rownames = F)
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
      tmp = peak_files()
      tmp[[names(peak_gr)[i]]] = peak_gr[[i]]
      peak_files(tmp)
    }
    
    sets = input$ChooseIntervalSets
    output$SetChooser = renderUI({
      pfl = sets$left
      #add to active
      pfr = unique(c(sets$right, names(peak_gr)))
      
      # print(pf)
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    # peak_prev_file(NULL)
    # peak_prev_name(NULL)
  })
  
  gr_to_plot = reactiveVal(label = "gr_to_plot")
  
  # observeEvent(input$BtnAnalyze, handlerExpr = {
  #   print(input$StrategyRadio)
  #   to_anlayze = input$ChooseIntervalSets$right
  #   names(to_anlayze) = to_anlayze
  #   peak_df = peak_files()
  #   peak_df = peak_df[to_anlayze]
  #   peak_gr = lapply(peak_df, GRanges)
  #   print(names(peak_gr))
  #   gr_to_plot(intersectR(peak_gr, use_first = input$StrategyRadio == "serial", ext = input$NumericMergeExtension))
  # })
  # 
  
  #insulate other elements from needless updates when ChooseIntervalSets not really updated
  observeEvent(input$ChooseIntervalSets, {
    new_anlayze = input$ChooseIntervalSets$right
    names(new_anlayze) = new_anlayze
    was_analyze = to_analyze()
    if(length(was_analyze) != length(new_anlayze)){
      to_analyze(new_anlayze)
    }else if(!all(was_analyze == new_anlayze)){
      to_analyze(new_anlayze)
    }
    
  })
  
  to_analyze = reactiveVal(character())
  
  output$AnalysisPlot = renderPlot({
    if(is.null(input$NumericMergeExtension) || 
       is.null(isolate(input$ChooseIntervalSets))) return(NULL)
    if(length(isolate(input$ChooseIntervalSets$right)) == 0) return(NULL)
    # print(input$StrategyRadio)
    to_analyze_ = to_analyze()
    names(to_analyze_) = to_analyze_
    peak_df = peak_files()
    peak_df = peak_df[to_analyze_]
    peak_gr = lapply(peak_df, GRanges)
    print(names(peak_gr))
    gr_to_plot(intersectR(peak_gr, use_first = input$StrategyRadio == "serial", ext = input$NumericMergeExtension))
    if(is.null(gr_to_plot())) return(NULL)
    gr = gr_to_plot()
    nam = names(peak_files())
    df = as.data.frame(elementMetadata(gr))
    tp = which(colnames(df) != "group")
    input$StrategyRadio
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
